import numpy as np
np.seterr(divide='ignore', invalid='ignore')
import scipy.sparse as sp 
from scipy import signal
from scipy.interpolate import interp1d
import pandas as pd 
from copy import copy 
import cooler 
from cooltools.lib.numutils import logbins
import sys 
from mirnylib.numutils import observedOverExpected, logbinsnew, zoomArray
from multiprocessing import Pool


def chr_cooler_format(c,df):
    c_chr = bool(np.all(np.array(['chr' in name for name in c.chromnames])))
    df_chr = bool(np.all(np.array(['chr' in name for name in df['chr1'].unique()])))

    df['chr1'] = df['chr1'].apply(lambda x: 'chr'*c_chr+x)
    df['chr2'] = df['chr2'].apply(lambda x: 'chr'*c_chr+x)
    return df


def fetchCooler(c, regions, coolerFetch = lambda coo, ext:coo.matrix(balance=True, sparse=True).fetch(ext),
               mask=True,  force=False, ):
    """    
    Takes cooler, a set of regions (chromosomes, chromosomal arms, regions) 
    For each region, it fetches a sparse cooler matrix.  Returns (as a generator) a list of (matrix, bins) pairs 
    
    Parameters
    ----------    
    c: cooler  -  cooler to use 
    regions: list of (chrom, start, end)  
        None is interpreted as start or end of a chromosome
        Start must be a multiple of resolution 
    coolerFetch: a function to fetch an region from cooler (optional) 
        If you want to pileup raw matrices rather than balanced, change this function accordingly
    mask: True, False,None, array, or a function 
        A mask of bins that should be included in the pileup 
        if a function, it is called as a function of a matrix for a given region 
        If an array, ones should correspond to the bins that are "valid", and the array should have the same length as bins
        if True use the default mask, which is the bins with nonzero marginals in a given region (not genome-wide!) 
        if None, assumes that cooler.bins() already has column named "mask" 
        if False, do not use the mask at all (include all bins) 
    
    Returns 
    -------
        A generator returning (matrix, bins) pairs 
        
    """
    regions = [list(i) for i in regions]
    resolution = c.binsize

    for i in regions:
        if i[1] == None:
            i[1] = 0 
        if i[2] == None:
            i[2] = c.chromsizes[i[0]]

  
    for a in regions: 
        if str(a[0]) not in c.chromnames:
            raise ValueError("Chromosome {0} from regions not found in cooler".format(a))
        if (a[1] % resolution) != 0:
            raise ValueError("Start of an region should be a multiple fo resolution")
    
#     bins = c.bins()[:]
    
#     # managing masks 
#     if mask is False: 
#         bins["mask"] = 1 
#     elif mask is None:
#         assert "mask" in bins.columns
#     elif mask is True: 
#         pass 
#     elif callable(mask):
#         pass 
#     else:
#         bins["mask"] = mask 
   
    
    for region in regions:
        matrix = coolerFetch(c, region)
        try: # setting matrix nans to zeros.
            matrix.data = np.nan_to_num(matrix.data, copy=False)
        except TypeError:  #workaround for old numpy versions
            matrix.data = np.nan_to_num(matrix.data)
#         st,end = c.extent(region)
#         subbins = bins[st:end].copy()
        if mask is True: 
            newmask = np.array((matrix.sum(axis=0) > 0 ))[0]
#         if callable(mask):
#             new_mask = mask(matrix)
#         subbins["mask"] = newmask 

        assert len(newmask) == matrix.shape[0]

        yield matrix, newmask

def chunkDataFrame(c, regions, positionDataFrame, 
               columns=[("chrom1","pos1"), ("chrom2","pos2")], force=False, ): 
    """
    Takes cooler, a set of regions (chromosomes, chromosomal arms, regions), and a dataframe. 
    For each region, it converts (chr1, pos1) columns of the dataframe 
    (specified in the "columns" argument) to a column named "ind1", which is an index of this position 
    in the final matrix. Returns a dict of {region:dataframe} 
    
    Parameters
    ----------    
    c: cooler  -  cooler to use 
    regions: list of (chrom, start, end)  
        None is interpreted as start or end of a chromosome
        Start must be a multiple of resolution 
    positionDataFrame: pd.DataFrame containing positions of genomic elements as pairs of columns (chromosome end)
        This dataframe can contain any number of positions. For example, for averaging loops, you may need 
        two positions: start and end of the loop. For each position, you need a chromosome and a genomic 
        coordinate columns in the dataframe. 
        These columns are specified in the next argument 
    columns: list of length-2 lists/tuples 
        These are columns to use for pileups. Columns may repeat (e.g. [('chr','pos1'),('chr','pos2')])
        For the first entry of this list, a column "ind1" will be created. For the second - "ind2". 

    Returns 
    -------
        A dict {region:dataframe} for each region  
        
    """
    originalRegions = [tuple(i) for i in regions]
    regions = [list(i) for i in regions]
    resolution = c.binsize
    
    positionDataFrame = positionDataFrame.copy()
    for i in regions:
        if i[1] == None:
            i[1] = 0 
        if i[2] == None:
            i[2] = c.chromsizes[i[0]]
            
    columns = list(columns)
    firstChr = columns[0][0]  # first column specifying a chromosomes 
    for a in columns: # check that all columns are in the dataframe
        for b in a:
            if b not in positionDataFrame:
                raise ValueError("Column {0} not found in dataframe".format(b))
    for a in columns[1:]:   # check that all pairs are cis 
        if not (positionDataFrame[a[0]] == positionDataFrame[firstChr]).all():
            raise ValueError("Only cis pileups are currently supported")
    
    if not force: # check that dataframe has chromosomes that are in cooler ('not' added by Sameer - 3/18)
        for a in positionDataFrame[firstChr].astype(str).unique():  
            if a not in c.chromnames:
                raise ValueError("Chrom. {0} from Dataframe notin cooler (force=True to override)".format(a))    
    for a in regions: 
        if str(a[0]) not in c.chromnames:
            raise ValueError("Chromosome {0} from regions not found in cooler".format(a))
        if (a[1] % resolution) != 0:
            raise ValueError("Start of an region should be a multiple fo resolution")
    for a in columns: 
        positionDataFrame[a[0]] = positionDataFrame[a[0]].astype(str)
    grouped = positionDataFrame.groupby(firstChr)
    bins = c.bins()[:]
    
    result = {}
    for orig_region, region in zip(originalRegions, regions):
        frame = grouped.get_group(region[0])
        
        framemask = np.ones(len(frame), dtype = np.bool)
        for column in columns: # Selects for parts of the dataframe that are in the region. First selects chr1 then chr2.
            pos = frame[column[1]].values
            framemask[pos < region[1]] = False
            framemask[pos >= region[2]] = False
        frame = frame.iloc[framemask]
        for j,column in enumerate(columns): # Creates associated indices for the matrices. 
            indColumn = (frame[column[1]] - region[1] )// c.binsize
            frame["ind{0}".format(j+1)] = indColumn
        result[orig_region] = frame 
    return result


def sparseSymmetricOOE(matrix, mask):
    '''Quick OOE operation for sparse symmetric matrices. This will be used by the LocalObsOverExp object to 
    compute OOE on support regions.'''
    
    if matrix.shape[0] == 1:
        return matrix
    
    #Finding number of valid bins per diagonal using FFT convolve
    count_per_diag = signal.fftconvolve(mask, mask[::-1], mode='full')
    count_per_diag = np.round(count_per_diag[len(count_per_diag)//2:])
    count_per_diag = count_per_diag.astype(int)
    
    row, col, data = matrix.row, matrix.col, matrix.data
    nan_indices = ~np.isfinite(data)
    data[nan_indices]=0
    diff = abs(row-col)
    
    #Summing by diagonal
    scaling = np.bincount(diff, weights=data, minlength=len(count_per_diag))/2
    assert len(scaling)==len(count_per_diag)
    
    
    bins = logbins(1,len(scaling),ratio=1.2)
    bin_mids = np.sqrt(bins[1:]*bins[0:-1])
    lab = np.concatenate(tuple((i+1)*np.ones(bins[i+1]-bins[i], dtype=int) for i in range(len(bins)-1)))
    log_scaling = np.bincount(lab,weights=scaling[1:])
    log_count = np.bincount(lab, weights=count_per_diag[1:])
    coarse_expected = log_scaling[1:]/log_count[1:]
    
    f = interp1d(np.log10(bin_mids), np.log10(coarse_expected), kind='linear')
    y = f(np.log10(np.arange(2,np.floor(bin_mids[-1]))))
    x = np.log10(np.arange(2,np.floor(bin_mids[-1])))
    
    xremaining = np.log10(np.arange(np.round(10**x[-1]+1),len(scaling)))
    yremaining = y[-1] + ((y[-1]-y[-2])/(x[-1]-x[-2]))*(xremaining - x[-1])
    x = np.append(x,xremaining)
    y = np.append(y,yremaining)

    fine_expected = 10**y
    fine_bins = np.round(10**x)
    
    for i in range(1,-1,-1):
        fine_expected = np.insert(fine_expected,0,scaling[i]/count_per_diag[i])
        fine_bins = np.insert(fine_bins,0,i).astype(int)
        
    assert np.all((fine_bins[1:]-fine_bins[0:-1])==1)
    
    matrix.data = data/fine_expected[diff]
#     matrix.data[nan_indices] = np.nan
    return matrix


def loopPileup(chunked_feature_df, tup, columns=["ind1", "ind2"], pad=40):
    matrixAndBins, region = tup
    frame = chunked_feature_df[region]
    print(region)
    matrix, mask = matrixAndBins[0], matrixAndBins[1]
    m, n = matrix.shape
    matrix = sparseSymmetricOOE(matrix, mask)
    matrix = sp.csc_matrix(matrix, copy=False)
    locations1 = frame[columns[0]].values
    locations2 = frame[columns[1]].values

    total_PU = []
    mask_PU = []

    for loc1, loc2 in zip(locations1, locations2): 
        hi1 = loc1 + pad + 1
        lo1 = loc1 - pad
        hi2 = loc2 + pad + 1
        lo2 = loc2 - pad
        dm, dn = hi1 - lo1, hi2 - lo2
        
        out_of_bounds = False
        pad_left = pad_right = pad_bottom = pad_top = None
        if lo1 < 0:
            pad_bottom = -lo1
            out_of_bounds = True
        if lo2 < 0:
            pad_left = -lo2
            out_of_bounds = True
        if hi1 > m:
            pad_top = dm - (hi1 - m)
            out_of_bounds = True
        if hi2 > n:
            pad_right = dn - (hi2 - n)
            out_of_bounds = True
            
        if out_of_bounds:        
            i0 = max(lo1, 0)
            i1 = min(hi1, m)
            j0 = max(lo2, 0)
            j1 = min(hi2, n)
            
            submatrix = matrix[i0:i1, j0:j1].toarray().astype(float)
            submask1 = mask[i0:i1]
            submask2 = mask[j0:j1]
            outer = submask1[:,None] * submask2[None,:]
            submatrix[~outer] = 0 
            
            snippet =  np.full((dm, dn), 0.0)
            snip_mask = np.full((dm, dn), 0.0)
            
            snippet[pad_bottom:pad_top, 
                    pad_left:pad_right] = submatrix
            snip_mask[pad_bottom:pad_top, 
                    pad_left:pad_right] = outer
        else:
            submatrix = matrix[lo1:hi1, lo2:hi2].toarray().astype(float)
            submask1 = mask[lo1:hi1]
            submask2 = mask[lo2:hi2]
            outer = submask1[:,None] * submask2[None,:]
            submatrix[~outer] = 0
            
            snippet = submatrix
            snip_mask = outer

        nan_rows = np.sum(snippet, axis=0) == 0
        nan_cols = np.sum(snippet, axis=1) == 0
        snippet[nan_rows, :] = np.nan
        snippet[:, nan_cols] = np.nan
        total_PU.append(snippet)
        mask_PU.append(snip_mask)
    
    if not len(total_PU):
        total_PU = np.nan
        mask_PU = np.nan
    else:
        total_PU = np.dstack(tuple(total_PU))
        mask_PU = np.dstack(tuple(mask_PU))
    
    return total_PU, mask_PU


# def TADPileup(chunked_feature_df, pad, tup, columns=["ind1", "ind2"]):
#     matrixAndBins, region = ziplist
#     frame = chunked_feature_df[region]
#     print(region)
#     matrix, mask = matrixAndBins[0], matrixAndBins[1]
#     m, n = matrix.shape
#     matrix = sparseObservedOverExpected(matrix, mask)
#     matrix = sp.csc_matrix(matrix, copy=False)
#     starts = frame[columns[0]]
#     ends = frame[columns[1]]
#     centers = (starts+ends)//2

    
#         pad = int(1.5*wd)
#         if mid - pad < 0:
#             continue
#         if mid + pad + 1 > M:
#             continue
        
#         submatrix = matrix[mid-pad:mid+pad+1, mid-pad:mid+pad+1].todense()
#         submask = mask[mid-pad:mid+pad+1]
#         outer = submask[:,None] * submask[None,:]
#         submatrix[~outer] = 0 
#         assert submatrix[~outer].sum() == 0
#         zoomedMatrix = zoomArray(submatrix, (81,81), sameSum=True, order=1)
#         zoomedMask = zoomArray(outer, (81,81), sameSum=True, order=1)
#         matrices = matrices + zoomedMatrix
#         masks = masks + zoomedMask
    
#     return matrices, masks
#     total_PU = []
#     mask_PU = []

#     for mid in centers: 
#         hi = mid + pad + 1
#         lo = mid - pad
        
#         dm = hi - lo
#         out_of_bounds = False
#         pad_left = pad_right = pad_bottom = pad_top = None
#         if lo < 0:
#             pad_bottom = -lo
#             pad_left = -lo
#             out_of_bounds = True
#         if hi > m:
#             pad_top = dm - (hi - m)
#             out_of_bounds = True
#         if hi > n:
#             pad_right = dm - (hi - n)
#             out_of_bounds = True
            
#         if out_of_bounds:        
#             i0 = max(lo, 0)
#             i1 = min(hi, m)
#             j0 = max(lo, 0)
#             j1 = min(hi, n)
            
#             submatrix = matrix[i0:i1, j0:j1].toarray().astype(float)
#             submask1 = mask[i0:i1]
#             submask2 = mask[j0:j1]
#             outer = submask1[:,None] * submask2[None,:]
#             submatrix[~outer] = 0 
            
#             snippet =  np.full((dm, dn), 0.0)
#             snip_mask = np.full((dm, dn), 0.0)
            
#             snippet[pad_bottom:pad_top, 
#                     pad_left:pad_right] = submatrix
#             snip_mask[pad_bottom:pad_top, 
#                     pad_left:pad_right] = outer
#         else:
#             submatrix = matrix[lo1:hi1, lo2:hi2].toarray().astype(float)
#             submask1 = mask[lo1:hi1]
#             submask2 = mask[lo2:hi2]
#             outer = submask1[:,None] * submask2[None,:]
#             submatrix[~outer] = 0
            
#             snippet = submatrix
#             snip_mask = outer

#         nan_rows = np.sum(snippet, axis=0) == 0
#         nan_cols = np.sum(snippet, axis=1) == 0
#         snippet[nan_rows, :] = np.nan
#         snippet[:, nan_cols] = np.nan
#         total_PU.append(snippet)
#         mask_PU.append(snip_mask)
    
#     if not len(total_PU):
#         total_PU = np.nan
#         mask_PU = np.nan
#     else:
#         total_PU = np.dstack(tuple(total_PU))
#         mask_PU = np.dstack(tuple(mask_PU))
    
#     return total_PU, mask_PU
        