import os
import sys

sout = sys.stdout 
sys.stdout = open(os.devnull, 'w')

path = '/net/levsha/share/sameer/github/mirnylab-experimental/sameer/'
subfolders = [item for item in os.listdir(path) if item[0]!='.']
for item in subfolders:
    sys.path.insert(0,path+item)
    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mirnylib
import mirnylib.plotting
from mirnylib.genome import Genome

from mirnylib.numutils import logbinsnew
import scipy.sparse as sp
from scipy import signal
import cooler
from cooltools.lib.numutils import _logbins_numba as logbins
from functools import partial
# from eigendecomposition import condense_eigenvector
from bioframe import bedslice
from bioframe.tools import tsv, bedtools
import bioframe
import multiprocess as mp

import DNA_info
sys.stdout = sout

def cooler_matrix_generator(cool, header='weight'):
    '''Creates matrix fetcher function for a balanced cooler.'''
    
    def _matrix_fetcher(region1, region2):
        matrix = cool.matrix(balance=header, sparse=True).fetch(region1, region2)
        return matrix
    
    return _matrix_fetcher


def cooler_mask(cool, header='weight', thres=None):
    '''Creates mask fetcher function for a balanced cooler based on NaN values column in cool.bins() defined by header 
       variable. If header is False, the mask is determined by marginals of the matrix having a value greater than thres.'''
    
    def _mask_fetcher(region):
        if not header:
            if thres is None:
                raise ValueError('If header is False, then a value must be provided to thres')
                
            matrix = cool.matrix(balance=False, sparse=True).fetch(region, region)
            margs = np.bincount(matrix.row, weights=matrix.data, minlength=matrix.shape[0])
            mask = (margs>=thres).astype(int)

        else:
            bins = cool.bins().fetch(region)
            if header not in bins.columns:
                raise ValueError("If header not found in cooler bins' columns")
                
            mask = bins[header].values
            mask = np.logical_not(np.isnan(mask)).astype(int)
            
        return mask
    
    return _mask_fetcher


def cooler_compartment_mask(cool, vector, val='A'):
    assert all(x in vector.columns for x in ['chrom', 'start', 'end', 'label'])
    assert val in vector.label.unique()
    
    vector = vector.groupby('chrom')
    
    def _mask_fetcher(region):
        bins = cool.bins().fetch(region)
        length = len(bins)
        bins['label'] = np.nan
        if region[0] not in vector.chrom.unique():
            return np.zeros(len(bins))
        
        region_eig = bedslice(vector, region[0], region[1], region[2])
        
        if len(region_eig) == 0:
            return np.zeros(len(bins))
        
        for i, row in region_eig.iterrows():
            cond = np.logical_and(bins.start >= row.start, bins.end <= row.end)
            bins.loc[bins[cond].index, 'label'] = row.label
            
        mask = np.logical_and(bins.label==val, ~np.isnan(bins.weight)).values.astype(int)
        
        assert not np.any(np.isnan(bins[mask.astype(bool)]['weight']))
        assert np.all(bins[mask.astype(bool)]['label'] == val)
        
        return mask
    
    return _mask_fetcher

def check_region_format(region, cis=True):
    row_region, col_region = region
    
    if cis and row_region[0] != col_region[0]:
        raise Exception(f'Region: {row_region} and Region:{col_region} come from different chromosomes. \nCannot perform cis scaling. None value will be returned.')
              
    if row_region[1] < 0 or row_region[2] < 0:
        raise Exception(f'Region: {row_region} is incorrectly formatted - Contains negative values. None value will be returned')
    
    if col_region[1] < 0 or col_region[2] < 0:
        raise Exception(f'Region: {col_region} is incorrectly formatted - Contains negative values. None value will be returned')

def orient_and_find_diagonal(row_region, col_region, res):
    '''Analyzes and modifies the regions for use by function '_cis_binning'. Does two things: 
    1) Checks if regions provided are in cis.
    2) Asserts upper-triangularity so that the start of the row region is always less than start 
    of the column region. Since cis matrices are symmetric this matrix doesn't change the 
    information contained in the regions, just the location of the diagonal in the matrix.
    3) Finds that location of the diagonal in the matrix. This is then used by '_cis_binning'
    
    Inputs:
        row_region:           3-tuple of format (chr, start, end) in bp/monomers that represent the row extent.
        col_region:           3-tuple of format (chr, start, end) in bp/monomers that represent the column extent.
        res:                  Resolution of the matrix that represents relation between regions and matrix bins.

    Outputs:
        (chr1, start1, end1): New row extent that is upper triangular.
        (chr2, start2, end2): New column extent that is upper triangular.
        diagonal_index:       Index of the diagonal under the schema that diagonal indices are shifted
                              so that no diagonal has a negative index.'''
    
    chr1, start1, end1 = row_region
    chr2, start2, end2 = col_region
    
    if start1 > start2: # Condition under which matrix isn't upper triangular.
        temp = start1
        start1 = start2
        start2 = temp
        
        temp = end1
        end1 = end2
        end2 = temp
        
    row_width = np.ceil(end1/res) - np.floor(start1/res) # This is just the number of rows in the matrix. This represents the 
                                                         # shift of diagonals under the diagonal indexing schema described 
                                                         # above.
    diag_diff = np.floor(start2/res) - np.floor(start1/res) # This represents the position of the diagonal if there were no 
                                                            # shift. 
    diagonal_index = int(row_width - diag_diff - 1) 
# To get a sense of why the -1 is required, think of the limiting case of when the diagonal is on the left most diagonal i.e. start2 = end1. Then row_width - diag_diff = np.ceil(end1/res) - np.floor(start2/res) = 1. But we need this (the main diagonal) to have index 0 so we subtract 1.
    
    return (chr1, start1, end1), (chr2, start2, end2), diagonal_index


def _cis_binning(matrix_fetcher, row_fetcher, col_fetcher, res, ignore_diags, region):
    ''' Does cis-binning by diagonal of matrix associated with a particular region with separate masking
    capabilities for rows and columns. If the regions overlap i.e. binning over main diagonal then it only 
    bins upper triangle/quadrilateral of the matrix. Otherwise bins over entire region.
        
    Input:
        matrix_fetcher:     Function that takes in two region arguments and returns a Sparse Coordinate matrix 
                            that is to be binned. Matrix does not have to be square.
        row_fetcher:        Function that takes in one region argument and returns an integer array associated with
                            with the mask for that region. Invalid bins must be 0. Specific to row region.
        col_fetcher:        Same as row fetcher but for columns.
        res:                Resolution of matrix in bp/monomers.
        ignore_diags:       Number of diagonals to ignore when binning.
        region:             2-tuple containing row and column regions over which binning is to take place.
                            Row and column regions must be a 3-tuple of a format (chr, start, end) and must 
                            acceptable by both matrix_fetcher, mask_fetcher.

    Outputs:
        df:                 DataFrame with columns 'counts' and 'valid_bins'. 'counts' correspond the total count on 
                            a diagonal. 'valid_bins' corresponds to number of valid_bins on a diagonal. If the region
                            has no non-zero values, then None is returned.

    TODO: 2) Implement 2D mask'''
    
    check_region_format(region)
    
    row_region, col_region = region
    row_region, col_region, d = orient_and_find_diagonal(row_region, col_region, res)

    if row_region[2] > col_region[2]:
        
        subregion1 = ((row_region[0], row_region[1], col_region[2]), col_region)
        subregion2 = ((row_region[0], col_region[2], row_region[2]), col_region)
        print(f'Breaking up ({row_region},{col_region}) into {subregion1} and {subregion2}')
            
        subbinned1 = _cis_binning(matrix_fetcher, row_fetcher, col_fetcher, res, ignore_diags, subregion1)
        subbinned2 = _cis_binning(matrix_fetcher, row_fetcher, col_fetcher, res, ignore_diags, subregion2)
        
        result = merge_regional_results([subbinned1, subbinned2])
        return result
    
    print(f'Binning region: {row_region} and {col_region}')
    
    matrix = matrix_fetcher(row_region, col_region)
    
    if matrix.nnz == 0:
        size = matrix.shape[0]+matrix.shape[1]
        df = pd.DataFrame({'diag':np.arange(size) - d,'total': np.zeros(size), 'valid_bins':np.zeros(size)})
    else:
        row_mask = row_fetcher(row_region)
        col_mask = col_fetcher(col_region)

        assert sp.isspmatrix_coo(matrix), 'Matrix must be a sparse matrix of format COO'
        assert matrix.shape[0] == row_mask.shape[0]
        assert matrix.shape[1] == col_mask.shape[0]


        # Calculates number of valid bins in a diagonal in an efficient manner. The valid bins in a diagonal 
        # is the cross-correlation of the column mask with the row mask. This basically the 
        # complement of the what Nezar does with the bad bins.
    #     np.fft.restore_all()
        cross_corr = signal.fftconvolve(row_mask, col_mask[::-1], mode='full')
        num_valid_bins_per_diag = np.round(cross_corr).astype(int)

        # Masking out data.
        row_ind_to_mask = np.where(row_mask==0)[0] # Tells me the indices of the masked out bins
        col_ind_to_mask = np.where(col_mask==0)[0] 
        masked_rows = np.isin(matrix.row, row_ind_to_mask) # Boolean array that tells us whether the row entry is to be masked
        masked_cols = np.isin(matrix.col, col_ind_to_mask) 
        masked_ind = np.logical_or(masked_rows, masked_cols) # The data is to be masked if the row OR the col are masked
        matrix.data[masked_ind] = 0

        # Binning counts by diagonal
        diff = matrix.col - matrix.row
        diff = diff + matrix.shape[0] - 1 # This ensures that the difference is always non-negative. Elements on the furthest 
                                          # negative diagonal should have diff values of zero.

    # The next line is the workhorse of this function. Bincount counts number of instances of a particular values of 'diff' starting from zero to the value set by parameter 'minlength'. The 'weights' parameter weights each value of diff by a corresponding weight (when weights are not specified, each count has an implicit wieght of 1). The net effect of this line of code is adding up all the matrix.data elements that have a particular value of diff. Since each unique value of diff specifies a diagonal, this line bins by diagonal.             
        total_diag_sum = np.bincount(diff, weights=matrix.data, minlength=len(num_valid_bins_per_diag)) 


        genomic_distance = (np.arange(len(total_diag_sum)).astype(int) - d) # Genomic locations relative to the diagonal in units of resolution.
        df = pd.DataFrame({'total': total_diag_sum,
                           'valid_bins': num_valid_bins_per_diag,
                           'diag': genomic_distance})

        
    if ignore_diags > 0:
        df = df[df['diag'] >= ignore_diags] # Only need this region. Negative regions represents a 
                                                        # subset of the same information.
    else:
        df = df[df['diag'] >= 0]
    
    df['diag'] = (df['diag']*res).astype(int)
    df['region1'] = f'{row_region[0]}:{row_region[1]}-{row_region[2]}'
    df['region2'] = f'{col_region[0]}:{col_region[1]}-{col_region[2]}'
    df.set_index(['region1','region2', 'diag'], inplace=True)

    return df

def cis_binning(region_list, matrix_fetch, row_fetch, col_fetch, res, ignore_diags, mapper=map):
    ''' Iterates _cis_binning over all regions. Supports parallelization e.g. multiprocess.Pool
    
        Input:
            region_list:        List of 2-tuples containing row and column regions over which 
                                binning is to take place. Row and column regions must of a format acceptable 
                                by both matrix_fetch, mask_fetch.
            matrix_fetch:       Function that takes in two region arguments and returns a Sparse Coordinate matrix 
                                that is to be binned. Matrix does not have to be square.
            row_fetch:          Function that takes in one region argument and returns an integer array associated with
                                with the mask for that region. Invalid bins must be 0. Specific to rows.
            col_fetch:          Same as row_fetch but for columns.
            res:                Resolution of matrix in bp/monomers.
            
        
        Outputs:
            binned_results:     List of DataFrames for each entry in 'regions'. Regions with no counts have value None.
                                Each DataFrame has columns 'counts' and 'valid_bins'. 'counts' correspond the 
                                total count on a diagonal. 'valid_bins' corresponds to number of valid_bins on a diagonal.'''
    
    
    binned_results = list(mapper(partial(_cis_binning, matrix_fetch, row_fetch, col_fetch, res, ignore_diags), region_list))

    return binned_results


def _trans_binning(matrix_fetcher, row_fetcher, col_fetcher, res, region):
    ''' Calculates mean of matrix associated with a particular region with separate masking
    capabilities for rows and columns. Since regions are in trans, there is not distance dependence.
        
    Input:
        matrix_fetcher:     Function that takes in two region arguments and returns a Sparse Coordinate matrix 
                            that is to be binned. Matrix does not have to be square.
        row_fetcher:        Function that takes in one region argument and returns an integer array associated with
                            with the mask for that region. Invalid bins must be 0. Specific to row region.
        col_fetcher:        Same as row fetcher but for columns.
        res:                Resolution of matrix in bp/monomers.
        region:             2-tuple containing row and column regions over which binning is to take place.
                            Row and column regions must be a 3-tuple of a format (chr, start, end) and must 
                            acceptable by both matrix_fetcher, mask_fetcher.

    Outputs:
        df:                 DataFrame with columns 'counts' and 'valid_bins'. 'counts' correspond the total count on 
                            a diagonal. 'valid_bins' corresponds to number of valid_bins on a diagonal. If the region
                            has no non-zero values, then None is returned.'''
    
    check_region_format(region, cis=False)
    
    row_region, col_region = region
    
    print('Binning region:', row_region, 'and', col_region)
    
    matrix = matrix_fetcher(row_region, col_region)
    row_mask = row_fetcher(row_region)
    col_mask = col_fetcher(col_region)
    
    assert sp.isspmatrix_coo(matrix), 'Matrix must be a sparse matrix of format COO'
    assert matrix.shape[0] == row_mask.shape[0]
    assert matrix.shape[1] == col_mask.shape[0]

    
    row_ind_to_mask = np.where(row_mask==0)[0] # Tells me the indices of the masked out bins
    col_ind_to_mask = np.where(col_mask==0)[0] 
    masked_rows = np.isin(matrix.row, row_ind_to_mask) # Boolean array that tells us whether the row entry is to be masked
    masked_cols = np.isin(matrix.col, col_ind_to_mask) 
    masked_ind = np.logical_or(masked_rows, masked_cols) # The data is to be masked if the row OR the col are masked
    matrix.data[masked_ind] = 0

    total = np.nansum(matrix.data)
    valid_bins = np.sum(row_mask)*np.sum(col_mask)

    df = pd.DataFrame({'total':[total], 'valid_bins':[valid_bins], 'region1':[row_region], 'region2':[col_region]})
    df.set_index(['region1', 'region2'], inplace=True)
    return df

def trans_binning(region_list, matrix_fetch, row_fetch, col_fetch, res, mapper=map):
    ''' Iterates _cis_binning over all regions. Supports parallelization e.g. multiprocess.Pool
    
        Input:
            region_list:        List of 2-tuples containing row and column regions over which 
                                binning is to take place. Row and column regions must of a format acceptable 
                                by both matrix_fetch, mask_fetch.
            matrix_fetch:       Function that takes in two region arguments and returns a Sparse Coordinate matrix 
                                that is to be binned. Matrix does not have to be square.
            row_fetch:          Function that takes in one region argument and returns an integer array associated with
                                with the mask for that region. Invalid bins must be 0. Specific to rows.
            col_fetch:          Same as row_fetch but for columns.
            res:                Resolution of matrix in bp/monomers.
            
        
        Outputs:
            binned_results:     List of DataFrames for each entry in 'regions'. Regions with no counts have value None.
                                Each DataFrame has columns 'counts' and 'valid_bins'. 'counts' correspond the 
                                total count on a diagonal. 'valid_bins' corresponds to number of valid_bins on a diagonal.'''
    
    
    binned_results = list(mapper(partial(_trans_binning, matrix_fetch, row_fetch, col_fetch, res), region_list))

    return binned_results

def exclude_regions(df, regions_to_keep=[], genome=None, print_final=False):
    if len(regions_to_keep):
        assert genome is not None, 'Please provide valid genome'
        chromsizes = bioframe.fetch_chromsizes(genome)
    else:
        if print_final:
            print(np.asarray(df.region.unique()))
        return df
    
    regions_to_keep = [bioframe.parse_region(reg, chromsizes) 
                  for reg in regions_to_keep]
    
    assert 'region' in df.columns
    
    regions = df['region'].apply(lambda x: bioframe.parse_region(x, chromsizes)).values
    chrom, start, end = list(zip(*regions))
    df['chrom'] = chrom
    df['start'] = start
    df['end'] = end

    new_df = []
    for chrom, start, end in regions_to_keep:
        sub_df = bioframe.bedslice(df, (chrom, start, end))
        new_df.append(sub_df)
    new_df = pd.concat(new_df)
    
    if print_final:
        print(np.asarray(new_df.region.unique()))
        
    del new_df['chrom']
    del new_df['start']
    del new_df['end']
    
    return new_df

def merge_regional_results(df, regions_to_keep=[], genome=None, print_final=False):
    ''' Merges a list of DataFrames into a single sums up entries that have the same index (index should ideally 
        represent genomic distance). The purpose of this function is to merge scalings obtained from individual regions
        into a single scaling. 
        
        Input:
            result_list: List of DataFrames. MAKE SURE LIST DOESN'T CONTAIN None. DataFrames must have 
                         columns named 'sum' and 'valid_bins'. Groupby summing 
                         only done on these columns.

        Outputs:
            combined_df: DataFrame with columns named 'sum' and 'valid_bins' that
                         contains aggregated values of the input DataFrames.'''
    
    df = df.reset_index().set_index('diag')
    df = exclude_regions(df, regions_to_keep=regions_to_keep, genome=genome, print_final=print_final)
    
    drop_cols = [col for col in df.columns if not np.issubdtype(df[col].dtype, np.number)]
    for col in drop_cols:
        del df[col]
        
    gb = df.groupby(by= lambda x:x)
    collapsed_list = []
    for col in df.columns:
        collapsed = gb[col].sum()
        collapsed_list.append(collapsed)
    
    combined_df = pd.concat(collapsed_list, axis=1)
    return combined_df

def log_binning(df, ratio=1.2):
    ''' Log bins a DataFrame in logarithmically spaced bins based on it's index.
        
        Input:
            df:     DataFrame to be log-binned. Must NOT contain a column named 'label'.

        Outputs:
            df_log: Log-binned DataFrame with index containing the mid-points of the log bins.'''
    
    assert 'label' not in df.columns
    bins = df.index.values
    
    if bins[0] == 0:
        df = df.drop(bins[0])
        bins = df.index.values

    log_bins = logbins(np.min(bins), np.max(bins), ratio=ratio)
    num_bins = len(log_bins)-1
    bin_limits = list(zip(log_bins[0:-1],log_bins[1:]))
    bin_mids = pd.DataFrame({'genomic_dist': np.sqrt(log_bins[0:-1] * log_bins[1:])})
    bin_mids['label'] = bin_mids.index
    bin_mids.set_index('label', inplace=True)
    df['label'] = -1
    
    for i in range(num_bins):
        start, end = bin_limits[i]
        where = np.logical_and(df.index>=start, df.index<end)
        df.loc[where, 'label'] = i
    
    df.loc[df.index[-1], 'label'] = i
    
    log_binned_df = df.groupby('label').sum()
    log_binned_df = log_binned_df.join(bin_mids, how='left')
    log_binned_df.set_index('genomic_dist', inplace=True)
    
    return log_binned_df

def cooler_global_scaling(cool, genome, trans=True, mapper=map, balance='weight', thres=None, ignore_diags=2):
    
    row_masker = col_masker = cooler_mask(cool, header=balance, thres=thres)
    matrix_fetcher = cooler_matrix_generator(cool, header=balance)
    resolution = cool.info['bin-size']
    
    chrom_arms = DNA_info.get_chromosome_arms(genome)
    cis_regions = [(arm, arm) for arm in chrom_arms]
    
    cis_results = cis_binning(cis_regions, matrix_fetcher, row_masker, col_masker, resolution, ignore_diags, mapper=mapper)
    cis_results = pd.concat(cis_results)
    cis_results = cis_results.reset_index().rename(columns={'region1':'region'})
    del cis_results['region2']
    cis_results.set_index(['region','diag'], inplace=True, drop=True)

    
    if trans:
        print('Computing trans expected')
        chromsizes = bioframe.fetch_chromsizes(genome)
        trans_regions = [(bioframe.parse_region(cool.chromnames[i], chromsizes=chromsizes), 
                          bioframe.parse_region(cool.chromnames[j], chromsizes=chromsizes))
                         for i in range(len(cool.chromnames))
                          for j in range(i+1, len(cool.chromnames))]

        trans_results = trans_binning(trans_regions, matrix_fetcher, row_masker, col_masker, resolution, mapper=mapper)
        trans_results = [result for result in trans_results if result is not None]
        trans_results = pd.concat(trans_results)
        trans_results['chrom1'] = trans_results.index.map(lambda x: x[0][0]).values
        trans_results['chrom2'] = trans_results.index.map(lambda x: x[1][0]).values
        trans_results.set_index(['chrom1', 'chrom2'], inplace=True)
        
        return cis_results, trans_results
    return cis_results

def create_dir(datasave):
    data_dir = os.path.dirname(datasave)
    if not os.path.exists(data_dir):
        try:
            os.makedirs(data_dir)
        except FileExistsError:
            print(datasave+' already exists')

if __name__ == '__main__':
    genome = sys.argv[1]
    res = int(sys.argv[2])
    do_trans = bool(int(sys.argv[3]))
    balance_col = sys.argv[4]
    savepath = sys.argv[5]
    n_cores = int(sys.argv[6])
    filepaths = sys.argv[7:]
    filepaths = [f'{file}::/resolutions/{res}' if '.mcool' in file else file for file in filepaths]
    
    for file in filepaths:
        print(f'\nComputing scalings for {file}\n')
        cool = cooler.Cooler(file)
        filename = file
        if '.mcool' in file:
            filename = file.split('::')[0]
        filename = filename.split('/')[-1].split('.')[0]
        save = savepath+filename+'.hdf5'
        print(f'Saving to: {save}\n')
        
        with mp.Pool(n_cores) as p:
            if do_trans:
                cis, trans = cooler_global_scaling(cool, genome, trans=True, mapper=p.map, balance=balance_col, ignore_diags=0)
            else:
                cis = cooler_global_scaling(cool, genome, trans=False, mapper=p.map, balance=balance_col, ignore_diags=0)
                
        create_dir(save)
        store = pd.HDFStore(save)
        store.put('scaling', cis, format='table', data_columns=True)
        if do_trans:
            store.put('trans_lvl', trans, format='table', data_columns=True)
        store.close()

