import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.interpolate import interp1d
import scipy.signal as signal
import multiprocess as mp

import cooler
import cooltools.snipping as snipping
from cooltools.lib.numutils import logbins
import bioframe

from mirnylib.numutils import zoomArray

import DNA_info

def clean_up_loops(loop_list, arms):
    '''Removes items for loop_list that are either not in any of the regions contained in arms or is contained
    in a region between two of the regions in arms. Use this to clean up a new list of loops and save the 
    modified list for future use.'''
    
    for each in ['chrom1', 'chrom2', 'pos1', 'pos2']:
        assert each in loop_list.columns
    
    features = loop_list.copy(deep=True)
    features['index1'] = -1
    features['index2'] = -1
    
    for i, arm in enumerate(arms):
        chrom = arm[0]
        start = arm[1]
        end = arm[2]
        
        features['index1'] = features.apply(lambda x: i 
                                            if (x['chrom1']==chrom and x['pos1'] > start)
                                            and x['pos1'] < end else x['index1'], axis=1)
        
        features['index2'] = features.apply(lambda x: i 
                                            if (x['chrom2']==chrom and x['pos2'] > start)
                                            and x['pos2'] < end else x['index2'], axis=1)
        
    features = features[np.logical_or(features.index1 != -1, features.index2 != -1)]
    features = features[features.index1 == features.index2]
    features = features[['chrom1', 'pos1', 'chrom2', 'pos2']]

    return features

def sparseSymmetricOOE(matrix, mask, log_binning=True):
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
    
    if log_binning:
        
        hi = len(scaling)
        lo = 1
        ratio = 1.2
        N = int(np.log(hi / lo) / np.log(ratio))
        
        bins = logbins(1, len(scaling), N=N)
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
    else:
        fine_expected = scaling/count_per_diag

    
    matrix.data = data/fine_expected[diff]
#     matrix.data[nan_indices] = np.nan
    return matrix

class LocalObsExpSnipper:
    '''Object whose methods are fed to cooltools.snipping.pileup function. Only works if regions that
    are fed to the select method are the same i.e. region1 MUST BE SAME AS region2.'''
    
    def __init__(self, clr, cooler_opts=None, log_binning=True):
        self.clr = clr
        self.log_binning = log_binning
        self.binsize = self.clr.binsize
        self.offsets = {}
        self.pad = True
        self.cooler_opts = {} if cooler_opts is None else cooler_opts
        self.cooler_opts.setdefault('sparse', True)
    
    def select(self, region1, region2):
        print(region1, region2)
        self.offsets[region1] = self.clr.offset(region1) - self.clr.offset(region1[0])
        self.offsets[region2] = self.clr.offset(region2) - self.clr.offset(region2[0])
        matrix = (self.clr.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        
        mask = self.clr.bins().fetch(region1)
        mask = np.isfinite(mask['weight'].values).astype(int)
        matrix = sparseSymmetricOOE(matrix, mask, log_binning=self.log_binning)
        
        if self.cooler_opts['sparse']:
            matrix = matrix.tocsr()

        return matrix

    def snip(self, matrix, region1, region2, tup):
        s1, e1, s2, e2 = tup
        offset1 = self.offsets[region1]
        offset2 = self.offsets[region2]
        binsize = self.binsize
        
        lo1, hi1 = (s1 // binsize) - offset1, (e1 // binsize) - offset1
        lo2, hi2 = (s2 // binsize) - offset2, (e2 // binsize) - offset2
        
        if hi1 < 0 or hi2 < 0:
            print(region1, s1, e1, region2, s2, e2)
            print(offset1, offset2)
            print(lo1, hi1, lo2, hi2)
        assert hi1 >= 0
        assert hi2 >= 0
        
        m, n = matrix.shape
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
            snippet =  np.full((dm, dn), 0.0)
            snippet[pad_bottom:pad_top, 
                    pad_left:pad_right] = matrix[i0:i1, j0:j1].toarray().astype(float)
#             print(m,n)

#             print(i0, i1, j0, j1)
#             print(matrix[i0:i1, j0:j1].toarray().astype(float).shape)
#             print(snippet[pad_bottom:pad_top, pad_left:pad_right].shape)
        else:
            snippet = matrix[lo1:hi1, lo2:hi2].toarray().astype(float)

        nan_rows = np.sum(snippet, axis=0) == 0
        nan_cols = np.sum(snippet, axis=1) == 0
        snippet[nan_rows, :] = np.nan
        snippet[:, nan_cols] = np.nan
            
        return snippet


class DifferenceSnipper:
    '''Object whose methods are fed to cooltools.snipping.pileup function. Only works if regions that
    are fed to the select method are the same i.e. region1 MUST BE SAME AS region2.'''
    
    def __init__(self, clr1, clr2, cooler_opts=None, log_binning=True):
        self.clr1 = clr1
        self.clr2 = clr2
        self.log_binning = log_binning
        assert clr1.binsize == clr2.binsize
        self.binsize = self.clr1.binsize
        self.offsets = {}
        self.pad = True
        self.cooler_opts = {} if cooler_opts is None else cooler_opts
        self.cooler_opts.setdefault('sparse', True)
    
    def select(self, region1, region2):
        print(region1, region2)
        self.offsets[region1] = self.clr1.offset(region1) - self.clr1.offset(region1[0])
        self.offsets[region2] = self.clr1.offset(region2) - self.clr1.offset(region2[0])
        
        matrix1 = (self.clr1.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        matrix2 = (self.clr2.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        
        mask1 = self.clr1.bins().fetch(region1)
        mask1 = np.isfinite(mask1['weight'].values).astype(int)
        mask2 = self.clr2.bins().fetch(region1)
        mask2 = np.isfinite(mask2['weight'].values).astype(int)
        
        matrix1 = sparseSymmetricOOE(matrix1, mask1, log_binning=self.log_binning)
        matrix2 = sparseSymmetricOOE(matrix2, mask2, log_binning=self.log_binning)
        
        matrix  = sp.coo_matrix(matrix1.todense() - matrix2.todense())
        
        if self.cooler_opts['sparse']:
            matrix = matrix.tocsr()

        return matrix

    def snip(self, matrix, region1, region2, tup):
        s1, e1, s2, e2 = tup
        offset1 = self.offsets[region1]
        offset2 = self.offsets[region2]
        binsize = self.binsize
        
        lo1, hi1 = (s1 // binsize) - offset1, (e1 // binsize) - offset1
        lo2, hi2 = (s2 // binsize) - offset2, (e2 // binsize) - offset2
        
        if hi1 < 0 or hi2 < 0:
            print(region1, s1, e1, region2, s2, e2)
            print(offset1, offset2)
            print(lo1, hi1, lo2, hi2)
        assert hi1 >= 0
        assert hi2 >= 0
        
        m, n = matrix.shape
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
            snippet =  np.full((dm, dn), 0.0)
            snippet[pad_bottom:pad_top, 
                    pad_left:pad_right] = matrix[i0:i1, j0:j1].toarray().astype(float)
#             print(m,n)

#             print(i0, i1, j0, j1)
#             print(matrix[i0:i1, j0:j1].toarray().astype(float).shape)
#             print(snippet[pad_bottom:pad_top, pad_left:pad_right].shape)
        else:
            snippet = matrix[lo1:hi1, lo2:hi2].toarray().astype(float)

        nan_rows = np.sum(snippet, axis=0) == 0
        nan_cols = np.sum(snippet, axis=1) == 0
        snippet[nan_rows, :] = np.nan
        snippet[:, nan_cols] = np.nan
            
        return snippet

class RatioSnipper:
    '''Object whose methods are fed to cooltools.snipping.pileup function. Only works if regions that
    are fed to the select method are the same i.e. region1 MUST BE SAME AS region2.'''
    
    def __init__(self, clr1, clr2, cooler_opts=None, log_binning=True):
        self.clr1 = clr1
        self.clr2 = clr2
        self.log_binning = log_binning
        assert clr1.binsize == clr2.binsize
        self.binsize = self.clr1.binsize
        self.offsets = {}
        self.pad = True
        self.cooler_opts = {} if cooler_opts is None else cooler_opts
        self.cooler_opts.setdefault('sparse', True)
    
    def select(self, region1, region2):
        print(region1, region2)
        self.offsets[region1] = self.clr1.offset(region1) - self.clr1.offset(region1[0])
        self.offsets[region2] = self.clr1.offset(region2) - self.clr1.offset(region2[0])
        
        matrix1 = (self.clr1.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        matrix2 = (self.clr2.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        
        mask1 = self.clr1.bins().fetch(region1)
        mask1 = np.isfinite(mask1['weight'].values).astype(int)
        mask2 = self.clr2.bins().fetch(region1)
        mask2 = np.isfinite(mask2['weight'].values).astype(int)
        
        matrix1 = sparseSymmetricOOE(matrix1, mask1, log_binning=self.log_binning)
        matrix2 = sparseSymmetricOOE(matrix2, mask2, log_binning=self.log_binning)
        
        matrix1 = matrix1.todense()
        matrix1[matrix1==0] = 1
        
        matrix2 = matrix2.todense()
        matrix2[matrix2==0] = 1
        
        matrix  = matrix1/matrix2
        
#         if self.cooler_opts['sparse']:
#             matrix = matrix.tocsr()

        return matrix

    def snip(self, matrix, region1, region2, tup):
        s1, e1, s2, e2 = tup
        offset1 = self.offsets[region1]
        offset2 = self.offsets[region2]
        binsize = self.binsize
        
        lo1, hi1 = (s1 // binsize) - offset1, (e1 // binsize) - offset1
        lo2, hi2 = (s2 // binsize) - offset2, (e2 // binsize) - offset2
        
        if hi1 < 0 or hi2 < 0:
            print(region1, s1, e1, region2, s2, e2)
            print(offset1, offset2)
            print(lo1, hi1, lo2, hi2)
        assert hi1 >= 0
        assert hi2 >= 0
        
        m, n = matrix.shape
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
            snippet =  np.full((dm, dn), 0.0)
            snippet[pad_bottom:pad_top, 
                    pad_left:pad_right] = np.asarray(matrix[i0:i1, j0:j1]).astype(float)
#             print(m,n)

#             print(i0, i1, j0, j1)
#             print(matrix[i0:i1, j0:j1].toarray().astype(float).shape)
#             print(snippet[pad_bottom:pad_top, pad_left:pad_right].shape)
        else:
            snippet = np.asarray(matrix[lo1:hi1, lo2:hi2]).astype(float)

        nan_rows = np.sum(snippet, axis=0) == 0
        nan_cols = np.sum(snippet, axis=1) == 0
        snippet[nan_rows, :] = np.nan
        snippet[:, nan_cols] = np.nan
            
        return snippet
    
class DifferenceSnipper:
    '''Object whose methods are fed to cooltools.snipping.pileup function. Only works if regions that
    are fed to the select method are the same i.e. region1 MUST BE SAME AS region2.'''
    
    def __init__(self, clr1, clr2, cooler_opts=None, log_binning=True):
        self.clr1 = clr1
        self.clr2 = clr2
        self.log_binning = log_binning
        assert clr1.binsize == clr2.binsize
        self.binsize = self.clr1.binsize
        self.offsets = {}
        self.pad = True
        self.cooler_opts = {} if cooler_opts is None else cooler_opts
        self.cooler_opts.setdefault('sparse', True)
    
    def select(self, region1, region2):
        print(region1, region2)
        self.offsets[region1] = self.clr1.offset(region1) - self.clr1.offset(region1[0])
        self.offsets[region2] = self.clr1.offset(region2) - self.clr1.offset(region2[0])
        
        matrix1 = (self.clr1.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        matrix2 = (self.clr2.matrix(**self.cooler_opts)
                          .fetch(region1, region2))
        
        mask1 = self.clr1.bins().fetch(region1)
        mask1 = np.isfinite(mask1['weight'].values).astype(int)
        mask2 = self.clr2.bins().fetch(region1)
        mask2 = np.isfinite(mask2['weight'].values).astype(int)
        
        matrix1 = sparseSymmetricOOE(matrix1, mask1, log_binning=self.log_binning)
        matrix2 = sparseSymmetricOOE(matrix2, mask2, log_binning=self.log_binning)
        
        matrix1 = matrix1.todense()
        matrix1[matrix1==0] = 1
        
        matrix2 = matrix2.todense()
        matrix2[matrix2==0] = 1
        
        matrix  = matrix1 - matrix2
        
#         if self.cooler_opts['sparse']:
#             matrix = matrix.tocsr()

        return matrix

    def snip(self, matrix, region1, region2, tup):
        s1, e1, s2, e2 = tup
        offset1 = self.offsets[region1]
        offset2 = self.offsets[region2]
        binsize = self.binsize
        
        lo1, hi1 = (s1 // binsize) - offset1, (e1 // binsize) - offset1
        lo2, hi2 = (s2 // binsize) - offset2, (e2 // binsize) - offset2
        
        if hi1 < 0 or hi2 < 0:
            print(region1, s1, e1, region2, s2, e2)
            print(offset1, offset2)
            print(lo1, hi1, lo2, hi2)
        assert hi1 >= 0
        assert hi2 >= 0
        
        m, n = matrix.shape
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
            snippet =  np.full((dm, dn), 0.0)
            snippet[pad_bottom:pad_top, 
                    pad_left:pad_right] = np.asarray(matrix[i0:i1, j0:j1]).astype(float)
#             print(m,n)

#             print(i0, i1, j0, j1)
#             print(matrix[i0:i1, j0:j1].toarray().astype(float).shape)
#             print(snippet[pad_bottom:pad_top, pad_left:pad_right].shape)
        else:
            snippet = np.asarray(matrix[lo1:hi1, lo2:hi2]).astype(float)

        nan_rows = np.sum(snippet, axis=0) == 0
        nan_cols = np.sum(snippet, axis=1) == 0
        snippet[nan_rows, :] = np.nan
        snippet[:, nan_cols] = np.nan
            
        return snippet