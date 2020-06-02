import os
import sys
import numpy as np
import pandas as pd
import cooler
import matplotlib.pyplot as plt
from cooltools.lib.numutils import observed_over_expected

def construct_cis_saddleplot(mat, vector, num_percentile=10, min_percent=2.5, max_percent = 97.5, dist_lims=None):
    
    mat[~np.isfinite(mat)] = 0
    mask = np.isfinite(vector)
    mat[~(mask[:,None]*mask[None,:])] = np.nan
    mat, _, _, _ = observed_over_expected(mat, (mask[:,None]*mask[None,:]))
    assert np.all(np.where(np.isnan(vector))[0] == np.where(~mask)[0])
    
    pureVector = vector[mask]
    percentiles = np.percentile(pureVector, np.linspace(min_percent, max_percent, num_percentile+1))
    percent_vector = np.searchsorted(percentiles, vector).astype(np.float64) #Ascertains which percentile bracket each element of the vector belongs to. 
    percent_vector = percent_vector - 1.0
    percent_vector[np.logical_or(percent_vector < 0, percent_vector >= num_percentile)] = np.nan # Eliminates values that are beyond limits set by min_percent and max_percent. 

    ar = np.arange(0,len(mat))
    diag_ind = np.abs(ar[:,None] - ar[None,:])
    diag_mask = np.ones_like(diag_ind).astype(bool)
    if dist_lims is None:
        udist_lim = len(mat)
        ldist_lim = 2
    else:
        ldist_lim = dist_lims[0]
        udist_lim = dist_lims[1]
        
    diag_mask[(diag_ind>ldist_lim) & (diag_ind<=udist_lim)] = False
    
    bin_index_saddleplot = percent_vector[:,None] * num_percentile + percent_vector[None,:] #Creates a two digit index matrix based on the percentile of the row and columns.
    to_remove = ~np.isfinite(bin_index_saddleplot)
    bin_index_saddleplot[to_remove] = num_percentile ** 2 # Dummy value associated with bad bins 
    bin_index_saddleplot[diag_mask] = num_percentile ** 2
    bin_index_saddleplot = bin_index_saddleplot.astype(int)
    
    saddle_sum = np.bincount(bin_index_saddleplot.flatten(), weights = mat.flatten(), minlength = num_percentile**2+1)
    saddle_count = np.bincount(bin_index_saddleplot.flatten())
    saddle_sum = saddle_sum[0:-1].reshape(num_percentile,num_percentile)
    saddle_count = saddle_count[0:-1].reshape(num_percentile,num_percentile)
    
    return saddle_sum, saddle_count


def construct_trans_saddleplot(mat, vector1, vector2, num_percentile=10, min_percent=2.5, max_percent = 97.5):
    
    mat[~np.isfinite(mat)] = 0
    mask1 = np.isfinite(vector1)
    mask2 = np.isfinite(vector2)
    mat[~(mask1[:,None]*mask2[None,:])]=np.nan
    exp = np.nanmean(mat)
    mat = mat/exp
    
#     assert np.all(np.where(np.isnan(vector))[0] == np.where(~mask)[0])
    
    pureVector1 = vector1[mask1]
    percentiles1 = np.percentile(pureVector1, np.linspace(min_percent, max_percent, num_percentile+1))
    percent_vector1 = np.searchsorted(percentiles1, vector1).astype(np.float64) #Ascertains which percentile bracket each element of the vector belongs to. 
    percent_vector1 = percent_vector1 - 1.0
    percent_vector1[np.logical_or(percent_vector1 < 0, percent_vector1 >= num_percentile)] = np.nan # Eliminates values that are beyond limits set by min_percent and max_percent. 

    pureVector2 = vector2[mask2]
    percentiles2 = np.percentile(pureVector2, np.linspace(min_percent, max_percent, num_percentile+1))
    percent_vector2 = np.searchsorted(percentiles2, vector2).astype(np.float64) #Ascertains which percentile bracket each element of the vector belongs to. 
    percent_vector2 = percent_vector2 - 1.0
    percent_vector2[np.logical_or(percent_vector2 < 0, percent_vector2 >= num_percentile)] = np.nan # Eliminates values that are beyond limits set by min_percent and max_percent. 

    
    bin_index_saddleplot = percent_vector1[:,None] * num_percentile + percent_vector2[None,:] #Creates a two digit index matrix based on the percentile of the row and columns.
    to_remove = ~np.isfinite(bin_index_saddleplot)
    bin_index_saddleplot[to_remove] = num_percentile ** 2 # Dummy value associated with bad bins 
    bin_index_saddleplot = bin_index_saddleplot.astype(int)
    
    saddle_sum = np.bincount(bin_index_saddleplot.flatten(), weights = mat.flatten(), minlength = num_percentile**2+1)
    saddle_count = np.bincount(bin_index_saddleplot.flatten())
    saddle_sum = saddle_sum[0:-1].reshape(num_percentile,num_percentile)
    saddle_sum = (saddle_sum + saddle_sum.T)
    saddle_count = saddle_count[0:-1].reshape(num_percentile,num_percentile)
    saddle_count= (saddle_count + saddle_count.T)
    
    return saddle_sum, saddle_count 

def comp_strength(saddle):
    side = len(saddle)-1
    return np.array([saddle[0,0], saddle[side,side], saddle[0, side]])


def genome_comp_strength(cooler, vector, func=np.nanmean, num_percentile=10, min_percent=2.5, max_percent = 97.5, dist_lim=None):
    chroms = vector['chrom'].unique()
    try:
        balanced=True
        mat = cooler.matrix(balance=True).fetch(chroms[0])
    except ValueError:
        print('Cooler not balanced... Working with unbalanced cooler')
        balanced=False
            
    cmp_str = [] 
    for ind ,chrom in enumerate(chroms):
        vec = vector[(vector.chrom==chrom)]['cis'].values
        mat = cooler.matrix(balance=balanced).fetch(chrom)

        assert len(mat) == len(vec)
        
        saddle = construct_saddleplot(mat, vec, num_percentile=num_percentile, min_percent=min_percent, max_percent = max_percent, dist_lim=dist_lim)
        cmp_str.append(comp_strength(saddle))
    
    cmp_str = np.array(cmp_str)
    return func(cmp_str)
