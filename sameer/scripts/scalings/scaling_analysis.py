import numpy as np
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d
from cooltools.lib.numutils import logbins
from new_scalings import log_binning

def format_scaling(scaling, ratio=1.2):
    bins = scaling.index.values
    if np.issubdtype(bins.dtype, np.number):
        diff = bins[1:] - bins[0:-1]
        if np.all(diff == diff[0]):
            scaling = log_binning(scaling, ratio=ratio)
            
    if 'prob' not in scaling.columns: 
        assert 'total' in scaling.columns
        assert 'valid_bins' in scaling.columns
        scaling['prob'] = scaling['total']/scaling['valid_bins']
    
    return scaling
        
def normalize(scaling, trans=None, pos=20000):
    '''Normalizes the the scalings & trans levels by making the value of the scaling at the 
    point 'pos' equal to 1.
    
        scale_df:   DataFrame containing scalings. The index must be genomic distance.
                    
        trans_df:   DataFrame containing trans expected. indexed by MultiIndex: 'chrom1', 'chrom2'.
                    
        pos:        Integer variable the contains the genomic distance at which 
                    scaling at which normalization occurs.'''
    
    scale_df = scaling.copy(deep=True)
    bins = scale_df.index.values
    loc = np.argmin(abs(bins-pos))
    entry = scale_df.iloc[loc]['prob']
    scale_df['prob'] = scale_df['prob']/entry

    if trans is not None:
        trans_df = trans.copy(deep=True)
        trans_df['prob'] = trans_df['prob']/entry
        return scale_df, trans_df
    else:
        return scale_df
    
def finite_difference(scaling_df, smooth):

    prob = np.log10(scaling_df['prob'].values)
    bins = np.log10(scaling_df.index.values)

    slope = (prob[1:]-prob[0:-1])/(bins[1:]-bins[0:-1])
    if smooth:
        slope = gaussian_filter1d(slope, smooth)
    new_bins = (bins[1:]+bins[0:-1])/2
    new_bins = 10**new_bins
    
    slope_df = pd.DataFrame({'der':slope, 'bins':new_bins})
    slope_df.set_index('bins', inplace=True)
    return slope_df

#def SG_derivative(scaling_df, window_size, poly_order):

#    prob = np.log10(scaling_df['prob'].values)
#    bins = np.log10(scaling_df.index.values)
    
    # Need to ensure that bins are equally separated. Usually log binning results in some empty bins 
    # at small distances that are dropped. We will artifically reinsert them by linear interpolation.
#    new_bins = []
#    new_prob = []
#    min_diff = np.min(np.round(bins[1:]-bins[0:-1], 4))
#    for i in range(len(prob)-1):

#        new_prob.append(prob[i])
#        new_bins.append(bins[i])
#        delta = bins[i+1] - bins[i]
#        slope = (prob[i+1] - prob[i])/delta
#        if  delta > min_diff:
#            n_fill = np.round(delta/min_diff).astype(int)
#            for j in range(1, n_fill):
#                new_prob.append(prob[i]+ slope*j*min_diff)
#                new_bins.append(bins[i]+ j*min_diff)
#    new_prob = np.asarray(new_prob)
#    new_bins = np.asarray(new_bins)
    
#    slope = savgol_filter(new_prob, window_size, poly_order, deriv=1, delta=min_diff, mode='nearest')
    
#    slope_df = pd.DataFrame({'der':slope, 'bins':10**new_bins})
#    return slope_df