import os
import importlib as imp

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, interp2d
import scipy.signal as signal
from scipy.ndimage.filters import gaussian_filter1d, gaussian_filter
from scipy.spatial import cKDTree, KDTree
from scipy.spatial.distance import minkowski
import multiprocess as mp
import h5py
import itertools


import cooler
import cooltools
import cooltools.snipping as snipping
# import cooltools.expected as expected
from cooltools.lib.numutils import logbins
import bioframe
from bioframe import fetch_chromsizes
#import mirnylib.plotting
from bioframe.tools import bedtools, tsv# import intersect

#import DNA_info
#import new_scalings
#import microc
#import cooltools_pileups

from PIL import Image
import io
import seaborn as sns

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
#from matplotlib_venn import venn2, venn3, venn2_circles
from mpl_toolkits.mplot3d import Axes3D
# import mirnylib.plotting
# %matplotlib notebook
import sys

def insert_margin(df, margin):
    df['start'] = df['pos'].apply(lambda x: x - margin if x - margin >= 0 else 0)
    df['end'] = df['pos'] + margin
    del df['pos']
    df = df.assign(dot_id=df.index.values).sort_values(['chrom', 'start'])
    return df

def overlap_2d_bedtools(target, reference, margin, return_ref=False):
    
    l_target = target[['chrom1', 'pos1']].rename(columns=lambda x: x.replace('1',''))
    l_target = insert_margin(l_target, margin)

    l_ref = reference[['chrom1', 'pos1']].rename(columns=lambda x: x.replace('1',''))
    l_ref = insert_margin(l_ref, margin)
    
    with tsv(l_ref) as a, tsv(l_target) as b: 
        l_intersect = bedtools.intersect(a=a.name, b=b.name, wa=True, wb=True)
        l_intersect.columns = [col+'_r' for col in l_ref.columns] + [col+'_t' for col in l_target.columns]
        l_intersect.set_index(['dot_id_r','dot_id_t'], inplace=True)
    
    
    r_target = target[['chrom2', 'pos2']].rename(columns=lambda x: x.replace('2',''))
    r_target = insert_margin(r_target, margin)
    
    r_ref = reference[['chrom2', 'pos2']].rename(columns=lambda x: x.replace('2',''))
    r_ref = insert_margin(r_ref, margin)
    
    with tsv(r_ref) as a, tsv(r_target) as b: 
        r_intersect = bedtools.intersect(a=a.name, b=b.name, wa=True, wb=True)
        r_intersect.columns = [col+'_r' for col in r_ref.columns] + [col+'_t' for col in r_target.columns]
        r_intersect.set_index(['dot_id_r','dot_id_t'], inplace=True)
    
    merged_df = l_intersect.merge(r_intersect, how='inner', left_index=True, right_index=True).reset_index()
    
    target_inds = merged_df.dot_id_t.values
    target_result = target.loc[target_inds].copy().sort_index().drop_duplicates()
    
    if return_ref:
        ref_inds = merged_df.dot_id_r.values
        reference_result = reference.loc[ref_inds].copy().sort_index().drop_duplicates()
        
        return target_result, reference_result
    
    return target_result

###
print(sys.argv)

loop_path = sys.argv[1]
cooler_path = sys.argv[2]
savepath = sys.argv[3]

f1=sys.argv[4]
f2=sys.argv[5]

f1_mod=f1.replace("-", "_")
f2_mod=f2.replace("-", "_")

loop_files = {
        f1_mod:f1+'/combineddots/cloops_'+f1+'.mapq_30.1000.mcool.combined.bedpe.postproc',
        f2_mod:f2+'/combineddots/cloops_'+f2+'.mapq_30.1000.mcool.combined.bedpe.postproc'
            }
coolers = {
        f1_mod:cooler.Cooler(f'{cooler_path}{f1}.mapq_30.1000.mcool::/resolutions/5000'),
        f2_mod:cooler.Cooler(f'{cooler_path}{f2}.mapq_30.1000.mcool::/resolutions/5000')
          }


all_dots = {}
for key, file in loop_files.items():
    print(key)

    dots = pd.read_csv(loop_path+file, sep='\t')#merge_proximal_entries(dots_5k, dots_10k, 10000)
    dots['pos1'] = (dots['start1'] + dots['end1'])//2
    dots['pos2'] = (dots['start2'] + dots['end2'])//2
    all_dots[key] = dots
    print('Number of dots:', dots.shape[0], '\n')



A, B = overlap_2d_bedtools(all_dots[f1_mod], all_dots[f2_mod], 10000, return_ref=True)
common_dots = {f1_mod: A, f2_mod:B}
print('# Common entries for in ',f1, A.shape[0], '\n')
print('# Common entries for in ', f2,B.shape[0], '\n')
A.to_csv(f1_mod+"_"+f2_mod+".txt",sep="\t")
B.to_csv(f2_mod+"_"+f1_mod+".txt",sep="\t")

unique_dots = {}
s=0
for key, all_df in all_dots.items():
    print(key, '\n')
    file1_name = list(all_dots.keys())[0]
    file2_name = list(all_dots.keys())[1]
    common_df = common_dots[key]
    common_inds = np.unique(common_df.index.values)
    all_inds = all_df.index.values
    
    unique_df = all_df.loc[~np.isin(all_inds, common_inds)]
    print('Number of unique dots:', unique_df.shape[0], '\n')
    unique_dots[key] = unique_df

    if(s>0):
        unique_df.to_csv(file2_name+"_uniq_comp_to_"+file1_name+".txt",sep="\t")

    else:
        unique_df.to_csv(file1_name+"_uniq_comp_to_"+file2_name+".txt",sep="\t")
    s=s+1
    assert common_df.shape[0] + unique_df.shape[0] == all_df.shape[0]




 



