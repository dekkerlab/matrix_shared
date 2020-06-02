import os
import sys
path = '/net/levsha/share/sameer/github/mirnylab-experimental/sameer/'
subfolders = [item for item in os.listdir(path) if item[0]!='.']
for item in subfolders:
    sys.path.insert(0,path+item)
    
import multiprocess as mp
import numpy as np
import pandas as pd

import mirnylib.plotting
import cooler
import cooltools
from cooltools import snipping
import bioframe
from sklearn.mixture import GaussianMixture

import DNA_info
import matrix_manager as mm
import cooltools_pileups


resolution = 5000
arms = DNA_info.get_chromosome_arms('hg38')
data = mm.Dataset()
table = data.get_tables()
table = mm.filter_data(table, filter_dict={'enzyme':'double'})
table = mm.get_coolers(table, res=resolution)

dot_path = f'./pileups/dots/rao2014/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist_hg38.txt'
dot_list = pd.read_csv(dot_path, sep='\t')
dot_list['pos1'] = (dot_list['start1'] + dot_list['end1'])//2
dot_list['pos2'] = (dot_list['start2'] + dot_list['end2'])//2
        
binned_list1 = snipping.make_bin_aligned_windows(resolution, dot_list['chrom1'].values, 
                                       dot_list['pos1'].values, flank_bp=40*resolution)
binned_list1 = snipping.assign_regions(binned_list1, arms)

binned_list2 = snipping.make_bin_aligned_windows(resolution, dot_list['chrom2'].values, 
                                       dot_list['pos2'].values, flank_bp=40*resolution)
binned_list2 = snipping.assign_regions(binned_list2, arms)
binned_list1 = binned_list1[['chrom', 'start', 'end', 'region']].merge(
                                               binned_list2[['chrom', 'start', 'end', 'region']], 
                                               left_index=True, right_index=True, suffixes=('1', '2'))

features = binned_list1[binned_list1.region1 == binned_list1.region2]
del features['region2']
features = features.rename(columns={'region1':'region'})
assert np.all(features['region'].apply(lambda x: isinstance(x, str)))

if len(features) != len(dot_list):
    dot_list.loc[features.index].to_csv(dot_path, sep='\t')

for _, row in table.iterrows():
    name = row['lib_name']
    print(name)

    cool = row[f'cooler_{resolution}']           
    snipper = cooltools_pileups.LocalObsExpSnipper(cool, cooler_opts={'balance':True})
    with mp.Pool(20) as p:
        pileup = snipping.pileup(features, snipper.select, snipper.snip, map=p.map)

    savepath = f'./pileups/dots/rao2014/{resolution}'
    os.makedirs(savepath, exist_ok=True)
    savepath = f'{savepath}/{name}.npy'
    np.save(savepath, pileup)
    print(f'Saved to: {savepath}')