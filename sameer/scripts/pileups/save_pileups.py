import sys
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mirnylib.genome import Genome
from pileups import *
import cooler
from multiprocessing import Pool
import h5py
from itertools import chain


def genome_pileup(c, regions, feature_locations_df, func):
   
    matrix_gen = fetchCooler(c,regions)
    chunked_features_dict = chunkDataFrame(c, regions, feature_locations_df, columns=[("chr1","pos1"), ("chr2","pos2")])

#    pileup_list = []
#    for i,region in zip(a, regions):
#        pileup_list.append(func((i, chunkedLoops[region],region)))
    with Pool(20) as p:
        pileup_list = p.map(func, [(mat, chunked_features_dict[region], region) for mat,region in zip(matrix_gen, regions)])

    totals, masks = list(zip(*pileup_list))
    totals = np.concatenate(tuple(totals), axis=0)
    masks = np.concatenate(tuple(masks), axis=0)
    
    retrieved_df = pd.concat(list(chunked_features_dict.values()))
    totals, masks = fill_missing_entries(feature_locations_df, retrieved_df, totals, masks)
    
    return totals, masks
    
    
def fill_missing_entries(original_features_df, retrieved_features_df, totals, masks):
    original_ind = original_features_df.index.sort_values()
    retrieved_ind = retrieved_features_df.index.sort_values()
    
    miss_ind = np.setdiff1d(original_ind, retrieved_ind)
    print('The following entries from the feature list were not captured
    i = 0
    while len(miss_ind) > 0:
        to_insert = miss_ind[miss_ind <= totals.shape[0]-1]
        miss_ind = miss_ind[miss_ind > totals.shape[0]-1]
        totals = np.insert(totals, to_insert, 0, axis=0)
        masks = np.insert(masks, to_insert, 0, axis=0)
        i += 1
        if i%100 == 0:
            print(str(i)+' iterations completed in trying to fill in place missing features.')
            print(miss_ind)
    
    return totals, masks

                  
def chr_cooler_format(c,df):
    c_chr = bool(np.all(['chr' in name for name in c.chromnames]))

    df['chr1'] = df['chr1'].apply(lambda x: 'chr'*c_chr+x[x.find('r')+1:])
    df['chr2'] = df['chr2'].apply(lambda x: 'chr'*c_chr+x[x.find('r')+1:])
    return df
    
def get_chrom_arms(c, gen_name):
    genome = Genome('/net/levsha/share/lab/genomes/'+gen_name)
    has_chr = 'chr' in c.chromnames[0]
    res = c.info['bin-size']
    labels = ['chr'*has_chr+label for label in genome.chrmLabels[0:-2]]
    chrmLens = genome.chrmLens[0:-2] 
    centSt = genome.cntrStarts[0:-2]
    centEnd = genome.cntrEnds[0:-2]
    
    arm1 = list(zip(labels, np.zeros(len(labels),dtype=int), (centSt//res + 1)*res))
    arm2 = list(zip(labels, centEnd//res*res, chrmLens//res*res))
    
    return list(chain.from_iterable(zip(arm1, arm2)))


#def chrom_corr(c):
#    if 'chr' in c.chromnames[0]:
#        regions = [('chr'+str(i),None,None) for i in range(1,23)]
#        regions.append(('chrX', None, None))
#        return regions
#    else:
#        regions = [(str(i),None,None) for i in range(1,23)]
#        regions.append(('X', None, None))
#        return regions
    
def drop_unbalanced(filepaths):
    new_fp = []
    bad_fp = []
    for file in filepaths:
        c = cooler.Cooler(file)
        if 'weight' in c.bins().columns:
            new_fp.append(file)
        else:
            bad_fp.append(file)
            
    if len(new_fp)==0:
        raise Exception('All coolers are unbalanced and hence cooler calculation cannot be performed')
    elif len(bad_fp)==0:
        print('All coolers are balanced :)')
    else:
        print('The following coolers were not balanced and hence were dropped:\n')
        print(bad_fp)
    
    return new_fp


def create_dir(datasave):
    data_dir = os.path.dirname(datasave)
    if not os.path.exists(data_dir):
        try:
            os.makedirs(data_dir)
        except FileExistsError:
            print(datasave+' already exists')
            

def batch_process(filepaths, genome, feature_type, savepath=None):

    if feature_type.lower() in ['loops','loop']:
        feature_type = 'loops'
    elif feature_type.lower() in ['tads','tad']:
        feature_type = 'TADs'
    else:
        raise ValueError("Please enter feature type as first argument. So far 'TADs' and 'loops' are implemented")
        
    if len(filepaths) == 0:
        raise Exception('Empty file list')
    elif isinstance(filepaths, str):
        filepaths = [filepaths]
    n = len(filepaths)
  
    filepaths = drop_unbalanced(filepaths)
#     file_df = pd.DataFrame({'filepath':filepaths})
#     file_df = extract_info(file_df, genomes, savepaths, feature_type)
    
    if feature_type == 'loops':
        if genomes == 'hg19':
            feature_list = pd.read_csv("/net/levsha/hic/DrosophilaSingleCell2015/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.txt", sep="\t")
        else:
            feature_list = pd.read_csv("/net/levsha/hic/DrosophilaSingleCell2015/GSE63525_mouse_lymphoblasts_HiCCUPS_looplist.txt", sep="\t")
        
        feature_list["pos1"] = (feature_list["x1"] + feature_list["x2"]) // 2 
        feature_list["pos2"] = (feature_list["y1"] + feature_list["y2"]) // 2
        func = loopPileup
        
    elif feature_type == 'TADs':
        if genomes == 'hg19':
            feature_list = pd.read_csv("/net/levsha/hic/DrosophilaSingleCell2015/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", sep="\t")
        else:
            feature_list = pd.read_csv("/net/levsha/hic/DrosophilaSingleCell2015/GSE63525_mouse_lymphoblasts_Arrowhead_domainlist.txt", sep="\t")
        feature_list["pos1"] = feature_list["x1"]
        feature_list["pos2"] = feature_list["x2"]
        func = TADPileup
        
    totals_list = []
    masks_list = []
    for i, filepath in enumerate(filepaths):
                  
        print(str(i)+" of "+str(n)+": "+filepath)
        c = cooler.Cooler(filepath)
        res = c.info['bin-size']
        feature_list = chr_cooler_format(c, feature_list)
        regions = get_chrom_arms(c,genome)
        
        totals, masks = genome_pileup(c, regions, feature_list, func)
        totals_list.append(totals)
        masks_list.append(masks)
        if savepath is not None:
            if filepath.find('levsha/share/lab') != -1:
                subfolder = filepath.replace('/net/levsha/share/lab/','')
            elif filepath.find('levsha/share/sameer/scalings_project/') != -1:
                subfolder = filepath.replace('levsha/share/sameer/scalings_project/','')
            else:
                subfolder = filepath.replace('/net/levsha/hic/','')
            subfolder = subfolder[0:subfolder.find('/')]+ '/'+feature_type+'/'
            
            filename = filepath[filepath.rfind('/')+1:filepath.find('.'+str(res)+'.cool')]+'.'+str(res)+'.hdf5'
            save = savepath+subfolder+filename
            print('Saved to: '+save)
            create_dir(save)
                  
            HDFStore = h5py.File(save, 'w')
            HDFStore.create_dataset('masks', data=masks)
            HDFStore.create_dataset('totals', data=totals)
            HDFStore.close()
                  
    return totals_list, masks_list
    
if __name__ == "__main__":
    feature_type = sys.argv[1]
    genomes = sys.argv[2]
    filepaths = sys.argv[3:]
    savepaths = '/net/levsha/share/sameer/scalings_project/'
    
    batch_process(filepaths, genomes, feature_type, savepaths)