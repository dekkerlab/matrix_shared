import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mirnylib
import mirnylib.plotting
from mirnylib.genome import Genome
import cooler
from hiclib.fragmentHiC import HiCdataset
import json
from multiprocessing import Pool


def which_enzyme(name):
    enz_list = ['MboI','NcoI','DpnII','MflI','NheI','DdeI','HindIII']
    for enz in enz_list:
        if name.find(enz)!=-1:
            return enz
    return 'HindIII'


def save_scaling(series):
    filepath = series['filepath']
    gen_name = series['gen_name']
    savepath = series['savepath']
    fragsave = savepath.replace('.hdf5','.frag')
    print(filepath)
    
    genome = Genome('/net/levsha/share/lab/genomes/'+gen_name)
        
    name = filepath[(filepath.rfind('/')+1):filepath.find('.frag')]
    enz = which_enzyme(name)
    
    fragdata = HiCdataset(fragsave, genome, enzymeName=enz, inMemory=False)
    fragdata.load(filepath)
    
    try:
        b,v = fragdata.plotScaling(plot=False)
        scaling_df = pd.DataFrame({'bins': b.tolist(), 'prob': v.tolist()})
    except:
        print('Error in calculating scaling for '+filepath+' ...Continuing')
        return None
    
    if savepath is not None:
        store = pd.HDFStore(savepath)
        store.put('scaling', scaling_df, format='table', data_columns=True)
        store.close()
    os.remove(fragsave)
    
    return scaling_df
        
    
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
    else:
        print('The following coolers were not balanced and hence were dropped:\n')
        print(bad_fp)
    
    return new_fp


def extract_info(file_df, genomes, savepaths):
    
    file_df['subfolder'] = file_df['filepath'].map(lambda x: x.replace(x[(x.rfind('/')+1):],'').replace('/net/levsha/share/lab/','') if x.find('levsha/share/lab/')!=-1 else x.replace(x[(x.rfind('/')+1):],'').replace('/net/levsha/hic/',''))

    file_df['name'] = file_df['filepath'].apply(lambda x: x[x.rfind('/')+1:x.find('.frag')])

    if genomes is None:
        file_df['gen_name'] = file_df['filepath'].map(lambda x: 'mm10' if x.find('mm10') != -1 else 'mm9' if x.find('mm9') != -1 else 'hg19')
    elif isinstance(genomes, str):
        genomes = [genomes]*len(file_df)
        file_df['gen_name'] = pd.Series(genomes)
        
    elif len(genomes) != len(filepaths):
        raise Exception("Number of genomes don't match number of filepaths")
        
    else:
        file_df['gen_name'] = pd.Series(genomes)

        
    if savepaths is None:
        savepaths = [None]*len(file_df)
        file_df['savepath'] = pd.Series(savepaths)
        
    if isinstance(savepaths, str):
        if savepaths[-1] == '/':
            savepaths = savepaths[0:-1]
        file_df['savepath'] = file_df[['subfolder','name']].apply(lambda x:savepaths+'/'+x['subfolder']+'scalings/'+x['name']+'.json', axis=1)
            
    elif len(savepaths) != len(filepaths):
        raise Exception("Number of savepaths don't match number of filepaths")
         
    else:
        file_df['savepath'] = pd.Series(savepaths)
            
    del file_df['name'] 
    del file_df['subfolder']
    
    
    return file_df


def create_dir(datasave):
    data_dir = os.path.dirname(datasave)
    if not os.path.exists(data_dir):
        try:
            os.makedirs(data_dir)
        except FileExistsError:
            print(datasave+' already exists')
            
            
def generate_scalings(filepaths, genomes=None, savepaths=None):
    if len(filepaths) == 0:
        raise Exception('Empty file list')
    elif isinstance(filepaths, str):
        filepaths = [filepaths]
  
    file_df = pd.DataFrame({'filepath':filepaths})
    file_df = extract_info(file_df, genomes, savepaths)
    if savepaths is not None:
        file_df['savepath'].apply(create_dir)
#    data = []
    args = [row for ind, row in list(file_df.iterrows())]
#    for arg in args:
#        data.append(save_scaling(arg))
    with Pool(20) as p:
        data = list(p.map(save_scaling,args))
        
    file_df['data'] = pd.Series(data)
    file_df = file_df[(file_df['data'] != None)]
    
    return file_df