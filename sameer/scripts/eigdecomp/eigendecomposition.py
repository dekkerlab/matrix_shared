import os
import sys

path = '/net/levsha/share/sameer/github/mirnylab-experimental/sameer/'
subfolders = [item for item in os.listdir(path) if item[0]!='.']
for item in subfolders:
    sys.path.insert(0,path+item)

sout = sys.stdout 
with open(os.devnull, 'w') as devnull:
    sys.stdout = devnull

import numpy as np
import pandas as pd
import mirnylib
from mirnylib.genome import Genome
import cooler
import cooltools.eigdecomp as eigdecomp
import fnmatch
from bioframe import bedslice, fetch_chromsizes, to_bigwig

import DNA_info
sys.stdout = sout


def condense_eigenvector(vector):
    assert 'label' in vector.columns
    assert not np.any(vector.label.isnull())

    cluster_cond = np.logical_or(vector.chrom != vector.chrom.shift(1), vector.label != vector.label.shift(1)) 
    vector['cluster'] = np.cumsum(cluster_cond)
        
    condensed_vector = vector.groupby('cluster').agg({'chrom': lambda x: x.values[0], 'start':'min', 
                                         'end':'max', 'label': lambda x: x.values[0]})

    assert np.all(condensed_vector['start'] <= condensed_vector['end'])
    
    return condensed_vector

def sort_by_eigenvalue(lams, vectors):
    lam_list = [] 
    vector_list = []
    for reg, lambdas in lams.iterrows():
        if fnmatch.fnmatch(reg, '*:*-*'):
            chrom = reg[0:reg.find(':')]
            start = int(reg[reg.find(':')+1:reg.find('-')])
            end = int(reg[reg.find('-')+1:])
        else:
            chrom = reg
            start, end = None, None

        if start is None and end is None:
            region_vector = vectors[vectors.chrom == chrom].copy(deep=True)
        else:
            region_vector = bedslice(vectors.groupby('chrom'), chrom, start, end)


        if np.any(np.isnan(lambdas.values)):
            srtd_idx = np.array([0,1,2])
        else:
            srtd_idx = np.argsort(-np.abs(lambdas.values))
            
        region_vector[['E1', 'E2', 'E3']] = region_vector[['E1', 'E2', 'E3']].values[:, srtd_idx]

        lam_list.append(lambdas.values[srtd_idx])
        vector_list.append(region_vector)

    sorted_vectors = pd.concat(vector_list)
    missing = [ch for ch in vectors.chrom.unique() if ch not in sorted_vectors.chrom.unique()]

    for item in missing:
        vector_list.append(vectors[vectors.chrom == item].copy(deep=True))

    sorted_lams = pd.DataFrame(data=np.concatenate(tuple(lam_list)).reshape(-1,3), columns=lams.columns)
    sorted_lams['region'] = lams.index
    sorted_lams.set_index('region', inplace=True)
    sorted_vectors = pd.concat(vector_list).drop_duplicates()

    return sorted_lams, sorted_vectors

def find_lam_discrepencies(lams):
    discrepencies = []
    for reg, lambdas in lams.iterrows():
        
        srtd_idx = np.array([0,1,2])
        if not np.any(np.isnan(lambdas.values)):
            srtd_idx = np.argsort(-np.abs(lambdas.values))
        
        if srtd_idx[0] != 0:
            discrepencies.append(reg)

    return discrepencies


def drop_unbalanced(filepaths, res):
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
            
def save_bigwig(vectors, savepath, genome, columns=['E1', 'E2', 'E3']):
    chroms = fetch_chromsizes(genome)
    for item in columns:
        save = savepath+'.{}.bw'.format(item)
        create_dir(save)
        to_bigwig(vectors, chroms, save, value_field=item)
        
def batch_process(filepaths, genome, resolution, balance_col='weight', do_trans=False, savepath=None, bigwig=False):
    if isinstance(filepaths, str):
        filepaths = [filepaths]
    assert len(filepaths) != 0, 'Empty file list.'
  
    filepaths = drop_unbalanced(filepaths, resolution)
    n = len(filepaths)
    
    cis_regions = DNA_info.get_chromosome_arms(genome)
    genes = DNA_info.gene_content(genome, resolution, gc=True)
        
    lam_list = []
    vector_list = []
    discrep_list = []
    for i, filepath in enumerate(filepaths):
        print(str(i+1)+" of "+str(n)+": "+filepath)
        cool = cooler.Cooler(filepath)
#         bins = cool.bins()[:]
        
        lams, vectors = eigdecomp.cooler_cis_eig(cool, genes, regions=cis_regions, balance=balance_col,
                                                 phasing_track_col='frac_gc')
        print('Decomposition in cis is done')
        
        if do_trans:
            trans_partition = np.r_[[cool.offset(chrom) for chrom in cool.chromnames[0:23]]] #Here we ignore chrX, chrY, chrM
            trans_lam, trans_vecs = eigdecomp.cooler_trans_eig(cool, genes, partition=trans_partition, balance=balance_col,
                                                               phasing_track_col='frac_gc')
            print('Decomposition in trans is done')

            trans_lam['region'] = 'trans'
            trans_lam.set_index('region',drop=True, inplace=True)
            lams = pd.concat([lams, trans_lam])

            trans_vecs = trans_vecs[['chrom','start','end','E1','E2','E3']]
            vectors = vectors.merge(trans_vecs, on=['chrom','start','end'], how='outer',  suffixes=['_cis','_trans'])
        
        discrep = find_lam_discrepencies(lams)
        
        lam_list.append(lams)
        vector_list.append(vectors)
        discrep_list.append(discrep)
        
        if savepath is not None:
            if '.mcool' in filepath:
                filename = filepath.split(':')[0]
            filename = '.'.join(filename.split('/')[-1].split('.')[0:-2])
            save = savepath+filename+'.hdf5'
            create_dir(save)
            store = pd.HDFStore(save)
            store.put('vectors', vectors, format='table', data_columns=True)
            store.put('lams', lams, format='table', data_columns=True)
            store.close()
            
            if not os.path.isfile(savepath+'discrepencies.txt'):
                with open(savepath+'discrepencies.txt','w+') as f:
                    f.write(f"{filename}\t{', '.join(discrep)}\n")
            else:
                with open(savepath+'discrepencies.txt','a+') as f:
                    f.write(f"{filename}\t{', '.join(discrep)}\n")
                            
            print('Saved to: '+save, '\n')
            if bigwig:
                wig_save = savepath+'bigwigs/{}'.format(filename)
                save_bigwig(vectors, wig_save, genome)
        
        print('Eigendecomposition for '+filepath+': DONE!\n')
        
    return lam_list, vector_list, discrep_list
    
if __name__ == "__main__":
    genome = sys.argv[1]
    res = int(sys.argv[2])
    savepath = sys.argv[3]
    balance_col = sys.argv[4]
    bw = bool(int(sys.argv[5]))
    trans = bool(int(sys.argv[6]))
    filepaths = sys.argv[7:]
    filepaths = [f'{file}::/resolutions/{res}' if '.mcool' in file else file for file in filepaths]
    
    batch_process(filepaths, genome, res, balance_col=balance_col, do_trans=trans, savepath=savepath, bigwig=bw)
    