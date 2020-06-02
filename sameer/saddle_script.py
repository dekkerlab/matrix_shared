import sys
import os

path = '/net/levsha/share/sameer/github/mirnylab-experimental/sameer/'
subfolders = [item for item in os.listdir(path) if item[0]!='.']
for item in subfolders:
    sys.path.insert(0,path+item)
    
import pandas as pd
import numpy as np

import cooler
import bioframe
import cooltools

import saddleplot
import eigendecomposition
import DNA_info
import matrix_manager as mm

res = 100000
genome='hg38'

qlo, qhi = 0.02, 0.98
n_bins = 10
q_edges = np.linspace(qlo, qhi, n_bins)

db = mm.Dataset()
data = db.get_tables()
data = mm.filter_data(data, filter_dict={'seq':'deep'})
# data = mm.filter_data(data, filter_dict={'celltype':'HelaS3', 'xlink':'FA', 'enzyme':'MNase'})
# data = pd.concat([df1, df2])

data = mm.get_coolers(data, res=res)
data = mm.get_compartments(data, res=res)

savepath = './saddles/'
for i, row in data.iterrows():
    name = row['lib_name']
    print(name)
    os.makedirs(f'{savepath}{name}/100000/cis', exist_ok=True)
    os.makedirs(f'{savepath}{name}/100000/trans', exist_ok=True)
    c = row['cooler_100000']
    lams = row['lams_100000']
    vector = row['vectors_100000']
    
#     supports = {'cis': DNA_info.get_chromosome_arms(genome),
#                 'trans': [(chrom, 0, c.chromsizes[chrom]) 
#                                    for chrom in c.chromnames[0:22]]}
    
    for region in DNA_info.get_chromosome_arms(genome):
        print(region)
        chrom, start, end = region
        mat = c.matrix(balance=True).fetch(region)
        vec = bioframe.bedslice(vector.groupby('chrom'), chrom, start, end)
        vec = vec['E1_cis'].values
        if np.all(np.isnan(vec)):
            continue
        if len(vec) == 0:
            continue
        S, C = saddleplot.construct_cis_saddleplot(mat, vec, num_percentile=20)
        
        np.save(f'{savepath}{name}/100000/cis/{chrom}:{start}-{end}.npy', np.dstack((S,C)))
        
    for i in np.arange(len(c.chromnames[0:22])):
        for j in np.arange(i+1, len(c.chromnames[0:22])):
            
            chrom1 = c.chromnames[i]
            start1, end1 = 0, c.chromsizes[chrom1]
            chrom2 = c.chromnames[j]
            start2, end2 = 0, c.chromsizes[chrom2]
            print((chrom1,start1,end1), (chrom2,start2,end2))
            
            mat = c.matrix(balance=True).fetch((chrom1, start1, end1), (chrom2, start2, end2))
            
            vec1 = bioframe.bedslice(vector.groupby('chrom'), chrom1, start1, end1)
            vec1 = vec1['E1_cis'].values
            vec2 = bioframe.bedslice(vector.groupby('chrom'), chrom2, start2, end2)
            vec2 = vec2['E1_cis'].values

            if np.all(np.isnan(vec1)) or np.all(np.isnan(vec2)):
                continue
            if (len(vec1) == 0) or  (len(vec2) == 0):
                continue
            S, C = saddleplot.construct_trans_saddleplot(mat, vec1, vec2, num_percentile=20)
              
            np.save(f'{savepath}{name}/100000/trans/{chrom1}-{chrom2}.npy', np.dstack((S,C)))
            
        