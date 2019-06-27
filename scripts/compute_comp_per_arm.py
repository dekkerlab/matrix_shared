import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocess as mp

#import mirnylib
#from mirnylib.genome import Genome
import cooler
import cooltools.eigdecomp as eigdecomp
import cooltools.saddle as saddle
import cooltools.expected as expected

import bioframe
from bioframe.tools import binnify, frac_gene_coverage, frac_gc
from bioframe.io.formats import load_fasta
from bioframe.io.formats import to_bigwig

# %matplotlib notebook
import cooltools
##
# Function that gets genomic information to phase eigenvectors
def gene_content(genome, binsize, gc=True, fasta_path=None):
    
    chrom_sizes = bioframe.fetch_chromsizes(genome)
    chrom_table = binnify(chrom_sizes, binsize)

    gene_count = frac_gene_coverage(chrom_table, genome)
    if gc:
        assert fasta_path is not None, 'Please provide valid fasta file path if you want GC content'
        fasta_records = load_fasta(fasta_path) 
        gene_count['frac_gc'] = frac_gc(chrom_table, fasta_records)
    
    return gene_count
##

# Function to return list of extents of chromosomes arms
def get_chromosome_arms(genome, exclude=None):
    # Uses bioframe to get chromosomal regions
    
    if exclude is not None:
        if isinstance(exclude, str):
            exclude = [exclude]
        exclude = [str(item) for item in exclude]
    else:
        exclude = []
    
    try:
        chromlengths = bioframe.fetch_chromsizes(genome)
        centromeres = bioframe.fetch_centromeres(genome).set_index('chrom')
    except:
        print(f'Information for genome {genome} could not be found.')
        return None
    
    arms = []
    for chrom, length in chromlengths.iteritems():
        if chrom in exclude:
            continue
            
        if chrom in centromeres.index:
            mid = centromeres.loc[chrom, 'mid']
            arms.append((chrom, 0, mid))
            arms.append((chrom, mid, length))
        else:
            arms.append((chrom, 0, length))

    return arms

##

# Initializing library vairables
genome = 'hg38'
# fasta_path = f'/net/levsha/share/lab/genomes/{genome}/{genome}.fa'
fasta_path='/nl/umw_job_dekker/cshare/reference/fasta/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz'
# /nl/umw_job_dekker/cshare/reference/sorted_chromsizes/hg38.reduced.chrom.sizes
res=sys.argv[2]
# res=50000
cool_path = sys.argv[1]+'::/resolutions/'+res


file_path = sys.argv[1].split('/')
file_name_tmp = file_path[-1]
file_name=file_name_tmp.split('_')[0]
print(file_name)
out_prefix= file_name


arms = get_chromosome_arms(genome, exclude=['chrX','chrY','chrM'])


# Getting parameters for cooler_cis_eig
cool = cooler.Cooler(cool_path)
resolution = cool.info['bin-size']
genes = gene_content(genome, resolution, gc=False, fasta_path=None)

supports = {'cis': get_chromosome_arms(genome, exclude=['chrX','chrY','chrM']),
#             lams.dropna().index.map(lambda x: (x[0:x.find(':')], 
#                                                   x[x.find(':')+1:x.find('-')],
#                                                   x[x.find('-')+1:])).unique().values,
            'trans': [(chrom, 0, cool.chromsizes[chrom]) 
                               for chrom in cool.chromnames[0:22]]}
#supports

# Computing eigenvealues and eigenvectors
lams, vectors = eigdecomp.cooler_cis_eig(cool, genes, regions=supports['cis'], 
                                                   phasing_track_col='gene_count', 
                                                   sort_metric='spearmanr')
# cooler_cis_eig sorts eigenvectors by decreasing Spearman correlation with the phasing track 
# (gene count in this case). In the past, people have used the eigenvector associated with the max 
# eigenvalue. It is worth considering situations where these two sorting process differ.

# Output
contact_type='cis'
lams.to_csv(out_prefix + '.' + contact_type +'.'+  res+ '.lam.txt', sep='\t', index=False)
vectors.to_csv(out_prefix + '.' + contact_type +'.'+ res+'.vecs.tsv', sep='\t', index=False)
bioframe.to_bigwig(vectors,cool.chromsizes,out_prefix +'.' + contact_type + '.'+ res+ '.bw',value_field='E1')


exp = {}
#Cis Expected
with mp.Pool(10) as p:
    result = expected.diagsum(cool, supports['cis'],
                    transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']},
                    ignore_diags=2, map=p.map)
result = pd.concat([result[arm] for arm in supports['cis']], 
                   keys=[(arm[0], f'{arm[0]}:{arm[1]}-{arm[2]}') for arm in supports['cis']], 
                   names=['chrom', 'region'])
result = result.groupby(['chrom','diag']).sum()
result['balanced.avg'] = result['balanced.sum'] / result['n_valid']
result = result.reset_index()
exp['cis'] = result

#Trans Expected
with mp.Pool(10) as p:
    result = expected.blocksum_pairwise(cool, supports['trans'],
                transforms={'balanced': lambda p: p['count']*p['weight1']*p['weight2']},
                map=p.map)
result = pd.DataFrame([{'chrom1': s1[0], 'chrom2': s2[0], **rec} for (s1, s2), rec in result.items()],
                            columns=['chrom1', 'chrom2', 'n_valid', 'count.sum', 'balanced.sum'])
result['balanced.avg'] = result['balanced.sum'] / result['n_valid']
exp['trans'] = result


## saddle plots based on cooltools CLI

qlo, qhi = 0.02, 0.98
n_bins = 10
q_edges = np.linspace(qlo, qhi, n_bins)
binedges = saddle.quantile(vectors['E1'], q_edges)
getmatrix = {'cis':saddle.make_cis_obsexp_fetcher(cool, (exp['cis'], 'balanced.avg')),
             'trans':saddle.make_trans_obsexp_fetcher(cool, (exp['trans'], 'balanced.avg'))}

saddleplots = {}
for contact_type in ['cis','trans']:
    digitized, hist = saddle.digitize_track(binedges, track=(vectors, 'E1'), regions=supports[contact_type])
    # Digitize track in the current cooltools introduces some duplicates that causes errors later on.
    digitized = digitized.drop_duplicates(subset=['chrom','start','end'])

    S, C = saddle.make_saddle(getmatrix[contact_type], binedges, (digitized, 'E1.d'), 
                              regions=supports[contact_type], contact_type=contact_type, verbose=False)
    
    saddleplots[contact_type] = np.log2(S/C)

## vmin and vmax 
# per chromosome 


fig, ax = plt.subplots(ncols=2)
for i, contact_type in enumerate(['cis','trans']):
    ax[i].imshow(saddleplots[contact_type], cmap='coolwarm',vmin=-0.5, vmax=0.5)
    ax[i].set_title(f'Saddleplot for {contact_type} data \n using eigenvectors computed \n in cis per arm')
#plt.colorbar()
plt.savefig(out_prefix+res+".png")

def find_lam_discrepencies(lams):
    discrepencies = []
    for reg, lambdas in lams.iterrows():
        
        srtd_idx = np.array([0,1,2])
        if not np.any(np.isnan(lambdas.values)):
            srtd_idx = np.argsort(-np.abs(lambdas.values))
        
        if np.all(srtd_idx != [0,1,2]):
            discrepencies.append(reg)

    return discrepencies

discrepencies = find_lam_discrepencies(lams)
assert len(discrepencies) == 0, f'Discrepencies detections in following regions: {discrepencies}'













