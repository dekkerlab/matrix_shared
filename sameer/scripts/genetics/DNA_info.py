import numpy as np
import pandas as pd
import bioframe
from bioframe import binnify, frac_gene_coverage, frac_gc
from bioframe import load_fasta
from mirnylib.genome import Genome



def gene_content(genome, binsize, gc=True):
    
    chrom_sizes = bioframe.fetch_chromsizes(genome)
    chrom_table = binnify(chrom_sizes, binsize)

    gene_count = frac_gene_coverage(chrom_table, genome)
    if gc:
        fasta_path = f'/net/levsha/share/lab/genomes/{genome}/{genome}.fa'    
        fasta_records = load_fasta(fasta_path) 
        gene_count['frac_gc'] = frac_gc(chrom_table, fasta_records)
    
    return gene_count

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

def process_centromere_file(genome):
    df = {'chrom':[],'start':[],'end':[]}
    with open(f'/net/levsha/share/lab/genomes/{genome}/centromeres.txt', 'r') as f:
        while True:
            line = f.readline()
            if line == '':
                break
            line = line.split('\t')
            df['chrom'].append(line[1])
            df['start'].append(line[2])
            df['end'].append(line[3])
    df = pd.DataFrame(df)

    gb = df.groupby('chrom')
    cent_dict = {}
    for chrom, group in gb:
        start = np.min(group['start'].values.astype(int))
        end = np.max(group['end'].values.astype(int))
        cent_dict[chrom] = (start+end)//2
    cent_dict['chrM'] = -1
    
    return pd.Series(cent_dict)