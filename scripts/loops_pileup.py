# bsub -q short -W 4:00 -R "rusage[mem=50000]" -oo multiple_dot_lists.out -eo multiple_dot_lists.err 'python multiple_dot_lists.py'


# %matplotlib inline
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')
import multiprocess as mp
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler
#import bbi
from cooltools import snipping
import sys

def pileup_multiple_dot_lists(cool_file,dot_file_list, exp_cool,resolution,flank,anchor_dist,anchor_flank,pileup_name):
    i=0
    filename1=cool_file[0].split("/")[-2].split("_hg38")[0]
    filename2=cool_file[1].split("/")[-2].split("_hg38")[0]
    filename3=cool_file[2].split("/")[-2].split("_hg38")[0]
    
    cool = [filename1,filename2,filename3]
    exp_cool = [exp_cool[0], exp_cool[1], exp_cool[2]]
    conditions = ['HiC-FA-DpnII', 'HiC-DSG-DpnII','MicroC-DSG-MNase']

    print(filename1)
    print(filename2)
    print(filename3)
        
    resolution=resolution
    flank = flank
    #resolution=sys.argv[4]
    hg38 = bioframe.fetch_chromsizes('hg38')
    chromsizes = bioframe.fetch_chromsizes('hg38')
    chromosomes = list(chromsizes.index)
    binsize = resolution

    cooler_paths = {    
    'HiC-FA-DpnII' : cool_file[0],
    'HiC-DSG-DpnII' : cool_file[1],
    'MicroC-DSG-MNase' : cool_file[2],
    }
    
    exp_paths = {    
    'HiC-FA-DpnII' : exp_cool[0],
    'HiC-DSG-DpnII'  : exp_cool[1],
    'MicroC-DSG-MNase' : exp_cool[2],
    }

    
    
    long_names = {
    'HiC-FA-DpnII': 'HiC-FA-DpnII',
    'HiC-DSG-DpnII': 'HiC-DSG-DpnII',
    'MicroC-DSG-MNase': 'MicroC-DSG-MNase',
    }

    pal = sns.color_palette('colorblind')
    colors = {
        filename1: pal[0],
        filename2 : '#333333',
        filename3: pal[2],
    }

    clrs = {
    cond: cooler.Cooler(cooler_paths[cond]) for cond in conditions
    }

    anchor_dist = anchor_dist
    anchor_flank = flank
    # dot file list
    gs = plt.GridSpec(nrows=len(conditions), ncols=len(dot_file_list) + 1)
    plt.figure(figsize=(6 * len(conditions)+1, 7))
            
    for dot_file in dot_file_list:
        print(dot_file)
        sites = pd.read_table(dot_file)
        mid1=(sites['start1']+sites['end1'])/2
        mid2=(sites['start2']+sites['end2'])/2
        new_file=pd.DataFrame()
        new_file = pd.concat([sites['chrom1'],mid1,sites['chrom2'],mid2],axis=1)

        # "convergent" orientation of paired CTCF motifs
        # sites = sites[(sites['strand1'] == '+') & (sites['strand2'] == '-')] ## not working 
        new_file.columns=['chrom1','mid1','chrom2','mid2']
        print(len(new_file))
        new_file.head()
        supports = [(chrom, 0, chromsizes[chrom]) for chrom in chromosomes]

        snippet_flank = flank
        windows1 = snipping.make_bin_aligned_windows(
            binsize, 
            new_file['chrom1'], 
            new_file['mid1'],
            flank_bp=snippet_flank)
        # windows1['strand'] = sites['strand1']
        windows2 = snipping.make_bin_aligned_windows(
            binsize, 
            new_file['chrom2'], 
            new_file['mid2'],
            flank_bp=snippet_flank)
        windows = pd.merge(windows1, windows2, left_index=True, right_index=True, suffixes=('1', '2'))
        windows = snipping.assign_regions(windows, supports)
        windows = windows.dropna()
        windows.head()
        #stacks = {}
        piles = {}
        k=0        
        for cond in conditions:
            expected = pd.read_table(exp_paths[cond])
            snipper = snipping.ObsExpSnipper(clrs[cond], expected)
            print(snipper)
            stack = snipping.pileup(windows, snipper.select, snipper.snip)
            #stacks[cond] = stack
            piles[cond] = np.nanmean(stack, axis=2)

            opts = dict(
                vmin=-0.25,
                vmax=0.25,
                extent=[-flank//1000, flank//1000, -flank//1000, flank//1000],
                cmap='coolwarm'
            )

            ax = plt.subplot(gs[k,i])
            img = ax.matshow(
                np.log10(piles[cond]), #piles[cond]),  
                **opts)
            #plt.title(dot_name_list[i],fontsize=7)
            #ax.xaxis.tick_bottom()
            k=k+1
            ax.yaxis.set_visible(True)
            ax.xaxis.set_visible(False)
            if k > 0:
                ax.yaxis.set_visible(True)
                ax.xaxis.set_visible(False)
                ax = plt.subplot(gs[len(conditions)])
                #plt.suptitle(f'Dot calls ({anchor_dist//1000} +/- {anchor_flank//1000})kb apart\n'
                #                 f'Hi-C resolution = {binsize//1000}kb; # of pairs = {len(windows)}')
        #plt.title(dot_name_list[i])
        i=i+1
    #plt.colorbar(img, cax=ax)
    plt.savefig(pileup_name)

out_path="/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/"
import glob, os
esc_coolers_path=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-ESC4DN-FA-DpnII-R1-R2_hg38.mapq_30.1000.mcool::resolutions/5000",
                 "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-ESC4DN-DSG-DpnII-R1-R2_hg38.mapq_30.1000.mcool::resolutions/5000",
                 "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38.mapq_30.1000.mcool::resolutions/5000"]

esc_exp_path=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-ESC4DN-FA-DpnII-R1-R2_hg38_5000.cis.expected",
              "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-ESC4DN-DSG-DpnII-R1-R2_hg38_5000.cis.expected",
              "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38_5000.cis.expected"]

dot_file_list_esc=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-ESC4DN-DSG-DpnII-R1-R2_hg38_U54_ESC4DN_FA_DpnII_R1_R2_hg38_uniq_comp_to_U54_H1ESC4DN_FA_DSG_MNase_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54_H1ESC4DN_FA_DSG_MNase_R1_R2_hg38_U54_ESC4DN_FA_DpnII_R1_R2_hg38_uniq_comp_to_U54-ESC4DN-DSG-DpnII-R1-R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-ESC4DN-FA-DpnII-R1-R2_hg38_uniq_comp_to_U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38_union_cloops_U54-ESC4DN-DSG-DpnII-R1-R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-ESC4DN-DSG-DpnII-R1-R2_hg38_uniq_comp_to_U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38_union_cloops_U54-ESC4DN-FA-DpnII-R1-R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38_U54_ESC4DN_FA_DpnII_R1_R2_hg38_U54_ESC4DN_DSG_DpnII_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-ESC4DN-DSG-DpnII-R1-R2_hg38_U54_H1ESC4DN_FA_DSG_MNase_R1_R2_hg38_uniq_comp_to_U54_ESC4DN_FA_DpnII_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-H1ESC4DN-FA-DSG-MNase-R1-R2_hg38_uniq_comp_to_U54-ESC4DN-DSG-DpnII-R1-R2_hg38_union_cloops_U54-ESC4DN-FA-DpnII-R1-R2_hg38.txt"]


#pileup_multiple_dot_lists(esc_coolers_path,dot_file_list_esc,esc_exp_path,5000,300000,15000,10000,out_path+"dot_pileup_all_ESC.pdf")
pileup_multiple_dot_lists(esc_coolers_path,dot_file_list_esc,esc_exp_path,5000,200000,75000,200000,out_path+"dot_pileup_all_ESC_200kb.pdf")

#pileup_multiple_dot_lists(cool_file,dot_file_list,exp_cool,resolution,flank,anchor_dist,anchor_flank,pileup_name):



out_path="/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/"
import glob, os
hff_coolers_path=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-HFFc6-FA-DpnII-R1-R2_hg38.mapq_30.1000.mcool::resolutions/5000",
                 "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-HFFc6-DSG-DpnII-R1-R2_hg38.mapq_30.1000.mcool::resolutions/5000",
                 "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/U54-HFFc6-FA-DSG-MNase-R1-R3.hg38.mapq_30.1000.mcool::resolutions/5000"]

hff_exp_path=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-HFFc6-FA-DpnII-R1-R2_hg38_5000.cis.expected",
              "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-HFFc6-DSG-DpnII-R1-R2_hg38_5000.cis.expected",
              "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/expected/U54-HFFc6-FA-DSG-MNase-R1-R3_5000.cis.expected"]

dot_name_list=[ "MicroC-DSG-MNase ∩ (HiC-DSG-DpnII - HiC-FA-DpnII)",
                "HiC-FA-DpnII - (MicroC-DSG-MNase ∪ HiC-DSG-DpnII)",
                "MicroC-DSG-MNase ∩ (HiC-FA-DpnII - HiC-DSG-DpnII)",
                "HiC-DSG-DpnII - (MicroC-DSG-MNase ∪ HiC-FA-DpnII)",
                "HiC-FA-DpnII ∩ (HiC-DSG-DpnII - MicroC-DSG-MNase)",
                "MicroC-DSG-MNase ∩ (HiC-FA-DpnII ∩ HiC-DSG-DpnII)",
                "MicroC-DSG-MNase-(HiC-DSG-DpnII ∪ HiC-FA-DpnII)"]

dot_file_list_hff=["/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DpnII-R1-R2_hg38_U54_HFFc6_DSG_DpnII_R1_R2_hg38_uniq_comp_to_U54_HFFc6_FA_DSG_MNase_R1_R3.hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DpnII-R1-R2_hg38_uniq_comp_to_U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_union_cloops_U54-HFFc6-DSG-DpnII-R1-R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_U54_HFFc6_FA_DpnII_R1_R2_hg38_uniq_comp_to_U54_HFFc6_DSG_DpnII_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-DSG-DpnII-R1-R2_hg38_uniq_comp_to_U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_union_cloops_U54-HFFc6-FA-DpnII-R1-R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_U54_HFFc6_DSG_DpnII_R1_R2_hg38_uniq_comp_to_U54_HFFc6_FA_DpnII_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_U54_HFFc6_FA_DpnII_R1_R2_hg38_U54_HFFc6_DSG_DpnII_R1_R2_hg38.txt",
                "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots/for_pileup/U54-HFFc6-FA-DSG-MNase-R1-R3.hg38_uniq_comp_to_U54-HFFc6-DSG-DpnII-R1-R2_hg38_union_cloops_U54-HFFc6-FA-DpnII-R1-R2_hg38.txt"]




#pileup_multiple_dot_lists(hff_coolers_path,dot_file_list_hff,hff_exp_path,5000,300000,150000,10000,out_path+"dot_pileup_all_HFF.pdf")
#pileup_multiple_dot_lists(cool_file,dot_file_list,exp_cool,resolution,flank,anchor_dist,anchor_flank,pileup_name):
pileup_multiple_dot_lists(hff_coolers_path,dot_file_list_hff,hff_exp_path,5000,200000,75000,200000,out_path+"dot_pileup_all_HFF_200kb.pdf")



