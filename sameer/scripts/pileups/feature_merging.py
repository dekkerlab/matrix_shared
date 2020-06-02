import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from itertools import chain
import cooler
import bioframe

def bedpeslice(df, chrom, start, end):

    index = df.index.values
    df_l = df[['chrom1', 'start1', 'end1']].rename(columns=lambda x: x[0:-1])
    gb_l = df_l.groupby('chrom')
    subset_l = bioframe.bedslice(gb_l, chrom, start, end)
    index_l = subset_l.index.values
    mask_l = np.isin(index, index_l)

    df_r = df[['chrom2', 'start2', 'end2']].rename(columns=lambda x: x[0:-1])
    gb_r = df_r.groupby('chrom')
    subset_r = bioframe.bedslice(gb_r, chrom, start, end)
    index_r = subset_r.index.values
    mask_r = np.isin(index, index_r)

    sliced_df = df.loc[(mask_l & mask_r)].copy()
    
    return sliced_df


def union(df1, df2, radius=10000):

    gb1 = df1.groupby('chrom1')
    gb2 = df2.groupby('chrom1')

    unique_df2 = []
    for chrom, group1 in gb1:
        
        group1 = group1.reset_index(drop=True)
        points1 = group1[['pos1','pos2']].values
        kdt1 = cKDTree(points1, copy_data=True)

        group2 = gb2.get_group(chrom).reset_index(drop=True)
        idx2 = group2.index.values
        points2 = group2[['pos1','pos2']].values
        kdt2 = cKDTree(points2, copy_data=True)

        common_idx2 = kdt1.query_ball_tree(kdt2, r=radius, p=np.inf)
        common_idx2 = list(chain.from_iterable(common_idx2))
        common_idx2 = np.unique(common_idx2)

        unique_mask = ~np.isin(idx2, common_idx2)
        unique_df2.append(group2.loc[unique_mask])

    unique_df2 = pd.concat(unique_df2)

    union_df = pd.concat([df1, unique_df2]).reset_index(drop=True)
    union_df = union_df.sort_values(['chrom1','start1','end1','chrom2','start2','end2'])
    
    return union_df


def intersection(df1, df2, radius=10000, return_unique=False):

    gb1 = df1.groupby('chrom1')
    gb2 = df2.groupby('chrom1')

    common_df1 = []
    unique_df1 = []
    for chrom, group1 in gb1:

        group1 = group1.reset_index(drop=True)
        idx1 = group1.index.values
        points1 = group1[['pos1','pos2']].values
        kdt1 = cKDTree(points1, copy_data=True)

        group2 = gb2.get_group(chrom).reset_index(drop=True)
        points2 = group2[['pos1','pos2']].values
        kdt2 = cKDTree(points2, copy_data=True)

        common_idx1 = kdt2.query_ball_tree(kdt1, r=radius, p=np.inf)
        common_idx1 = list(chain.from_iterable(common_idx1))
        common_idx1 = np.unique(common_idx1)

        common_mask = np.isin(idx1, common_idx1)
        common_df1.append(group1.loc[common_mask])
        unique_df1.append(group1.loc[~common_mask])

    common_df1 = pd.concat(common_df1).reset_index(drop=True)
    common_df1 = common_df1.sort_values(['chrom1','start1','end1','chrom2','start2','end2'])
    
    if return_unique:
        unique_df1 = pd.concat(unique_df1).reset_index(drop=True)
        unique_df1 = unique_df1.sort_values(['chrom1','start1','end1','chrom2','start2','end2'])
        
        return common_df1, unique_df1
    
    return common_df1