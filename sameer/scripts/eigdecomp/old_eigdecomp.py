# -*- coding: utf-8 -*-
"""
Class-free trans eigenvector decomposition on whole genome Hi-C numpy array
"""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from mirnylib import numutils


def _filter_heatmap(A, transmask, perc_top, perc_bottom):
    # Truncate trans blowouts
    lim = np.percentile(A[transmask], perc_top)
    tdata = A[transmask]
    tdata[tdata > lim] = lim
    A[transmask] = tdata

    # Remove bins with poor coverage in trans
    marg = np.sum(A, axis=0)
    marg_nz = marg[np.sum(A, axis=0) > 0]
    min_cutoff = np.percentile(marg_nz, perc_bottom)
    dropmask = (marg > 0) & (marg < min_cutoff)
    idx = np.flatnonzero(dropmask)
    A[idx, :] = 0
    A[:, idx] = 0
    return A


def _fake_cis(A, cismask):
    cismask = cismask.astype(np.int64)
    s = np.abs(np.sum(A, axis=0)) <= 1e-10
    cismask[:, s] = 2
    cismask[s, :] = 2
    numutils.fakeCisImpl(A, cismask)
    return A


def _orient_eigs(lam, vecs, gc):
    # If GC is provided reorder and change signs of E1, E2, etc. by GC 
    # correlation. If not, reorder by descending eigenvalue magnitude and 
    # don't change signs.
    if gc is not None:
        corrs = [spearmanr(gc, vec, nan_policy='omit')[0] for vec in vecs]
        signs = np.sign(corrs)
        idx = np.argsort(-np.abs(corrs))
        lam, vecs, signs = lam[idx], vecs[idx], signs[idx]
        # change signs
        for i in range(len(vecs)):
            vecs[i] = signs[i] * vecs[i]
    else:
        idx = np.argsort(-np.abs(lam))
        lam, vecs = lam[idx], vecs[idx]
    return lam, vecs


def _eig(A, k):
    # Compute eigs
    vecs, lam = numutils.zeroEIG(A, numPCs=k)

    # eigsh returns unit vectors, but just making sure
    # then rescale by sqrt(eigval)
    for j in range(len(lam)):
        vecs[j] /= np.sqrt(np.sum(vecs[j]**2))
        vecs[j] *= np.sqrt(np.abs(lam[j]))

    return lam, np.array(vecs)


def trans_eig(A, partition, k=3, perc_top=99.95, perc_bottom=1, gc=None):
    """
    Compute compartmentalization eigenvectors on trans contact data
    Parameters
    ----------
    A : 2D array
        balanced whole genome contact matrix
    partition : sequence of int
        bin offset of each contiguous region to treat separately (e.g., 
        chromosomes or chromosome arms)
    k : int
        number of eigenvectors to compute; default = 3
    perc_top : float (percentile)
        filter - clip trans blowout contacts above this cutoff; default = 99.95
    perc_bottom : float (percentile)
        filter - remove bins with trans coverage below this cutoff; default=1
    gc : 1D array, optional
        GC content per bin for reordering and orienting the primary compartment 
        eigenvector; not performed if no array is provided
    Returns
    -------
    eigenvalues, eigenvectors
    """
    if A.shape[0] != A.shape[1]:
        raise ValueError("A is not symmetric")
    
    A = np.array(A)
    A[np.isnan(A)] = 0
    n_bins = A.shape[0]
    if not (partition[0] == 0 and 
            partition[-1] == n_bins and 
            np.all(np.diff(partition) > 0)):
        raise ValueError("Not a valid partition. Must be a monotonic sequence "
                         "from 0 to {}.".format(n_bins))

    # Delete cis data and create trans mask
    extents = zip(partition[:-1], partition[1:])
    part_ids = []
    for n, (i0, i1) in enumerate(extents):
        A[i0:i1, i0:i1] = 0
        part_ids.extend([n] * (i1 - i0))
    part_ids = np.array(part_ids)
    transmask = (part_ids[:, None] != part_ids[None, :])

    # Filter heatmap
    A = _filter_heatmap(A, transmask, perc_top, perc_bottom)
    
    # Fake cis and re-balance
    A = _fake_cis(A, ~transmask)
    A = numutils.iterativeCorrection(A)[0]
    A = _fake_cis(A, ~transmask)
    A = numutils.iterativeCorrection(A)[0]
    
    # Compute eig
    Abar = A.mean()
    O = (A - Abar) / Abar
    lam, vecs = _eig(O, k)
    lam, vecs = _orient_eigs(lam, vecs, gc)
    
    return lam, vecs


def cis_eig(A, k=3, robust=True, gc=None, classic=False):
    """
    Compute compartment eigenvector on a cis matrix
    Parameters
    ----------
    A : 2D array
        balanced whole genome contact matrix
    k : int
        number of eigenvectors to compute; default = 3
    robust : bool
        Clip top 0.1 percentile and smooth first two diagonals
    gc : 1D array, optional
        GC content per bin for choosing and orienting the primary compartment 
        eigenvector; not performed if no array is provided
    classic : bool
        Do it old-school
    Returns
    -------
    eigenvalues, eigenvectors
    """
    A = np.array(A)
    A[~np.isfinite(A)] = 0
    
    mask = A.sum(axis=0) > 0

    if A.shape[0] <= 5 or mask.sum() <= 5:
        return (
            np.array([np.nan for i in range(k)]), 
            np.array([np.ones(A.shape[0]) * np.nan for i in range(k)])
        )

    if robust:
        A = np.clip(A, 0, np.percentile(A, 99.9))
        fill_value = np.mean(np.diag(A, 2) * 2) 
        for d in [-1, 0, 1]:
            numutils.fillDiagonal(A, fill_value, d) 
            A[~mask, :] = 0
            A[:, ~mask] = 0

    OE = numutils.observedOverExpected(A[mask, :][:, mask])

    if robust:
        OE = np.clip(OE, 0, np.percentile(OE, 99.9))

    if classic:
        OE = numutils.iterativeCorrection(OE)[0]
        if (~np.isfinite(OE)).sum() > 0:
            return (
                np.array([np.ones(A.shape[0]) * np.nan for i in range(k)]),
                np.array([np.nan for i in range(k)]),
            )
        # mean-centered (subtract mean)
        eigvecs_compressed, eigvals = numutils.EIG(OE, k)
    else:
        eigvecs_compressed, eigvals = numutils.EIG((OE - 1.0), k, 
            subtractMean=False, divideByMean=False)
    
    # Restore full eigs
    eigvecs = []
    for i in range(k):
        v = np.ones(mask.shape[0]) * np.nan
        v[mask] = eigvecs_compressed[i]
        eigvecs.append(v)
    eigvecs = np.array(eigvecs)

    # Orient and reorder
    eigvals, eigvecs = _orient_eigs(eigvals, eigvecs, gc)

    return eigvals, eigvecs


def cooler_cis_eigs(clr, bins, gc_col='GC', **kwargs):
    bins_grouped = bins.groupby('chrom')
    def _each(chrom):
        A = clr.matrix(balance=True).fetch(chrom)
        
        gc = None
        if gc_col in bins:
            gc = bins_grouped.get_group(chrom)[gc_col].values
        
        lam, vec = cis_eig(A, robust=True, gc=gc, **kwargs)
        
        eig = bins_grouped.get_group(chrom).copy()
        for i in range(len(vec)):
            eig['E{}'.format(i+1)] = vec[i]
        
        return lam, eig
    
    bins_chroms = list(bins_grouped.groups.keys())
    map_chroms = [chrom for chrom in clr.chromnames if chrom in bins_chroms]

    lams, vecs = zip(*map(_each, map_chroms))

    lams = pandas.DataFrame(
        index=map_chroms,
        data=np.vstack(lams),
        columns=['lam1', 'lam2', 'lam3'],
    )
    lams.index.name = 'chrom'
    vecs = pandas.concat(vecs, axis=0, ignore_index=True)
    return lams, vecs


def cooler_trans_eig(clr, bins, partition=None, gc_col='GC', **kwargs):

    if partition is None:
        partition = np.r_[ 
            [clr.offset(chrom) for chrom in clr.chromnames], len(clr.bins())]
    
    lo = partition[0]
    hi = partition[-1]
    A = clr.matrix(balance=True)[lo:hi, lo:hi]

    gc = None
    if gc_col in bins:
        gc = bins[gc_col].values
    lam, PCs = trans_eig(A, partition, gc=gc, **kwargs)

    eigs = bins.copy()
    for i in range(len(PCs)):
        eigs['E{}'.format(i+1)] = PCs[i]

    return lam, eigs

def max_eigenvector(lams, vectors):
    ''' Used to pull out eigenvector associated with largest eigenvalue. This is now redundant as 
    new cooltools code sorts lams, vectors by decreasig order of some phasing track. New lams, vectors
    generated been sorted by Pearson correlation with gene counts as a proxy for activity.'''
    
    max_lam = pd.DataFrame({'ind':[],'vals':[]})
    max_lam['ind'] = lams.dropna().apply(lambda x: str(np.argmax(abs(x)))[-1], axis=1)
    max_lam['vals'] = lams.dropna().apply(lambda x: x[np.argmax(abs(x))], axis=1)
    
    for chrom, row in max_lam.iterrows():
        v = vectors[(vectors['chrom']==chrom)]
        v['val'] = v['E'+row['ind']]
        if 'max_vec' not in locals():
            max_vec = v[['chrom', 'start', 'end', 'val']]
        else:
            max_vec = pd.concat([max_vec, v[['chrom', 'start', 'end', 'val']]])
        
    return max_lam, max_vec

def condense_eigenvector(vector):
    assert 'label' in vector.columns
    assert not np.any(vector.label.isnull())

    cluster_cond = np.logical_or(vector.chrom != vector.chrom.shift(1), vector.label != vector.label.shift(1)) 
    vector['cluster'] = np.cumsum(cluster_cond)
        
    condensed_vector = vector.groupby('cluster').agg({'chrom': lambda x: x.values[0], 'start':'min', 
                                         'end':'max', 'label': lambda x: x.values[0]})

    assert np.all(condensed_vector['start'] <= condensed_vector['end'])
    
    return condensed_vector