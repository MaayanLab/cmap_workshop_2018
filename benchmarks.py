import os, sys, json
from itertools import combinations
import h5py
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import distance
from sklearn import metrics
from joblib import delayed, Parallel

from plots import *

def _gesa_enrichment_score(ranks_s):
    '''Calculate enrichment score from a rank ordered boolean array.
    ranks_s: np.array([0., 1., 0., 0.])
        - 1.: hits
        - 0.: misses
    '''
    n_hits = ranks_s.sum()
    n_misses = ranks_s.shape[0] - n_hits
    
    p_hit = np.cumsum(ranks_s) / n_hits
    p_miss = np.cumsum(1 - ranks_s) / n_misses
    p_diff = np.absolute(p_hit - p_miss)
    idx = np.argmax(p_diff)
    es = p_hit[idx] - p_miss[idx]
    return es
    
def gsea_score(sig1, sig2, n_sig=100):
    '''GSEA-based Kolmogorov-Smirnov statsitics.
    n_sig: number of top ranked genes to be treated as significant
    '''
    # number of genes
    n = len(sig1)
    # rank genes in sig1 (0: most down gene, 977: most up genes)
    ranks1 = stats.rankdata(sig1) - 1 
    # identify top up/down genes in sig1
    sig1_down = ranks1 < n_sig
    sig1_up = ranks1 > (n-1-n_sig)
    # argsort sig2
    sig2_srt_idx = sig2.argsort()
    # Compute ES: sig1 as query, sig2 as ref rank
    es_up1 = _gesa_enrichment_score( sig1_up[sig2_srt_idx].astype(float) )
    es_down1 = _gesa_enrichment_score( sig1_down[sig2_srt_idx].astype(float) )
    
    # rank genes in sig2
    ranks2 = stats.rankdata(sig2) - 1
    # identify top up/down genes in sig2
    sig2_down = ranks2 < n_sig
    sig2_up = ranks2 > (n-1-n_sig)
    # argsort sig1
    sig1_srt_idx = sig1.argsort()
    # Compute ES: sig2 as query, sig1 as ref rank
    es_up2 = _gesa_enrichment_score( sig2_up[sig1_srt_idx].astype(float) )
    es_down2 = _gesa_enrichment_score( sig2_down[sig1_srt_idx].astype(float) )
        
    # es_up is using up gene set to find hits in a list ascending ordered, 
    # therefore, the desirable sign should be negative
    score = (es_down1 - es_up1 + es_down2 - es_up2) / 4. 
    return score

def cosine_sim(sig1, sig2):
    '''Cosine similarity'''
    return 1 - distance.cosine(sig1, sig2)

def correlation(sig1, sig2):
    '''Pearson correlation'''
    return 1 - distance.correlation(sig1, sig2)

def pscore(mat, func, n_jobs=1, **kwargs):
    '''mat is a signature by gene matrix, apply func to all pairwise signatures.
    Similar to pdist
    '''
    n = mat.shape[0]
    n_scores = n * (n-1) / 2
    scores = np.zeros(n_scores)
    c = 0
    if n_jobs == 1:
        for i, j in combinations(range(n), 2):
            scores[c] = func(mat[i], mat[j])
            c += 1
    else:
        scores = Parallel(n_jobs=n_jobs, backend='multiprocessing', verbose=10)(
            delayed(func)(mat[i], mat[j], **kwargs) for i, j in combinations(range(n), 2))
        scores = np.array(scores)
    return scores


def compute_pairwise_connectivity_scores(sig_mat, meta_df, n_jobs=1):
    '''Given a meta_df of signatures, compute the pairwise scores and return a df.'''
    sig_ids = meta_df.index.tolist()

    # scores_es50 = pscore(sig_mat, gsea_score, n_jobs=n_jobs, n_sig=50)

    scores_cosine = pscore(sig_mat, cosine_sim, n_jobs=n_jobs)
    scores_corr = pscore(sig_mat, correlation, n_jobs=n_jobs)
    
    n_scores = len(scores_cosine)
    sigs_i = [None] * n_scores
    sigs_j = [None] * n_scores
    c = 0
    for sig_i, sig_j in combinations(sig_ids, 2):
        sigs_i[c] = sig_i
        sigs_j[c] = sig_j
        c += 1

    df = pd.DataFrame({
        'sig_i': sigs_i, 'sig_j': sigs_j, 
        'cosine': scores_cosine,
        'corr': scores_corr,
        # 'ES50': scores_es50
    })
    return df


def group_scores(res_scores, meta_df, same_cell=False, batch='batch_prefix'):
    '''Group scores to get the follow 4 groups:
    - same drug same batch
    - same drug diff batch
    - diff drug same batch
    - diff drug diff batch
    '''
    d_sig_pert = dict(zip(meta_df.index, meta_df['pert_id']))
    drugs_i = np.array([d_sig_pert[s] for s in res_scores['sig_i']])
    drugs_j = np.array([d_sig_pert[s] for s in res_scores['sig_j']])
    
    if batch == 'batch_prefix':
        # Batch: CPC004
        batches_i = np.array(map(lambda x:x.split('_')[0], res_scores['sig_i']))
        batches_j = np.array(map(lambda x:x.split('_')[0], res_scores['sig_j']))
    else:
        # Batch: CPC004_MCF7_6H
        batches_i = np.array(map(lambda x:x.split(':')[0], res_scores['sig_i']))
        batches_j = np.array(map(lambda x:x.split(':')[0], res_scores['sig_j']))
        
    m_same_batch = batches_i == batches_j
    m_diff_batch = batches_i != batches_j
    
    m_same_drug = drugs_i == drugs_j
    m_diff_drug = drugs_i != drugs_j
        
    masks = (m_same_drug & m_same_batch), (m_same_drug & m_diff_batch), \
        (m_diff_drug & m_same_batch), (m_diff_drug & m_diff_batch)
    if same_cell:
        cells_i = np.array(map(lambda x:x.split('_')[1], res_scores['sig_i']))
        cells_j = np.array(map(lambda x:x.split('_')[1], res_scores['sig_j']))
        m_same_cell = cells_i == cells_j
        masks = [m & m_same_cell for m in masks]
    return masks
    

def density_plot_scores(res_scores, meta_df, same_cell=False, batch='batch_prefix'):
    fig = plt.figure(figsize=(12, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    m1, m2, m3, m4 = group_scores(res_scores, meta_df, same_cell=same_cell, batch=batch)
    print(m1.sum(), m2.sum(), m3.sum(), m4.sum())
    
    for metric, ax in zip(['cosine', 'corr'], [ax1, ax2]):
        ax = sns.distplot(res_scores.loc[m1][metric], 
                           hist=False, kde=True, ax=ax,
                           label='same drug same batch'
                          )
        ax = sns.distplot(res_scores.loc[m2][metric],
                           hist=False, kde=True, ax=ax,
                           label='same drug diff batch'
                          )
        ax = sns.distplot(res_scores.loc[m3][metric], 
                           hist=False, kde=True, ax=ax,
                           label='diff drug same batch'
                          )
        ax = sns.distplot(res_scores.loc[m4][metric],
                           hist=False, kde=True, ax=ax,
                           label='diff drug diff batch'
                          )        
        
        ax.legend(loc='upper left', prop={'size':10})
        ax.set_xlabel(metric)
        ax.set_ylabel('Density')
        ax.set_xlim([-1, 1])

    fig.tight_layout()
    
    return fig

