from __future__ import division
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats.mstats import zscore
from joblib import delayed, Parallel

import gctx_utils

def chdir(data, sampleclass):
    m1 = sampleclass == 1
    m2 = sampleclass == 2
    
    gamma = 0.5
    
    data = zscore(data)
    
    ## start to compute
    n1 = m1.sum() # number of controls
    n2 = m2.sum() # number of experiments

    ## the difference between experiment mean vector and control mean vector.
    meanvec = data[:,m2].mean(axis=1) - data[:,m1].mean(axis=1) 

    ## initialize the pca object
    pca = PCA(n_components=None)
    pca.fit(data.T)

    ## compute the number of PCs to keep
    cumsum = pca.explained_variance_ratio_ # explained variance of each PC
    keepPC = len(cumsum[cumsum > 0.001]) # number of PCs to keep
    v = pca.components_[0:keepPC].T # rotated data 
    r = pca.transform(data.T)[:,0:keepPC] # transformed data

    dd = ( np.dot(r[m1].T,r[m1]) + np.dot(r[m2].T,r[m2]) ) / float(n1+n2-2) # covariance
    sigma = np.mean(np.diag(dd)) # the scalar covariance

    shrunkMats = np.linalg.inv(gamma*dd + sigma*(1-gamma)*np.eye(keepPC))

    b = np.dot(v, np.dot(np.dot(v.T, meanvec), shrunkMats))

    b /= np.linalg.norm(b) # normalize b to unit vector
    
    return b


def compute_chdir_signature(sig_id, distil_ids_pert, distil_ids_sub, mat, mat_centered):
    # distil_ids_pert = row['distil_id']
    # Make the sample_class
    mask_pert = np.in1d(distil_ids_sub, distil_ids_pert)
    sample_class = mask_pert.astype(int) + 1
    # Apply CD on the original mat
    cd_coefs = chdir(mat.T, sample_class)
    # Apply CD on the mean centered mat
    cd_coefs_centered = chdir(mat_centered.T, sample_class)
    # Averaging profiles after mean centering
    avg_vals = mat_centered[mask_pert].mean(axis=0)

    doc = {
        'sig_id': sig_id,
        'CD_noncenter_LM': cd_coefs,
        'CD_center_LM': cd_coefs_centered,
        'avg_center_LM': avg_vals,
    }
    return doc




