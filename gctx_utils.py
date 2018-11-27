import h5py
import numpy as np

def slice_matrix(gctx, cids, rids):
    '''Slice the mat by cids and rids and ensure the mat 
    is ordered by cids and rids.'''    
    all_cids = gctx['/0/META/COL/id']
    c_mask = np.in1d(all_cids, cids)
    cids_subset = all_cids[c_mask].tolist()
    c_indices = np.array([cids_subset.index(id_) 
                          for id_ in cids])

    mat = gctx['/0/DATA/0/matrix']
    submat = mat[c_mask, :][c_indices, :]
    
    all_rids = gctx['/0/META/ROW/id']
    r_mask = np.in1d(all_rids, rids)
    rids_subset = all_rids[r_mask].tolist()
    r_indices = np.array([rids_subset.index(id_) 
                          for id_ in rids])
    submat = submat[:, r_mask][:, r_indices]

    return submat

def mean_center(mat, centerby):
    '''Mean center a mat based on centerby. mat is a samples x genes matrix'''
    mat_centered = np.zeros_like(mat)
    
    for group in set(centerby):
        mask = np.in1d(centerby, [group])
        mat_centered[mask] = mat[mask] - mat[mask].mean(axis=0)
    
    return mat_centered

    