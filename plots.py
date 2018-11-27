import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk', style='whitegrid')

from scipy.stats.mstats import zscore
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D

def PCA_plot(mat, meta_df, standardize=3, hue=None, style=None, size=None):
    '''
    mat: genes by samples np.array
    meta_df: pd.DataFrame of the metadata

    '''
    if standardize == 2: # standardize along rows/genes
        mat = zscore(mat, axis=1)
    elif standardize == 1: # standardize along cols/samples
        mat = zscore(mat, axis=0)
    
    pca = PCA(n_components=3)
    ## get variance captured
    pca.fit(mat.T)
    variance_explained = pca.explained_variance_ratio_[0:3]
    variance_explained *= 100
    ## compute PCA and plot
    pca_transformed = pca.transform(mat.T)
    pca_transformed_df = pd.DataFrame(pca_transformed[:, :2], 
        columns=['PC1', 'PC2'], 
        index=meta_df.index)
    pca_transformed_df = pca_transformed_df.merge(meta_df, 
        left_index=True,
        right_index=True)

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax = sns.scatterplot('PC1', 'PC2', data=pca_transformed_df, 
                         hue=hue, size=size, style=style,
        alpha=0.7, s=20, marker='o'
        )

    ax.set_xlabel('PC1 (%.2f'%variance_explained[0] + '%' + ' variance captured)', fontsize=20)
    ax.set_ylabel('PC2 (%.2f'%variance_explained[1] + '%' + ' variance captured)', fontsize=20)
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    return fig