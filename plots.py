import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk', style='whitegrid')

from scipy.stats.mstats import zscore
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D

import plotly
import plotly.plotly as py
import plotly.graph_objs as go
plotly.offline.init_notebook_mode(connected=True)

COLORS20 = [
    '#1f77b4',
    '#aec7e8',
    '#ff7f0e',
    '#ffbb78',
    '#2ca02c',
    '#98df8a',
    '#d62728',
    '#ff9896',
    '#9467bd',
    '#c5b0d5',
    '#8c564b',
    '#c49c94',
    '#e377c2',
    '#f7b6d2',
    '#7f7f7f',
    '#c7c7c7',
    '#bcbd22',
    '#dbdb8d',
    '#17becf',
    '#9edae5',
    ]

def scatter_plot(coords, meta_df, hue=None, style=None, size=None):
    df = pd.DataFrame(coords, columns=['x', 'y'], 
        index=meta_df.index)
    df = df.merge(meta_df, 
        left_index=True,
        right_index=True)

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111)
    ax = sns.scatterplot('x', 'y', data=df, 
                         hue=hue, size=size, style=style,
        alpha=0.7, s=20, marker='o'
        )    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # fig.tight_layout()
    return fig

def plotly_webgl_scatter(coords, meta_df, hue=None, label_cols=[]):
    
    n_categories = meta_df[hue].nunique()
    
    if n_categories < 11:
        colors = sns.color_palette(n_colors=n_categories).as_hex()
    else:
        # colors = sns.color_palette('husl', n_colors=n_categories).as_hex()
        colors = COLORS20
    df = pd.DataFrame(coords, columns=['x', 'y'], 
        index=meta_df.index)
    df = df.merge(meta_df, 
        left_index=True,
        right_index=True)
    
    data = []
    i = 0
    
    for category, df_grouped in df.groupby(hue):
        labels = [] 
        for index, row in df_grouped.iterrows():
            label = '%s<br>' % index
            subrow = row.loc[label_cols]
            label += '<br>'.join(map(lambda x: ': '.join(x), zip(subrow.index, subrow.astype(str))))
            labels.append(label)

        trace = go.Scattergl(
            x = df_grouped['x'].values,
            y = df_grouped['y'].values,
            name = category,
            text = labels,
            mode = 'markers',
            hoverinfo = 'text',
            marker = dict(
                line = dict(
                    width = 1,
                    color = colors[i]
                )
            )
        )
        data.append(trace)
        i += 1
        
    layout = go.Layout(
             title="",
             width=1000,
             height=1000,
             showlegend=True,
            hovermode='closest'
    )

    fig=go.Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)
    return 


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

    fig = plt.figure(figsize=(7, 7))
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
    # fig.tight_layout()
    return fig


def plotly_network(G, layout, meta_df, hue=None, label_cols=[]):
    N = len(layout.coords)
    edges = [e.tuple for e in list(G.es)]
    
    n_categories = meta_df[hue].nunique()
    
    if n_categories < 11:
        colors = sns.color_palette(n_colors=n_categories).as_hex()
    else:
        colors = COLORS20

    Xe=[]
    Ye=[]
    Ze=[]
    for e in edges:
        Xe+=[layout[e[0]][0], layout[e[1]][0], None]# x-coordinates of edge ends
        Ye+=[layout[e[0]][1], layout[e[1]][1], None]  
        Ze+=[layout[e[0]][2], layout[e[1]][2], None]  

    # edges
    trace1=go.Scatter3d(x=Xe,
                   y=Ye,
                   z=Ze,
                   mode='lines',
                   line=dict(color='rgb(125,125,125)', width=1),
                   hoverinfo='none',
                   name='edges'
                   )
    # nodes
    df = pd.DataFrame(layout.coords, columns=['x', 'y', 'z'], 
        index=meta_df.index)
    df = df.merge(meta_df, 
        left_index=True,
        right_index=True)

    data = [trace1]
    i = 0
    
    for category, df_grouped in df.groupby(hue):
        labels = [] 
        for index, row in df_grouped.iterrows():
            label = '%s<br>' % index
            subrow = row.loc[label_cols]
            label += '<br>'.join(map(lambda x: ': '.join(x), zip(subrow.index, subrow.astype(str))))
            labels.append(label)

        trace = go.Scatter3d(
            x = df_grouped['x'].values,
            y = df_grouped['y'].values,
            z = df_grouped['z'].values,
            name = category,
            text = labels,
            mode = 'markers',
            hoverinfo = 'text',
            marker = dict(
                line = dict(
                    width = 1,
                    color = colors[i]
                )
            )
        )
        data.append(trace)
        i += 1

    axis=dict(showbackground=False,
              showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title=''
              )

    layout = go.Layout(
             title="",
             width=1000,
             height=1000,
             showlegend=True,
             scene=dict(
                 xaxis=dict(axis),
                 yaxis=dict(axis),
                 zaxis=dict(axis),
            ),
         margin=dict(
            t=100
        ),
        hovermode='closest',
        annotations=[
               dict(
               showarrow=False,
                text="",
                xref='paper',
                yref='paper',
                x=0,
                y=0.1,
                xanchor='left',
                yanchor='bottom',
                font=dict(
                size=14
                )
                )
            ],    
        )

    fig=go.Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)
    return



