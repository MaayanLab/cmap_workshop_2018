from collections import Counter, OrderedDict

import numpy as np
import networkx as nx
import igraph as ig
from sklearn import preprocessing, neighbors
import matplotlib.pyplot as plt

def compute_adjcency_mat(X, metric='euclidean'):
    pdist = pairwise_distances(X, metric=metric)
    adj_mat = 1 - pdist / pdist.max()
    # remove 1's on the diagnal
    adj_mat -= np.eye(X.shape[0])
    return adj_mat


def plot_degree_distribution(G):
    """G: a nx.Graph object
    """
    fig, ax = plt.subplots(figsize=(5,5))
    
    degrees = [k for n, k in G.degree()]
    degrees = dict(Counter(degrees))
    x = degrees.keys()
    y = degrees.values()

    ax.scatter(x, y, s=10, alpha=.6)
    ax.set_xlabel('Degree', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)    
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.tight_layout()
    return ax
    

def create_knn_graph(X, k=30, metric='euclidean'):
    '''Create a graph from a data matrix (sample x features).
    '''
    adj_mat = neighbors.kneighbors_graph(X, k, mode='connectivity', metric=metric)
    G = nx.from_scipy_sparse_matrix(adj_mat)
    return G

# thresholding-Firework
def create_graph_by_threshold(adj_mat, percentile):
    triu_idx = np.tril_indices(adj_mat.shape[0], 1)
    thresold = np.percentile(adj_mat[triu_idx], percentile)
    adj_mat_ = adj_mat.copy()
    adj_mat_[adj_mat<thresold] = 0
    G = nx.from_numpy_matrix(adj_mat_)
    return G

def create_graph_by_threshold_knn(adj_mat, percentile, k=1, X=None):
    '''combine the graph from `create_graph_by_threshold` with a kNN graph.
    '''
    G_thres = create_graph_by_threshold(adj_mat, percentile)
    G_knn = create_knn_graph(X, k=k)
    return nx.compose(G_thres, G_knn)

def nx_graph_to_igraph(G):
    # convert via adjacency matrix
    g = ig.Graph.Adjacency((nx.to_numpy_matrix(G) > 0).tolist(), mode=ig.ADJ_UNDIRECTED)
    return g

def plotly_network(G, layout, meta_df, hue=None):
    N = len(layout.coords)
    edges = [e.tuple for e in list(G.es)]

    Xn=[layout[k][0] for k in range(N)]# x-coordinates of nodes
    Yn=[layout[k][1] for k in range(N)]# y-coordinates
    Zn=[layout[k][2] for k in range(N)]# z-coordinates
    Xe=[]
    Ye=[]
    Ze=[]
    for e in edges:
        Xe+=[layout[e[0]][0], layout[e[1]][0], None]# x-coordinates of edge ends
        Ye+=[layout[e[0]][1], layout[e[1]][1], None]  
        Ze+=[layout[e[0]][2], layout[e[1]][2], None]  



