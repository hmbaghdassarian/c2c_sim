#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tempfile
import os
import warnings

import networkx as nx
import igraph

import scipy
import pandas as pd
import numpy as np
from collections import OrderedDict
import random

# notebook
import sys
sys.path.append("../..")
from definitions import ROOT_DIR

# # script version
# path_ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__), '../'))
# import sys
# sys.path.insert(1, path_)
# from definitions import ROOT_DIR


# In[3]:


class graph_generator():
    def __init__(self, seed = None):
        '''Initializer for generating and analyzing various networks
        
        Parameters
        ----------
        seed: int, optional
            value to set seed for graph output
        '''
        self.seed = seed
    
    def bipartite_sf(self, nodes, degrees, alpha = 2, edges = None, seed = None, 
                     check_properties = True, compare_barabasi = True):
        '''Wraps bigraph_r and retrieves properties from the resultant random, undirected, scale-free, bipartite network

        Parameters
        ----------
        nodes: int
            2*n is the number of nodes in graph. Restrictions: number of nodes in group 1 = number of nodes in group 2. 
        degrees: int
            average degree of nodes 
        alpha: float
            scale-free exponent for network degree distribution (recommended 2<alpha<3)
        edges: int, optional
            number of (directed) edges
        seed: int, optional
            set a seed for random graph output
        check_properties: bool
            checks network properties (scale-free power distribution, undirected)
        compare_barabasi: bool
            checks how closely the bipartite network matched a Barabasi alg generated network with similar input parameters

        Returns
        ----------
        B: dict
            keys: 'nx' and/or 'ig'
            values: networkx and/or igraph objects of the generated random, undirected, scale-free, bipartite network 
        node_groups: dict
            grouping of nodes into ligands (key = 'L') or receptors (key = 'R')
        fit: igraph.FittedPowerLaw
            scale-free network parameters for B_ig (p-value from Kolmogrov-Smirnov test)
        comp: pd.DataFrame or None
            summary of differences in network properties between  bipartite network and similar Barabasi network
        '''

        curdir = os.path.dirname(os.path.abspath(os.curdir)) ##### set to bigraph_r directory
        os.chdir(ROOT_DIR + '/core/')
        output = tempfile.mkstemp(suffix = '_bipartite_sf.csv', dir = os.getcwd())[1]
        beta = alpha
        cmd = 'Rscript bigraph_r.r ' + '--beta=' + str(beta) + ' --nodes=' + str(nodes) + ' --degrees=' + str(degrees)
        cmd += ' --output=' + str(output) 
        if edges is not None:
            cmd += ' --edges=' + str(edges)
        if self.seed is not None:
            cmd += ' --seed=' + str(self.seed)

        print('Generate undirected, bipartite, scale-free graph')
        os.system(cmd)

        adj_matrix = pd.read_csv(output, index_col = 0)
        try:
            os.remove(output)
        except:
            warnings.warn('Graph output file not found')
        os.chdir(curdir)
       
        adj_matrix.index = adj_matrix.index -1 # 0 indexing
        adj_matrix.columns = adj_matrix.index

        # extract information
        node_groups = {'L': list(range(nodes)), 'R': list(range(nodes, 2*nodes))}

        # igraph
        A = adj_matrix.replace(0, float('nan'), inplace = False).stack().reset_index()
        
        # note, need to double check the following:
        msg = 'Need to make sure the symmetric matrix (due to bipartite) removes the bidirectionality "
        msg += '(only take the upper triangle of A rather than the full matrix)'
        # only taking the upper triangle will ensure one column of the adjacency list is only ligands and the other receptors
        # otherwise will have both ligands and receptors present in both columns
        warnings.warn(msg)
        
        B_ig = igraph.Graph.Bipartite(edges = list(zip(A.level_0, A.level_1)), 
                              types = [True]*nodes + [False]*nodes)

        fit = self.power_fit(B_ig)
        if check_properties:
            print('Check network properties')
            perfect = True
            if fit.p < 0.05:
                perfect = False
                warnings.warn('Network is not scale-free')

            change = alpha/fit.alpha 
            if (change > 2) or (change < 0.5):
                perfect = False
                warnings.warn('Power law alpha fit is off by more than 2-fold')

            if igraph.Graph.is_directed(B_ig):
                perfect = False
                warnings.warn('Netork is directed, expected undirected')

            if perfect:
                print('All properties are as expected')

        comp = None
        if compare_barabasi:
            if isinstance(degrees,float):
                degrees = int(round(degrees))
            Bt = igraph.Graph.Barabasi(n = B_ig.vcount(), m = degrees, directed = False)
            fit_t = self.power_fit(Bt) 

            comp = pd.DataFrame(data = {'Bipartite': [B_ig.vcount(), B_ig.ecount(), fit.alpha, fit.p ],
                         'Barabasi': [Bt.vcount(), Bt.ecount(), fit_t.alpha, fit_t.p ]}, 
                index = ['nodes', 'edges', 'fitted alpha', 'p-val'])
            comp['change (Barabasi/Bipartite)'] = comp.Barabasi/comp.Bipartite
        
        B = dict({'ig': B_ig})
        B_nx = nx.bipartite.from_biadjacency_matrix(scipy.sparse.bsr_matrix(adj_matrix.loc[node_groups['L'], node_groups['R']]))
        for (n1, n2, d) in B_nx.edges(data=True): # convert to unweighted
            d.clear()
        B['nx'] = B_nx
        
        return B, node_groups, fit, comp
            
    def nx_to_edgelist(self, B, n):
        '''Convert an undirected, binary bipartite network to an edge list.
        
        Parameters
        ----------
        B: nx.Graph
            an undirected, unweighted, bipartite network. Group 1 nodes are ordered as described in nx.bipartite.random_graph
        n: int
            the number of nodes in group 1
        
        Returns
        ----------
        B: nx.Graph
            same as input, with disconnected nodes removed from network
        edge_list: list
            each entry is a tuple representing an edge between a nodes in group 1 and a node in group 2
        node_groups: dict
            keys: '1' for group 1, '2' for group 2
            values: lists corresponding to node labels for each group

        '''
        
        if not nx.bipartite.is_bipartite(B) or nx.is_directed(B):
            raise ValueError('Network must be bipartite and undirected')
        
        group1 = list(B.nodes)[:n]
        group2 = list(B.nodes)[n:]

        disconnected = list(nx.isolates(B))
        if len(disconnected) > 0:
            mssg = '{} nodes are disconnected, removing from network'.format(len(disconnected))
            warnings.warn(mssg)

            group1 = sorted(set(group1).difference(disconnected))
            group2 = sorted(set(group2).difference(disconnected))

            B.remove_nodes_from(disconnected)

        node_groups = {'1': group1, '2': group2}
        edge_list = list(B.edges)
        
        return B, edge_list, node_groups
    
    def drop_disconnected_nodes(self, G):
        '''Remove disconnected nodes from network
        
        Parameters
        ----------
        B: nx.Graph 
        
        Returns
        ---------
        B: nx.Graph
            same as input but with nodes that have degree = 0 removed
        
        '''
        G.remove_nodes_from(list(nx.isolates(G)))

    def power_fit(self, G):
        '''Gets the power law fitted parameters for a graph

        Parameters
        ----------
        G: nx.Graph or igraph.Graph
            Don't have an options for weighted bipartite networkx graphs yet.

        Returns
        ----------
        fit: igraph.FittedPowerLaw
            scale-free network parameters for B_ig (p-value from Kolmogrov-Smirnov test)
        '''
        if isinstance(G, nx.Graph):
            G_ = self.ig_from_nx(G)
        else:
            G_ = G.copy()
        fit = igraph.power_law_fit(G_.degree())
        return fit

    def ig_from_nx(self, G):
        '''Conversion from networkx object to igraph object.

        Parameters
        ----------
        G: nx.Graph 
            Don't have an options for weighted bipartite networkx graphs yet.

        Returns
        ----------
        G_ig: igraph.Graph
        '''

        if nx.bipartite.is_bipartite(G):
            if nx.is_weighted(G):
                raise ValueError('No code for converting weighted bipartite networkx graphs to igraph')
            # must be continuously sorted for ig
            G_ = nx.relabel_nodes(G, mapping = {v:k for k,v in dict(enumerate(G.nodes)).items()}, copy = True)
            groups = nx.bipartite.sets(G_)
            mapper = {k: True for k in groups[0]}
            mapper.update({k: False for k in groups[1]})
            mapper = OrderedDict({k: mapper[k] for k in sorted(mapper)})
            G_ig = igraph.Graph.Bipartite(types = list(mapper.values()), edges = sorted(G_.edges), 
                                         directed = nx.is_directed(G_))
            for new, old in dict(enumerate(G.nodes)).items():
                G_ig.vs[new]['id'] = old
        else:
            mode = 'DIRECTED' if nx.is_directed(G) else 'UNDIRECTED'
            G_ig = igraph.Graph.Adjacency((nx.to_numpy_matrix(G) > 0).tolist(), mode = mode)    
        return G_ig
    
    def summarize(self, G):
        '''Generate network level summary statistics
        
        Parameters
        ----------
        G: nx.Graph
        
        Returns
        ----------
        summary: pd.DataFrame
            alpha is fitted the power law exponent
            p-value is from the Kolmogorov-Smignrov test for power-law
            remainder are network level statistics
        
        '''
        summary = pd.DataFrame(index = ['alpha', 'p-val', 'clustering_coefficient', 
                                  'median_degree', 'median_eigenvector_centrality'])
        try:
            fit = self.power_fit(G)
            alpha, p = fit.alpha, fit.p
        except:
            alpha, p = float('nan'), float('nan')
        
        C = nx.average_clustering(G)
        
        try:
            centrality = np.median(list(nx.centrality.eigenvector.eigenvector_centrality(G).values()))
        except:
            centrality = float('nan')

        summary[''] = [alpha,p, C, np.median([i[1] for i in G.degree]), centrality]
        return summary
    
    def subset_nodes(self, G, nodes = None, subset_size = None, drop = True):
        '''Generate a network subset by randomly removing nodes

        Parameters
        ----------
        G: nx.Graph 
        nodes: list
            Each entry is a node in G to be removed, takes precedence of subset_size
        subset_size: float (0,1)
            removes a random subset, leaving a network with size proprtional to subset_size*G
        drop: bool
            whether to drop disconnected nodes

        Returns
        ---------
        _G: nx.Graph
            same as input but with specified nodes removed

        '''
        _G = G.copy()
        if nodes is None:
            if subset_size is None or subset_size >= 1 or subset_size <= 0:
                raise ValueError('Must specify an appropriate subset size when nodes is not specified')
            
            random.seed(self.seed)
            nodes = random.sample(G.nodes, round(len(G.nodes)*(1-subset_size)))
        
        _G.remove_nodes_from(nodes)
        if drop:
            self.drop_disconnected_nodes(_G)
        return _G
    def subset_edges(self, G, edges = None, subset_size = None, drop = True):
        '''Generate a network subset by randomly removing edges

        Parameters
        ----------
        G: nx.Graph 
        edges: list
            Each entry is an edge in G to be removed, takes precedence of subset_size
        subset_size: float (0,1)
            removes a random subset, leaving a network with size proprtional to subset_size*G
        drop: bool
            whether to drop disconnected nodes

        Returns
        ---------
        _G: nx.Graph
            same as input but with specified edges removed

        '''
        _G = G.copy()
        if edges is None:
            if subset_size is None or subset_size >= 1 or subset_size <= 0:
                raise ValueError('Must specify an appropriate subset size when nodes is not specified')
            random.seed(self.seed)
            edges = random.sample(G.edges, round(len(G.edges)*(1-subset_size)))

        _G.remove_edges_from(edges)
        if drop:
            self.drop_disconnected_nodes(_G)
        return _G


# In[ ]:




