#!/usr/bin/env python
# coding: utf-8

# In[12]:


import tempfile
import os
import warnings

import networkx as nx
import igraph

import scipy
import pandas as pd


# In[ ]:


# curdir = os.getcwd()
# abspath = os.path.dirname(os.path.abspath(__file__))


# In[2]:


class graph_generator():
    def __init__(self, obj_type):
        '''Initializer for generating various networks
        
        Parameters
        ----------
        obj_type: str
            "nx" to return a networkx object, "ig" to return an igraph object, "all" to return both
        
        '''
        
        if obj_type not in ['nx', 'ig', 'all']:
            raise ValueError('Must specify appropriate network object type')
            
        self.obj_type = obj_type
        
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
    
    def bipartite_sf(self, nodes, degrees, alpha = 2, edges = None, check_properties = True, compare_barabasi = True):
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
    #     beta: float
    #         fitness of node, alpha = 1 / beta + 1, alpha is scale-free exponent
    #     beta = (1/alpha)-1

        output = tempfile.mkstemp(suffix = '_bipartite_sf.csv', dir = './')[1]
#         output = tempfile.mkstemp(suffix = '_bipartite_sf.csv', 
#                                   dir = 'abspath')[1]
        print(output)
        beta = alpha
        cmd = 'Rscript bigraph_r.r ' + '--beta=' + str(beta) + ' --nodes=' + str(nodes) + ' --degrees=' + str(degrees)
        cmd += ' --output=' + str(output) 
        if edges is not None:
            cmd += ' --edges=' + str(edges)

        print('Generate undirected, bipartite, scale-free graph')
        os.system(cmd)

        # format adjacency matrix
#         os.chdir(abspath)
#         adj_matrix = pd.read_csv(os.path.join('../../scripts/simulation', os.path.basename(output)), index_col = 0)
        adj_matrix = pd.read_csv(output, index_col = 0)
#         os.chdir(curdir)
       
        os.remove(output)
        adj_matrix.index = adj_matrix.index -1 # 0 indexing
        adj_matrix.columns = adj_matrix.index

        # extract information
        node_groups = {'L': list(range(nodes)), 'R': list(range(nodes, 2*nodes))}

        # igraph
        A = adj_matrix.replace(0, float('nan'), inplace = False).stack().reset_index()
        B_ig = igraph.Graph.Bipartite(edges = list(zip(A.level_0, A.level_1)), 
                              types = [True]*nodes + [False]*nodes)

        fit = igraph.power_law_fit(B_ig.degree())

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
            Bt = igraph.Graph.Barabasi(n = B_ig.vcount(), m = degrees, directed = False)
            fit_t = igraph.power_law_fit(Bt.degree())

            comp = pd.DataFrame(data = {'Bipartite': [B_ig.vcount(), B_ig.ecount(), fit.alpha, fit.p ],
                         'Barabasi': [Bt.vcount(), Bt.ecount(), fit_t.alpha, fit_t.p ]}, 
                index = ['nodes', 'edges', 'fitted alpha', 'p-val'])
            comp['change (Barabasi/Bipartite)'] = comp.Barabasi/comp.Bipartite
        
        B = dict()
        if self.obj_type == 'nx' or self.obj_type == 'all':
            B_nx = nx.bipartite.from_biadjacency_matrix(scipy.sparse.bsr_matrix(adj_matrix.loc[node_groups['L'], node_groups['R']]))
            B['nx'] = B_nx
        if self.obj_type == 'ig' or self.obj_type == 'all':
            B['ig'] = B_ig
        

        return B, node_groups, fit, comp


# In[5]:



# alpha = 2
# nodes = 100
# degrees = 5

# gg = graph_generator('all')
# B, node_groups, fit, comp = gg.bipartite_graph(alpha = alpha, nodes = nodes, degrees = degrees, 
#                                                  check_properties = True, compare_barabasi = True)

