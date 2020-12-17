#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import warnings
import copy

import itertools
import math
import pandas as pd
import networkx as nx

import sys
sys.path.insert(1, '../../scripts/')
from simulation.graphs import graph_generator as gg_


# In[2]:


class Simulate():
    def __init__(self):
        '''Initialize self

        '''
    
    def LR_network(self, network_type = None, B = None, **params):
        '''
        Simulates a PPI network of *potential* ligand-receptor interactions. \
        Defines one tensor dimension
        Caveats: for a scale-free network, the number of ligands = the number of receptors \
                 for either network, there may be disconnected edges depending on "p"
        
        Parameters
        ----------
        network_type: str
             "scale-free" to generate a scale-free network or "normal" to generate a network with a normal degree distribution
        B: nx.Graph
            a user provided undirected bipartite network. Assumes in nx.Graph.nodes, ligands are listed \
            before receptors. Takes precedence over network_type.
        **params: dict
            the required parameters for generating a bipartite, undirected normal network either scale-free or not. \
            network_type = scale-free: keys - (see graphs.graph_generator.bipartite_sf for description) - nodes, degrees, alpha, edges
            network_type = normal: keys - n_ligands, n_receptors, p analogous to n,m,p in nx.bipartite.gnmk_random_graph
            B != None: keys - n_ligands as described above
        
        Returns
        ----------
        self.B: nx.Graph
            undirected bipartite graph with specified degree distribution (power or normal), or user specified B \
            disconnected nodes are removed
        self.edge_list: list
            each entry is a tuple representing a potential interaction between a ligand-receptor pair
        self.node_groups: dict
            keys: 'L' for ligands, 'R' for receptors
            values: lists corresponding to protein IDs for each category
        self.fit, self.comp: 
            see graphs.graph_generator.bipartite_sf

        '''
        gg = gg_('nx') # return networkx object for graphs
        self.comp, self.fit = None, None
        if B is not None:
            # properties checked when calling gg.nx_to_edgelist
            if network_type is not None:
                warnings.warn('You have specified a network type and provided a network, B will take priority over network type')
            if 'n_ligands' not in params:
                raise ValueError('For a provided B, you must specify n_ligands in params')
            self.B = B
        elif network_type == 'scale-free': 
            if 'degrees' not in params or 'nodes' not in params:
                raise ValueError('Must specify degrees and nodes in **params')
            if 'alpha' not in params: 
                params['alpha'] = 2 # also default in gg obj, didn't make it a **kwrag
            if 'edges' not in params:
                B, node_groups, self.fit, self.comp = gg.bipartite_sf(nodes = params['nodes'], degrees = params['degrees'], 
                                                            alpha = params['alpha'])
            else:
                B, node_groups, self.fit, self.comp = gg.bipartite_sf(nodes = params['nodes'], degrees = params['degrees'], 
                                                        alpha = params['alpha'], edges = ['edges'])  
            self.B = B['nx']
            params['n_ligands'] = params['nodes'] # same no. of ligands and receptors
        elif network_type == 'normal':
            if sorted(params) != ['n_ligands', 'n_receptors', 'p']:
                raise ValueError('Must specify n_ligands, n_receptors in **params')
            else:
                self.B = nx.bipartite.gnmk_random_graph(params['n_ligands'],params['n_receptors'], params['p'])
        else:
            raise ValueError('Must specify an appropriate network_type or provide a network B')
        
        self.B, self.edge_list, ng = gg.nx_to_edgelist(self.B, params['n_ligands'])
        
        self.node_groups = {'L': ng['1'], 'R': ng['2']}
        
    def tensor_slice(self, n_cells, binary = True):

        '''Simulates a static time point tensor slice
        
        Parameters
        ----------
        n_cells: int
            the total number of cells to simulate 
        binary: bool
            whether L-R scores are binary or continuous b/w [0,1]
        
        Returns
        -------
        self.ts: pd.DataFrame
            matrix with cell network_type pairs as columns, ligand-receptor pairs as rows, scores as entries
        '''
        # TO DO: generate biased distributions of LR scores
        self.ts = pd.DataFrame(columns = list(itertools.combinations(range(n_cells), 2)), index = self.edge_list)
        
        
    def copy(self):
        return copy.deepcopy(self)


# In[3]:


# output = tempfile.mkstemp(suffix = '_bipartite_sf.csv', dir = './')[1]
# __location__ = os.path.realpath(
#     os.path.join(os.getcwd(), os.path.dirname(os.path.basename(output))))


# In[4]:


# init
sim_sf = Simulate() 
sim_norm = Simulate()

# simulate a randomly connected ligand-receptor network
sim_sf.LR_network(network_type = 'scale-free', **{'nodes': 10, 'degrees': 1, 'alpha': 2}) #scale-free
sim_norm.LR_network(network_type = 'normal', **{'n_ligands': 10, 'n_receptors': 10, 'p': 0.5}) # normally distributed


# In[ ]:




