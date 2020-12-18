#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import warnings
import copy

import itertools
import math
import pandas as pd
import numpy as np
import networkx as nx

import sys
sys.path.insert(1, '../../scripts/')
from simulation.graphs import graph_generator as gg_


# In[2]:


class LR():
    '''object to store metadata and relevant information for the ligan-recept dimension of tensor
    for internal use
    
    '''
    def __init__(self, B, ligands, receptors, edge_list, network_type = None, alpha = None, 
                 fit = None, comp = None, p = None):
        '''Initialize
        
        Parameters
        ----------

        B: nx.Graph
            a undirected bipartite network representing PPI between ligands and receptors (direction would always be L-->R)
        ligands: list 
            ligand IDs for each protein
        receptors: list
            receptor IDs for each protein
        self.edge_list: list
            each entry is a tuple representing a potential interaction between a ligand-receptor pair, ligands on 0 index of each tuple
        network_type: str
            "scale-free" indicates scale-free network, "normal" indicates a normal degree distribution
        alpha: float
            scale-free exponent for network degree distribution (recommended 2<alpha<3)
        p: float
            probability of adding an edge when using network_type option = 'normal'
        fit: igraph.FittedPowerLaw
            scale-free network parameters for B_ig (p-value from Kolmogrov-Smirnov test)
        comp: pd.DataFrame or None
            summary of differences in network properties between  bipartite network and similar Barabasi network
        '''
        if network_type is None:
            self.network_type = 'user-speficied'
        else:
            self.network_type = network_type
        self.B = B
        self.ligands = ligands
        self.receptors = receptors
        self.edge_list = edge_list
        self.alpha = alpha
        self.fit = fit
        self.comp = comp
        self.p = p


# In[3]:


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
        self.LR: 
            populates LR object, key parameters outlined here
        self.B: nx.Graph
            undirected bipartite graph with specified degree distribution (power or normal), or user specified B \
            disconnected nodes are removed
        self.edge_list: list
            each entry is a tuple representing a potential interaction between a ligand-receptor pair, ligands on 0 index of each tuple


        '''
        gg = gg_() # return networkx object for graphs
        user = False
        if B is not None:
            user = True
            # properties checked when calling gg.nx_to_edgelist
            if network_type is not None:
                warnings.warn('You have specified a network type and provided a network, B will take priority over network type')
            if 'n_ligands' not in params:
                raise ValueError('For a provided B, you must specify n_ligands in params')
        elif network_type == 'scale-free': 
            if 'degrees' not in params or 'nodes' not in params:
                raise ValueError('Must specify degrees and nodes in **params')
            if 'alpha' not in params: 
                params['alpha'] = 2 # also default in gg obj, didn't make it a **kwrag
            if 'edges' not in params:
                B, node_groups, fit, comp = gg.bipartite_sf(nodes = params['nodes'], degrees = params['degrees'], 
                                                            alpha = params['alpha'])
            else:
                B, node_groups, fit, comp = gg.bipartite_sf(nodes = params['nodes'], degrees = params['degrees'], 
                                                        alpha = params['alpha'], edges = ['edges'])  
            B = B['nx']
            params['n_ligands'] = params['nodes'] # same no. of ligands and receptors
        elif network_type == 'normal':
            if sorted(params) != ['n_ligands', 'n_receptors', 'p']:
                raise ValueError('Must specify n_ligands, n_receptors in **params')
            else:
                B = nx.bipartite.random_graph(params['n_ligands'],params['n_receptors'], params['p'])
        else:
            raise ValueError('Must specify an appropriate network_type or provide a network B')
        
        B, edge_list, ng = gg.nx_to_edgelist(B, params['n_ligands']) # formatting/extract info
        
        # store PPI information in LR()
        if user:
            self.LR = LR(B, ng['1'], ng['2'], edge_list)
        elif network_type == 'scale-free':
            self.LR = LR(B, ng['1'], ng['2'], edge_list, network_type = network_type, 
                         alpha = params['alpha'], fit = fit, comp = comp)
        elif network_type == 'normal':
            self.LR = LR(B, ng['1'], ng['2'], edge_list, network_type = network_type, p = params['p'])

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
        
        # TO DO: 
        # first separate into groupings/metadata covariates?
        # generate biased distributions of LR scores
        self.ts = pd.DataFrame(columns = list(itertools.combinations(range(n_cells), 2)), index = self.LR.edge_list)
        
        
    def emulate_sf_network(self, G):
        '''Emulate a user-provided network
        
        Parameters
        ----------
        G: nx.Graph
            user-provided network (recommended scale-free)
        
        Returns
        ---------
        G2: nx.Graph
            random bipartite scale-free network built using G's properties
        
        '''
        gg = gg_()
        fit = gg.power_fit(G)
        if fit.p < 0.05:
            warnings.warn('Input network is not scale-free')
        print('----Simulated network------')
        
        G2, node_groups, fit2, comp = gg.bipartite_sf(nodes = round(len(G.nodes)), # should be 1/2 the nodes, but many are disconnected 
                                 degrees = np.mean([i[1] for i in G.degree]),
                                 alpha = fit.alpha, edges = len(G.edges),
                                 check_properties = True, compare_barabasi = False)
        G2 = G2['nx']
        gg.drop_disconnected_nodes(G2)
        return G2
    def copy(self):
        return copy.deepcopy(self)


# In[4]:


# output = tempfile.mkstemp(suffix = '_bipartite_sf.csv', dir = './')[1]
# __location__ = os.path.realpath(
#     os.path.join(os.getcwd(), os.path.dirname(os.path.basename(output))))


# In[5]:


# init
sim_sf = Simulate() 
sim_norm = Simulate()

# simulate a randomly connected ligand-receptor network
sim_sf.LR_network(network_type = 'scale-free', **{'nodes': 10, 'degrees': 1, 'alpha': 2}) #scale-free
sim_norm.LR_network(network_type = 'normal', **{'n_ligands': 10, 'n_receptors': 5, 'p': 0.5}) # normally distributed

