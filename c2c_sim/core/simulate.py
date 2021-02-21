#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import warnings
import copy
import pickle

import uuid
import random
from tqdm import tqdm
from tqdm import trange

import networkx as nx

import pandas as pd
import numpy as np
import tensorly

import math
import itertools
from scipy.stats import skewnorm

from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

# script version: 
path_ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__), '../'))
import sys
sys.path.insert(1, path_)

from core.graphs import graph_generator as gg_
from core import utils


# In[2]:


def weight_bias(n, skew):
    '''
    n: int
        Number of entries to bin
    skew: extent to which skew binning
    '''
    unbiased = [1/n]*n
    if skew == 0:
        biased = np.array(unbiased)
    elif skew == 1:
        biased = [0]*len(unbiased)
        biased[-1] = 1
        biased = np.array(biased)
    else:
        X = np.linspace(0, len(unbiased), len(unbiased))
        biased = skewnorm.pdf(X, a = 1, loc = len(unbiased), scale = (1-skew/1)*len(unbiased))
        biased = biased/sum(biased)
    return biased*100

class LR():
    '''object to store metadata and relevant information for the ligan-receptor dimension of tensor
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
    def generate_metadata(self, n_LR_cats = {2: 0}, cat_skew = 0):
        '''Generate metadata groupings for the L-R pairs. Categories are defined as distinct types of 
        metadata associated with the LR pair, e.g. "signaling pathway". Subcategories are
        the associated labels within a category, e.g. "growth" and "inflammation" within the "signaling pathway" category.

        Note: For skew, 0 means evenly distributed, 1 means all LR pairs fall into the first category/subcategory. 

        n_LR_cats: dict
            The length of the dictionary represents the total number of categories associated with the LR
            Each key is an integer representing the number of subcategories for the particular category. 
            Each value is a float [0,1] indicating the skew of distribution of LRs across 
            subcategories within each category. 
        cat_skew: float [0,1]
            Skew of distribution of LRs across categories

        '''
        if len(n_LR_cats) > 1:
            raise ValueError('Currently, only one metadata category can be considered')
        # group each LR into the categories above
        # generate categories
        LR_categories = [str(uuid.uuid4()).split('-')[-1] for i in range(len(n_LR_cats))]
        cat_bias = weight_bias(n = len(LR_categories), skew = cat_skew)
        self.LR_metadata = pd.DataFrame(data = {'LR_id': self.edge_list, 
                        'category': random.choices(population = LR_categories, weights=cat_bias, 
                                                  k=len(self.edge_list))})

        # generate subcategories
        self.LR_metadata['subcategory'] = float('nan')
        i = 0
        for n_subcat, subcat_skew in n_LR_cats.items():
            sub = self.LR_metadata[self.LR_metadata.category == LR_categories[i]]
            subcat_bias = weight_bias(n = n_subcat, skew = subcat_skew)
            self.LR_metadata.loc[sub.index, 'subcategory'] = random.choices(population = [str(uuid.uuid4()).split('-')[-1] for i in range(n_subcat)], 
                           weights=subcat_bias, k=sub.shape[0])
            i += 1

class CCI_MD():
    '''Generate the CCI network for the tensor slice at time point 0'''
    
    def cci_network(self, n_cells, directional = True, autocrine = True):
            '''Initialize the cell-cell interaction network.

            n_cells: int
                the total number of cells to simulate; all cell-cell pairs will have a potential interaction, but only
                those that actually interact will have a score > 0 in the tensor slice
            directional: bool
                whether cell-cell interactions are directional (tuple of cell (A,B) indicates interaction from A-->B) or 
                not
            autocrine: bool
                whether cells can interact with themselves
            '''
            # generate random cell ids
            self._cell_ids = [str(uuid.uuid4()).split('-')[-1] for i in range(n_cells)]
            if directional:
                self.cell_interactions = list(itertools.permutations(self._cell_ids, 2))
            else:
                self.cell_interactions = list(itertools.combinations(self._cell_ids, 2))
            
            if autocrine:
                self.cell_interactions += [(id_, id_) for id_ in self._cell_ids]
    def generate_metadata(self, n_cell_cats = {2: 0}, cat_skew = 0, 
                         remove_homotypic = None):
        '''Generate metadata groupings for the cells (individual). Categories are defined as distinct types of 
        metadata associated with the cell or protein, e.g. "cell type" and "cell cycle phase". Subcategories are
        the associated labels within a category, e.g. "T-cell" and "dendritic cell" within the "cell type" category.
        
        Note: For skew, 0 means evenly distributed, 1 means all cells fall into the first category/subcategory. 
        
        n_cell_cats: dict
            The length of the dictionary represents the total number of categories associated with the cell
            Each key is an integer representing the number of subcategories for the particular category. 
            Each value is a float [0,1] indicating the skew of distribution of cells across 
            subcategories within each category. 
        cat_skew: float [0,1]
            Skew of distribution of cells across categories
        remove_homotypic: int
            whether to remove homotypic ineractions between cells by cell category; how many categories to consider; 
            must be <= the number of categories present in the metadata. Recommended to leave as default and 
            instead toggle consider_homotypic option in Simulate().generate_tensor_md()

        '''
        if len(n_cell_cats) > 1:
            raise ValueError('Currently, only one metadata category can be considered')
        if remove_homotypic > len(n_cell_cats):
            raise ValueError('The value for "remove_homotypic" cannot be larger than the total number of categories associated with the cells')

        # group each cell into the categories above
        # generate categories
        cell_categories = [str(uuid.uuid4()).split('-')[-1] for i in range(len(n_cell_cats))]
        cat_bias = weight_bias(n = len(cell_categories), skew = cat_skew)
        self.cell_ids = pd.DataFrame(data = {'cell_id': self._cell_ids, 
                        'category': random.choices(population = cell_categories, weights=cat_bias, 
                                                  k=len(self._cell_ids))})
        
        # generate subcategories
        self.cell_ids['subcategory'] = float('nan')
        i = 0
        for n_subcat, subcat_skew in n_cell_cats.items():
            sub = self.cell_ids[self.cell_ids.category == cell_categories[i]]
            subcat_bias = weight_bias(n = n_subcat, skew = subcat_skew)
            self.cell_ids.loc[sub.index, 'subcategory'] = random.choices(population = [str(uuid.uuid4()).split('-')[-1] for i in range(n_subcat)], 
                           weights=subcat_bias, k=sub.shape[0])
            i += 1
        self.cell_metadata = self.cell_ids
        del self.cell_ids
        

        if remove_homotypic is not None and remove_homotypic > 0: # remove homotypic interactions of a given category
            print('Remove homotypic cell interactions for {} categories'.format(remove_homotypic))
            i = 0
            to_remove = list()
            while i < remove_homotypic:
                cat = cell_categories[i]
                sub = self.cell_metadata[self.cell_metadata.category == cell_categories[i]]
                cell_ids = sub.cell_id.tolist()

                for ccp in self.cell_interactions:
                    sub_ = sub[(sub.cell_id == ccp[0]) | (sub.cell_id == ccp[1])]
                    if sub_.shape[0] == 2 and sub_.subcategory.unique().shape[0] == 1:
                        to_remove.append(ccp)

                i += 1

            # remove homotypic interactions as identified above for categories 1-i
            self.cell_interactions = list(set(self.cell_interactions).difference(to_remove))

            # filter out any cells that no longer are present 
            cell_ids = list(set(sum(list(zip(*self.cell_interactions)), ())))
            self.cell_metadata = self.cell_metadata[self.cell_metadata.cell_id.isin(cell_ids)]
            self.cell_metadata.reset_index(inplace = True, drop = True)

class sim_tensor():
    '''Class to consolidate and store all label information for the input of tensor-cell2cell'''
    def __init__(self, tensor_cci, conditions, lr_labels, senders, receivers, 
                lr_metadata, cell_metadata, cell_lr_metadata):
        '''
        Note: object stores all input arguments + self.rank (see below)
        
        Paramters
        -------
        tensor_cci: tensorly.tensor
            A tensor formatted in the method needed to run the tensor-cell2cell pipeline
        conditions: list
            label for the conditions, defines first dimension of tensor
        lr_labels: list 
            label for the LR labels, defines the second dimension of tensor
        sender: list
            label for sender cells, defines the third dimension of the tensor
        receiver: list
            label for sender cells, defines the fourth dimension of the tensor
        lr_metadata: pd.DataFrame
            metadata associated with the ligand-receptor pairs
        cell_metadata: pd.DataFrame
            metadata associated with the cells
        cell_lr_metadata: pd.DataFrame
            metadata associated with the expected groups of CC-LR categories that are expected to have patterns of 
            interaction
        
        Returns
        -------
        self.rank: int
            the expected rank after decomposition (= number of patterns in cell_lr_metadata)
        '''
        self.tensor_cci = tensor_cci
        self.conditions = conditions
        self.lr_labels = lr_labels
        self.senders = senders
        self.receivers = receivers
        
        self.lr_metadata = lr_metadata
        self.cell_metadata = cell_metadata
        self.cell_lr_metadata = cell_lr_metadata
        self.rank = self.cell_lr_metadata.shape[0]


# In[ ]:


def fold_change_pattern(initial_value, score_change = 'max'):
    '''The maximum change in the average LR score given the starting value'''
    decrease = False
    if initial_value > 0.5:
        initial_value = 0.5 - (initial_value - 0.5)
        decrease = True
    
    if score_change == 'max':
        change = 1 - initial_value
    elif score_change == 'scale':
        if initial_value >= 0.2:
            change = 2*initial_value
        else:
            change = initial_value + 0.2
        change = change - initial_value
    else:
        raise ValueError('Argument score_change can only be "scale" or "max"')
    
    if decrease:
        change = - change
    
    return change

def linear(x, n_conditions):
    return list(np.linspace(x[1], x[1] + x[0], n_conditions))

def fit_increasing_power(x,adj1,adj2,pw):
    return ((x+adj1) ** pw) * adj2

def power(x, n_conditions):
    change, initial_val = x[0], x[1]
    if change > 0: 
        # https://stackoverflow.com/questions/33186740/fitting-exponential-function-through-two-data-points-with-scipy-curve-fit
        p1 = [1, n_conditions] 
        p2 = [initial_val, initial_val + change]

        pw = 0.2
        A = np.exp(np.log(p2[0]/p2[1])/pw)
        a = (p1[0] - p1[1]*A)/(A-1)
        b = p2[0]/(p1[0]+a)**pw
        xf=np.linspace(1,n_conditions, n_conditions)
        vector = fit_increasing_power(xf, a, b, pw)
    else:
        if initial_val + change == 0:
            vector = np.geomspace(initial_val, initial_val + change + 1e-9, num=n_conditions)
        else:
            vector = np.geomspace(initial_val, initial_val + change, num=n_conditions)
    return vector

def pulse(x, n_conditions):
    change = x[0]
    initial_val = x[1]
    
    vector = [initial_val] * n_conditions # initialize
    
    if n_conditions % 2 == 1:
        mid_point = [math.floor(n_conditions/2)]
    else:
        mid_point = [n_conditions/2 - 1, n_conditions/2]

    periph = None
    if n_conditions >= 5: 
        periph = [min(mid_point)-1, max(mid_point)+1]
    
    for m in mid_point:
        vector[int(m)] = initial_val + change
    if periph is not None:
        for p in periph:
            vector[int(p)] = initial_val + (change*0.5)
    return vector

def oscillate(x, n_conditions):
    osc_period = 3
    if n_conditions > 3:
        iter_vals = list(np.linspace(x[1], x[1] + x[0], osc_period))
        iter_vals += [iter_vals[1]]#iter_vals[1:-1][::-1]

        vector = list()
        for i,j in enumerate(itertools.cycle(iter_vals)):
            vector.append(j)
            if i >= n_conditions - 1:
                break
        return vector
    else:
        return pulse(x, n_conditions)

pattern_mapper = {'linear': linear, 'pulse': pulse, 'oscillate': oscillate, 
                 'power': power}

def generate_pattern(x, n_conditions):
    vector = pattern_mapper[x[0]](x[1:], n_conditions)
    vector[0] = x[2]
    return vector

class Simulate():
    def __init__(self):
        '''Initialize self

        '''
        self.cci = None
        self._convert_bulk = None
    
    def LR_network(self, network_type = None, B = None, subset = False, **params):
        '''
        Simulates a PPI network of *potential* ligand-receptor interactions, or extracts information. \
        from a use provided network.Defines one tensor dimension
        Caveats: for a scale-free network, the number of ligands = the number of receptors \
                 for either network, there may be disconnected edges depending on "p"
        
        Parameters
        ----------
        network_type: str
             "scale-free" to generate a scale-free network or "normal" to generate a network with a normal degree distribution
        B: nx.Graph
            a user provided undirected, unweighted bipartite network. Assumes in B.nodes, ligands are listed \
            before receptors. Takes precedence over network_type.
        subset: bool
            if B is provided and subset is true, this will take a random subset of the network, dropping disconnected nodes \
            (of a specified size, specfied in params)
        **params: dict (keys for each option specified below)
            the required parameters for generating a bipartite, undirected random network either scale-free or not. \
            
            network_type = scale-free: keys - nodes, degrees, alpha, edges (see graphs.graph_generator.bipartite_sf for description) 
            network_type = normal: keys - n_ligands, n_receptors, p analogous to n,m,p in nx.bipartite.gnmk_random_graph
            B != None: keys - n_ligands as described above
            subset = True: keys - 
                n_ligands as described above
                'subset_size' a value between (0,1) indicating the proportional \
                size of the subset (by nodes) compared to the network
                'subset_type' either 'edges' or 'nodes' indicating whether to subset by removing nodes or edges \
                (edges recommended because they maintain the scale-free property)
        
        Returns
        ----------
        self.LR: 
            populates LR object, key outputs outlined here
        self.LR.B: nx.Graph
            undirected bipartite graph with specified degree distribution (power or normal), or user specified B \
            disconnected nodes are removed
        self.LR.edge_list: list
            each entry is a tuple representing a potential interaction between a ligand-receptor pair, ligands on 0 index of each tuple


        '''
        gg = gg_() # return networkx object for graphs
        user = False
        if B is not None: #untested
            user = True
            # properties checked when calling gg.nx_to_edgelist
            if network_type is not None:
                warnings.warn('You have specified a network type and provided a network, B will take priority over network type')
            if 'n_ligands' not in params:
                raise ValueError('For a provided B, you must specify n_ligands in params')
            
            if subset:
                if 'subset_size' not in params or 'subset_type' not in params:
                    raise ValueError('To subset B, you must provide a desired subset_size and subset_type')
                if params['subset_type'] == 'edges':
                    B = gg.subset_edges(B, subset_size = params['subset_size'], drop = True)
                elif params['subset_type'] == 'nodes': 
                    B = gg.subset_nodes(B, subset_size = params['subset_size'], drop = True)
                else:
                    raise ValueError("The subset_type param must be either 'edges' or 'nodes'")
            
            
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

    def emulate_sf_network(self, G):
        '''Emulate a user-provided (recommended scale-free) network for L-R pair tensors dimension
        
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
                                 degrees = np.median([i[1] for i in G.degree]),
                                 alpha = fit.alpha, edges = len(G.edges),
                                 check_properties = True, compare_barabasi = False)
        G2 = G2['nx']
        gg.drop_disconnected_nodes(G2)
        return G2
    
    def generate_tensor_md(self, n_patterns, n_conditions, patterns = ['pulse', 'linear', 'oscillate', 'power'], 
                          consider_homotypic = True, score_change = 'max'):
        '''Generates cell-LR metadata pairs for tensor slices.
        
        Parameters
        ----------
        n_patterns: int (>0)
            the number of CC - LR metadata pairs for which to form distinct interactions 
            the remaining background will default to 0, with noise increasing this value
            the groups with distinct interactions will each have distinct average values spaced b/w (0,1]
        n_conditions: int (>=0)
            the number of conditions across which to generate tensor slices 
        patterns: list
            list of strings, each of which should be included as a potential pattern for a given cell 
            metadata - LR metadata pair. Options: ['pulse', 'linear', 'power', 'oscillate']
        consider_homotypic: bool
            whether to allow homotypic interactions to have patterns. If False, homotypic interactions are guaranteed
            to be a part of the background. 
        score_change: str
            one of ['max', 'scale']. max will make the change in scores as large as possible, scale will scale 
            the change relative to the distance of the comunication from 0.5 at condition '0'
        
        Returns
        ---------
        self.clrm: pd.DataFrame
            a list of metadata CC-LR pairs for which patterns of scores will change across conditions
            alongside the expected average score for each condition
        '''
        #checks------------------------------------------------------------------------------------------------
        allowed_patterns = ['pulse', 'linear', 'oscillate', 'power']
        if patterns is not None:
            if len(set(patterns).difference(allowed_patterns)) > 0:
                raise ValueError('Patterns can only include: ' + ', '.join(allowed_patterns))
        else:
            patterns = allowed_patterns 
        
        if n_conditions == 2 and set(['oscillate', 'pulse']).difference(patterns): 
            warnings.warn('At least 3 conditions are required for oscillations or pulses, only linear patterns will be considered')
            patterns = ['linear']
        if n_conditions <= 1:
            warnings.warn('No conditions specified, only a tensor slice will be generated')
            if n_conditions == 0:
                n_conditions = 1
        
        self.n_conditions = n_conditions
        
       
        #------------------------------------------------------------------------------------------------
        
        n_lr_cat = len(self.LR.LR_metadata.subcategory.unique())
        n_cc_cat = len(self.cci.cell_metadata.subcategory.unique())
    
        # all possible groups that have patterns of expression

        lr_group = list()
        for i in range(1, n_lr_cat + 3):
            lr_group += list(itertools.combinations(self.LR.LR_metadata.subcategory.unique(), i))

        # all possible cell groups that have patterns of expression
        ccat_map = dict(zip(self.cci.cell_metadata.cell_id, self.cci.cell_metadata.subcategory))
        ccati = list()
        for ci in self.cci.cell_interactions:
            ccati.append((ccat_map[ci[0]], ccat_map[ci[1]]))
        ccati = pd.Series(ccati).unique().tolist()
        if not consider_homotypic:
            ccati = [cp for cp in ccati if cp[0] != cp[1]]

        # all possible groups of cell metadata - LR metadata pairs
        clrm = pd.DataFrame(columns = ['cell_subcat', 'LR_subcat'])
        counter = 0
        for i in list(itertools.product(ccati, lr_group)):
            clrm.loc[counter, : ]= [i[0], i[1]]
            counter += 1

        # no all-all combinations
        # clrm.drop(index = [clrm.shape[0] - 1], inplace = True)

        if n_patterns > clrm.shape[0]:
            warnings.warn('More patterns than possible specificed, setting to maximum possible: {}'.format(clrm.shape[0]))
            n_patterns = clrm.shape[0]

        # chose random subset to assign patterns to 
        clrm = clrm.loc[sorted(random.sample(clrm.index.tolist(), k = n_patterns)),]
        clrm.reset_index(inplace = True, drop = True)
        clrm.LR_subcat = clrm.LR_subcat.apply(lambda x: x[0])
        
        self.clrm = clrm
        self.n_patterns = n_patterns
        self.ts_frame = pd.DataFrame(columns = self.cci.cell_interactions, index = self.LR.edge_list)

        # sort metadata categories
        ccat_map = dict(zip(self.cci.cell_metadata.cell_id, self.cci.cell_metadata.subcategory))
        LR_map = dict(zip(self.LR.LR_metadata.LR_id, self.LR.LR_metadata.subcategory))
        self._lrcats = [LR_map[lri] for lri in self.ts_frame.index]
        ccats = [(ccat_map[ci[0]], ccat_map[ci[1]]) for ci in self.ts_frame.columns]
        
        # get tensor slice coordinates for CC-LR pairs with expected patterns
        # get tensor slice coordinates for CC-LR pairs with expected patterns
        def get_coords(x):
            '''Get the coords for each CC-LR metadata pair'''
            col_coord = [j for j,cmd in enumerate(ccats) if cmd == x[0]]
            row_coord = [i for i,lrmd in enumerate(self._lrcats) if lrmd == x[1]]
            coords = list(itertools.product(row_coord, col_coord))
            coords = [tuple(i[0] for i in coords), tuple([i[1] for i in coords])]
            return coords
        self.clrm['ts_coordinates'] = self.clrm[['cell_subcat', 'LR_subcat']].apply(lambda x: get_coords(x), axis = 1).tolist()
        
        # initial value
        self.clrm[['0']] = list(np.arange(1/self.n_patterns, 1+1/self.n_patterns, 1/self.n_patterns))
        
        # patterns over time
        ap = list()
        for i in range(math.ceil(self.n_patterns/len(patterns))):
            random.shuffle(patterns)
            ap += patterns
        self.clrm.insert(3, 'pattern', ap[:self.n_patterns])

        self.clrm = pd.concat([self.clrm,
                  pd.DataFrame(index = self.clrm.index, columns = [str(i) for i in range(1,self.n_conditions)])], axis = 1)
        self.clrm.insert(3, 'change', self.clrm['0'].apply(fold_change_pattern, args = (score_change,)))

        # apply patterns to get averages across conditions
        if self.n_conditions > 1:
            self.clrm[[str(i) for i in range(self.n_conditions)]] = self.clrm[['pattern', 'change', '0']].apply(generate_pattern, args = (self.n_conditions,), axis = 1).tolist()

    def generate_tensor(self, noise, binary = False, bulk = False, noise_max = None):
        '''Generates the tensor.
        
        Parameters
        ----------
        noise: float [0,1]
            extent from which to perturb scores from the expected average value, including background
        binary: bool
            whether to have scores be continuous b/w [0,1] (False) or binary (True). Binary scoring not currently 
            implemented and must be set to False
        bulk: bool
            whether single-cells should be simulated (False) or cells should be grouped into metadata category (True)
        noise_max: float [0,1]
            the maximum average value of the background noise (@ noise = 1). If None, defaults to the minimum
            average communication score value at condition '0'. Will likely fail at values > 0.1. 
        
        Returns
        ---------
        self.ts: dictionary
            keys are labels for each condition (0 through n_conditions-1). Values are tensor slices with
            columns as cell-cell pairs and rows as ligand-receptor pairs
        '''
        if binary:
            binary = False
            warnings.warn('Only continuous scoring is currently implemented')
        if noise > 0 and (noise_max == 0):
            raise ValueError('Noise set to > 0, yet maximum noise is specified at 0. Specify a value larger than 0 or None')
        

        self.ts = dict() # intitialize
        c_labels = [str(i) for i in range(self.n_conditions)]
        
        if bulk: # consolidate cells into the metadata category
            if self._convert_bulk is None: 
                ccat_map = dict(zip(self.cci.cell_metadata.cell_id, self.cci.cell_metadata.subcategory))
                ccat_map = {cpi: (ccat_map[cpi[0]], ccat_map[cpi[1]]) for cpi in self.ts_frame.columns}

                # change matrix to just include cell metadata groupings rather than individual cell IDs
                self.ts_frame = pd.DataFrame(index = self.ts_frame.index, columns = sorted(set(ccat_map.values())))


                # rewrite coordinates according to bulk groupings
                def get_coords(x):
                    '''Get the coords for each CC-LR metadata pair'''
                    col_coord = [j for j,cmd in enumerate(self.ts_frame.columns.tolist()) if cmd == x[0]]
                    row_coord = [i for i,lrmd in enumerate(self._lrcats) if lrmd == x[1]]
                    coords = list(itertools.product(row_coord, col_coord))
                    coords = [tuple(i[0] for i in coords), tuple([i[1] for i in coords])]
                    return coords
                self.clrm['ts_coordinates'] = self.clrm[['cell_subcat', 'LR_subcat']].apply(lambda x: get_coords(x), axis = 1).tolist()
                self._convert_bulk = True
                # _convert_bulk records if .generate_tensor() method with bulk = True being run for first time
                # allows .generate_tensor() method to be run multiple times (e.g., to test different levels of noise)

        if noise == 0:
            for cond in c_labels:
                df = self.ts_frame.copy()
                for idx in self.clrm.index:
                    avg_val = self.clrm.loc[idx, cond]
                    coords = self.clrm.loc[idx, 'ts_coordinates']
                    df.values[coords] = avg_val # non-background 
                df.fillna(0, inplace = True) # background
                self.ts[cond] = df
        else:
            if noise_max is None:
                noise_max = self.clrm['0'].min()
            scale = noise_max/np.array([utils.piecewise_fit(noise_max, *utils.fit_params)])[0]
            for cond in c_labels:
                df = self.ts_frame.copy()
                # background
                vals = utils.get_truncated_normal(n = self.ts_frame.shape[0]*self.ts_frame.shape[1], 
                                                  sd = noise*noise_max, mean = 0)*scale
                vals[vals>1] = 1
                df[:] = vals.reshape(self.ts_frame.shape)
                for idx in self.clrm.index: # non-background
                    vc = self.clrm.loc[idx,c_labels]
                    min_val = vc.min()
                    avg_val = self.clrm.loc[idx, cond]
                    coords = self.clrm.loc[idx, 'ts_coordinates']
                    
                    # boundaries generated like background
                    if min_val == 0 and abs(min_val - avg_val) < 1e-3:
                        noise_max_ = vc[vc>0].min()
                        scale_ = noise_max_/np.array([utils.piecewise_fit(noise_max_, *utils.fit_params)])[0]
                        vals = utils.get_truncated_normal(n = len(coords[0]), 
                                                  sd = noise*noise_max_, mean = 0)*scale_
                        vals[vals>1] = 1
                    else: # sd function of the average value
                        vals = utils.get_truncated_normal(n = len(coords[0]), 
                                                                   sd = noise*avg_val*1.1, mean = avg_val)
                    df.values[coords] = vals
                self.ts[cond] = df

    def reshape(self):
        '''Creates tensor in the format necessary for running with tensor-cell2cell and stores 
        as a sim_tensor object under self.sim_tensor
        
        Formatting:
        - communication matrix: a matrix of sender by receiver cells for a given LR pair
        - 3D tensor: each slice is a communication matrix, the full tensor is all slices for all LR pairs
        - final tensor: each slice is the 3D tensor, the full tensor is all slices for all conditions
        '''
        senders = sorted(set([i[0] for i in self.ts_frame.columns]))
        receivers = sorted(set([i[1] for i in self.ts_frame.columns]))


        # map CC coordinates in CC-LR matrix to .iloc coordinates for fast filling of "communication" matrix
        tcs_coords = dict()
        for coord_i, cell_i in enumerate(senders):
            for coord_j, cell_j in enumerate(receivers):
                tcs_coords[(cell_i,cell_j)] = (coord_i, coord_j)

        coords = [tcs_coords[col] for col in self.ts_frame.columns]
        coords = [tuple([i[0] for i in coords]), tuple([i[1] for i in coords])]

        tcs_ = np.full([len(senders),len(receivers)], np.nan)

        tensor_list = list()
        print('Generate reshaped tensor')
        for cond in tqdm(self.ts):
            tcs_list = list()
            ts = self.ts_frame.copy()
            ts[:] = self.ts[cond].values # assign values for specific condition in a CC-LR matrix
            for idx in range(self.ts_frame.shape[0]):
                tcs = tcs_.copy() # generate a "communication" matrix and fill with each row of the CC-LR matrix
                tcs[coords] = ts.iloc[idx,:].values
                tcs_list.append(tcs)
            tensor_list.append(tensorly.tensor(tcs_list))
        tensor_cci = tensorly.tensor(tensor_list)
        del tensor_list

        # store in a sim_tensor object
        self.sim_tensor = sim_tensor(tensor_cci, conditions = list(self.ts.keys()), 
                                     lr_labels = self.ts_frame.index.tolist(), senders = senders, receivers = receivers, 
                                    lr_metadata = self.LR.LR_metadata, cell_metadata = self.cci.cell_metadata, 
                                     cell_lr_metadata = self.clrm.drop(columns=['ts_coordinates']))
    
    def get_background(self):
        '''Get the coordinates for the tensor background (no patterns)
        
        Returns
        ---------
        background_coords: list
            List of tuples for background coordinates (non-pattern). 
            Formatted to allow easy subsetting of numpy matrix
            
        '''
        
        coords = list()
        for x in self.clrm.ts_coordinates:
            i,j = x[0], x[1]
            coords += [(i[k], j[k]) for k in range(len(i))]
        all_coords = list(itertools.product(range(self.ts_frame.shape[0]), range(self.ts_frame.shape[1])))
        background_coords = list(set(all_coords).difference(coords))
        background_coords = [tuple([i[0] for i in background_coords]), tuple([i[1] for i in background_coords])]
        
        return background_coords
    
    def visualize(self, subset = None, include_background = True, background_coords = None, 
                 file_name = None, colors = None):
        '''Generate lineplot of communication score across conditions for each pattern
        Parameters
        ----------
        subset: float [0,1]
            proportion by which to subset tensor for faster graphing (not subsetted if None or 1)
        include_background: whether to include the background coordinates in the visualization
            whether to have scores be continuous b/w [0,1] (False) or binary (True). Binary scoring not currently 
            implemented and must be set to False
        background_coords: list
            List of tuples for background coordinates (non-pattern). 
            Formatted to allow easy subsetting of numpy matrix
        file_name: str
            File name to save figure. "full/path/to/filename.ext". If None, will not save
        colors: list
            Length = n_patterns. Each element is a color to assign to one of patterns. 
        
        '''
        if include_background and background_coords is None:
            raise ValueError('Must provide background coordinates, see .get_background()')
        res_all = pd.DataFrame(columns = ['score', 'cat', 'condition'])

        for cond in self.ts:
            for cat in self.clrm.index:
                res_ = pd.DataFrame(columns = ['score', 'cat', 'condition'])
                res_['score'] = self.ts[cond].values[self.clrm.loc[cat,'ts_coordinates']]
                res_['cat'] = str(cat+1)
                res_['condition'] = cond
                res_all = pd.concat([res_all, res_], axis = 0)
            if include_background:
                res_ = pd.DataFrame(columns = ['score', 'cat', 'condition'])
                res_['score'] = self.ts['0'].values[background_coords]
                res_['cat'] = 'background'
                res_['condition'] = cond
                res_all = pd.concat([res_all, res_], axis = 0)

        res_all.reset_index(inplace = True, drop = True)
        res_all['score'] = res_all['score'].astype(float)

        if subset is not None:
            res_all = res_all.loc[random.sample(res_all.index.tolist(), k = round(subset*res_all.shape[0])),]
            res_all.reset_index(inplace = True, drop = True)

        cats = res_all.cat.unique().tolist()

        if include_background:
            cats.remove('background')
        cat_ordering = [str(int(j)) for j in sorted([int(i) for i in cats])]
        if include_background:
            cat_ordering = ['background'] + cat_ordering
        res_all["cat"] = pd.Categorical(res_all["cat"], categories=cat_ordering) 

        cond_ordering = [str(int(j)) for j in sorted([int(i) for i in res_all.condition.unique()])]
        res_all["condition"] = pd.Categorical(res_all["condition"], categories=cond_ordering) 

        print('Generate graph')
        if colors is None:
            if len(cats) > 10:
                raise ValueError('You have assigned more than 10 colors and must provide a palette to colors')
            colors = sns.color_palette("tab10")[:len(cat_ordering)]
        else:
            colors = colors

        fig, ax = plt.subplots(figsize = (10,5))
        sns.lineplot(data = res_all, y = 'score', x = 'condition', hue =  'cat', palette = colors,
                     ax = ax)

        handles, labels = ax.get_legend_handles_labels()
        labels = [labels[0]] + ['{:.0f}'.format(float(i)) for i in labels[1:]]
        ax.legend(labels=labels, handles = handles,bbox_to_anchor=(-0.1, 1), title = 'CC-LR subcategory pair')

        if file_name is not None:
            plt.savefig(file_name, bbox_to_inches = 'tight')
        ("")

    def copy(self):
        return copy.deepcopy(self)
    
    def pickle(self, filename):
        '''Store Simulate() object as a pickled file. 
        
        Parameters
        ----------
        filename: str
            full/path/to/filename.pickle
        '''
        
        with open(filename, 'wb') as f:
            pickle.dump(self, f)


# In[ ]:




