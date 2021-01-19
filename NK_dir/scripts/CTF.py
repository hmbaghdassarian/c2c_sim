#!/usr/bin/env python
# coding: utf-8

# In[92]:


import os
import warnings
import qiime2 as q2

from qiime2.plugins.emperor.actions import (plot, biplot)
from qiime2.plugins.diversity.actions import (beta_phylogenetic, pcoa)
from deicode import rpca #from qiime2.plugins.deicode.actions import rpca

import pandas as pd
from pandas import DataFrame
import random
import numpy as np
import h5py
import biom

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib

from gemelli.ctf import ctf #from gemelli.actions import ctf
from gemelli.ctf import ctf_helper
from gemelli.factorization import TensorFactorization
from gemelli.preprocessing import build, rclr
from gemelli._ctf_defaults import (DEFAULT_COMP, DEFAULT_MSC,
                                   DEFAULT_MFC, DEFAULT_MAXITER,
                                   DEFAULT_FMETA as DEFFM)
from qiime2.plugins.longitudinal.actions import (volatility, linear_mixed_effects)

# import sys
# sys.path.insert(1, '../scripts/')
# from load_environmental_variables import *
local_data_path = '/data2/hratch/immune_CCI_pseudotime/'


# In[465]:


# load files
print('load files')
# read in t0 cci timepoint
colnames = open(local_data_path + 'interim/velocyto_analyses/column_names.txt').read().splitlines() # cell barcodes in order of CCI distance matrices
# cell_ids = pd.read_csv(local_data_path + 'processed/5k_pbmc_celltypes_velocytoformatted.csv', index_col = 0) # map cell barcodes to cell type
# cell_id_map = dict(zip(cell_ids.SampleID, cell_ids.Cell_Type)) # 

cci_dt = h5py.File(local_data_path + 'interim/CCI_dt.h5', "r")
# cci_t0 = pd.DataFrame(cci_dt[sorted(cci_dt.keys())[0]], columns = colnames, index = colnames)
cci_t0 = biom.Table(np.array(cci_dt[sorted(cci_dt.keys())[0]]), 
                    sample_ids = colnames, observation_ids = colnames)

# # FIX THIS
# # duplicate time poitns not causing error, drop for now
# metadata.drop_duplicates(subset = ['velocity_pseudotime'], inplace = True)
# cci_t0 = cci_t0.loc[metadata.index, metadata.index]
# cci_t0 = biom.Table(cci_t0.values, sample_ids = metadata.index, observation_ids = metadata.index)

# load me`tadata (including pseudotime) from velocyto analysis
metadata = pd.read_csv(local_data_path + 'interim/velocyto_analyses/velocyto_attributes.csv', index_col = 0)
metadata['cell_ids'] = metadata.index

n_bins = 100
metadata['time'] = pd.qcut(x = metadata.velocity_pseudotime, q = n_bins, 
                              labels = list(range(n_bins)), duplicates = 'drop')


# In[ ]:


print('get ctf')
ctf_results = ctf(table = cci_t0, sample_metadata = metadata, 
                                           individual_id_column = 'cell_ids', 
                                           state_column = 'time', 
                 feature_metadata = None)
subject_biplot = ctf_results[0]
state_biplot = ctf_results[1]
distance_matrix = ctf_results[2]
state_subject_ordination = ctf_results[3]
state_feature_ordination = ctf_results[4]

state_feature_ordination.index = state_feature_ordination.index.astype(int)
print('save ctf')

with open(local_data_path + 'interim/ctf_results.pickle', 'wb') as handle:
    pickle.dump(ctf_results)
    
print('complete')


