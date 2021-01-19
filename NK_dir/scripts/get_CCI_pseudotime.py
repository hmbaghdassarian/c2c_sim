#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.insert(1, '../scripts/') # comment out in python script
from load_environmental_variables import *


# In[41]:


import cell2cell as c2c
import numpy as np
import pandas as pd

from tqdm.auto import tqdm

import glob
import os
from multiprocessing import Pool


# In[7]:


print('Load params')
files = {'ppi': local_data_path + 'raw/Human-2020-Cabello-Aguilar-LR-pairs.csv', 
        'output_folder': local_data_path + 'interim/get_CCI_psuedotime/'}

if not os.path.isdir(files['output_folder']):
    os.mkdir(files['output_folder'])

rnaseq_setup, ppi_setup, cutoff_setup, analysis_setup = dict(), dict(), dict(), dict()

rnaseq_setup['gene_col'] = None
rnaseq_setup['drop_nangenes'] = True
rnaseq_setup['log_transform'] = False

ppi_setup['protein_cols'] = ['ligand', 'receptor']

cutoff_setup['type'] = 'constant_value'
cutoff_setup['parameter'] = 0.1

analysis_setup['communication_score'] = 'expression_thresholding'
analysis_setup['cci_score'] = 'bray_curtis'
analysis_setup['cci_type'] = 'undirected'


# In[9]:


cell_types = c2c.io.load_table(local_data_path + 'processed/5k_pbmc_celltypes_velocytoformatted.csv')
celltype_mapper = cell_types[['SampleID', 'Cell_Type']].set_index('SampleID').to_dict()['Cell_Type']


# In[33]:


filenames = glob.glob(local_data_path + 'interim/velocyto_analyses/projected_gene_expression_csvs/T-*.csv')


def get_CCI(file, counter):
    basename = os.path.basename(file).split('.csv')[0]
    print('Filename:' + basename)
    print('File number {} of {}'.format(counter, len(filenames)))
    try:
        print('Format rnaseq data')
        rnaseq_data = c2c.io.load_rnaseq(rnaseq_file=file,
                                                 gene_column=rnaseq_setup['gene_col'],
                                                 drop_nangenes=rnaseq_setup['drop_nangenes'],
                                                 log_transformation=rnaseq_setup['log_transform'],
                                                 format='auto',
                                                 **{'index_col' : 0})
        print('Format ppi data')
        ppi_data = c2c.io.load_ppi(ppi_file=files['ppi'],
                                   interaction_columns=ppi_setup['protein_cols'],
                                   rnaseq_genes=list(rnaseq_data.index),
                                   format='auto')
        
        print('Set up parameters')
        if analysis_setup['cci_type'] == 'undirected':
            bi_ppi_data = c2c.preprocessing.bidirectional_ppi_for_cci(ppi_data=ppi_data, verbose=False)
            ref_ppi = ppi_data
        else:
            bi_ppi_data = ppi_data.copy()
            ref_ppi = None
            
        print('Setup interaction space')
        interaction_space = c2c.core.InteractionSpace(rnaseq_data=rnaseq_data,
                                                              ppi_data=bi_ppi_data,
                                                              gene_cutoffs=cutoff_setup,
                                                              communication_score=analysis_setup['communication_score'],
                                                              cci_score=analysis_setup['cci_score'],
                                                              cci_type=analysis_setup['cci_type'],
                                                              verbose=False)
        # compute interactions
        print('Compute interactions')
        interaction_space.compute_pairwise_cci_scores(use_ppi_score=False, verbose=False)

        # # untested - don't need communication for now
        # # compute communication
        # interaction_space.compute_pairwise_communication_scores(ref_ppi_data=ref_ppi, cci_type='directed', verbose=False)

        print('Save CCI dataframe')
        interaction_space.interaction_elements['cci_matrix'].to_csv(files['output_folder'] + basename + '_CCI.csv')
        # interaction_space.interaction_elements['communication_matrix'].to_csv(files['output_folder'] + basename + '_CCC.csv')
    except:
        print('CCI failed on ' + basename)
    print('-----------------------------------------------------------------------------')


# In[ ]:


print('Begin parallelization')
pool = Pool(processes = 20)
pool.starmap(get_CCI, zip(filenames, list(range(1, len(filenames) + 1))))
pool.close()


# In[60]:


##run the following command in terminal
# cmd = 'python ' + root_path + 'scripts/get_CCI_psuedotime.py > ' + local_data_path 
# cmd += 'interim/get_CCI_pseudotime_terminal_output.txt'
# print(cmd)

