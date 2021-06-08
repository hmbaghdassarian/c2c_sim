#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import random
import numpy as np
import scanpy as sc
# import scipy
from tqdm import tqdm

data_path = '/data2/hratch/immune_CCI/covid/covid_atlas/'
load_h5 = True


# In[3]:


pbmc_covid = sc.read_mtx(data_path + 'raw/GSE158055_covid19_counts.mtx.gz') # raw counts
if load_h5: 
    pch5 = sc.read_h5ad(data_path + 'raw/COVID19_ALL.h5ad') # load dataset
# pbmc_covid = sc.read_10x_mtx(data_path + 'raw/counts/')

print('Finished loading covid datasets')


# In[ ]:


# exclude samples with fewer than 2000 cells
md_cell = pd.read_csv(data_path + 'raw/GSE158055_cell_annotation.csv.gz')
n_samples = md_cell.sampleID.value_counts()
samples_to_keep = n_samples[n_samples > 2000].index.tolist()



# In[ ]:


md = pd.read_excel(data_path + 'raw/GSE158055_sample_metadata.xlsx', sheet_name = 0, skiprows=20)
md = md.iloc[range(304 - 20), range(25)]
md = md[md['Sample name'].isin(samples_to_keep)]


contexts = md['characteristics: CoVID-19 severity'].unique()
n_contexts = contexts.shape[0]

context_counts = md['characteristics: CoVID-19 severity'].value_counts() 
min_context_type = context_counts[context_counts == context_counts.min()].index.tolist()[0]
min_context_count = len(md[md['characteristics: CoVID-19 severity'] == min_context_type]['Patients'].unique())

max_samples = min_context_count*n_contexts

context_map = {context: md[md['characteristics: CoVID-19 severity'] ==                            context][['Sample name', 'Patients']].reset_index(drop = True) for context in contexts}


# In[ ]:


# randomly select samples subsetted from the entire dataset
# make sure to choose an even number of each context
# make sure not to repeat patients within a context

n_iter = 1 # number of times to run subsetting
seed = 0

Samples = pd.DataFrame(columns = ['iteration', 'n_samples', 'sample_names'])
idx = 0

sample_iters = [3, 6, 12, 24, 36, 48, 60, 75]#list(range(n_contexts, max_samples + 1, n_contexts))

for iteration in range(n_iter):
    for n_samples in sample_iters:
        cmap_temp = context_map.copy()
        n_sample_per_context = int(n_samples/n_contexts)
        samples = list()
        for context in contexts:
            df = context_map[context].sample(frac=1, random_state = seed).drop_duplicates(subset = 'Patients') # shuffle rows to randomly drop duplicates
            samples += df.sample(n = n_sample_per_context, random_state = seed)['Sample name'].tolist()
            seed += 1
        Samples.loc[idx,:] = [iteration, n_samples, samples]
        idx += 1

Samples['sample_names'] = Samples.sample_names.apply(lambda x: '; '.join(x))            
Samples.to_csv(data_path + 'interim/timing_inputs/samples_for_timing.csv')            


# In[ ]:


#cell_ids = pd.read_csv(data_path + 'raw/GSE158055_covid19_barcodes.tsv.gz', header = None)
gene_ids = pd.read_csv(data_path + 'raw/GSE158055_covid19_features.tsv.gz', header = None)
md_cell.set_index('cellName', drop = True, inplace = True)

pbmc_covid = pbmc_covid.transpose() 


# In[ ]:


if load_h5:
    pbmc_covid.obs = pch5.obs
    pbmc_covid.var = pch5.var
else:
    pbmc_covid.obs = md_cell
    pbmc_covid.var = gene_ids.set_index(0, drop = True)


# In[ ]:


# split by sample id 

def flatten_list(t):
    return [item for sublist in t for item in sublist]

def create_raw_counts(sample_id):
    df = pbmc_covid[pbmc_covid.obs.sampleID == sample_id]
    sc.pp.filter_cells(df, min_genes=50) 
#     sc.pp.filter_genes(df, min_cells = 3) # avoid filtering genes, will need intersection of remaining genes, which filters to many out when subsequently filtering for LR pairs
    return df

sample_ids = list(set(flatten_list([sn.split('; ') for sn in Samples.sample_names.tolist()])))
sample_counts = {sample_id: create_raw_counts(sample_id) for sample_id in sample_ids}


min_cells = min([df.n_obs for sample_id, df in sample_counts.items()])


# In[ ]:


# seed = 24
seed += 1

# subset to min cells_to_keep and write to csv
cells_to_keep = list()
for sample_id in tqdm(sample_counts):
    random.seed(seed)
    df = sample_counts[sample_id]
    df = df[df.obs.index.isin(random.sample(df.obs.index.tolist(), min_cells_to_keep))] # subset
    df.to_df().to_csv(data_path + 'interim/umi_for_timing/' + sample_id + '.csv') # write
    cells_to_keep += df.obs.index.tolist()
    seed += 1

pbmc_covid.obs[pbmc_covid.obs.index.isin(cells_to_keep)].to_csv(data_path + 'interim/timing_inputs/metadata_for_timing.csv')


# In[ ]:


print('Complete')

