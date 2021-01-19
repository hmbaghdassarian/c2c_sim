#!/usr/bin/env python
# coding: utf-8

# # **Analysis**

# run all code manually: click runtime>run all
# 
# 
# 
# 
# 

# In[3]:


# !pip install velocyto
# !pip install -U loompy
# # !pip install sparse


# In[8]:


# import loompy 
import velocyto as vcy
import numpy as np
import pandas as pd
import scipy as sc
# from google.colab import files
import h5py
import os
# import matplotlib.pyplot as plt


# In[2]:


local_data_path = '/data2/hratch/immune_CCI/'


# In[5]:


# take the output of velocyto run10x (5k_pbmc_v3_count.loom) and proceed with velocyto.py analyses 
vlm = vcy.VelocytoLoom(local_data_path + 'NK_files/velocyto/5k_pbmc_v3_count_celltypes.loom')


# In[6]:


vlm.normalize("S", size=True, log=True)
vlm.S_norm 


# In[7]:


vlm.plot_fractions()


# In[8]:


vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))
vlm.set_clusters(vlm.ca["Cell_Type"])
vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)


# In[9]:


vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)


# In[10]:


vlm._normalize_S(relative_size=vlm.S.sum(0),
             target_size=vlm.S.sum(0).mean())
vlm._normalize_U(relative_size=vlm.U.sum(0),
             target_size=vlm.U.sum(0).mean())
vlm.perform_PCA()
vlm.knn_imputation(n_pca_dims=20, k=500, balanced=True, b_sight=3000, b_maxl=1500, n_jobs=16)


# In[7]:


vlm.fit_gammas()


# In[ ]:


vlm.plot_phase_portraits(['HES4'])


# In[ ]:


vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
# vlm.extrapolate_cell_at_t(delta_t=100)


# In[ ]:


# pd.DataFrame(vlm.Sx_sz).equals(pd.DataFrame(vlm.Sx))
# since this returns True, assume vlm.Sx_sz is the expression matrix at t = 0, should double check this is correct
# pd.DataFrame(vlm.Sx_sz).equals(pd.DataFrame(vlm.Sx_sz_t)) # this should not return True (should be changing with dt)


# In[6]:


if not os.path.isdir(local_data_path + 'NK_files/velocyto_analyses/'):
    os.mkdir(local_data_path + 'NK_files/velocyto_analyses/')


# In[ ]:


# assume vlm.Sx_sz is at t = 0, create h5 object
expression_dt = h5py.File(local_data_path + 'NK_files/velocyto_analyses/expression_dt.h5', 'w')
expression_dt.create_dataset(str(0), data = vlm.Sx_sz) # expression_over_time = {0: vlm.Sx_sz}

min_dt, max_dt, interval = 1,100,1 
counter = 1
for dt in np.arange(min_dt, max_dt, interval):
    print(counter)
    counter += 1
    vlm.extrapolate_cell_at_t(delta_t=dt)
    expression_dt.create_dataset(str(dt), data=vlm.Sx_sz_t)
  # expression_over_time[dt] = sc.sparse.csc_matrix(vlm.Sx_sz_t, dtype = np.uint16) 

expression_dt.close()


# In[ ]:


# write file with cell names
with open(local_data_path + 'NK_files/velocyto_analyses/column_names.txt', 'w') as f:
    for name in vlm.ca['CellID']:
        f.write(name + '\n')


# repeat for rownames
with open(local_data_path + 'NK_files/velocyto_analyses/row_names.txt', 'w') as f:
    for name in vlm.ra['Gene']:
        f.write(name + '\n')

