#!/usr/bin/env python
# coding: utf-8

# In[1]:


from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd
import numpy as np
import itertools
import random

from tqdm import tqdm

import gc

import sys
sys.path.insert(1, '../../scripts/')
from simulation.utils import get_truncated_normal, piecewise_fit, fit_params
from simulation.utils import relat,x
from simulation.simulate import Simulate, CCI_MD

dp = '/data2/hratch/cci_dt/figures/'


# # What does a tensor slice look like across noise?

# In[2]:


# generate simulation object
sim = Simulate() 
sim.LR_network(network_type = 'scale-free', **{'nodes': 1000, 'degrees': 3, 'alpha': 2}) #scale-free
sim.LR.generate_metadata(n_LR_cats = {3: 0}, cat_skew = 0)

cci = CCI_MD()
cci.cci_network(n_cells = 50, directional = False)
cci.generate_metadata(n_cell_cats = {3: 0}, cat_skew = 0, remove_homotypic = 1)

sim.cci = cci 
sim.generate_tensor_md(n_patterns = 4, n_conditions = 0, patterns = ['pulse', 'linear', 'oscillate'])


# In[ ]:


# get background coordinates
coords = list()
for x in sim.clrm.ts_coordinates:
    i,j = x[0], x[1]
    coords += [(i[k], j[k]) for k in range(len(i))]
all_coords = list(itertools.product(range(sim.ts_frame.shape[0]), range(sim.ts_frame.shape[1])))
background_coords = list(set(all_coords).difference(coords))
background_coords = [tuple([i[0] for i in background_coords]), tuple([i[1] for i in background_coords])]


# In[3]:


print('Generate melted df for one tensor slice')
res_all = pd.DataFrame(columns = ['score', 'noise', 'cat'])
for noise in tqdm(np.arange(0,1.01, 0.1)):
    viz = sim.copy()
    viz.generate_tensor(noise = noise, binary = False)
    for cat in viz.clrm.index:
        res_ = pd.DataFrame(columns = ['score', 'cat', 'noise'])
        res_['score'] = viz.ts['0'].values[viz.clrm.loc[cat,'ts_coordinates']]
        res_['cat'] = str(cat)
        res_['noise'] = noise
        res_all = pd.concat([res_all, res_], axis = 0)
    
    res_ = pd.DataFrame(columns = ['score', 'cat', 'noise'])
    res_['score'] = viz.ts['0'].values[background_coords]
    res_['cat'] = 'background'
    res_['noise'] = noise
    res_all = pd.concat([res_all, res_], axis = 0)
    
res_all.reset_index(inplace = True, drop = True)
res_all['score'] = res_all['score'].astype(float)
list_ordering = ['background','0', '1', '2', '3']
res_all["cat"] = pd.Categorical(res_all["cat"], categories=list_ordering) 


# In[9]:


# res_sub = res_all.loc[random.sample(res_all.index.tolist(), k = round(res_all.shape[0]*0.2))]
# del res_all


# In[ ]:


print('Graph one tensor slice')
fig,ax = plt.subplots(figsize = (10,10))
sns.boxplot(data = res_all, x = 'cat', y = 'score', hue = 'noise')
ax.legend(labels=['{:.1f}'.format(n) for n in sorted(res_all.noise.unique())], 
         bbox_to_anchor=(-0.1, 1))

plt.savefig(dp + 'ts0_comparison.pdf', bbox_to_inches = 'tight')
plt.savefig(dp + 'ts0_comparison.png', bbox_to_inches = 'tight')
print('Complete graph one tensor slice')
("")


# In[1]:


gc.collect()


# # What does a CC-LR metadata pair look like across conditions?

# In[2]:


# generate simulation object
sim = Simulate() 
sim.LR_network(network_type = 'scale-free', **{'nodes': 1000, 'degrees': 3, 'alpha': 2}) #scale-free
sim.LR.generate_metadata(n_LR_cats = {3: 0}, cat_skew = 0)

cci = CCI_MD()
cci.cci_network(n_cells = 50, directional = False)
cci.generate_metadata(n_cell_cats = {3: 0}, cat_skew = 0, remove_homotypic = 1)

sim.cci = cci 
sim.generate_tensor_md(n_patterns = 3, n_conditions = 5, patterns = ['pulse', 'linear', 'oscillate'])



# In[14]:


# get background coordinates
coords = list()
for x in sim.clrm.ts_coordinates:
    i,j = x[0], x[1]
    coords += [(i[k], j[k]) for k in range(len(i))]
all_coords = list(itertools.product(range(sim.ts_frame.shape[0]), range(sim.ts_frame.shape[1])))
background_coords = list(set(all_coords).difference(coords))
background_coords = [tuple([i[0] for i in background_coords]), tuple([i[1] for i in background_coords])]


# In[17]:


print('Generate melted df for multiple conditions')
res_all = pd.DataFrame(columns = ['score', 'noise', 'cat', 'condition'])

for noise in tqdm(np.arange(0,1.01, 0.25)):
    print('---------------------------------------')
    viz = sim.copy()
    viz.generate_tensor(noise = noise, binary = False)
    for cond in viz.ts:
        for cat in viz.clrm.index:
            res_ = pd.DataFrame(columns = ['score', 'cat', 'noise', 'condition'])
            res_['score'] = viz.ts[cond].values[viz.clrm.loc[cat,'ts_coordinates']]
            res_['cat'] = str(cat)
            res_['noise'] = noise
            res_['condition'] = cond
            res_all = pd.concat([res_all, res_], axis = 0)

        res_ = pd.DataFrame(columns = ['score', 'cat', 'noise', 'condition'])
        res_['score'] = viz.ts['0'].values[background_coords]
        res_['cat'] = 'background'
        res_['noise'] = noise
        res_['condition'] = cond
        res_all = pd.concat([res_all, res_], axis = 0)
    
res_all.reset_index(inplace = True, drop = True)
res_all['score'] = res_all['score'].astype(float)

cat_ordering = ['background','0', '1', '2']
res_all["cat"] = pd.Categorical(res_all["cat"], categories=cat_ordering) 

cond_ordering = sorted(res_all.condition.unique())
res_all["condition"] = pd.Categorical(res_all["condition"], categories=cond_ordering) 

noise_ordering = sorted(res_all.noise.unique())
res_all["noise"] = pd.Categorical(res_all["noise"], categories=noise_ordering) 


# In[63]:


# res_sub = res_all.loc[random.sample(res_all.index.tolist(), k = round(res_all.shape[0]*0.005))]


# In[ ]:


print('Generate graph multiple conditions')
cats = sorted(res_al.cat.unique().tolist())
cats.remove('background')
labels_ = ['{:.2f}'.format(n) for n in sorted(res_al.noise.unique())]
colors = sns.color_palette("tab10")[:len(labels_)]

fig, ax = plt.subplots(nrows = len(cats), figsize = (15,12))
counter = 0
for cat in tqdm(cats):
    viz = res_al[res_al.cat == cat]
    sns.lineplot(data = viz, y = 'score', x = 'condition', hue =  'noise', palette = colors,
                 ax = ax[counter])
    if counter == 0:
        ax[counter].legend(labels=labels_, bbox_to_anchor=(-0.1, 1))
    else:
        ax[counter].get_legend().remove()

    ax[counter].set_title('Condition: ' + cat)
    counter += 1

plt.savefig(dp + 'pattern_delta_condition.pdf', bbox_to_inches = 'tight')
plt.savefig(dp + 'pattern_delta_condition.pdf', bbox_to_inches = 'tight')
print('Complete graph multiple conditions')
("")


# In[ ]:


del res_all
gc.collect()

