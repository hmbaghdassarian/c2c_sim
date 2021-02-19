#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import pickle
import numpy as np
import multiprocessing

import sys
sys.path.insert(1, '../../../cell2cell/')
import cell2cell as c2c
from cell2cell.tensor.tensor import BaseTensor
from cell2cell.tensor.factorization import _compute_norm_error

sys.path.insert(1, '../../../scripts/')
from simulation.simulate import Simulate, CCI_MD


# specify path to save figures
fig_path = ''#'/data2/hratch/cci_dt/figures/'
data_path = ''#'/data2/hratch/cci_dt/noise_v_error/'
version = 1

# # Error vs Noise

# In[32]:


n_iter = 1000 # of iterations per noise value
n_cores = 20 # number of cores for parallelization

with open(data_path + 'sim_obj_v' + str(version) + '.pickle', 'rb') as f:
    sim = pickle.load(f)

fn = data_path + 'error_vs_noise_' + str(version) + '.tab' # save results here
if not os.path.isfile(fn):
    with open(fn, 'w') as f:
        f.write('Noise' + '\t' + 'Error' + '\n')
else:
    raise ValueError('Unexpected results file already exists')


# In[33]:


def get_error(noise):
    try:
        sim.generate_tensor(noise = noise, binary = False, bulk = True, noise_max = None)
        tensor_ = BaseTensor()
        tensor_.tensor = sim.sim_tensor.tensor_cci

        r, err = _compute_norm_error(rank = sim.sim_tensor.rank, interaction_tensor = tensor_)
        with open(fn, 'a') as f:
            f.write(str(noise) + '\t' + str(err) + '\n')
    except:
        with open(fn, 'a') as f:
            f.write('fail' + '\t' + 'fail' + '\n')


# In[ ]:


for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
    print('------------------------------------------')
    print('Noise: {:.2f}'.format(noise))
    
    pool = multiprocessing.Pool(processes = n_cores)
    try:
        pool.map(get_error, [noise]*n_iter)
        pool.close()
        pool.join()
    except:
        pool.close()
        pool.join()
        raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))

