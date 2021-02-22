#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pickle
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt

import sys
sys.path.insert(1, '../cell2cell/')
import cell2cell as c2c
from cell2cell.tensor.tensor import BaseTensor
from cell2cell.tensor.factorization import _compute_norm_error, _compute_tensor_factorization

sys.path.insert(1, '../c2c_sim/')
from core.simulate import Simulate, CCI_MD


# specify path to save figures
fig_path = ''#'/data2/hratch/cci_dt/figures/'
data_path = ''#'/data2/hratch/cci_dt/noise_v_error/'
version = 1

n_iter = 100 # of iterations per noise value
n_cores = 20 # number of cores for parallelization


# In[5]:


with open(data_path + 'sim_obj_v' + str(version) + '.pickle', 'rb') as f:
    sim = pickle.load(f)

fn_1 = data_path + 'error_vs_noise_' + str(version) + '_noise0_05.tab' # save results here
fn_2 = data_path + 'error_vs_noise_' + str(version) + '_noise0_25.tab' # save results here

for fn in [fn_1, fn_2]:
    if not os.path.isfile(fn):
        with open(fn, 'w') as f:
            f.write('Noise' + '\t' + 'Error' + '\n')
    else:
        raise ValueError('Unexpected results file already exists')


# In[33]:


def get_error(noise, noise_max, fn):
    try:
        with open(data_path + 'sim_obj_v' + str(version) + '.pickle', 'rb') as f:
            sim = pickle.load(f)
        
        sim.generate_tensor(noise = noise, binary = False, bulk = True, noise_max = noise_max)
        tensor_ = BaseTensor()
        tensor_.tensor = sim.sim_tensor.tensor_cci

        tl_object_ = _compute_tensor_factorization(tensor=tensor_.tensor,
                                           rank=sim.sim_tensor.rank,
                                           tf_type='non_negative_cp',
                                           init='svd',
                                           random_state=None,
                                           mask=tensor_.mask,
                                           verbose=False)
        err = _compute_norm_error(tensor = tensor_.tensor, tl_object = tl_object_)
        
        with open(fn, 'a') as f:
            f.write(str(noise) + '\t' + str(err) + '\n')
    except:
        with open(fn, 'a') as f:
            f.write('fail' + '\t' + 'fail' + '\n')


# In[ ]:


print('Cap background noise at 0.05')
for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
    print('------------------------------------------')
    print('Noise: {:.2f}'.format(noise))
    
    pool = multiprocessing.Pool(processes = n_cores)
    try:
        pool.starmap(get_error, zip([noise]*n_iter, [0.05]*n_iter, [fn_1]*n_iter))
        pool.close()
        pool.join()
    except:
        pool.close()
        pool.join()
        raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))


# In[ ]:


print('Cap background noise at 0.25')
for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
    print('------------------------------------------')
    print('Noise: {:.2f}'.format(noise))
    
    pool = multiprocessing.Pool(processes = n_cores)
    try:
        pool.starmap(get_error, zip([noise]*n_iter, [None]*n_iter, [fn_2]*n_iter))
        pool.close()
        pool.join()
    except:
        pool.close()
        pool.join()
        raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))
print('Complete')

