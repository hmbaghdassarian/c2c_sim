#!/usr/bin/env python
# coding: utf-8

# In[1]:


#---prevent tensorly non_negative_parafac from parallelizing
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import tensorly as tl
tl.set_backend(tl.get_backend(), local_threadsafe=True)

#------
import pickle
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from tqdm import trange
import copy

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

n_iter = 100#0 # of iterations per noise value
n_cores = 15 # number of cores for parallelization
par = True


# In[2]:


with open(data_path + 'sim_obj_v' + str(version) + '.pickle', 'rb') as f:
    sim_0 = pickle.load(f)

fn_1 = data_path + 'error_vs_noise_' + str(version) + '_noise0_05.tab' # save results here
fn_2 = data_path + 'error_vs_noise_' + str(version) + '_noise0_25.tab' # save results here

for fn in [fn_1, fn_2]:
    if not os.path.isfile(fn):
        with open(fn, 'w') as f:
            f.write('Noise' + '\t' + 'Error' + '\t' + 'Rank' + '\n')
    else:
        raise ValueError('Unexpected results file already exists')


# In[4]:


def get_error(noise, noise_max, fn):
    try:
        #-------new patterns using original LR and CC networks (otherwise sim = sim_0)
        sim = Simulate()
        sim.LR = copy.deepcopy(sim_0.LR)
        sim.cci = copy.deepcopy(sim_0.cci)
        sim.generate_tensor_md(n_patterns = 4, n_conditions = 12, patterns = ['pulse', 'linear', 'oscillate', 'power'], 
                              consider_homotypic = True, score_change = 'max')

        #-------new communication scores
        sim.generate_tensor(noise = noise, binary = False, bulk = True, noise_max = noise_max)
        sim.reshape()
        #-------decompose
        tensor_ = BaseTensor()
        tensor_.tensor = sim.sim_tensor.tensor_cci

        error = list()
        for r_add in range(3):
            tensor_.compute_tensor_factorization(rank=sim.sim_tensor.rank + r_add,
                                               tf_type='non_negative_cp',
                                               init='svd',
                                               random_state=None,
                                               verbose=False)
            err = _compute_norm_error(tensor_.tensor, tensor_.tl_object)
            error.append(err)

#         err = min(error)
        for err, rank in zip(error, list(np.array(range(3)) + sim.sim_tensor.rank)):
            with open(fn, 'a') as f:
                f.write(str(noise) + '\t' + str(err) + '\t' + str(rank) + '\n')
    except:
        with open(fn, 'a') as f:
            f.write('fail' + '\t' + 'fail' + '\n')


# In[ ]:


print('Cap background noise at 0.05')
for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
    print('------------------------------------------')
    print('Noise: {:.2f}'.format(noise))
    
    if par:
        pool = multiprocessing.Pool(processes = n_cores)
        try:
            pool.starmap(get_error, zip([noise]*n_iter, [0.05]*n_iter, [fn_1]*n_iter))
            pool.close()
            pool.join()
        except:
            pool.close()
            pool.join()
            raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))
    else:
        for i in trange(n_iter):
            get_error(noise, 0.05, fn_1)


# In[ ]:


print('Cap background noise at 0.25')
for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
    print('------------------------------------------')
    print('Noise: {:.2f}'.format(noise))
    
    if par:
        pool = multiprocessing.Pool(processes = n_cores)
        try:
            pool.starmap(get_error, zip([noise]*n_iter, [None]*n_iter, [fn_2]*n_iter))
            pool.close()
            pool.join()
        except:
            pool.close()
            pool.join()
            raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))
    for i in trange(n_iter):
        get_error(noise, None, fn_1)
print('Complete')


# In[ ]:





# In[ ]:




