#!/usr/bin/env python
# coding: utf-8

# In[1]:


#---prevent tensorly non_negative_parafac from parallelizing
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import tensorly as tl
tl.set_backend(tl.get_backend(), local_threadsafe=True)

#---prevent tqdm called from other functions when parallelizing
n_iter = 1000 # of iterations per noise value
n_cores = 15 # number of cores for parallelization
par = True # whether to parallelize

if par:
    import tqdm
    def nop(it, *a, **k):
        return it
    tqdm.tqdm = nop
else:
    from tqdm import trange

#------
import pickle
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
import copy

import warnings
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

import sys
sys.path.insert(1, '../cell2cell/')
import cell2cell as c2c
from cell2cell.tensor.tensor import BaseTensor
from cell2cell.tensor.factorization import _compute_norm_error, _compute_tensor_factorization

sys.path.insert(1, '../c2c_sim/')
from core.simulate import Simulate, CCI_MD


# specify path to save figures
fig_path = 'figures/'#'/data2/hratch/cci_dt/figures/'
data_path = 'data/'#'/data2/hratch/cci_dt/noise_v_error/'
version = 1


# In[2]:


with open(data_path + 'sim_obj_v' + str(version) + '.pickle', 'rb') as f:
    sim_0 = pickle.load(f) # from sim_tensor_decompose_vis script

fn_0 = data_path + 'error_vs_noise_' + str(version) + '_noise0_1.tab' # save results here
fn_1 = data_path + 'error_vs_noise_' + str(version) + '_noise0_15.tab' # save results here
fn_2 = data_path + 'error_vs_noise_' + str(version) + '_noise0_20.tab' # save results here
fns = [fn_0,fn_1,fn_2]

for fn in fns:
    if not os.path.isfile(fn):
        with open(fn, 'w') as f:
            f.write('Noise' + '\t' + 'Error' + '\t' + 'Rank' + '\n')
    else:
        raise ValueError('Unexpected results file already exists')


# In[4]:


def get_error(noise, noise_max, fn):
    with warnings.catch_warnings(): # ignore warnings
        warnings.simplefilter("ignore")
        with HiddenPrints(): #ignore prints of imported functions
            try:
                #-------new patterns/simulations
                sim = Simulate() 
                # simulate a scale_free randomly connected ligand-receptor network (potential interactions)
                sim.LR_network(network_type = 'scale-free', **{'nodes': 100, 'degrees': 3, 'alpha': 2}) #scale-free
                print('{} ligands, {} receptors, and {} edges'.format(len(sim.LR.ligands), len(sim.LR.receptors), 
                                                                      len(sim.LR.edge_list)))
                sim.LR.generate_metadata(n_LR_cats = {3: 0}, cat_skew = 0)
                cci = CCI_MD()
                cci.cci_network(n_cells = 50, directional = True, autocrine = True)
                cci.generate_metadata(n_cell_cats = {3: 0}, cat_skew = 0, remove_homotypic = 0)
                # add cell metadata to simulation object
                sim.cci = cci
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

                for err, rank in zip(error, list(np.array(range(3)) + sim.sim_tensor.rank)):
                    with open(fn, 'a') as f:
                        f.write(str(noise) + '\t' + str(err) + '\t' + str(rank) + '\n')
            except:
                with open(fn, 'a') as f:
                    f.write('fail' + '\t' + 'fail' + '\n')


# In[ ]:


counter = 0
for nm in [0.1, 0.15, 0.20]:
    print('Cap background noise at {:.2f}'.format(nm))
    for noise in list(np.arange(0.01,0.1, 0.01)) + list(np.arange(0.1,1.01,0.1)):
        print('------------------------------------------')
        print('Noise: {:.2f}'.format(noise))

        if par:
            pool = multiprocessing.Pool(processes = n_cores)
            try:
                pool.starmap(get_error, zip([noise]*n_iter, [nm]*n_iter, [fns[counter]]*n_iter))
                pool.close()
                pool.join()
            except:
                pool.close()
                pool.join()
                raise ValueError('Parallelization failed at noise: {:.2f}'.format(noise))
        else:
            for i in trange(n_iter):
                get_error(noise, nm, fn[counter])
    counter += 1
print('Complete')

