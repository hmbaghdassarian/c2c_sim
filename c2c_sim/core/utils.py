#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import random
import numpy as np
from scipy import optimize
import numpy as np
from scipy.stats import truncnorm
import itertools


# In[2]:


def get_truncated_normal(n, sd, mean=0, low=0, upp=1, seed = None):
    np.random.seed(seed) # reproducibility
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs(n)


# In[18]:


relat = pd.DataFrame(columns = ['mean', 'noise_mean'])
counter = 0
iter_ = np.arange(0.01,1,0.001)
seeds = list(range(len(iter_))) # reproducibility

for mean in iter_:
    relat.loc[counter,:] = [mean, np.mean(get_truncated_normal(n = 10**5, sd = mean, 
                                                               seed = seeds[counter]))] # noise = 1
    counter += 1

y = relat['noise_mean'].astype(float).values
x = relat['mean'].astype(float).values
             
def piecewise_fit(x, x0, m, b, a,c, d):
    return np.piecewise(x, [x < x0], [lambda x: m*x+b, lambda x: ((x-d)/c)**(1/a)])

fit_params , e = optimize.curve_fit(piecewise_fit, x, y, 
                          p0 = [0.37, 7.93278070e-01, 5.11093685e-04, 
                               6.21557667, 78.41240029, 0.34450856])


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

