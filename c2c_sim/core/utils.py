#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import numpy as np
from scipy import optimize
import numpy as np
from scipy.stats import truncnorm


# In[1]:


def get_truncated_normal(n, sd, mean=0, low=0, upp=1):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs(n)


# In[4]:


relat = pd.DataFrame(columns = ['mean', 'noise_mean'])
counter = 0
for mean in np.arange(0.01,1,0.001):
    relat.loc[counter,:] = [mean, np.mean(get_truncated_normal(n = 10**5, sd = mean))] # noise = 1
    counter += 1

y = relat['noise_mean'].astype(float).values
x = relat['mean'].astype(float).values
             
def piecewise_fit(x, x0, m, b, a,c, d):
    return np.piecewise(x, [x < x0], [lambda x: m*x+b, lambda x: ((x-d)/c)**(1/a)])

fit_params , e = optimize.curve_fit(piecewise_fit, x, y, 
                          p0 = [0.37, 7.93278070e-01, 5.11093685e-04, 
                               6.21557667, 78.41240029, 0.34450856])

