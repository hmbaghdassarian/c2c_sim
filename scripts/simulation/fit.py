#!/usr/bin/env python
# coding: utf-8

# In[11]:


import numpy as np
import scipy.optimize as op
import scipy.special as sp
import time


# In[ ]:


""" Contains functions used in fitting the power-law, exponential, log-normal,
Weibull (stretched exponential), and power-law with exponential cutoff, as well
as plpval() to find a p-value for the power-law fit. All should be called
directly. All distributions are discrete.

This code is directly copied from https://github.com/adbroido/SFAnalysis/blob/master/code/fit.py
Citation: Anna D. Broido & Aaron Clauset, "Scale-free networks are rare", Nature Communications 10, 1017 (2019).
"""


# In[22]:


def pl(x):
    """ Fits a tail-conditional power-law to a data set. This implements brute
    force optimization (grid search) instead of using a built in optimizer. The
    grid on alpha runs from alstart to alstart+shift. This is based on Aaron's
    plfit.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).
    
    This code is directly copied from https://github.com/adbroido/SFAnalysis/blob/master/code/fit.py
    Citation: Anna D. Broido & Aaron Clauset, "Scale-free networks are rare", Nature Communications 10, 1017 (2019).

    Input:
        x           ndarray, ndim = 1, dtype = integer
    Output:
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        ntail        int, number of datapoints above (and including) xmin
        L            float, log likelihood of the returned fit
        ks           float, goodness of fit statistic (Kolmogorov-Smirnov)
    
    """
    # find the and sort unique possible xmin values
    xminV = np.trim_zeros(np.unique(x))
    # initialize array of the fits for every xmin
    fitV = np.zeros([len(xminV),2])

    start_time = time.time()
    # initialize vector of constants
    # where the xmins start
    xminprev = min(xminV) - 1
    # initialize array of possible alpha values
    alstart = 1.01
    shift = 9.50
    alphaV = np.arange(alstart,alstart+shift,0.01)
    zetaV = sp.zeta(alphaV)
    constV = zetaV
    # shift up to start at the smallest xmin
    for j in range(xminprev):
        constV += -(1+j)**(-alphaV)

    # loop over the xmin values at find the best fit at each
    for i in range(len(xminV)):
        xmin = xminV[i]
        xtail = x[x>=xmin]
        ntail = len(xtail)
        # optimize over alpha
        # find the corresponding array of conditional log likelihoods
        Ls = -alphaV*np.sum(np.log(xtail)) - ntail*np.log(constV)
        # pull out the location of the best fit alpha
        aind = Ls.argmax()
        # find what alpha value is at this index
        alpha = alphaV[aind]
        # compute the KS statistic
        # theoretical cdf
        cdf = np.cumsum(range(np.min(xtail), np.max(xtail)+1)**(-alpha)/constV[aind])
        #  binned data
        xhist = np.histogram(xtail,range(np.min(xtail), np.max(xtail)+2))
        # empirical cdf
        edf = np.cumsum(xhist[0])/float(ntail)
        # KS stat
        ks = np.max(np.abs(cdf-edf))
        # add this KS stat and alpha to the array of fits
        fitV[i] = np.array([ks, alpha])
        # update the constants
        for j in range(xmin-xminprev):
            constV += -(xmin+j)**(-alphaV)
        xminprev = xmin


    # pull out the index of the smallest KS stat
    ksind = fitV[:,0].argmin()
    ks = fitV[ksind,0]
    # find the corresponding xmin
    xmin = xminV[ksind]
    # find the corresponding alpha
    alpha = fitV[ksind,1]
    # evaluate the likelihood here
    xtail = x[x>=xmin]
    ntail = len(xtail)
    start_time = time.time()
    const = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    L = -alpha * np.sum(np.log(xtail)) - ntail*np.log(const)
    # print "-------%s seconds -----------" %(time.time()-start_time)
    # print "alpha = %s" %alpha
    # print "xmin = %s" %xmin
    return [alpha,xmin, ntail, L, ks]

def plpval(x, alpha, xmin, gof):
    """ Finds p-value for the power-law fit using a KS test. This is based on
    Aaron's plpva.m Matlab code (http://tuvalu.santafe.edu/~aaronc/powerlaws/).
    
    This code is directly copied from https://github.com/adbroido/SFAnalysis/blob/master/code/fit.py
    Citation: Anna D. Broido & Aaron Clauset, "Scale-free networks are rare", Nature Communications 10, 1017 (2019).
    
    Input:
        x            ndarray, ndim = 1, dtype = integer
        alpha        float, exponent on x, must be > 1
        xmin         int, starting point for power law tail, must be >= 1
        gof           float, goodness of fit statistic (Kolmogorov-Smirnov)
    Output:
        p            p-value of the returned fit (reject PL hypothesis for p<0.1)
    """
    # set desired precision level in p-value
    eps = 0.01
    #num_resamps = int(np.ceil((1./4)*eps**(-2)))
    num_resamps = 1000
    bootstraps = np.zeros(num_resamps)
    n = len(x)
    xmax = np.max(x)
    tailinds = x>=xmin
    xtail = x[tailinds]
    xhead = x[~tailinds]
    ntail = len(xtail)
    nhead = len(xhead)
    ptail = float(ntail)/n
    mmax = 20*xmax
    # set the tail of the pdf
    #const_tail = ic.plconst(np.array(alpha),xmin)
    const_tail = sp.zeta(alpha) - np.sum(np.arange(1,xmin)**(-alpha))
    pdf_tail = np.arange(xmin,mmax+1)**(-alpha)/const_tail # i.e.; end at mmax
    # pad this with zeros (rewrite if we don't need to do this)
    pdf = np.zeros(mmax+1)
    pdf[xmin:] = pdf_tail
    # clean up in case this is a huge array
    del pdf_tail
    # set the cdf. rows are x-val or cdf(xval). So cdf(x=10) is cdf[1,10]
    cdf = np.array( [ np.arange(mmax+1), np.cumsum(pdf) ] )
    # tack on a last entry
    cdf = np.concatenate( (cdf , np.array([[mmax+1,1]]).T) , axis = 1 )

    # semi-parametric bootstrap
    starttime = time.time()
    for resamp_ind in range(num_resamps):
        # non-parametric bootstrap from the head of x
        # count how many of n random numbers are in the head, based on the probability of being in the head of x
        nnewhead = n
        while nnewhead >= n:
            nnewhead = np.sum(np.random.rand(n)>ptail)
        headinds = np.array([np.floor(nhead*np.random.rand(nnewhead))],dtype=int)
        newhead = xhead[headinds][0]
        nnewtail  = n-nnewhead

        # parametric bootstrap for the powerlaw tail
        rtail = np.sort(np.random.rand(nnewtail))
        newtail = np.zeros(nnewtail, dtype = int)
        indrtail = 0
        indnewtail = 0
        for xval in range(xmin, mmax+2):
            while (indrtail < len(rtail)) and (rtail[indrtail] <= cdf[1, xval]):
                indrtail += 1
            newtail[indnewtail:indrtail] = xval
            indnewtail = indrtail
            if indnewtail > nnewtail:
                break
        # combine into new sample
        newx = np.concatenate((newhead, newtail))
        if (newx == np.zeros_like(newx)).all():
            import pdb; pdb.set_trace()
        # fit this new sample
        [newalpha, newxmin, newntail, newLpl, newgof] = pl(newx)
        # print where we are
        current_p = np.sum(bootstraps[0:resamp_ind]>=gof)/(float(resamp_ind+1))
        # print "[%s]    p = %f" %(resamp_ind, current_p)
        # store gof stat
        bootstraps[resamp_ind] = newgof
        # if it's taking forever and we can end, do it
        if time.time() - starttime > 500:
            if resamp_ind > num_resamps/20.:
                if current_p<0.05 or current_p>0.5:
                    print("current p = %s   elapsed time = %s" %(current_p, time.time()-starttime))
                    return current_p
    p = np.sum(bootstraps>=gof)/float(num_resamps)
    print ("p = %.3f   elapsed time = %s" %(p, time.time()-starttime))
    return p


# In[43]:


# x = np.random.power(2,1000)
# x *= 1000
# x = np.array([round(i) for i in x])
x = np.array([8, 10, 18, 8, 12, 34, 15, 25, 31, 48, 46, 48, 40, 32, 34, 35, 19, 21, 29, 20, 28, 35, 27, 21, 17, 20, 21, 24, 22, 20, 20, 19, 19, 21, 13, 29, 31, 17, 12, 13, 16, 13, 16, 11, 14, 12, 14, 14, 10, 14, 12, 11, 11, 11, 14, 13, 14, 11, 14, 12, 15, 13, 15, 10, 14, 10, 11, 13, 10, 12, 12, 10, 10, 11, 12, 10, 9, 13, 10, 9, 10, 10, 9, 9, 9, 9, 10, 10, 9, 10, 11, 9, 9, 11, 10, 9, 9, 9, 9, 9])
[alpha, xmin, ntail,  L, ks] = pl(x)
p = plpval(x,alpha, xmin, ks)

