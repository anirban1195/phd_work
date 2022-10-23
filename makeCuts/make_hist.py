#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 22:47:36 2021

@author: dutta26
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab

store = np.load('/home/dutta26/test14_4.npy')
store[np.isnan(store)] = 0
temp1=[]
temp =store[:,3,0]
psf = store[:,3,6]
temp2 = store[:,3,6]
for j in range(len(temp)):
    if(temp[j] != 0 and abs(temp[j])< 200 ):
        temp1.append(temp2[j])

# =============================================================================
# temp =np.sqrt(store[:,2,1]**2 + store[:,2,2]**2)
# temp1 =[]
# for j in range(len(temp)):
#     if(temp[j] != 0 ):
#         temp1.append(temp[j])
# =============================================================================

# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = plt.hist(x=temp1,  color='black',
                            alpha=0.3, rwidth=0.95, normed=1)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('PSF')
plt.ylabel('Frequency')
#plt.title('Faint Object')
#plt.savefig('/home/dutta26/sine.png')
mu,sigma =norm.fit(temp1)
best_fit_line = norm.pdf(bins,mu,sigma)
plt.plot(bins, best_fit_line, 'k--')
#plt.plot([0.0,0.00, [0,10.851], 'y--')