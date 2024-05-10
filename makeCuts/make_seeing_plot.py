#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 14:02:51 2024

@author: dutta26
"""

import os
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/r_sf.npy')
i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/i_sf.npy')
g_sf_df = np.load('/scratch/bell/dutta26/abell_2390/g_sf.npy')
seeing =[]
for j in range(len(r_sf_df)):
    loc = np.where(r_sf_df[j,:,16] > 0)
    if(len(loc[0]) == 0):
        continue
    seeing.append(np.median(r_sf_df[j,loc,16]))
    
print (len(seeing))
plt.hist(x=seeing, bins=15, density=True,histtype=u'step', color='r', label='r')
seeing=[]
for j in range(len(i_sf_df)):
    loc = np.where(i_sf_df[j,:,16] > 0)
    if(len(loc[0]) == 0):
        continue
    seeing.append(np.median(i_sf_df[j,loc,16]))
print (len(seeing))
plt.hist(x=seeing, bins=15, density=True,histtype=u'step', color='k', label='i')
seeing=[]


for j in range(len(g_sf_df)):
    loc = np.where(g_sf_df[j,:,16] > 0)
    if(len(loc[0]) == 0):
        continue
    seeing.append(np.median(g_sf_df[j,loc,16]))
print (len(seeing))
plt.hist(x=seeing, bins=15, density=True,histtype=u'step', color='g', label='g')
seeing=[]  

plt.xlabel('Seeing in arcseconds')
plt.ylabel('Frequency')
plt.legend()