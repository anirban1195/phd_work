#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  9 07:34:05 2023

@author: dutta26
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
from scipy.ndimage import rotate
import pandas as pd
import os
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess

ir_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1'))
loc = np.where((ir_df[:,3]>3) & (ir_df[:,2]==0) )[0]
e1Arr=[]
e2Arr=[]
thetaArr=[]

for j in loc:
    if(ir_df[j, 7] == 0 or ir_df[j, 7]<0 or ir_df[j, 7]== None or np.isnan(ir_df[j, 7])):
        continue
    if(ir_df[j, 8] == 0 or ir_df[j, 8]<0 or ir_df[j, 8]== None or np.isnan(ir_df[j, 8])):
        continue
    if(ir_df[j, 59] == 1):
        continue
    xx = ir_df[j, 35] - ir_df[j, 38]
    yy = ir_df[j, 36] - ir_df[j, 39]
    xy = ir_df[j, 37] - ir_df[j, 40]
    e1 = (xx-yy)/(xx+yy)
    e2 = 2*xy/(xx+yy)
    if(xx<0.1 or yy<0.1):
        continue
    e1Arr.append(e1)
    
    e2Arr.append(e2)
    theta = 0.5*np.arctan2(e2,e1)*180/np.pi
    thetaArr.append(theta)

e1Arr = np.array(e1Arr)
e2Arr = np.array(e2Arr)

eArr = np.sqrt(e1Arr**2 + e2Arr**2)

# =============================================================================
# n, bins, patches = plt.hist(x=e1Arr, bins=40, color='r',
#                             alpha=0.7)
# n, bins, patches = plt.hist(x=e2Arr, bins=40, color='b',
#                             alpha=0.7)
# 
# n, bins, patches = plt.hist(x=eArr, bins=20, color='g',
#                             alpha=0.7)
# =============================================================================

n, bins, patches = plt.hist(x=thetaArr, bins=40, color='r',
                            alpha=0.7)
print (sigma_clipped_stats(eArr), len(eArr))