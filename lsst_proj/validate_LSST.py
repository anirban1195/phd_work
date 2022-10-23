#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:09:13 2022

@author: dutta26
"""
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess
plt.style.use('seaborn-whitegrid')

band='5'

ir_df = np.array(pd.read_pickle('/scratch/halstead/d/dutta26/m38/df_'+band+'.pk1'))
store = np.load('/scratch/halstead/d/dutta26/m38/singleFrame_'+band+'.npy')
fluxArr=[]
meanArr=[]
std_devArr=[]
for j in range(20):
    loc = (np.where((store[j,:,2] == 1) & (ir_df[:,3] < 50000) & (store[j,:,12] == 99)  & 
                    (ir_df[:,3] > 2000) & 
                    (store[j,:,14] == 0)))[0]

    mean,median, std = sigma_clipped_stats(store[j,loc, 7])
    #print (sigma_clipped_stats(store[j,loc, 7]))
    
    #a = store[j,loc, 31]/ ir_df[loc, 3]
    e1= (store[j,loc, 7] - store[j,loc, 8])/(store[j,loc, 7] + store[j,loc, 8])
    e2 = 2 * store[j,loc, 9]/(store[j,loc, 7] + store[j,loc, 8])
    a = np.sqrt(e1**2 +e2**2)
    
    print (sigma_clipped_stats(a))
    meanArr.append(sigma_clipped_stats(a)[0])
    std_devArr.append(sigma_clipped_stats(a)[2])


    
plt.errorbar(np.arange(0,20), meanArr, yerr=std_devArr, fmt='.k');    

# =============================================================================
# highArr=np.array([5000, 2000, 1000, 750, 500, 250, 100, 75])
# lowArr =np.array([2000, 1000, 750, 500, 250, 100, 75, 50])
# median =[]
# err =[]
# for k in range(len(lowArr)):
#     loc = np.where((ir_df[:,3] < highArr[k]) & (ir_df[:,3] > lowArr[k]) ) [0]
#     totArr=[]
#     for index in loc:
#         temp =store[:,index,:]
#         fluxArr=[]
#         for j in range(20):
#             if(temp[j,12] == 99):
#                 fluxArr.append(temp[j,3])
#             elif(temp[j,12] == 1):
#                 fluxArr.append(temp[j,31])
#         fluxArr = np.array(fluxArr)
#         fluxArr = fluxArr[fluxArr>0]
#         totArr.append(sigma_clipped_stats(fluxArr)[1]/ ir_df[index, 3])
#         
#     print (sigma_clipped_stats( totArr))
#     median.append(sigma_clipped_stats( totArr)[1] )
#     err.append(sigma_clipped_stats( totArr) [2])
# avg = (highArr+lowArr)/2   
# plt.errorbar(np.log10(avg/np.sqrt(300*9)), median, yerr=err, fmt='.k')
# =============================================================================
