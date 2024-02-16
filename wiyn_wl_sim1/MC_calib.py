#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 15:22:43 2023

@author: dutta26
"""

#1764,3363, 5067

import helper 
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
ratioArr = np.arange(0.1, 2, 0.1)
x=[]
y=[]
for ratio in ratioArr:
    arr=[]
    
    temp1=[]
    ellip =[]
    
    for j in range(500):
        sigxx = np.random.uniform (1, 10)
        sigyy= np.random.uniform (sigxx+3, 3*sigxx)
        sigxy = 0
        e1 = (sigxx-sigyy)/(sigxx+sigyy)
        err= ratio*np.sqrt(sigxx+sigyy)
        #e2 = 2*sigxy/(sigxx+sigyy)
        corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                             0, 0, 0, 
                                                             err,err,err*0.707, 
                                                             0,0,0)
        corr_e1= (corr_xx-corr_yy)/(corr_xx+corr_yy)
        arr.append(e1/corr_e1)
        temp1.append(e1/corr_e1)
        ellip.append(abs(e1))
        #plt.plot(ratio, e1/corr_e1, 'b.', markersize=2)
    
    plt.scatter(ratio*np.ones(len(ellip)), temp1, s=2, c=ellip)    
    temp1=[]
    ellip =[]
    for j in range(500):
        sigyy = np.random.uniform (1, 20)
        sigxx= np.random.uniform (sigyy+3, 3*sigyy)
        sigxy = 0
        e1 = (sigxx-sigyy)/(sigxx+sigyy)
        err= ratio*np.sqrt(sigxx+sigyy)
        #e2 = 2*sigxy/(sigxx+sigyy)
        corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                             0, 0, 0, 
                                                             err,err,err*0.707, 
                                                             0,0,0)
        corr_e1= (corr_xx-corr_yy)/(corr_xx+corr_yy)
        arr.append(e1/corr_e1)
        temp1.append(e1/corr_e1)
        ellip.append(abs(e1))
        #plt.plot(ratio, e1/corr_e1, 'b.', markersize=2)
    plt.scatter(ratio*np.ones(len(ellip)), temp1, s=2, c=ellip)    
    
    mean, med, std = sigma_clipped_stats(arr)
    plt.plot(ratio, mean, 'k.')
    x.append(ratio)
    y.append(mean)
plt.plot(x,y, 'k-')
    
# =============================================================================
# ratioArr = np.arange(0.1, 2, 0.1)
# x=[]
# y=[]
# for ratio in ratioArr:
#     arr=[]
#     for j in range(500):
#         sigxx = np.random.uniform (1, 30)
#         sigyy= sigxx
#         sigxy = np.random.uniform(0.3, 0.25*(sigxx+sigyy))
#         #e1 = (sigxx-sigyy)/(sigxx+sigyy)
#         err= ratio*np.sqrt(sigxx+sigyy)
#         e2 = 2*sigxy/(sigxx+sigyy)
#         corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
#                                                              0, 0, 0, 
#                                                              err,err,err*0.707, 
#                                                              0,0,0)
#         corr_e2 = 2*corr_xy/(corr_xx+corr_yy)
#         arr.append(e2/corr_e2)
#         plt.plot(ratio, e2/corr_e2, 'b.', markersize=2)
#         
#     for j in range(500):
#         sigxx = np.random.uniform (1, 30)
#         sigyy= sigxx
#         sigxy = np.random.uniform(-0.25*(sigxx+sigyy), -0.5)
#         #e1 = (sigxx-sigyy)/(sigxx+sigyy)
#         err= ratio*np.sqrt(sigxx+sigyy)
#         e2 = 2*sigxy/(sigxx+sigyy)
#         corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
#                                                              0, 0, 0, 
#                                                              err,err,err*0.707, 
#                                                              0,0,0)
#         corr_e2 = 2*corr_xy/(corr_xx+corr_yy)
#         arr.append(e2/corr_e2)
#         plt.plot(ratio, e2/corr_e2, 'b.', markersize=2)
#         
#     mean, med, std = sigma_clipped_stats(arr)
#     plt.errorbar(ratio, mean, yerr= std, fmt='k.')
#     x.append(ratio)
#     y.append(mean)
# plt.plot(x,y, 'k-')   
#         #corr_e2 = 2*corr_xy/(corr_xx+corr_yy)
#         
#         #print (e1, corr_e1, e1/corr_e1)
#         #print (e2, corr_e2, e2/corr_e2)
# =============================================================================
