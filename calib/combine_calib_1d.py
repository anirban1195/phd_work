#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 12:33:23 2023

@author: dutta26
"""

import numpy as np
# =============================================================================
# a = np.load('/scratch/bell/dutta26/abell_2390/calib_1001.npy')
# b = np.load('/scratch/bell/dutta26/abell_2390/calib_1002.npy')
# c = np.load('/scratch/bell/dutta26/abell_2390/calib_1003.npy')
# 
# c= np.vstack((a[:, :],b[:,:], c[:,:]))
# 
# c[:,0] = np.log10(c[:,0]/c[:,1])
# c[:,2] = c[:,2]/c[:,1]
# c=c[c[:, 0].argsort()]
# 
# 
# count = 0
# part = 0.025
# l = len(np.arange(-3.3, 2.15, part))
# store = np.zeros((l, 3))
# for val in np.arange(-3.3, 2.15, part):
#     indices = np.where((c[:,0]>val-part )  & (c[:,0]<val+part))
#     print (len(indices[0]))
#     store[count, 0], store[count, 1], store[count, 2] = val, np.median(c[indices,2]), np.median(c[indices,3])
#     
#     count += 1
# np.save('/scratch/bell/dutta26/abell_2390/calib_final.npy', store)
# =============================================================================


a = np.load('/scratch/bell/dutta26/abell_2390/calib_1001_inf.npy')
b = np.load('/scratch/bell/dutta26/abell_2390/calib_1002_inf.npy')
c = np.load('/scratch/bell/dutta26/abell_2390/calib_1003_inf.npy')
c1 = np.load('/scratch/bell/dutta26/abell_2390/calib_1004_inf.npy')

c= np.vstack((a[:, :],b[:,:], c[:,:], c1[:,:]))

c[:,0] = np.log10(c[:,0]/c[:,1])
c[:,2] = c[:,2]/c[:,1]
c=c[c[:, 0].argsort()]


count = 0
part = 0.05
l = len(np.arange(1.175, 2.9, part))
store = np.zeros((l, 3))
for val in np.arange(1.175, 2.9, part):
    indices = np.where((c[:,0]>val-part )  & (c[:,0]<val+part))
    print (len(indices[0]))
    store[count, 0], store[count, 1], store[count, 2] = val, np.median(c[indices,2]), np.median(c[indices,3])
    
    count += 1
np.save('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy', store)

