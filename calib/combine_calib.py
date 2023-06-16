#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 21:48:26 2023

@author: dutta26
"""

import numpy as np
a = np.load('/scratch/bell/dutta26/abell_2390/calib_1inf.npy')
b = np.load('/scratch/bell/dutta26/abell_2390/calib_2inf.npy')
b1 = np.load('/scratch/bell/dutta26/abell_2390/calib_3inf.npy')

c= np.vstack((a[:, :],b[:,:], b1[:,:]))
np.save('/scratch/bell/dutta26/abell_2390/calib_combined1_inf.npy', c)    
# =============================================================================
# temp = np.arange(0, 6642, 41)
# ind = np.argsort(c[temp,1])
# 
# d = np.zeros((6642, 4))
# counter = 0
# for index in ind:
#     d[counter:counter+41, :] = c[index*41:index*41+41,:]
#     counter += 41
#     
# np.save('/scratch/bell/dutta26/abell_2390/calib_combined.npy', d)    
# =============================================================================
