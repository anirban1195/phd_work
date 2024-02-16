#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 12:36:37 2023

@author: dutta26
"""
import numpy as np
import matplotlib.pyplot as plt


c=  np.load('/scratch/bell/dutta26/abell_2390/calib_forcedDist1_err1001.npy')

c[:,0] = np.log10(c[:,0]/c[:,1])
c=c[c[:, 0].argsort()]


count = 0
part = 0.1
l = len(np.arange(-2.8, 2.1, part))
store = np.zeros((l, 4))
for val in np.arange(-2.8, 2.1, part):
    indices = np.where((c[:,0]>val-part )  & (c[:,0]<val+part))
    print (len(indices[0]))
    store[count, 0], store[count, 1], store[count, 2],store[count, 3] = val, np.median(c[indices,2]), np.median(c[indices,3])*2,np.median(c[indices,8])*2
    
    count += 1

np.save('/home/dutta26/codes/forced_truncated_calib.npy', store)

