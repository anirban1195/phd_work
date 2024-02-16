#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 12:13:15 2023

@author: dutta26
"""
import numpy as np
import matplotlib.pyplot as plt
errArr=[]
sizeArr=[]
ellipArr =[]
for j in range(10000):
    corr_sigxx = np.random.uniform(1, 8)
    corr_sigyy = np.random.uniform(1, 8)
    corr_sigxy = np.random.uniform(-(corr_sigxx+corr_sigyy)/3, (corr_sigxx+corr_sigyy)/3)
    e1_coadd = (corr_sigxx-corr_sigyy)/(corr_sigxx+corr_sigyy)
    e2_coadd= 2*corr_sigxy/(corr_sigxx+corr_sigyy)
    err_xx = err_yy = 0.45
    err_xy =0.2
    a= np.sqrt( err_xx**2 * 4* (corr_sigxx**2+corr_sigyy**2)/(corr_sigxx+corr_sigyy)**4)
    b = np.sqrt( (err_xx**2 * 8* (corr_sigxy**2)/(corr_sigxx+corr_sigyy)**4) + (err_xy**2 * 4/(corr_sigxx+corr_sigyy)**2))
    c = np.sqrt((e1_coadd**2 *a**2 / (e1_coadd**2 + e2_coadd**2)) + (e2_coadd**2 *b**2 / (e1_coadd**2 + e2_coadd**2)) )
    size = np.sqrt(corr_sigxx + corr_sigyy)
    
    ellipArr.append(e1_coadd**2 + e2_coadd**2)
    sizeArr.append(np.sqrt(corr_sigxx+corr_sigyy))
    errArr.append(c)
plt.plot(sizeArr, errArr, 'b.', markersize =1)