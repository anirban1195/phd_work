#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 11:49:59 2023

@author: dutta26
"""
import helper
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
th1Arr=[]
th2Arr=[]
e1Arr=[]
e2Arr =[]
e1DiffArr =[]
e2DiffArr =[]
cnt = 0
for  j in range(10000):
    x = np.random.uniform(0.2,1.5)
    #y = np.random.uniform(0.75*x,1.5*x)
    y = np.random.normal(1.2*x , 0.2*x)
    area =  2* np.pi*np.sqrt(6*6 )
    B =70000
    N = np.random.uniform(8,12)*np.sqrt(4*area*B)
    
    
    angle = np.random.uniform(0, np.pi)
    xx = 0.5*( (np.sin(angle)* y)**2  + (np.cos(angle)* x)**2 )
    yy = 0.5*( (np.sin(angle)* x)**2  + (np.cos(angle)* y)**2 )
    xy = 0.5*(x**2 - y**2)*np.cos(angle)*np.sin(angle)
    
    e1 = (xx-yy)/(xx+yy)
    e2 = 2*xy/(xx+yy)
    theta = 0.5*np.arctan2(e2,e1)*180/np.pi
    th1Arr.append(theta)
    e1Arr.append(np.sqrt(e1**2 + e2**2))
    
    
    s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * np.sqrt(6 + 6) 
    print (s)
    median_size = 3.5
    s_xx = np.sqrt(((1+e1)*s)**2 + ((median_size/30)*median_size*1.414)**2)
    s_yy = np.sqrt(((1-e1)*s)**2 + ((median_size/30)*median_size*1.414)**2)
    s_xy = np.sqrt((0.707*(1+abs(e2))*s)**2 + ((median_size/30)*median_size*1.414 * 0.707)**2)
    
# =============================================================================
#     s = np.sqrt( (1.3/(np.pi*N) + 2.2*4*area * B/(np.pi * N**2)) ) * 6 *2
#     print (s)
#     median_size = 3.5
#     s_xx = np.sqrt(((1+0)*s)**2 + ((median_size/30)*median_size*1.414)**2)
#     s_yy = np.sqrt(((1-0)*s)**2 + ((median_size/30)*median_size*1.414)**2)
#     s_xy = np.sqrt((0.707*(1+abs(0))*s)**2 + ((median_size/30)*median_size*1.414 * 0.707)**2)
#     
# =============================================================================
    
    
    
    xx += np.random.normal(0*s_xx, 0.75*s_xx)
    yy += np.random.normal(0*s_yy, 0.75*s_yy)
    xy += np.random.normal(0*s_xy, 0.75*s_xy)
    
    
    area = 2* np.pi*np.sqrt(6*6 )
    
    median_psf_flux = 500000
    s_psf = np.sqrt( (area/(np.pi*median_psf_flux) + 4*area**2 * B/(np.pi * median_psf_flux**2)) ) * np.sqrt(6 + 6)
    e1_psf = np.random.uniform(-0.1, 0.1)
    e2_psf = np.random.uniform(-0.1, 0.1)
    s_xx_psf = np.sqrt(((median_size/30)*median_size*1.414)**2 + ((1+e1_psf)*s_psf)**2)
    s_yy_psf = np.sqrt(((median_size/30)*median_size*1.414)**2 + ((1-e1_psf)*s_psf)**2)
    s_xy_psf = np.sqrt(((median_size/30)*median_size*1.414*0.707 )**2 + ((1+abs(e2_psf))*s_psf*0.707)**2)
    
    
    print (s_xx, s_xy)
    corr_xx, corr_yy, corr_xy , success = helper.correct(xx, yy, xy, 
                                                         0, 0, 0, 
                                                         s_xx,s_yy, s_xy, 
                                                         s_xx_psf, s_yy_psf, s_xy_psf)
    
    if(corr_xx == 0):
        corr_xx = corr_yy = 1
        corr_xy = 0
        cnt += 1
    e1_c = (corr_xx-corr_yy)/(corr_xx+corr_yy)
    e2_c = 2*corr_xy/(corr_xx+corr_yy)
    e2Arr.append(np.sqrt(e1_c**2 + e2_c**2))
    
    theta = 0.5*np.arctan2(e2_c,e1_c)*180/np.pi
    th2Arr.append(theta)
    
    e1DiffArr.append(e1-e1_c)
    e2DiffArr.append(e2-e2_c)

plt.subplot(221)

n, bins, patches = plt.hist(x=th1Arr, bins=40, color='r',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'Input')

n, bins, patches = plt.hist(x=th2Arr, bins=40, color='b',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'Output')
plt.xlabel('Theta')
plt.legend()
plt.subplot(222)
n, bins, patches = plt.hist(x=e1Arr, bins=40, color='r',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'Input')

n, bins, patches = plt.hist(x=e2Arr, bins=40, color='b',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'Output')

plt.xlabel('Elliptcity')
plt.legend()
plt.subplot(223)

n, bins, patches = plt.hist(x=e1DiffArr, bins=40, color='r',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'For e1')

n, bins, patches = plt.hist(x=e2DiffArr, bins=40, color='b',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0, label = 'For e2')
plt.xlabel(r'$e_n True - e_n corrected (n=1,2)$')
plt.legend()

plt.subplot(224)

n, bins, patches = plt.hist(x=np.array(th1Arr) - np.array(th2Arr), bins=40, color='r',density=True,
                            alpha=0.95, histtype='step', linewidth=3.0)
plt.xlabel('Actual Phi - Measured Phi')
plt.legend()