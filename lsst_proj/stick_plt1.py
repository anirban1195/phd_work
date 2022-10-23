#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 16:01:40 2022

@author: dutta26
"""

import matplotlib.pyplot as plt

import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import numpy as np
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord   
from astropy.coordinates import FK5 
import measure_pythonV
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import griddata
import sys,helper,os

ir_coadd_data = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df2_2.pk1')
store = np.array(ir_coadd_data)
#loc = np.where((store[:,2] == 1) & (store[:,3] > 500) & (store[:,7] < 2.5) &(store[:,7] > 0.5))[0]

loc = np.where((store[:,2] == 1) & (store[:,3] > 10000) & 
               (store[:,3] < 5000000) &(store[:,7] > 0.35) 
               &(store[:,7] < 2.85))[0]

# =============================================================================
# ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
# store = np.array(ir_coadd_data)
# loc = np.where((store[:,2] == 1) & (store[:,3] > 500) & (store[:,7] < 8) &(store[:,7] > 3))[0]
# =============================================================================
# =============================================================================
# store = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data_final.npy')
# loc= (np.where((store[:,2] == 1) & (store[:,3] > 1000)  &(store[:,3] < 40000) 
#                & (store[:,13] == 1) & (store[:,7] >3) & (store[:,7] <6) ))[0]
# =============================================================================
a=store[loc,7]
b=store[loc,7]- store[loc,38]
xList=store[loc,10]
yList=store[loc, 11]
sizeList = np.sqrt(store[loc,7] + store[loc,8])
#e1 = (store[loc,7] - store[loc,8])/(store[loc,7] + store[loc,8])
#e2 = 2* store[loc,9] /(store[loc,7] + store[loc,8])
sigmaxx = store[loc,7] -store[loc,38] + 0.5
sigmayy = store[loc,8] -store[loc,39] + 0.5
sigmaxy = store[loc,9] -store[loc,40] 

print (sigma_clipped_stats(sigmaxx))
print (sigma_clipped_stats(sigmayy))
print (sigma_clipped_stats(sigmaxy))

e1 = (sigmaxx-sigmayy)/(sigmaxx+sigmayy)
e2 = 2*sigmaxy / (sigmaxx+sigmayy)

print (sigma_clipped_stats(e1))
print (sigma_clipped_stats(e2))
    
# =============================================================================
# e1 = (store[loc,7] - store[loc,8] +store[loc,39])/(store[loc,7] + store[loc,8] )
# e2 =  2* (store[loc,9])  /(store[loc,7] + store[loc,8])
# e1 = e1+0.01
# e2 = e2+0.01
# =============================================================================

ellip = np.sqrt(e1**2+e2**2)
theta = 0.5*np.arctan2(e2,e1)
cosArr = np.cos(theta)
sinArr = np.sin(theta)

fact=5
#Get the rage of values of x and y
xShape = int(10000/fact) + 10
yShape = int(10000/fact) + 10
xList= xList/fact
yList = yList/fact


count = 0
for j in range(len(xList) - 1):
    dist = (ellip[j]/ 0.01) *5
    if(dist > 100):
        continue
    xNew = xList[j] + dist*cosArr[j]
    yNew = yList[j] + dist*sinArr[j]
    xArr = [int(xList[j]), int(xNew)]
    yArr = [int(yList[j]), int(yNew)]
    plt.plot(xArr, yArr, 'b')
    count += 1
print (count)        
sys.exit()


# =============================================================================
# points=[]
# values =[] 
# grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
# 
# for j in range(len(xList)- 1):
#     x = int(round(xList[j]))
#     y = int(round(yList[j]))
#     
#     points.append([y,x])
#     values.append( ellip[j])
# 
# 
# grid_z1 = griddata(points, values, (grid_y, grid_x), method='linear')    
# grid_z1 = np.array(grid_z1, dtype = np.float32).T
# hdu = fits.PrimaryHDU(grid_z1) 
# hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/residue.fits', clobber=True)     
# =============================================================================
