#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:59:54 2023

@author: dutta26
"""


import sys,os
from astropy.io import fits
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
import helper
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


# =============================================================================
# ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
# loc_coadd = np.where((ir_coadd_data[:,2] == 1) & (ir_coadd_data[:,3]>500) & (ir_coadd_data[:,3]<2000))[0]
# 
# temp = np.zeros((100, len(loc_coadd), 77))
# band = 'r'
# for sliceIndex in range(1,100):
#     #sliceIndex = 1
#     sizex =sizey = 100
#     tot = np.zeros((sizex, sizey))
#     
#     store = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')
#     
#     star_arr = store[sliceIndex,(np.where((store[sliceIndex,:,2] == 1) & (store[sliceIndex,:,12] == 99) 
#                                & (store[sliceIndex,:,13] == 0) & (store[sliceIndex,:,14] == 0)))[0],  : ]
#     
#     
#     #Now tuse k sigma clip to find usable stars. Just do for sigxx
#     mean,median, std = sigma_clipped_stats(star_arr[:, 7])
#     mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
#     mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
#     print (mean,median, std)
#     loc = np.where((store[sliceIndex,:,2] == 1) & 
#                                           (store[sliceIndex,:,7] >= mean-3*std) &
#                                           (store[sliceIndex,:,7] <= mean+3*std) &
#                                           (store[sliceIndex,:,8] >= mean1-3*std1) &
#                                           (store[sliceIndex,:,8] <= mean1+3*std1) &
#                                           (store[sliceIndex,:,3] < 1e9) &
#                                           (store[sliceIndex,:,12] == 99) & 
#                                           (store[sliceIndex,:,13] == 0) & 
#                                           (store[sliceIndex,:,14] == 0)&
#                                           (ir_coadd_data[:,2] == 1) & 
#                                           (ir_coadd_data[:,3]>500) & 
#                                           (ir_coadd_data[:,3]<2000))[0]
#     for index in loc:
#         index_loc = np.where(loc_coadd == index)[0]
#         #print (index_loc)
#         temp[sliceIndex, index_loc, :] = store[sliceIndex,index,:]
# 
# np.save('/scratch/bell/dutta26/abell_2390/stars_sf/bright_stars.npy', temp)      
# 
# 
# 
# 
# 
# 
# 
# 
# ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
# loc_coadd = np.where((ir_coadd_data[:,2] == 1) & (ir_coadd_data[:,3]>100) & (ir_coadd_data[:,3]<500))[0]
# 
# temp = np.zeros((100, len(loc_coadd), 77))
# band = 'r'
# for sliceIndex in range(1,100):
#     #sliceIndex = 1
#     sizex =sizey = 100
#     tot = np.zeros((sizex, sizey))
#     
#     store = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')
#     
#     star_arr = store[sliceIndex,(np.where((store[sliceIndex,:,2] == 1) & (store[sliceIndex,:,12] == 99) 
#                                & (store[sliceIndex,:,13] == 0) & (store[sliceIndex,:,14] == 0)))[0],  : ]
#     
#     
#     #Now tuse k sigma clip to find usable stars. Just do for sigxx
#     mean,median, std = sigma_clipped_stats(star_arr[:, 7])
#     mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
#     mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
#     print (mean,median, std)
#     loc = np.where((store[sliceIndex,:,2] == 1) & 
#                                           (store[sliceIndex,:,7] >= mean-3*std) &
#                                           (store[sliceIndex,:,7] <= mean+3*std) &
#                                           (store[sliceIndex,:,8] >= mean1-3*std1) &
#                                           (store[sliceIndex,:,8] <= mean1+3*std1) &
#                                           (store[sliceIndex,:,3] < 1e9) &
#                                           (store[sliceIndex,:,12] == 99) & 
#                                           (store[sliceIndex,:,13] == 0) & 
#                                           (store[sliceIndex,:,14] == 0)&
#                                           (ir_coadd_data[:,2] == 1) & 
#                                           (ir_coadd_data[:,3]>100) & 
#                                           (ir_coadd_data[:,3]<500))[0]
#     for index in loc:
#         index_loc = np.where(loc_coadd == index)[0]
#         #print (index_loc)
#         temp[sliceIndex, index_loc, :] = store[sliceIndex,index,:]
# 
# np.save('/scratch/bell/dutta26/abell_2390/stars_sf/medium_stars.npy', temp)      
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
# loc_coadd = np.where((ir_coadd_data[:,2] == 1) & (ir_coadd_data[:,3]>50) & (ir_coadd_data[:,3]<100))[0]
# 
# temp = np.zeros((100, len(loc_coadd), 77))
# band = 'r'
# for sliceIndex in range(1,100):
#     #sliceIndex = 1
#     sizex =sizey = 100
#     tot = np.zeros((sizex, sizey))
#     
#     store = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')
#     
#     star_arr = store[sliceIndex,(np.where((store[sliceIndex,:,2] == 1) & (store[sliceIndex,:,12] == 99) 
#                                & (store[sliceIndex,:,13] == 0) & (store[sliceIndex,:,14] == 0)))[0],  : ]
#     
#     
#     #Now tuse k sigma clip to find usable stars. Just do for sigxx
#     mean,median, std = sigma_clipped_stats(star_arr[:, 7])
#     mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
#     mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
#     print (mean,median, std)
#     loc = np.where((store[sliceIndex,:,2] == 1) & 
#                                           (store[sliceIndex,:,7] >= mean-3*std) &
#                                           (store[sliceIndex,:,7] <= mean+3*std) &
#                                           (store[sliceIndex,:,8] >= mean1-3*std1) &
#                                           (store[sliceIndex,:,8] <= mean1+3*std1) &
#                                           (store[sliceIndex,:,3] < 1e9) &
#                                           (store[sliceIndex,:,12] == 99) & 
#                                           (store[sliceIndex,:,13] == 0) & 
#                                           (store[sliceIndex,:,14] == 0)&
#                                           (ir_coadd_data[:,2] == 1) & 
#                                           (ir_coadd_data[:,3]>50) & 
#                                           (ir_coadd_data[:,3]<100))[0]
#     for index in loc:
#         index_loc = np.where(loc_coadd == index)[0]
#         #print (index_loc)
#         temp[sliceIndex, index_loc, :] = store[sliceIndex,index,:]
# 
# np.save('/scratch/bell/dutta26/abell_2390/stars_sf/faint_stars.npy', temp)      
# 
# 
# 
# 
# 
# 
# 
# =============================================================================












fileArr =['bright', 'medium', 'faint']
colorArr= ['r', 'g', 'b']
for j in range(len(fileArr)):
    temp = np.load('/scratch/bell/dutta26/abell_2390/'+fileArr[j]+'_stars.npy')   
    a,b,c = np.shape(temp)
    plotArr =[]    
    for index in range(b):
# =============================================================================
#             fluxArr = temp[0:100,index,3]
#             fluxArr= fluxArr[fluxArr>0]
#             fluxArr = fluxArr/np.median(fluxArr)
#             for val in fluxArr:
#                 plotArr.append(val)
# =============================================================================
    
# =============================================================================
#             offsetArr = np.sqrt(temp[0:100,index,4]**2 + temp[0:100,index,5]**2 )
#             offsetArr = offsetArr[offsetArr!=0]
#             for val in offsetArr:
#                 plotArr.append(val)
# =============================================================================
                
# =============================================================================
#             sizeArr = np.sqrt(temp[0:100,index,7] + temp[0:100,index,8] )
#             sizeArr = sizeArr[sizeArr!=0]
#             for val in sizeArr:
#                 plotArr.append(val)
# =============================================================================

            e1 =(temp[0:100,index,7] - temp[0:100,index,8] )/(temp[0:100,index,7] + temp[0:100,index,8] )
            e2 = 2*temp[0:100,index,9]/(temp[0:100,index,7] + temp[0:100,index,8] )
            loc = np.where((e1!=0) & (np.abs(e1)<1) & (np.abs(e2)<1))[0]
            e1 = e1[loc]
            e2 = e2[loc]
            
            for val1, val2 in zip(e1,e2):
                plotArr.append(np.sqrt(val1**2+val2**2))
            
    n, bins, patches = plt.hist(x=plotArr, bins='auto', histtype=u'step', color=colorArr[j], label='yy', density=True)    