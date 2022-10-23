#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 22 21:00:58 2022

@author: dutta26
"""
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess
from scipy.signal import convolve2d as conv2
import scipy.stats as stats


# =============================================================================
# 
# low = [71, 1130, 3129, 3714, 4340, 4853, 4915, 6378, 6484, 9182, 11466, 
#        11756, 12283, 13170, 13424, 13445, 13948, 15470, 16716, 17423, 
#        18744, 18797, 21231, 21529, 21663, 23015, 23488, 23876, 24338,
#        24962, 26647, 27633, 30153, 30508, 30699, 32003, 32224, 32355,
#        35333, 35751, 36680, 37094, 38237, 39252, 39304, 39355, 40729, 
#        40968, 41355, 41484, 42878, 43335, 43965, 44003, 45147, 45180,
#        45296, 45736, 46020, 46652, 46653, 46735, 47147, 48130, 49555, 
#        51223, 51866, 51893, 52763]
# 
# high = [816, 1114, 1715, 2532, 2788, 3270, 4693, 5023, 6112, 6713, 6883, 
#         6941, 8422, 9534, 9544, 9724, 10149, 11820, 13548, 13753, 15143, 
#         15626, 15665, 16062, 17092, 18289, 19246, 19495, 20774, 22774, 25733,
#         26646, 26818, 26954, 27351, 27998, 29052, 29111, 30071, 31135, 
#         31467, 32428, 32945, 33696, 34379, 34442, 35245, 35688, 36077, 
#         37757, 38455, 38666, 39839, 39910, 40130, 40345, 40509, 42285, 
#         42590, 42677, 42710, 42790, 43094, 43166, 43703, 44007, 44325, 
#         45503, 45755, 45891, 46264, 46452, 46573, 46922, 47159, 47778, 
#         47920, 48334, 48590, 49177, 49966, 50630, 51196, 51383]
# =============================================================================


band = 'i'
frame_data = np.load('/scratch/halstead/d/dutta26/abell_2390/test4_'+band+ '.npy')
#frame_data = np.load('/home/dutta26/codes/singleFrame_'+band+'.npy')
coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
coadd_data = np.array(coadd_data)
coadd_flux = coadd_data[:,3]
a,b,c = np.shape(frame_data)

arr = []
mjd = []
airmass =[]



maxFlux = 20000
minFlux = 0
valid = np.where((coadd_data[:,2] == 1) & (coadd_data[:,3] >= minFlux) & (coadd_data[:,3] <= maxFlux) )
arr_high =[]
for lowIndex in valid[0]:
    temp = frame_data[:,lowIndex,3]
    temp = temp[temp != 0]
    temp = temp [temp > 10 ]
    temp = temp[temp < 1e9]
    mean, med, std = sigma_clipped_stats(temp)
    temp = temp [temp > (med-2*std) ]
    temp = temp[temp < (med+2*std)]
    mean, med, std = sigma_clipped_stats(temp)
    mean, med, std = sigma_clipped_stats(temp)
    flux_dev = (temp)/med
    for j in range(len(temp)):
        if(flux_dev[j]> 2):
            continue
        arr_high.append(flux_dev[j])
# =============================================================================
#     valid_temp = np.where(( frame_data[:,lowIndex,3]> (med-5*std) ) & (frame_data[:,lowIndex,3]< (med+5*std)))
#     for j in valid_temp[0]:
#         #print (j)
#         arr_high.append(np.sqrt(frame_data[j,lowIndex,6]**2 + frame_data[j,lowIndex,7]**2))
# 
# =============================================================================

print (len(arr_high))        
arr_high=np.array(arr_high)   
arr_high= arr_high[arr_high>0]
arr_high= arr_high[arr_high<5]    
           
n, bins, patches = plt.hist(x=arr_high, color = 'b' , label = '2000-1000',
                             alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=100)



# =============================================================================
# maxFlux = 100
# minFlux = 50
# valid = np.where((coadd_data[:,2] == 1) & (coadd_data[:,3] >= minFlux) & (coadd_data[:,3] <= maxFlux) )
# arr_low =[]
# for lowIndex in valid[0]:
#     temp = frame_data[:,lowIndex,3]
#     #print ('aaaaaaaa')
#     temp = temp[temp != 0]
#     temp = temp [temp > 10 ]
#     temp = temp[temp < 1e8]
#     mean, med, std = sigma_clipped_stats(temp)
#     temp = temp [temp > (med-5*std) ]
#     temp = temp[temp < (med+5*std)]
#     
#     mean, med, std = sigma_clipped_stats(temp)
#     flux_dev = (temp)/med
#     for j in range(len(temp)):
#         if(flux_dev[j]> 2):
#             continue
#         arr_low.append(flux_dev[j])
# print (len(arr_low))        
# n, bins, patches = plt.hist(x=arr_low, color = 'r',  label = '50-100',
#                              alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=100)
# 
# plt.legend()
# plt.xlabel('Flux/ Median Flux')
# =============================================================================




# =============================================================================
# mu = np.mean(arr)
# variance = 1000
# sigma = math.sqrt(variance)
# x = np.linspace(mu - 5*sigma, mu + 5*sigma, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma))
# plt.show()
# =============================================================================
# =============================================================================
# for j in range (a):
#     flux_dev =  (frame_data[j,low,3]- coadd_data[low,3]) / np.sqrt(coadd_data[low,3])
#     valid = np.where((flux_dev>-40) & (flux_dev<40) & (frame_data[j,low,2] == 1) & 
#                      (frame_data[j,low,15] == 0) &(frame_data[j,low,16] == 0) &  (frame_data[j,low,3] > 0) &
#                      (frame_data[j,low,17] == 0) & (coadd_data[low,3] > minFlux) & (coadd_data[low,3] < maxFlux))
#     print (len(valid[0]))
#     for index in valid[0]:
#         arr.append(flux_dev[index])
# n, bins, patches = plt.hist(x=arr,  
#                             alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=100)
# plt.grid(axis='y', alpha=0.75)
# =============================================================================
# =============================================================================
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as stats
# import math
# 
# mu = np.mean(arr)
# variance = 1
# sigma = math.sqrt(variance)
# x = np.linspace(mu - 5*sigma, mu + 5*sigma, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma))
# plt.show()
# 
# 
# for j in range (a):
#     flux_dev =  (frame_data[j,high,3]- coadd_data[high,3]) / np.sqrt(coadd_data[high,3])
#     valid = np.where((flux_dev>-40) & (flux_dev<40) & (frame_data[j,high,2] == 1) & 
#                      (frame_data[j,high,15] == 0) &(frame_data[j,high,16] == 0) &  (frame_data[j,high,3] > 0) &
#                      (frame_data[j,high,17] == 0) & (coadd_data[high,3] > minFlux) & (coadd_data[high,3] < maxFlux))
#     print (len(valid[0]))
#     for index in valid[0]:
#         arr.append(flux_dev[index])
# n, bins, patches = plt.hist(x=arr,  
#                             alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=100)
# plt.grid(axis='y', alpha=0.75)
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as stats
# import math
# 
# mu = np.mean(arr)
# variance = 1
# sigma = math.sqrt(variance)
# x = np.linspace(mu - 5*sigma, mu + 5*sigma, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma))
# plt.show()
# =============================================================================
