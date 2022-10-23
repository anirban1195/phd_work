#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 18:08:10 2021

@author: dutta26
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats



band = 'r'
frame_data = np.load('/scratch/halstead/d/dutta26/abell_2390/test3_r.npy')
#frame_data_unsl = np.load('/scratch/halstead/d/dutta26/abell_2390/test3_'+band +'.npy')
#frame_data = np.load('/scratch/halstead/d/dutta26/abell_2390/test2_r.npy')
coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
coadd_data = np.array(coadd_data)
coadd_flux = coadd_data[:,3]
a,b,c = np.shape(frame_data)
arr = []
arr1=[]
mjd = []
airmass =[]
bkg = []
std =[]
cnt = 0
maxFlux = 20000 
minFlux = 300



meanArr=[]
medianArr=[]
stdArr=[]
for j in [6]:
    flux_dev =  (frame_data[j,:,3]*frame_data[j,:, 30]- coadd_data[:,3]) / np.sqrt(coadd_data[:,3])
    valid = np.where((flux_dev>-15) & (flux_dev<15) & (frame_data[j,:,12] == 0) &
                      (frame_data[j,:,13] == 0) & (frame_data[j,:,14] == 0) & (frame_data[j,:,2] == 1)  
                     & (frame_data[j,:,3]* frame_data[j,:,30] > minFlux)
                     & (frame_data[j,:,3]* frame_data[j,:,30] < maxFlux))
    

    
    e1 = (frame_data[j,valid[0], 7] - frame_data[j,valid[0], 8])/(frame_data[j,valid[0], 7] + frame_data[j,valid[0], 8])
    e1_pred = (frame_data[j,valid[0], 38] - frame_data[j,valid[0], 39])/(frame_data[j,valid[0], 38] + frame_data[j,valid[0], 39])
    
    e2 = (2*frame_data[j,valid[0], 9])/(frame_data[j,valid[0], 7] + frame_data[j,valid[0], 8])
    e2_pred = (2*frame_data[j,valid[0], 40])/(frame_data[j,valid[0], 38] + frame_data[j,valid[0], 39])
    
# =============================================================================
#     delEllip = np.sqrt((e1-e1_pred)**2 + (e2-e2_pred)**2)
#     mean, median, std = sigma_clipped_stats(delEllip)   
#     meanArr.append(mean)
#     stdArr.append(std)
#     medianArr.append(median)
#     print (len(valid[0]))
# =============================================================================
    a = (e1-e1_pred)**2 
    mean1,median1,std1 = sigma_clipped_stats(a)
    
    b = (e2-e2_pred)**2 
    mean2,median2,std2 = sigma_clipped_stats(b)
    
    rms = np.sqrt(np.nanmedian((e1-e1_pred)**2 )+ np.nanmedian((e2-e2_pred)**2))
    #rms = np.sqrt(median1+ median2)
    stdArr.append(rms)
    for index in valid[0]:
        meanArr.append(frame_data[j,index, 8]-frame_data[j,index, 39] )  
    print (rms, len(valid[0]))
    
    n, bins, patches = plt.hist(x=a[a<0.04],  color = 'b', label = '2000-1000',
                            alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=100)
# =============================================================================
#     for index in valid[0]:
#         arr1.append(frame_data[j,index, 7] -frame_data[j,index, 38]) 
# =============================================================================
# =============================================================================
#     diff = e1-e1_pred    
#     mean, med, std = sigma_clipped_stats(diff)
#     print (mean, std)
#     plt.plot(mean,std,'r.')
# =============================================================================
# =============================================================================
# stdArr=np.array(stdArr)
# n, bins, patches = plt.hist(x=stdArr[stdArr<0.05],  color = 'b', label = '2000-1000',
#                             alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=20)
# 
# =============================================================================
# =============================================================================
# n, bins, patches = plt.hist(x=arr1,  color = 'r',label = '2000-1000',
#                             alpha=0.3, rwidth=0.95, density=True, stacked=True, bins=500)
# 
# 
# mean, med, std = sigma_clipped_stats(arr)
# print (mean, med, std)
# 
# mean, med, std = sigma_clipped_stats(arr1)
# print (mean, med, std)
# 
# =============================================================================
# =============================================================================
# plt.grid(axis='y', alpha=0.75)
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as stats
# import math
# 
# mu = np.mean(arr)
# variance = 1/((maxFlux+minFlux)/2.0)
# sigma = math.sqrt(variance)
# x = np.linspace(mu - 5*sigma, mu + 5*sigma, 100)
# plt.plot(x, stats.norm.pdf(x, mu, sigma))
# plt.show()
# plt.xlabel('Variance in units of Sqrt N')
# plt.title('Flux ='+str(minFlux) + ' - '+str(maxFlux) + ' Mean = '+str(np.mean(arr))[0:4] + ' Variance = '+str(np.std(arr))[0:4] )
# =============================================================================
# =============================================================================
# arr1= []
# for j in range(107):
#     for k in range(len(arr[j])):
#         arr1.append(arr[j][k])
# =============================================================================
#plt.plot(std,airmass, 'b.')