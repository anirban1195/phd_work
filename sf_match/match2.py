#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 09:16:05 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
from astropy.stats import sigma_clipped_stats
from scipy.optimize import curve_fit
import os
def gaussian(x, mu, sig, A):
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


#Read the npy file
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
pArr =[]
sArr =[]
kArr =[]
eArr= []
yArr =[]
bArr=[]
fArr=[]
for slice_no in np.arange(1, 100, 1):
#for slice_no in [32]:
    #slice_no = 8
    
    
    #Find star loc 
    star_arr = r_sf_df[slice_no,(np.where((r_sf_df[slice_no,:,2] == 1) & (r_sf_df[slice_no,:,12] == 99) 
                               & (r_sf_df[slice_no,:,13] == 0) & (r_sf_df[slice_no,:,14] == 0)))[0],  : ]
    
    
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    #print (mean,median, std)
    
    star_arr = r_sf_df[slice_no, (np.where((r_sf_df[slice_no,:,2] == 1) & 
                                          (r_sf_df[slice_no,:,7] >= mean-3*std) &
                                          (r_sf_df[slice_no,:,7] <= mean+3*std) &
                                          (r_sf_df[slice_no,:,8] >= mean1-3*std1) &
                                          (r_sf_df[slice_no,:,8] <= mean1+3*std1) &
                                          (r_sf_df[slice_no,:,3] < 1e9) &
                                          (r_sf_df[slice_no,:,12] == 99) & 
                                          (r_sf_df[slice_no,:,13] == 0) & 
                                          (r_sf_df[slice_no,:,14] == 0)))[0],  : ]
    
    
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    if(median>30 or median1>30):
        #print ('aa')
        continue
    
    #Read scale and adjust
    scale_def = 0.0064953
    #find scale of current frame 
    temp = r_sf_df[slice_no,:,30]
    temp = temp[temp>0]
    scale = np.median(temp)
    
    
    size = np.sqrt(r_sf_df[slice_no,:,7] + r_sf_df[slice_no,:,8] )
    maxFlux = 10**6 * (scale_def/scale)
    minFlux = 10**4.8 * (scale_def/scale)
    loc = np.where((r_sf_df[slice_no,:,3]>minFlux) & (r_sf_df[slice_no,:,3]< maxFlux) &
                   (r_sf_df[slice_no,:,12]== 99) & (r_sf_df[slice_no,:,13]== 0) & (r_sf_df[slice_no,:,14]== 0)
                  & (r_sf_df[slice_no,:,7]< mean+6*std) & (r_sf_df[slice_no,:,8]< mean1+6*std1) & 
                  (r_sf_df[slice_no,:,7]> mean-6*std) & (r_sf_df[slice_no,:,8]> mean1-6*std1) & 
                  (r_sf_df[slice_no,:,2] == 1))[0]
    
    if(len(loc)< 100):
        #print ('bb')
        continue
    #data = r_sf_df[slice_no,loc,9]
    data = np.sqrt(r_sf_df[slice_no,loc,7] +r_sf_df[slice_no,loc,8])
    
    counts,bin_edges = np.histogram(data, bins=40)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.0
    err = []
    for j in range(1, len(bin_edges)):
        bin_start = bin_edges[j-1]
        bin_end = bin_edges[j]
        loc1 = np.where((data>=bin_start)& (data<bin_end))[0]
        if(len(loc1) <=5):
            err.append(0)
        else:
            err.append(np.nanstd(data[loc1]))
    #plt.errorbar(bin_centres, counts, xerr=None, yerr=np.sqrt(counts), fmt='o', markersize = 5)
    #plt.xlabel('Size')
    #plt.ylabel('Frequency')
    #plt.title('Log flux in range '+str(np.log10(minFlux))[0:4]+' - '+str(np.log10(maxFlux))[0:4])
    if(np.median( star_arr[:,38]) == -99):
       continue
            
    
    mean_xx = np.median( star_arr[:,7])
    mean_yy = np.median( star_arr[:,8])
    mean_xy = np.median( star_arr[:,9])
    mean_psf_size = np.median( np.sqrt( star_arr[:,7] + star_arr[:,8]) )
    
    e1 = (mean_xx - mean_yy)/(mean_xx + mean_yy)
    e2 = 2*mean_xy/(mean_xx + mean_yy)
    eArr.append(np.sqrt(e1**2 + e2**2))
    bArr.append(np.median( star_arr[:,6]))
    
    pArr.append(mean_psf_size)
    sArr.append(scale)
    fluxArr = (maxFlux+minFlux)/2
    area = 3.14 * np.sqrt( 2 * np.median(star_arr[:,7]) *2 * np.median(star_arr[:,8]))
    bkg = np.median(star_arr[:,15])
    
    err_size = np.sqrt( (area/(np.pi*fluxArr) + (4* area**2 * bkg)/(np.pi * fluxArr**2) ) )
    err_xx = err_yy = np.sqrt(( mean_psf_size* err_size)**2)
    err_xy = err_xx * 0.707
    k_guess_xx = np.median(star_arr[:,41])
    k_guess_yy = np.median(star_arr[:,42])
    k_guess_xy = np.median(star_arr[:,43])
    k_guess_size = (k_guess_xx+k_guess_yy)/(mean_psf_size * 1.414*2)
    
    xVal =[]
    yVal =[]
    tot = 0
    start = mean_psf_size - 10*np.sqrt(k_guess_size**2 + err_size**2)
    cutoff = mean_psf_size + 10*np.sqrt(k_guess_size**2 + err_size**2)
    #center = 13.6
    for j in range(len(bin_centres)):
        if(bin_centres[j]> cutoff or bin_centres[j]< start):
            continue
        xVal.append(bin_centres[j])
        yVal.append(counts[j])
    
    
    parameters, covariance = curve_fit(gaussian, xVal, yVal, [mean_psf_size, np.sqrt(k_guess_size**2 + err_size**2) , 80])    
    
    x_values = np.linspace(mean_psf_size - 10*np.sqrt(k_guess_size**2 + err_size**2),mean_psf_size + 10*np.sqrt(k_guess_size**2 + err_size**2), 25000)
    #plt.plot(x_values, gaussian(x_values,parameters[0], parameters[1], parameters[2]), 'r-',linewidth = 2)
    #print (parameters)
    kArr.append(np.sqrt(parameters[1]**2 - err_size**2))
    #if(np.sqrt(parameters[1]**2 - err_xx**2) > (0.2 *mean_psf_size )):
    #    print (slice_no)
    
    f_no = str(np.median(star_arr[:,44]))
    fArr.append(f_no)
    for files in os.listdir('/scratch/bell/dutta26/abell_2390/i/temp/'):
        if('weight' in files):
            continue
        if(f_no in files):
            if('2017' in files):
                yArr.append(2017)
            elif('2018' in files):
                yArr.append(2018)
            elif('2019' in files):
                yArr.append(2019)
            elif('2020' in files):
                yArr.append(2020)
            elif('2021' in files):
                yArr.append(2021)
    
    
eArr = np.array(eArr)
pArr = np.array(pArr)
kArr = np.array(kArr)
sArr = np.array(sArr)
yArr = np.array(yArr)
bArr =np.array(bArr)
fArr=np.array(fArr)
plt.plot(np.arange(1,20), np.arange(1,20)/30, 'b-')
plt.plot(pArr, kArr, 'r.')
plt.xlabel('Mean PSF')
plt.ylabel('Best Fit k (Size)')
# =============================================================================
# for j in range(len(eArr)):
#     if(yArr[j] == 2017):
#         plt.plot(pArr[j], kArr[j], 'r.')
#     if(yArr[j] == 2018):
#         plt.plot(pArr[j], kArr[j], 'b.')
#     if(yArr[j] == 2019):
#         plt.plot(pArr[j], kArr[j], 'k.')
#     if(yArr[j] == 2020):
#         plt.plot(pArr[j], kArr[j], 'r+')
#     if(yArr[j] == 2021):
#         plt.plot(pArr[j], kArr[j], 'b+')
# =============================================================================
