#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:03:37 2023

@author: dutta26
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
# =============================================================================
# r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
# size = np.sqrt(r_sf_df[1,:,7] + r_sf_df[1,:,8] )
# 
# loc = np.where((r_sf_df[1,:,12]== 99) & (r_sf_df[1,:,13]== 0) & (r_sf_df[1,:,14]== 0))[0]
# plt.plot(size[loc], np.log10(r_sf_df[1,loc,3]), 'b.', markersize = 2)
# 
# mean_psf_size = 5
# fluxArr = np.arange(100, 9000000, 100)
# area = 3.14 * np.sqrt(22*27)
# bkg = 400
# errSize = np.sqrt( (area/(np.pi*fluxArr) + 4*area**2 * bkg/(np.pi * fluxArr**2) )  + 0.043**2 )
# #errSize_sq = np.sqrt((2*mean_psf_size* errSize)**2 +0.13**2 +0.13**2)
# 
# 
# plt.plot(5-errSize, np.log10(fluxArr), 'r-')
# plt.plot(5+errSize, np.log10(fluxArr), 'r-', label= '1-Sigma')
# 
# plt.plot(5-2*errSize, np.log10(fluxArr), 'y-')
# plt.plot(5+2*errSize, np.log10(fluxArr), 'y-', label= '2-Sigma')
# 
# plt.plot(5-3*errSize, np.log10(fluxArr), 'k-')
# plt.plot(5+3*errSize, np.log10(fluxArr), 'k-', label= '3-Sigma')
# 
# plt.xlabel(r'Size')
# plt.ylabel('Log Flux')
# plt.legend()
# plt.title('Theoretical and Interpolation')
# 
# =============================================================================
# =============================================================================
# size =ir_coadd_df[:,7] 
# loc = np.where((ir_coadd_df[:,3]>0) & (ir_coadd_df[:,3]< 50000))[0]
# plt.plot(size[loc], np.log10(ir_coadd_df[loc,3]), 'b.', markersize = 2)
# 
# mean_psf_size = 3
# fluxArr = np.arange(0.001, 5000, 0.01)*30000
# area = 3.14 * np.sqrt(9*9)
# bkg = 70000
# errSize = np.sqrt( (area/(np.pi*fluxArr) + 4*area**2 * bkg/(np.pi * fluxArr**2)))
# err_sigxx = errSize * np.sqrt(2*4.5 + 2*4.5)
# err_sigxx = np.sqrt(err_sigxx ** 2 + 0.2**2)
# 
# plt.plot(4.4-err_sigxx, np.log10(fluxArr/30000), 'r-')
# plt.plot(4.4+err_sigxx, np.log10(fluxArr/30000), 'r-', label= '1-Sigma')
# 
# plt.plot(4.4-2*err_sigxx, np.log10(fluxArr/30000), 'y-')
# plt.plot(4.4+2*err_sigxx, np.log10(fluxArr/30000), 'y-', label= '2-Sigma')
# 
# plt.plot(4.4-3*err_sigxx, np.log10(fluxArr/30000), 'k-')
# plt.plot(4.4+3*err_sigxx, np.log10(fluxArr/30000), 'k-', label= '3-Sigma')
# 
# plt.xlabel('Size')
# plt.ylabel('Log Flux')
# plt.legend()
# plt.title('Theoretical and Interpolation of sigmaxx')
# =============================================================================


#Plot histograms of specific regions
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
slice_no = 1
size = np.sqrt(r_sf_df[slice_no,:,7] + r_sf_df[slice_no,:,8] )
maxFlux = 10**4
minFlux = 10**3.8
loc = np.where((r_sf_df[slice_no,:,3]>minFlux) & (r_sf_df[slice_no,:,3]< maxFlux) &
               (r_sf_df[slice_no,:,12]== 99) & (r_sf_df[slice_no,:,13]== 0) & (r_sf_df[slice_no,:,14]== 0)
              & (r_sf_df[slice_no,:,7]< 30) & (r_sf_df[slice_no,:,8]< 30))[0]

data = size[loc]




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
plt.errorbar(bin_centres, counts, xerr=None, yerr=np.sqrt(counts), fmt='o', markersize = 5)



plt.xlabel('Size')
plt.ylabel('Frequency')
plt.title('Log flux in range '+str(np.log10(minFlux))[0:4]+' - '+str(np.log10(maxFlux))[0:4])
mean_psf_size = 5
fluxArr = (maxFlux+minFlux)/2
area = 3.14 * np.sqrt(22*27)
bkg = 400   

errSize = np.sqrt( (area/(np.pi*fluxArr) + (4* area**2 * bkg*1.2)/(np.pi * fluxArr**2) ) )
# =============================================================================
# errSize1 = np.sqrt(( mean_psf_size* errSize*1.414)**2)
# errSize3 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.27**2) #for xx
# #errSize3 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.39**2) #for yy
# 
# =============================================================================
#For size only
errSize1 = np.sqrt((errSize)**2 )
errSize3 = np.sqrt((errSize)**2 + 0.048**2)
# =============================================================================
##Fo sigmaxy
# errSize1 = np.sqrt(( mean_psf_size* errSize/1.414)**2)
# errSize2 = np.sqrt(( mean_psf_size* errSize/1.414)**2 + 0.05**2 )
# errSize3 = np.sqrt(( mean_psf_size* errSize/1.414)**2 + 0.05**2 + 0.05**2 )
# =============================================================================

def gaussian(x, mu, sig, A):
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x_values = np.linspace(0,25, 25000)
n=counts
temp= np.sort( n )
loc1 = np.where( n == temp[len(temp)-1])[0]
loc2 = np.where( n == temp[len(temp)-2])[0]

center = 4.94
ht = 32.1
plt.plot(x_values, gaussian(x_values,center, errSize1, ht), 'r-',linewidth = 2)
plt.plot(x_values, gaussian(x_values,center, errSize3, ht), 'g-',linewidth = 2)

xVal =[]
yVal =[]

tot = 0
cutoff = 4.94 + errSize3
#center = 13.6
for j in range(len(bin_centres)):
    if(bin_centres[j]> cutoff):
        continue
    xVal.append(bin_centres[j])
    yVal.append(counts[j])
    
def gaussian1(x, A):
    mu = 4.94
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(errSize3, 2.)))

#xVal = np.arange(1,20, 0.01)
#yVal = gaussian1(xVal,10, 100) + np.random.normal(0,1,len(xVal))
parameters, covariance = curve_fit(gaussian1, xVal, yVal, [ 50])    
    
