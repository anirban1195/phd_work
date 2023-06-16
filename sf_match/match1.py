#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 19:05:32 2023

@author: dutta26
"""



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


#Plot histograms of specific regions
ir_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_temp3r_temp.pk1'))
size = np.sqrt(ir_coadd_df[:,7] + ir_coadd_df[:,8] )
maxFlux = 10**3.5
minFlux = 10**2.5
loc = np.where((ir_coadd_df[:,3]>minFlux) & (ir_coadd_df[:,3]< maxFlux) &
               (ir_coadd_df[:,12]== 99) & (ir_coadd_df[:,13]== 0) & (ir_coadd_df[:,14]== 0)
              & (ir_coadd_df[:,7]< 20) & (ir_coadd_df[:,8]< 20))[0]

data = ir_coadd_df[loc,8]




counts,bin_edges = np.histogram(data, bins=90)
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



plt.xlabel('Sigmayy')
plt.ylabel('Frequency')
plt.title('Log flux in range '+str(np.log10(minFlux))[0:4]+' - '+str(np.log10(maxFlux))[0:4])
mean_psf_size = 3.6
fluxArr = (maxFlux+minFlux)/2
area = 3.14 * np.sqrt(6*6.86*4)
bkg = 8000  

fluxArr = fluxArr*2507
errSize = np.sqrt( (area/(np.pi*fluxArr) + (4* area**2 * bkg)/(np.pi * fluxArr**2) ) )
errSize1 = np.sqrt(( mean_psf_size* errSize*1.414)**2)
#errSize3 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.8**2) #for xx
errSize3 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.6**2) #for yy

# =============================================================================
# #For size only
# errSize1 = np.sqrt((errSize)**2 )
# errSize3 = np.sqrt((errSize)**2 + 0.14**2)
# =============================================================================
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

center = 6.86
ht = 140
plt.plot(x_values, gaussian(x_values,center, errSize1, ht), 'r-',linewidth = 2)
plt.plot(x_values, gaussian(x_values,center, errSize3, ht), 'g-',linewidth = 2)

xVal =[]
yVal =[]

tot = 0
cutoff = 6.86 + errSize3
#center = 13.6
for j in range(len(bin_centres)):
    if(bin_centres[j]> cutoff):
        continue
    xVal.append(bin_centres[j])
    yVal.append(counts[j])
    
def gaussian1(x, mu, A):
    
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(errSize3, 2.)))

#xVal = np.arange(1,20, 0.01)
#yVal = gaussian1(xVal,10, 100) + np.random.normal(0,1,len(xVal))
parameters, covariance = curve_fit(gaussian1, xVal, yVal, [6.86,100])    
    
