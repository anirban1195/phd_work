#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 10:40:53 2023

@author: dutta26
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 

# =============================================================================
# ir_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1'))
# size = ir_coadd_df[:,8] - ir_coadd_df[:,39] + ir_coadd_df[:,7] - ir_coadd_df[:,38]
# 
# loc = np.where((ir_coadd_df[:,3]>0) & (ir_coadd_df[:,3]< 50000))[0]
# plt.plot(size[loc], np.log10(ir_coadd_df[loc,3]), 'b.', markersize = 2)
# 
# mean_psf_size = 3
# fluxArr = np.arange(0.001, 5000, 0.01)*30000
# area = 3.14 * np.sqrt(9*9)
# bkg = 70000
# errSize = np.sqrt( (area/(np.pi*fluxArr) + 4*area**2 * bkg/(np.pi * fluxArr**2) )  )
# errSize_sq = np.sqrt((2*mean_psf_size* errSize)**2 +0.13**2 +0.13**2)
# 
# 
# plt.plot(0-errSize_sq, np.log10(fluxArr/30000), 'r-')
# plt.plot(0+errSize_sq, np.log10(fluxArr/30000), 'r-', label= '1-Sigma')
# 
# plt.plot(0-2*errSize_sq, np.log10(fluxArr/30000), 'y-')
# plt.plot(0+2*errSize_sq, np.log10(fluxArr/30000), 'y-', label= '2-Sigma')
# 
# plt.plot(0-3*errSize_sq, np.log10(fluxArr/30000), 'k-')
# plt.plot(0+3*errSize_sq, np.log10(fluxArr/30000), 'k-', label= '3-Sigma')
# 
# plt.xlabel(r'Residual')
# plt.ylabel('Log Flux')
# plt.legend()
# plt.title('Theoretical and Interpolation')
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
ir_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1'))
size = np.sqrt(ir_coadd_df[:,7] + ir_coadd_df[:,8] )
maxFlux = 35
minFlux = 25
loc = np.where((ir_coadd_df[:,3]>minFlux) & (ir_coadd_df[:,3]< maxFlux) & (size >0) & (size<4.5))[0]
# =============================================================================
# n, bins, patches = plt.hist(x=size[loc], bins='auto', histtype='step', edgecolor='k', density = True,
#                             alpha=0.7, rwidth=1,linewidth = 2)
# plt.xlabel('Size')
# plt.ylabel('Frequency')
# =============================================================================
counts,bin_edges = np.histogram(size[loc], 20)
bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.0
err = []
for j in range(1, len(bin_edges)):
    bin_start = bin_edges[j-1]
    bin_end = bin_edges[j]
    loc1 = np.where((size[loc]>=bin_start)& (size[loc]<bin_end))[0]
    if(len(loc1) <=5):
        err.append(0)
    else:
        err.append(np.nanstd(size[loc][loc1]))
plt.errorbar(bin_centres, counts, xerr=err, fmt='o', markersize = 5)

plt.title('Log flux in range '+str(np.log10(minFlux))[0:4]+' - '+str(np.log10(maxFlux))[0:4])
mean_psf_size = 3
fluxArr = (maxFlux+minFlux)*15000
area = 3.14 * np.sqrt(9*9)
bkg = 70000

errSize = np.sqrt( (area/(np.pi*fluxArr) + 4*area**2 * bkg/(np.pi * fluxArr**2) ) )
#errSize = 0.2811
# =============================================================================
# errSize1 = np.sqrt(( mean_psf_size* errSize*1.414)**2)
# errSize2 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.13**2 )
# errSize3 = np.sqrt(( mean_psf_size* errSize*1.414)**2 + 0.13**2 + 0.13**2 )
# =============================================================================
errSize1 = np.sqrt((errSize)**2 )
errSize2 = np.sqrt((errSize)**2 + 0.031**2)
errSize3 = np.sqrt((errSize)**2 + 0.043**2)
# =============================================================================
# errSize1 = np.sqrt(( mean_psf_size* errSize/1.414)**2)
# errSize2 = np.sqrt(( mean_psf_size* errSize/1.414)**2 + 0.05**2 )
# errSize3 = np.sqrt(( mean_psf_size* errSize/1.414)**2 + 0.05**2 + 0.05**2 )
# =============================================================================

def gaussian(x, mu, sig, A):
    return (A)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

x_values = np.linspace(1, 6, 5000)
#temp= np.sort( n )
#loc1 = np.where( n == temp[len(temp)-1])[0]
#loc2 = np.where( n == temp[len(temp)-2])[0]

xShift = 0.0004
#plt.plot(x_values, gaussian(x_values,bins[loc1]+0.02, errSize1, np.max(n)), 'r-',linewidth = 2)
#plt.plot(x_values, gaussian(x_values,bins[loc1]+xShift, errSize, np.max(n)), 'r-',linewidth = 2)
#plt.plot(x_values, gaussian(x_values,bins[loc1]+xShift, errSize3, np.max(n)), 'g-',linewidth = 2)

plt.plot(x_values, gaussian(x_values,3, errSize, np.max(counts)), 'r-',linewidth = 2)
plt.plot(x_values, gaussian(x_values,3, errSize3, np.max(counts)), 'g-',linewidth = 2)

