#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 21:37:39 2021

@author: dutta26
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab
from astropy.stats import sigma_clipped_stats

store = np.load('/home/dutta26/test18_r.npy')
store[np.isnan(store)] = 0
store_star = store[:,0:200, :]
store_gal = store[:,300:15000, :]
j = 1000
fluxArr=[]
e1_err=[]
e2_err=[]

for j in np.arange(100):
    
    temp_sigxx = store_star[j,:, 3]
    temp_sigyy = store_star[j,:, 4]
    temp_sigxy = store_star[j,:, 6]
    arr = np.logical_and(temp_sigxx >0 , temp_sigyy >0)
    
    #print (temp_sigxy[arr], np.mean(temp_sigxy[arr]))
    print ('************************')
    print (j, np.mean(temp_sigxx[arr]), np.mean(temp_sigyy[arr]), np.std(temp_sigxx[arr]),np.std(temp_sigyy[arr]) )
    store_gal[j,:,3 ] = store_gal[j,:,3 ] - np.mean(temp_sigxx[arr])
    store_gal[j,:,4 ] = store_gal[j,:,4 ] - np.mean(temp_sigyy[arr])
    store_gal[j,:,6 ] = store_gal[j,:,6 ] - np.mean(temp_sigxy[arr])
    
    print (j, np.mean(temp_sigxx[arr]), np.mean(temp_sigyy[arr]), np.std(temp_sigxx[arr]),np.std(temp_sigyy[arr]) )
    print ('************************')
    store_star[j,arr,3 ] = store_star[j,arr,3 ] - np.mean(temp_sigxx[arr])
    store_star[j,arr,4 ] = store_star[j,arr,4 ] - np.mean(temp_sigyy[arr])
    store_star[j,arr,6 ] = store_star[j,arr,6 ] - np.mean(temp_sigxy[arr])
    
for j in range(4000):
    flux = store_gal[:,j, 0]
    mux = store_gal[:,j, 1]
    muy = store_gal[:,j, 2]
    sigxx = store_gal[:,j, 3]
    sigyy = store_gal[:,j, 4]
    sigxy = store_gal[:,j, 6]
    zp = store_gal[:,j, 10]
    seeing = store_gal[:,j, 9]
    bkg = store_gal[:,j, 8]
    
    arr= np.logical_and(flux>0, np.abs(flux) <100000)
    flux = flux[arr]
    mux = mux[arr]
    muy = muy[arr]
    sigxx = sigxx[arr]
    sigyy = sigyy[arr]
    sigxy = sigxy[arr]
    zp = zp[arr]
    seeing = seeing[arr]
    bkg = bkg[arr]
    
    
    
    wt = 100*np.power(10,(zp-25)/2.5)/((seeing)**2 *bkg)
    
    
    e1 = (sigxx-sigyy)/ (sigxx+sigyy)
    e2 = (2*sigxy)/ (sigxx+sigyy)
    
# =============================================================================
#     flux1 = flux.flatten()
#     mux1 = mux.flatten()
#     muy1 = muy.flatten()
# =============================================================================
    
# =============================================================================
#     arr= np.logical_and(flux>0, np.abs(flux) <100000)
#     flux = flux[arr]
#     mux = mux[arr]
#     muy = muy[arr]
#     e1=e1[arr]
#     e2= e2[arr]
# =============================================================================
    mean, median, stddev = sigma_clipped_stats(e1[1:], cenfunc=np.median)
    e1_err.append( np.median(e1[1:]) - e1[0])
    mean, median, stddev = sigma_clipped_stats(e2[1:], cenfunc=np.median)
    e2_err.append( np.median(e2[1:]) - e2[0])
    print (j,mean, e2[0], flux[0], stddev , len(e1), len(e2))
    fluxArr.append(flux[0])

plt.plot(fluxArr, e1_err, 'b.', markersize = 1)
plt.plot(fluxArr, e2_err, 'r.', markersize = 1)
plt.xlabel('Flux')
plt.ylabel('Ellipticity Error')
fluxArr =np.array(fluxArr)    
e1_err = np.array(e1_err)    
e2_err = np.array(e2_err)    
