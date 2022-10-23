#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 16:12:52 2022

@author: dutta26
"""


import numpy as np 
from astropy.io import fits
import helper, helper1, measure_pythonV,sys
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
plt.style.use('seaborn-whitegrid')
# =============================================================================
# raArr =np.arange(237.0083-(5*0.00555), 237.0083+(5*0.00555), 0.00555 )
# decArr = np.arange(-15.2308 -(5*0.00555) ,  -15.2308 + (5*0.00555) , 0.005551)  
# 
# raList=[]
# decList=[]
# mag = 16
# magArr=[]
# for j in range(len(raArr)):
#     mag = mag+ 1
#     for k in range(len(decArr)):
#         raList.append(raArr[j])
#         decList.append(decArr[k])
#         magArr.append(mag)
# 
# zp = 31.18
# magArr = np.array(magArr)
# ideal_fluxArr = 10**((zp - magArr)/2.5)
# 
# =============================================================================
raList = np.load ('/scratch/halstead/d/dutta26/lsst/raArr_wiyn.npy')
decList = np.load ('/scratch/halstead/d/dutta26/lsst/decArr_wiyn.npy')   
magArr = np.load('/scratch/halstead/d/dutta26/lsst/magArr_wiyn.npy')     
fileName = '/scratch/halstead/d/dutta26/lsst/3_bi.fits'
f=fits.open(fileName)
data = np.array(f[0].data)
f.close()

# =============================================================================
# xList, yList = helper.convertToXY(raList, decList, fileName)
# fluxArr= np.zeros((10, 10))
# sigxxArr=np.zeros((10, 10))
# sigyyArr=np.zeros((10, 10))
# for j in range(len(xList)):
#     y=int(yList[j])
#     x=int(xList[j])
#     cut = data[y-15:y+15, x-15:x+15]
#     flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#     index = magArr[j] - 17
#     for j in range(10):
#         if(fluxArr[index, j] == 0):
#             index1 = j
#             break
#     
#     if(flux == 0 or flux==None or np.isnan(flux)):
#         fluxArr[index, index1] = None
#         sigxxArr[index, index1] = None
#         sigyyArr[index, index1] = None
#         continue
#     fluxArr[index, index1]=flux
#     sigxxArr[index, index1] = sigxx
#     sigyyArr[index, index1]= sigyy
# 
# sizeArr = np.sqrt(sigxxArr**2 + sigyyArr**2)   
# for j in range(10):
#     if(j != 0):
#         plt.errorbar(np.log10(np.nanmean(fluxArr[j,:])), np.nanmean(sizeArr[j,:]), fmt='+r')
#     else:
#         plt.errorbar(np.log10(np.nanmean(fluxArr[j,:])), np.nanmean(sizeArr[j,:]), fmt='+r', label = '30x30 cut')
#         
# sys.exit()        
# =============================================================================
        
########################################################################################################
xList, yList = helper.convertToXY(raList, decList, fileName)
fluxArr= np.zeros((10, 10))
sigxxArr=np.zeros((10, 10))
sigyyArr=np.zeros((10, 10))
for j in range(len(xList)):
    y=int(yList[j])
    x=int(xList[j])
    cut = data[y-25:y+25, x-25:x+25]
    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
    print (flux)
    index = magArr[j] - 17
    for j in range(10):
        if(fluxArr[index, j] == 0):
            index1 = j
            break
    
    if(flux == 0 or flux==None or np.isnan(flux)):
        fluxArr[index, index1] = None
        sigxxArr[index, index1] = None
        sigyyArr[index, index1] = None
        continue
    fluxArr[index, index1]=flux
    sigxxArr[index, index1] = sigxx
    sigyyArr[index, index1]= sigyy

sizeArr = np.sqrt(sigxxArr + sigyyArr)   
for j in range(10):
    if(j != 110):
        a = sigma_clipped_stats(fluxArr[j,:])[0]
        b = sigma_clipped_stats(sizeArr[j,:])[0]
        plt.errorbar(np.log10(a), b, fmt='+g', markersize = 8)
    else:
        plt.errorbar(np.log10(np.nanmean(fluxArr[j,:])), np.nanmean(sizeArr[j,:]), fmt='+k', label = 'Combined')
        
        
        
   
###########################################################################################################
# =============================================================================
# xList, yList = helper.convertToXY(raList, decList, fileName)
# fluxArr= np.zeros((10, 10))
# sigxxArr=np.zeros((10, 10))
# sigyyArr=np.zeros((10, 10))
# for j in range(len(xList)):
#     y=int(yList[j])
#     x=int(xList[j])
#     cut = data[y-25:y+25, x-25:x+25]
#     flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
#     print (flux, sigxx)
#     index = magArr[j] - 17
#     for j in range(10):
#         if(fluxArr[index, j] == 0):
#             index1 = j
#             break
#     
#     if(flux == 0 or flux==None or np.isnan(flux)):
#         fluxArr[index, index1] = None
#         sigxxArr[index, index1] = None
#         sigyyArr[index, index1] = None
#         continue
#     fluxArr[index, index1]=flux
#     sigxxArr[index, index1] = sigxx
#     sigyyArr[index, index1]= sigyy
# 
# sizeArr = np.sqrt(sigxxArr + sigyyArr)   
# for j in range(10):
#     if(j != 0):
#         a= sigma_clipped_stats(fluxArr[j,:])[0]
#         b = sigma_clipped_stats(sizeArr[j,:])[0]
#         plt.errorbar(np.log10(a), np.nanmean(b), fmt='+b')
#     else:
#         plt.errorbar(np.log10(np.nanmean(fluxArr[j,:])), np.nanmean(sizeArr[j,:]), fmt='+b', label = '50x50 cut')
#         
# plt.legend()
# =============================================================================
# =============================================================================
# plt.plot(np.log10(fluxArr), sigxxArr, 'g+', label = '50x50 Cutout')
# plt.plot(np.log10(ideal_fluxArr), sigxxArr, 'r.', label= 'Theoretical Flux')
# =============================================================================
import shutil
shutil.copy(fileName, '/scratch/halstead/d/dutta26/lsst/test.fits')
f=fits.open('/scratch/halstead/d/dutta26/lsst/test.fits', mode='update')
for j in range(len(xList)):
    y=int(yList[j])
    x=int(xList[j])
    f[0].data[y-20, x-20:x+20] = 1000
    f[0].data[y+20, x-20:x+20] = 1000
    f[0].data[y-20:y+20, x-20] = 1000
    f[0].data[y-20:y+20, x+20] = 1000
f.flush()
plt.legend()
plt.xlabel('Log Flux')
plt.ylabel('Size')
