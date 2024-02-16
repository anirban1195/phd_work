#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 09:19:48 2023

@author: dutta26
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_forcedDist1_err1001.npy')
# =============================================================================
# lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_1002.npy')
# lut3 = np.load('/scratch/bell/dutta26/abell_2390/calib_1003.npy')
# =============================================================================
final = np.load('/home/dutta26/codes/forced_truncated_calib.npy')
for j in range(len(lut1)):
    flux = lut1[j,0]
    sqrtba= lut1[j,1]
    f = lut1[j,2]
    g = lut1[j,3]
    h = lut1[j,8]
    #plt.plot(np.log10(flux/sqrtba), f, 'b.', markersize = 1)
    #plt.plot(np.log10(flux/sqrtba), g, 'b.', markersize = 1)
    plt.plot(np.log10(flux/sqrtba), h, 'b.', markersize = 1)
    
# =============================================================================
# for j in range(len(lut2)):
#     flux = lut2[j,0]
#     sqrtba= lut2[j,1]
#     f = lut2[j,2]
#     g = lut2[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.',markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), g, 'b.', markersize = 1)
# 
# for j in range(len(lut3)):
#     flux = lut3[j,0]
#     sqrtba= lut3[j,1]
#     f = lut3[j,2]
#     g = lut3[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.',markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), g, 'b.', markersize = 1)
# =============================================================================

#plt.ylabel('f/sqrt(b)*A')
#plt.xlabel('Log (Flux/sqrt(b)*A)')    
#plt.ylabel(r'$ \frac{f}{A.\sqrt{B}}$ ')
#plt.ylabel('Ratio')v
#plt.title('Centroid x Bias')
plt.plot([0], [0], 'r.', markersize = 0)
plt.xlabel(r'$Log_{10}(\frac{N}{A.\sqrt{B}})$')   
plt.ylabel(r'$p_{\mu_{x}}$') 
plt.plot(final[:,0], final[:,3]/2, 'r-', markersize= 10, label = 'Fitted Function')
plt.plot(np.arange(-3.5, -2.6, 0.1), np.ones(9)*final[0,3]/2, 'r-', markersize= 10)
plt.plot(np.arange(2, 2.9, 0.1), np.ones(9)*final[48,3]/2, 'r-', markersize= 10)

plt.legend()
# =============================================================================
# xArr=[]
# yArr=[]
# for j in np.arange(-10, 10, 0.1):
#     xArr.append(-j+0.5)
#     yArr.append((erf(j)+1)*3.18/2)
# plt.plot(xArr, yArr, 'b--')
# 
# =============================================================================




# =============================================================================
# lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_1001_inf.npy')
# lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_1002_inf.npy')
# lut3 = np.load('/scratch/bell/dutta26/abell_2390/calib_1003_inf.npy')
# final = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
# 
# for j in range(len(lut1)):
#     flux = lut1[j,0]
#     sqrtba= lut1[j,1]
#     f = lut1[j,2]
#     g = lut1[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), g, 'b.')
#     
# for j in range(len(lut2)):
#     flux = lut2[j,0]
#     sqrtba= lut2[j,1]
#     f = lut2[j,2]
#     g = lut2[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), g, 'b.')
#     
# for j in range(len(lut3)):
#     flux = lut3[j,0]
#     sqrtba= lut3[j,1]
#     f = lut3[j,2]
#     g = lut3[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), g, 'b.')
#     
# plt.plot(final[:,0][:-6], final[:,2][:-6], 'r-', markersize= 10, label = 'Fitted Function')
# plt.ylabel('g')
# #plt.ylabel(r'$ \frac{f}{A.\sqrt{B}}$ ')
# plt.xlabel(r'$Log(\frac{N}{A.\sqrt{B}})$')    
# plt.legend()
# =============================================================================

