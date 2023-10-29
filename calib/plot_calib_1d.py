#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:02:45 2023

@author: dutta26
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

# =============================================================================
# lut1 = np.load('/home/dutta26/codes/calib/calib_1001.npy')
# lut2 = np.load('/home/dutta26/codes/calib/calib_1002.npy')
# lut3 = np.load('/home/dutta26/codes/calib/calib_1003.npy')
# final = np.load('/home/dutta26/codes/calib/calib_final.npy')
# for j in range(len(lut1)):
#     flux = lut1[j,0]
#     sqrtba= lut1[j,1]
#     f = lut1[j,2]
#     g = lut1[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
#     
# for j in range(len(lut2)):
#     flux = lut2[j,0]
#     sqrtba= lut2[j,1]
#     f = lut2[j,2]
#     g = lut2[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.',markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
# 
# for j in range(len(lut3)):
#     flux = lut3[j,0]
#     sqrtba= lut3[j,1]
#     f = lut3[j,2]
#     g = lut3[j,3]
#     #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.',markersize = 1)
#     plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
# 
# #plt.ylabel('f/sqrt(b)*A')
# #plt.xlabel('Log (Flux/sqrt(b)*A)')    
# #plt.ylabel(r'$ \frac{f}{A.\sqrt{B}}$ ')
# plt.ylabel('g')
# plt.xlabel(r'$Log(\frac{N}{A.\sqrt{B}})$')    
# plt.plot(final[:,0], 2*final[:,2], 'r-', markersize= 10, label = 'Fitted Function')
# plt.legend()
# =============================================================================
# =============================================================================
# xArr=[]
# yArr=[]
# for j in np.arange(-10, 10, 0.1):
#     xArr.append(-j+0.5)
#     yArr.append((erf(j)+1)*3.18/2)
# plt.plot(xArr, yArr, 'b--')
# =============================================================================





lut1 = np.load('/home/dutta26/codes/calib/calib_1001_inf.npy')
lut2 = np.load('/home/dutta26/codes/calib/calib_1002_inf.npy')
lut3 = np.load('/home/dutta26/codes/calib/calib_1003_inf.npy')
final = np.load('/home/dutta26/codes/calib/calib_final_inf.npy')

for j in range(len(lut1)):
    flux = lut1[j,0]
    sqrtba= lut1[j,1]
    f = lut1[j,2]
    g = lut1[j,3]
    #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
    plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
    
for j in range(len(lut2)):
    flux = lut2[j,0]
    sqrtba= lut2[j,1]
    f = lut2[j,2]
    g = lut2[j,3]
    #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
    plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
    
for j in range(len(lut3)):
    flux = lut3[j,0]
    sqrtba= lut3[j,1]
    f = lut3[j,2]
    g = lut3[j,3]
    #plt.plot(np.log10(flux/sqrtba), f/sqrtba, 'b.', markersize = 1)
    plt.plot(np.log10(flux/sqrtba), 2*g, 'b.', markersize = 1)
    
plt.plot(final[:,0][:-6], 2*final[:,2][:-6], 'r-', markersize= 10, label = 'Fitted Function')
plt.ylabel('g')
#plt.ylabel(r'$ \frac{f}{A.\sqrt{B}}$ ')
plt.xlabel(r'$Log(\frac{N}{A.\sqrt{B}})$')    
plt.legend()

