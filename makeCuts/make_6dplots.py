#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 07:27:35 2021

@author: dutta26
"""


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab

store = np.load('/home/dutta26/test14_uncorrected_i.npy')
store[np.isnan(store)] = 0

objIndex = 0

flux = []
mux = []
muy =[]
e1 =[]
e2= []
bkg =[]
psf =[]
mjd =[]
wt =[]
back_h =[]
seeing = [] 
zp =[] 
fwhm =[] 
airmass =[] 
mphase = [] 
mAngle =[]
mjd =[]
for objIndex in range(3):
    for j in range(len(store)):
        if(store[j,objIndex,0] != 0):
            flux.append(store[j,objIndex,0])
            mux.append(store[j,objIndex,1])
            muy.append(store[j,objIndex,2])
            e1.append(store[j,objIndex,3])
            e2.append(store[j,objIndex,4])
            bkg.append(store[j,objIndex,4])
            psf.append(store[j,objIndex,6])
            mjd.append(store[j,objIndex,7])
            #wt.append(store[j,objIndex,10])
            back_h.append(store[j,objIndex,8])
            seeing.append(store[j,objIndex,9])
            zp.append(store[j,objIndex,10])
            fwhm.append(store[j,objIndex,11])
            airmass.append(store[j,objIndex,12])
            mphase.append(store[j,objIndex,13])
            mAngle.append(store[j,objIndex,14])

        
flux = np.array(flux)
mux = np.array(mux)
muy = np.array(muy)
e1 = np.array(e1)
e2 = np.array(e2)
bkg = np.array(bkg)
psf = np.array(psf)
wt = np.array(wt)
mjd = np.array(mjd)
zp = np.array(zp)
sign =[]
for j in range(len(flux)):
    if(True):
        #flux[j] = 0
        sign.append('b+')
    else:
        sign.append('r.')

flux = flux * np.log10(zp/2.5)
       
for j in range(len(flux)):
# =============================================================================
#     if(zp[j]<24):
#         continue
# =============================================================================
    if(mjd[j] < 58200):
        plt.plot( airmass[j], flux[j], 'k.', markersize = 0)
    elif(mjd[j] > 58200 and mjd[j] <59000):
        plt.plot( airmass[j], flux[j], 'k+', markersize = 3)
    else:
        plt.plot( airmass[j], flux[j], 'k*', markersize = 3)
    
plt.xlabel('Airmass')
plt.ylabel(' PSF')
        
