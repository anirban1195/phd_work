#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 19:36:06 2023

@author: dutta26
"""
import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import helper
from numpy.random import multivariate_normal
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
from astropy.stats import sigma_clipped_stats


id_num = '1001'

import matplotlib.pyplot as plt
from scipy.special import erf
sizex = sizey= 100
#sigx = 6.25
#sigy = 6.25
sigxy = 0
#fluxArr = [10, 20, 50, 100, 250, 500, 750, 1000, 2500, 5000, 10000]
#bkgArr= [ 50, 100, 200]
fluxArr = np.hstack((np.arange(5, 10, 1), np.arange(10, 100, 10), np.arange(100, 1000, 100), 
                     np.arange(1000, 10000, 1000) , np.arange(10000, 110000, 10000)))

bkgArr = [50, 100, 300, 500, 750, 1000, 1500, 2000]
sigxArr = [2.5, 3, 3.5, 4, 4.5, 5, 5.5,6, 6.5]
totRows = len(fluxArr)*len(bkgArr)*len(sigxArr)
store = np.ones((totRows, 11))
store = store*np.nan


count = 0
for bkg in bkgArr:
    for sigx in sigxArr:
        if(id_num == '1001'):
            sigy = sigx
            sigxy = 0
        if(id_num == '1002'):
            sigy = 4
            sigxy = 0
        if(id_num == '1003'):
            sigy = 4
            sigxy = 6
        print (bkg, sigx, sigy)
        biasArr=[]
        biasArr2=[]
        biasArr3=[]
        area = np.pi*np.sqrt(sigx**2 * sigy**2 - sigxy**2)
        size= np.sqrt((sigx**2 + sigy**2)/2 )
        for flux in fluxArr:
            
            temp_measured =[]
            temp_measured2=[]
            temp_measured3=[]
            for j in range(500):
                muArr= [sizex/2.0-0.5, sizey/2.0-0.5]
                cov = [[(sigx**2), sigxy], [sigxy, (sigy**2)]]
                const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
                if(const < 1):
                    const = 1
                x, y = np.random.multivariate_normal(muArr, cov, const).T
                x = np.int32(np.round(x))
                y = np.int32(np.round(y))
                obj = np.zeros((sizex,sizey))
                np.add.at(obj, (y,x), 1)
                
                
                noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
                noise = np.round(noise)
                obj = obj + noise
                flux_calc, mux_calc, muy_calc, e1_calc, e2_calc, bkg_calc, psf_calc, sigxx_calc, sigyy_calc, sigxy_calc = helper.measure_new(obj, lut1, lut2, flux, 0, 0, sigx, sigy, sigxy, 1,1)#measure(obj, sigx, sigy, sigxy, 1,0)
                
                temp_measured.append(flux_calc)
                temp_measured2.append(sigxx_calc)
                temp_measured3.append(mux_calc)
            
            mean,median, std = sigma_clipped_stats(temp_measured)
            biasArr.append(std/np.sqrt(flux + 4*area*bkg))
            
            mean,median, std = sigma_clipped_stats(temp_measured2)
            biasArr2.append( std/np.sqrt(size**4/flux + 4*bkg*np.pi*size**6/flux**2 ) )
            
            mean,median, std = sigma_clipped_stats(temp_measured3)
            biasArr3.append(std/np.sqrt(size**2/flux + 4*bkg*np.pi*size**4/flux**2 ))
            
            store [count, 0], store [count, 1], store [count, 2], store [count, 3] = flux, np.sqrt(bkg)*area, biasArr[len(biasArr)-1], biasArr2[len(biasArr2)-1]
            store [count, 4], store [count, 5]=np.sqrt(flux + 4*area*bkg), np.sqrt(size**4/flux + 4*bkg*np.pi*size**6/flux**2 )
            store [count, 6], store [count, 7]=bkg, area
            store [count, 8], store [count, 9],  store [count, 10] = biasArr3[len(biasArr3)-1], np.sqrt(size**2/flux + 4*bkg*np.pi*size**4/flux**2 ), median
            #store2[count, 0], store2[count, 1], store2[count, 2] = flux, bkg*3.14*sigx*sigy, (np.nanmedian(temp_measured2) - sigx**2/2)/ sigx**2
            
            count += 1
        #biasArr= np.array(biasArr)
        #biasArr2= np.array(biasArr2)
        #plt.plot(np.log10(fluxArr), biasArr2, label = str(bkg*3.14*16)[0:5])
        
#plt.ylabel('Sigmaxx Bias')   
#plt.xlabel('Log True Flux')     
#plt.legend()
np.save('/scratch/bell/dutta26/abell_2390/calib_forcedDist1_err' +str(id_num)+'.npy', store)    
#np.save('/scratch/bell/dutta26/abell_2390/sigma_calib_66.npy', store2)    





