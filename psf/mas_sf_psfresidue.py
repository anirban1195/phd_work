#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:27:42 2023

@author: dutta26
"""


import sys,os
from astropy.io import fits
import numpy as np
import pandas as pd 
from astropy.stats import sigma_clipped_stats
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
import helper

for sliceIndex in range(19,40):
    #sliceIndex = 1
    sizex =sizey = 60
    tot = np.zeros((sizex, sizey))
    
    store = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
    
    
    
    
    star_arr = store[sliceIndex,(np.where((store[sliceIndex,:,2] == 1) & (store[sliceIndex,:,12] == 99) 
                               & (store[sliceIndex,:,13] == 0) & (store[sliceIndex,:,14] == 0)))[0],  : ]
    
    
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    print (mean,median, std)
    
    #Test my scheme
    x = np.linspace(0, sizex-1, sizex)
    y = np.linspace(0, sizey-1, sizey)
    x= x -sizex/2.0 + 0.5 
    y= y -sizey/2.0 + 0.5 
    x, y = np.meshgrid(x, y)
    mux = muy = 0
    flux = 1000
    alphax = np.sqrt(2*median)
    alphay = np.sqrt(2*median1)
    alphaxy = 2*median2
    arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
    A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
    k=(A * np.exp(-((x-mux)**2/(arbConst*alphax**2)+ (y-muy)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy)*(x-mux)/(arbConst* alphax**2 * alphay**2 ) )))
    k = k/np.sum(k)
    reference = k*flux
    hdu=fits.PrimaryHDU(reference)
    hdu.writeto('/scratch/bell/dutta26/abell_2390/psf/reference_'+str(sliceIndex)+'.fits', overwrite = True)
    
    star_arr = store[sliceIndex, (np.where((store[sliceIndex,:,2] == 1) & 
                                          (store[sliceIndex,:,7] >= mean-3*std) &
                                          (store[sliceIndex,:,7] <= mean+3*std) &
                                          (store[sliceIndex,:,8] >= mean1-3*std1) &
                                          (store[sliceIndex,:,8] <= mean1+3*std1) &
                                          (store[sliceIndex,:,3] < 1e9) &
                                          (store[sliceIndex,:,12] == 99) & 
                                          (store[sliceIndex,:,13] == 0) & 
                                          (store[sliceIndex,:,14] == 0)))[0],  : ]
    
    f_no = str(np.median(star_arr[:,44]))
    for files in os.listdir('/scratch/bell/dutta26/abell_2390/r/temp/'):
        if('weight' in files):
            continue
        if(f_no in files):
            
            f = fits.open('/scratch/bell/dutta26/abell_2390/r/temp/' + files)
            img_data = f[0].data
            f.close()
    
    
    #sys.exit()
    muxArr=[]
    muyArr=[]
    for j in range(len(star_arr)):
        print (j)
        x = int(round(star_arr[j, 10]))
        y = int(round(star_arr[j, 11]))
        size = int(sizex/2)
        cut = img_data[y-size:y+size, x-size:x+size]
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
        if(mux == None or muy == None):
            continue
        muxArr.append(mux)
        muyArr.append(muy)
        
        
    # =============================================================================
    #     #Test my scheme
    #     x = np.linspace(0, sizex-1, sizex)
    #     y = np.linspace(0, sizey-1, sizey)
    #     x= x -sizex/2.0 + 0.5 
    #     y= y -sizey/2.0 + 0.5 
    #     x, y = np.meshgrid(x, y)
    #     mux = muy = 1
    #     sigxx = sigyy= 4.5
    #     sigxy = 0
    #     flux = 1000
    #     alphax = np.sqrt(2*sigxx)
    #     alphay = np.sqrt(2*sigyy)
    #     alphaxy = 2*sigxy
    #     arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
    #     A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
    #     k=(A * np.exp(-((x-mux)**2/(arbConst*alphax**2)+ (y-muy)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy)*(x-mux)/(arbConst* alphax**2 * alphay**2 ) )))
    #     k = k/np.sum(k)
    #     k = k*flux
    # =============================================================================
        
        #Now create a perfect gaussian
        x = np.linspace(0, sizex-1, sizex)
        y = np.linspace(0, sizey-1, sizey)
        x= x -sizex/2.0 + 0.5 
        y= y -sizey/2.0 + 0.5 
        x, y = np.meshgrid(x, y)
        
        alphax = np.sqrt(2*sigxx)
        alphay = np.sqrt(2*sigyy)
        alphaxy = 2*sigxy
        arbConst = 2*(1- (alphaxy/(alphax*alphay))**2 )  
        A =1/(2*np.pi*alphax*alphay*np.sqrt(1- (alphaxy/(alphax*alphay))**2 ))
        k=(A * np.exp(-((x-mux)**2/(arbConst*alphax**2)+ (y-muy)**2/(arbConst*alphay**2) - 2*alphaxy*(y-muy)*(x-mux)/(arbConst* alphax**2 * alphay**2 ) )))
        k = k/np.sum(k)
        k = k*flux
        
        res = (cut - k)/flux
        res = np.roll(res, -int(round(muy)) , axis = 0)
        res = np.roll(res, -int(round(mux)) , axis = 1)
        tot += res
        
    tot = tot/len(star_arr)
    tot =tot - np.median(tot)
    hdu=fits.PrimaryHDU(tot)
    hdu.writeto('/scratch/bell/dutta26/abell_2390/psf/residue_'+str(sliceIndex)+'.fits', overwrite = True)
    
    # =============================================================================
    # tot = (tot/reference) *100 
    # tot[tot> 10000] = 0   
    # tot[tot< -10000] = 0   
    # hdu=fits.PrimaryHDU(tot)
    # hdu.writeto('/scratch/bell/dutta26/abell_2390/residue1.fits', overwrite = True)
    # =============================================================================
        
