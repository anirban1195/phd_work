#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 21:43:51 2023

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
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


band = 'g'
for sliceIndex in range(0,100):
    #sliceIndex = 1
    sizex =sizey = 100
    tot = np.zeros((sizex, sizey))
    
    store = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')
    
    star_arr = store[sliceIndex,(np.where((store[sliceIndex,:,2] == 1) & (store[sliceIndex,:,12] == 99) 
                               & (store[sliceIndex,:,13] == 0) & (store[sliceIndex,:,14] == 0)))[0],  : ]
    
    
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    print (mean,median, std)
    
    star_arr = store[sliceIndex, (np.where((store[sliceIndex,:,2] == 1) & 
                                          (store[sliceIndex,:,7] >= mean-3*std) &
                                          (store[sliceIndex,:,7] <= mean+3*std) &
                                          (store[sliceIndex,:,8] >= mean1-3*std1) &
                                          (store[sliceIndex,:,8] <= mean1+3*std1) &
                                          (store[sliceIndex,:,3] < 1e9) &
                                          (store[sliceIndex,:,12] == 99) & 
                                          (store[sliceIndex,:,13] == 0) & 
                                          (store[sliceIndex,:,14] == 0)))[0],  : ]
    
    index_loc = (np.where((store[sliceIndex,:,2] == 1) & 
                                          (store[sliceIndex,:,7] >= mean-3*std) &
                                          (store[sliceIndex,:,7] <= mean+3*std) &
                                          (store[sliceIndex,:,8] >= mean1-3*std1) &
                                          (store[sliceIndex,:,8] <= mean1+3*std1) &
                                          (store[sliceIndex,:,3] < 1e9) &
                                          (store[sliceIndex,:,12] == 99) & 
                                          (store[sliceIndex,:,13] == 0) & 
                                          (store[sliceIndex,:,14] == 0)))[0]
    
    f_no = str(np.median(star_arr[:,44]))
    for files in os.listdir('/scratch/bell/dutta26/abell_2390/'+band+'/temp/'):
        if('weight' in files):
            continue
        if(f_no in files):
            swarp_coadd_file = '/scratch/bell/dutta26/abell_2390/'+band+'/temp/' + files
            f = fits.open(swarp_coadd_file)
            shapey, shapex = np.shape(f[0].data)
            f.close()
            
    for files in os.listdir('/scratch/bell/dutta26/abell_2390/'+band+'/'):
        if('weight' in files or 'temp' in files):
            continue
        if(f_no in files):
            uncoadd_file = '/scratch/bell/dutta26/abell_2390/'+band+'/' + files
    
    f=fits.open(uncoadd_file)
    hdr = f[0].header
    f.close()
    
    
    
    os.mkdir('/scratch/bell/dutta26/abell_2390/stars_sf/'+band+'/'+str(sliceIndex)+'/')    
    for j in range(len(star_arr)):
        print (j)
        x= int(star_arr[j,10])
        y= int(star_arr[j,11])
        if(x<0 or y<0 or x>shapex or y>shapey):
            continue
        f = fits.open(swarp_coadd_file)[0]
        wcs = WCS(f.header)
        
        cutout = Cutout2D(f.data, position=(x,y), size=(100,100), wcs=wcs)
        f.data = cutout.data
        f.header.update(cutout.wcs.to_header())
        cnt = 0
        initial_len = len(f.header)
        cnt1 =0
        for key in hdr:
            cnt1 += 1
            if(key == ''):
                #f.header[key]=hdr[key][cnt]
                #f.header.append(('',hdr[key][cnt]), end= True)
                
                f.header.insert(len(f.header)-1, ('',hdr[key][cnt]))
                cnt += 1
                continue
            
            f.header[key] = (hdr[key], hdr.comments[key])
        f.header['TOTCNT'] =  (star_arr[j,3], 'Total Counts from source')
        f.header['SIGXX'] =  (2*star_arr[j,7], 'Sigmaxx of source')
        f.header['SIGYY'] =  (2*star_arr[j,8], 'Sigmayy of source')
        f.header['SIGXY'] =  (2*star_arr[j,9], 'Sigmaxy of source')
        f.header['INDEX'] =  index_loc[j]
        f.header['CHIPX'] =  star_arr[j,48]
        f.header['CHIPY'] =  star_arr[j,49]
        f.header['CCDNO'] =  star_arr[j,47]
        f.header['COMMENT'] = 'Created by Dutta,Peterson for MESSIER Project'
        f.header['COMMENT'] = 'From A2390 field'
        
        out_name = '/scratch/bell/dutta26/abell_2390/stars_sf/'+band+'/'+str(sliceIndex)+'/star_'+str(j)+'.fits'
        f.writeto(out_name, overwrite=True)
        

        
        
    