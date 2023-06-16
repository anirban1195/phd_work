#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 22:00:30 2023

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
band = 'i'
for sliceIndex in range(0,100):
    #sliceIndex = 1
    sizex =sizey = 100
    tot = np.zeros((sizex, sizey))
    
    store = np.load('/scratch/bell/dutta26/abell_2390/test2_'+band+'.npy')
    
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
    
    f_no = str(np.median(star_arr[:,44]))
    for files in os.listdir('/scratch/bell/dutta26/abell_2390/'+band+'/temp/'):
        if('weight' in files):
            continue
        if(f_no in files):
            
            f = fits.open('/scratch/bell/dutta26/abell_2390/'+band+'/temp/' + files)
            img_data = f[0].data
            f.close()
            
            print('/scratch/bell/dutta26/abell_2390/'+band+'/temp/' + files) 
            
    os.mkdir('/scratch/bell/dutta26/abell_2390/stars_sf/'+band+'/'+str(sliceIndex)+'/')    
    for j in range(len(star_arr)):
        print (j)
        x= int(star_arr[j,10])
        y= int(star_arr[j,11])
        
        cut = img_data[y-50:y+50, x-50:x+50]
        hdu=fits.PrimaryHDU(cut)
        out_name = '/scratch/bell/dutta26/abell_2390/stars_sf/'+band+'/'+str(sliceIndex)+'/star_'+str(j)+'.fits'
        hdu.writeto(out_name, overwrite = True)
        
        
    