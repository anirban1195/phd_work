#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 20:20:25 2024

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

avg_dev_xx = np.zeros((100, 31,8,8), dtype = np.float32)
avg_dev_yy = np.zeros((100, 31,8,8), dtype = np.float32)
avg_dev_xy = np.zeros((100, 31,8,8), dtype = np.float32)
avg_dev_size = np.zeros((100, 301,8,8), dtype = np.float32)
avg_dev_e1 = np.zeros((100, 31,8,8), dtype = np.float32)
avg_dev_e2 = np.zeros((100, 31,8,8), dtype = np.float32)
count = np.zeros((100, 31,8,8), dtype = np.float32)


bandList = ['g', 'r', 'i']
for band in bandList:
    for sliceIndex in range(0,100):
        print (sliceIndex)
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
        
        mean,median, std = sigma_clipped_stats(star_arr[:, 7])
        mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
        mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
        med_e1 = (median - median1)/(median+median1)
        med_e2 = 2*median2/(median+median1)
        med_size= np.sqrt(median+median1)
        for j in range(len(star_arr)):
            chip_no= int(star_arr[j,47])
            ccd_x = int(star_arr[j,48])
            ccd_y = int(star_arr[j,49])
            e1 = (star_arr[j,7] - star_arr[j,8])/(star_arr[j,7] + star_arr[j,8])
            e2 = 2*star_arr[j,9]/(star_arr[j,7] + star_arr[j,8])
            
            avg_dev_xx[sliceIndex, chip_no,ccd_x, ccd_y] = (median - star_arr[j,7])/median
            avg_dev_yy[sliceIndex, chip_no,ccd_x, ccd_y] = (median1 - star_arr[j,8])/median1
            avg_dev_xy[sliceIndex, chip_no,ccd_x, ccd_y] = (median2 - star_arr[j,9])/median2
            avg_dev_e1[sliceIndex, chip_no,ccd_x, ccd_y] = (med_e1 - e1)
            avg_dev_e2[sliceIndex, chip_no,ccd_x, ccd_y] = (med_e2 - e2)
            avg_dev_size[sliceIndex, chip_no,ccd_x, ccd_y] = (med_size - np.sqrt(star_arr[j,7]+star_arr[j,8]))/med_size
            count[sliceIndex, chip_no,ccd_x, ccd_y] += 1

my_dict={1: '33',
 2: '34',
 3: '43',
 4: '44',
 5: '32',
 6: '23',
 7: '24',
 8: '42',
 9: '35',
 10: '53',
 11: '45',
 12: '54',
 13: '22',
 14: '25',
 15: '52',
 16: '55',
 17: '31',
 18: '13',
 19: '41',
 20: '14',
 21: '36',
 22: '46',
 23: '21',
 24: '12',
 25: '15',
 26: '51',
 27: '26',
 28: '56',
 29: '11',
 30: '16',
 0: np.nan}  
stat_xx = np.zeros((2, 31,8,8), dtype = np.float32)
stat_yy = np.zeros((2, 31,8,8), dtype = np.float32)
stat_xy = np.zeros((2, 31,8,8), dtype = np.float32)
stat_size = np.zeros((2, 31,8,8), dtype = np.float32)
stat_e1 = np.zeros((2, 31,8,8), dtype = np.float32)
stat_e2 = np.zeros((2, 31,8,8), dtype = np.float32)
stat_size = np.zeros((2, 31,8,8), dtype = np.float32)


dev_xx_fits=np.zeros((44,53), dtype = np.float32)
err_xx_fits=np.zeros((44,53), dtype = np.float32)
dev_yy_fits=np.zeros((44,53), dtype = np.float32)
err_yy_fits=np.zeros((44,53), dtype = np.float32)
dev_xy_fits=np.zeros((44,53), dtype = np.float32)
err_xy_fits=np.zeros((44,53), dtype = np.float32)
dev_e1_fits=np.zeros((44,53), dtype = np.float32)
err_e1_fits=np.zeros((44,53), dtype = np.float32)
dev_e2_fits=np.zeros((44,53), dtype = np.float32)
err_e2_fits=np.zeros((44,53), dtype = np.float32)
dev_size_fits=np.zeros((44,53), dtype = np.float32)
err_size_fits=np.zeros((44,53), dtype = np.float32)
cnts_fits=np.zeros((44,53), dtype = np.float32)
for j in range(1,31):
    x_ccd =int(my_dict[j][0]) - 1
    y_ccd =int(my_dict[j][1]) - 1
    x_act = x_ccd*8 + x_ccd
    y_act = y_ccd*8 + y_ccd
    for k in range(8):
        for l in range(8):
            data = avg_dev_xx[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_xx[0,j,k,l] = median
            stat_xx[1,j,k,l] = std
            
            
            data = avg_dev_yy[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_yy[0,j,k,l] = median
            stat_yy[1,j,k,l] = std
            
            data = avg_dev_xy[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_xy[0,j,k,l] = median
            stat_xy[1,j,k,l] = std
            
            data = avg_dev_e1[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_e1[0,j,k,l] = median
            stat_e1[1,j,k,l] = std
            
            
            data = avg_dev_e2[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_e2[0,j,k,l] = median
            stat_e2[1,j,k,l] = std
            
            
            data = avg_dev_size[:, j,k,l]
            loc = np.where((data != 0) & ~np.isnan(data))[0]
            mean, median, std = sigma_clipped_stats(data[loc])
            stat_size[0,j,k,l] = median
            stat_size[1,j,k,l] = std
            
            
            dev_xx_fits[x_act+k, y_act+l] =stat_xx[0,j,k,l]
            err_xx_fits[x_act+k, y_act+l] =stat_xx[1,j,k,l]
            
            dev_yy_fits[x_act+k, y_act+l] =stat_yy[0,j,k,l]
            err_yy_fits[x_act+k, y_act+l] =stat_yy[1,j,k,l]
            
            dev_xy_fits[x_act+k, y_act+l] =stat_xy[0,j,k,l]
            err_xy_fits[x_act+k, y_act+l] =stat_xy[1,j,k,l]
            
            dev_e1_fits[x_act+k, y_act+l] =stat_e1[0,j,k,l]
            err_e1_fits[x_act+k, y_act+l] =stat_e1[1,j,k,l]
            
            dev_e2_fits[x_act+k, y_act+l] =stat_e2[0,j,k,l]
            err_e2_fits[x_act+k, y_act+l] =stat_e2[1,j,k,l]
            
            dev_size_fits[x_act+k, y_act+l] =stat_size[0,j,k,l]
            err_size_fits[x_act+k, y_act+l] =stat_size[1,j,k,l]
            
            cnts_fits[x_act+k, y_act+l] = len(data[loc])
            
hdu=fits.PrimaryHDU(np.rot90(dev_xx_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_xx.fits', overwrite=True)          
hdu=fits.PrimaryHDU(np.rot90(err_xx_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_xx.fits', overwrite=True)
                        
            
hdu=fits.PrimaryHDU(np.rot90(dev_yy_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_yy.fits', overwrite=True)           
hdu=fits.PrimaryHDU(np.rot90(err_yy_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_yy.fits', overwrite=True)         

hdu=fits.PrimaryHDU(np.rot90(dev_xy_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_xy.fits', overwrite=True)          
hdu=fits.PrimaryHDU(np.rot90(err_xy_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_xy.fits', overwrite=True)   


hdu=fits.PrimaryHDU(np.rot90(dev_e1_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_e1.fits', overwrite=True)          
hdu=fits.PrimaryHDU(np.rot90(err_e1_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_e1.fits', overwrite=True)   

hdu=fits.PrimaryHDU(np.rot90(dev_e2_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_e2.fits', overwrite=True)          
hdu=fits.PrimaryHDU(np.rot90(err_e2_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_e2.fits', overwrite=True)   


hdu=fits.PrimaryHDU(np.rot90(dev_size_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/dev_size.fits', overwrite=True)          
hdu=fits.PrimaryHDU(np.rot90(err_size_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/err_size.fits', overwrite=True)   

hdu=fits.PrimaryHDU(np.rot90(cnts_fits, k= 1, axes=(1,0)))
hdu.writeto('/scratch/bell/dutta26/abell_2390/stars_sf/counts.fits', overwrite=True)             