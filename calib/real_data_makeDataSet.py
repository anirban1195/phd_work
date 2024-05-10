#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 11:09:42 2024

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import sys,os
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
from astropy import wcs
import matplotlib.pyplot as plt
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')



band = 'i'
bandLoc = '/scratch/bell/dutta26/abell_2390/'
ir_coadd_df = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
band_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.npy')
band_sf_df = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')

# =============================================================================
# loc = np.where((i_coadd_data[:,2] != 1)  & (i_coadd_data[:,3] >0) & (i_coadd_data[:,12] == 99))[0]
# plt.figure(figsize=(13,10))
# plt.plot(np.sqrt(i_coadd_data[loc,7]+ i_coadd_data[loc,8]), 25-2.5*np.log10(i_coadd_data[loc,3]),'b.', markersize= 2)  
# 
# loc = np.where((i_coadd_data[:,2] == 1)  & (i_coadd_data[:,3] >0) & (i_coadd_data[:,12] == 99))[0]
# plt.plot(np.sqrt(i_coadd_data[loc,7]+ i_coadd_data[loc,8]), 25-2.5*np.log10(i_coadd_data[loc,3]),'r.', markersize= 2)  
# 
# plt.xlabel('Size')
# plt.ylabel('Magnitude')
# =============================================================================

loc = np.where((band_coadd_data[:,2] != 1)  & (band_coadd_data[:,3] >0) & (band_coadd_data[:,12] == 99))[0]
magArr = 25-2.5*np.log10(band_coadd_data[:,3])
sizeArr = np.sqrt(band_coadd_data[:,7]+ band_coadd_data[:,8])
loc = np.where((sizeArr>2.7) & (sizeArr<3.3) & (magArr>17) & (magArr<23.25))[0]
print (len(loc))
fileCount = -1
fileList = os.listdir(bandLoc +band+'/')
fluxDev_arr  = np.zeros( (160, len(loc)), dtype = np.float32)
dev_arr  = np.zeros( (160, len(loc), 10), dtype = np.float32)



for file in fileList:

    if('star' in file or '.weight' in file or 'temp' in file):
        continue
    fileCount += 1
    f=fits.open('/scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file))
    data = f[0].data
    f.close() 
    datay, datax = np.shape(data)
    sourceCount = -1
    print (file)
    for j in loc:
        sourceCount += 1
        #print(j)
        x = int(round(band_sf_df[fileCount,j,10]))
        y = int(round(band_sf_df[fileCount,j,11]))
        v_flag = band_sf_df[fileCount,j,13]
        b_flag = band_sf_df[fileCount,j,14] 
        if(v_flag != 0 or b_flag != 0 or band_sf_df[fileCount, j,62] == 1):
            continue
        
        cut = data[y-35: y+35, x-35: x+35]
        badLoc= np.where(cut <= 0)
        if(len(badLoc[0])> 0):
            continue
        
        sizey, sizex = np.shape(cut)
        if(sizex != 70 or sizey !=70):
            continue
        avgSigxx_frame = band_sf_df[fileCount, j,38]  
        avgSigyy_frame = band_sf_df[fileCount, j,39]  
        avgSigxy_frame = band_sf_df[fileCount, j,40]  
        #print (x,y, avgSigxx_frame, np.shape(cut))
        if(avgSigxx_frame> 25 or avgSigyy_frame>25 or x<0 or y<0 or x>datax or y>datay):
            continue
        
        guess_xx =  avgSigxx_frame
        guess_yy = avgSigyy_frame
        guess_xy = avgSigxy_frame
        guessmux = band_sf_df[fileCount,j,77]
        guessmuy = band_sf_df[fileCount,j,78]
        
        flux_expected = band_sf_df[fileCount,j,60]
        
        
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,1)
        
        flux1, mux1, muy1, e11, e21, bkg1, psf1, sigxx1, sigyy1, sigxy1 = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,0)
        
        back_h = band_sf_df[fileCount,j,15]
        dev_arr[fileCount, sourceCount, 0] = flux_expected/ np.sqrt(flux_expected + 4*3.14* (guess_xx+guess_yy)*back_h )
        
        
        if(band_sf_df[fileCount, j, 12] == 99):
            dev_arr[fileCount, sourceCount, 1] = (band_sf_df[fileCount, j, 3] - flux_expected)/flux_expected
            dev_arr[fileCount, sourceCount, 2] = (band_sf_df[fileCount, j, 7] - guess_xx) / guess_xx
            dev_arr[fileCount, sourceCount, 3] = (band_sf_df[fileCount, j, 9] - guess_xy) / guess_xy
        elif(band_sf_df[fileCount, j, 12] == 1):
            if(flux == None  or np.isnan(flux) or np.isinf(flux)):
                dev_arr[fileCount, sourceCount, 4] = -99
                dev_arr[fileCount, sourceCount, 5] = -99
                dev_arr[fileCount, sourceCount, 6] = -99
                
            else:
                dev_arr[fileCount, sourceCount, 4] = (flux - flux_expected)/flux_expected
                dev_arr[fileCount, sourceCount, 5] = (sigxx - guess_xx) / guess_xx
                dev_arr[fileCount, sourceCount, 6] = (sigxy - guess_xy) / guess_xy
               
                
                
            if(flux1 == None  or np.isnan(flux1) or np.isinf(flux1)):
                dev_arr[fileCount, sourceCount, 7] = -99
                dev_arr[fileCount, sourceCount, 8] = -99
                dev_arr[fileCount, sourceCount, 9] = -99
                hdu = fits.PrimaryHDU(cut)
                hdu.writeto('/scratch/bell/dutta26/abell_2390/bkg_test/'+str(fileCount) + '_'+str(sourceCount) +'.fits')
            else:
                dev_arr[fileCount, sourceCount, 7] = (flux1 - flux_expected)/flux_expected
                dev_arr[fileCount, sourceCount, 8] = (sigxx1 - guess_xx) / guess_xx
                dev_arr[fileCount, sourceCount, 9] = (sigxy1 - guess_xy) / guess_xy
        else:
            continue
            

#np.save('/scratch/bell/dutta26/test_data.npy', dev_arr)
#dev_arr = np.load('/scratch/bell/dutta26/test_data.npy')

# =============================================================================
# loc = np.where((dev_arr[:, :, 1] > -1) & (dev_arr[:, :, 1] < 1) & (dev_arr[:, :, 1] != 0)) 
# n, bins, patches = plt.hist(x=dev_arr[loc[0], loc[1], 1] , bins='auto',histtype=u'step', color='r', label='Tradional Moment Matching', density = True)
# print (np.std(dev_arr[loc[0], loc[1], 1]),len(loc[0]))
# 
# loc = np.where((dev_arr[:, :, 4] > -1) & (dev_arr[:, :, 4] < 1) & (dev_arr[:, :, 4] != 0))    
# n, bins, patches = plt.hist(x=dev_arr[loc[0], loc[1], 4], bins='auto',histtype=u'step', color='b', label='Forced Measurement', density = True)
# print (np.std(dev_arr[loc[0], loc[1], 4]), len(loc[0]), sigma_clipped_stats(dev_arr[loc[0], loc[1], 4]))
# 
# loc = np.where((dev_arr[:, :, 7] > -1) & (dev_arr[:, :, 7] < 1) & (dev_arr[:, :, 7] != 0))    
# n, bins, patches = plt.hist(x=dev_arr[loc[0], loc[1], 7], bins='auto',histtype=u'step', color='g', label='Forced Photometry', density = True)
# print (np.std(dev_arr[loc[0], loc[1], 7]),len(loc[0]), sigma_clipped_stats(dev_arr[loc[0], loc[1], 7]))
# plt.legend()
# plt.xlabel('Fractional Error in Flux')
# plt.savefig('/home/dutta26/codes/calib/flux_comp.png')
# plt.close()
# =============================================================================

