#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 13:27:56 2024

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


loc = np.where((band_coadd_data[:,2] != 1)  & (band_coadd_data[:,3] >0) & (band_coadd_data[:,12] == 99))[0]
magArr = 25-2.5*np.log10(band_coadd_data[:,3])
sizeArr = np.sqrt(band_coadd_data[:,7]+ band_coadd_data[:,8])
loc = np.where((sizeArr>2.8) & (sizeArr<3.4) & (magArr>17) & (magArr<22.5))[0]
print (len(loc))


bandList = ['i']

fileCount = -1
dev_arr  = np.zeros( (160, len(loc), 30), dtype = np.float32)

for band in bandList:
    band_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.npy')
    band_sf_df = np.load('/scratch/bell/dutta26/abell_2390/'+band+'_sf.npy')
    fileList = os.listdir(bandLoc +band+'/')
    for file in fileList:
        #Reject wt and temp files.
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
            if(guess_xx == 0 or guess_yy == 0):
                continue
            
            flux_c, mux_c, muy_c, e1_c, e2_c, bkg_c, psf_c, sigxx_c, sigyy_c, sigxy_c = helper.measure_new(cut, lut1, lut2)
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,1)
            flux_f, mux_f, muy_f, e1_f, e2_f, bkg_f, psf_f, sigxx_f, sigyy_f, sigxy_f = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,0)
            
            flux_t, mux_t, muy_t, e1_t, e2_t, bkg_t, psf_t, sigxx_t, sigyy_t, sigxy_t = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 100, 1)

            
            back_h = band_sf_df[fileCount,j,15]
            dev_arr[fileCount, sourceCount, 0] = flux_expected/ np.sqrt(flux_expected + 4*3.14* (guess_xx+guess_yy)*back_h )
            
            if(flux_c == None or flux_c <0 or psf_c == None or np.isnan(flux_c) or np.isnan(psf_c) or sigxx_c< 0 or np.abs(sigxx_c)> 70 or sigyy_c< 0 or np.abs(sigyy_c)> 70):
                dev_arr[fileCount, sourceCount, 1] = -99
                dev_arr[fileCount, sourceCount, 2] = -99
                dev_arr[fileCount, sourceCount, 3] = -99
                dev_arr[fileCount, sourceCount, 4] = -99
                dev_arr[fileCount, sourceCount, 18 ] = -99
            else:
                dev_arr[fileCount, sourceCount, 1] = flux_c
                dev_arr[fileCount, sourceCount, 2] = sigxx_c
                dev_arr[fileCount, sourceCount, 3] = sigyy_c
                dev_arr[fileCount, sourceCount, 4] = sigxy_c
                dev_arr[fileCount, sourceCount, 18 ] = np.sqrt(  (guessmux - mux_c)**2 +  (guessmuy - muy_c)**2 )
                
                
                
            if(flux == None or flux <0 or psf == None or np.isnan(flux) or np.isnan(psf) or sigxx< 0 or np.abs(sigxx)> 70 or sigyy< 0 or np.abs(sigyy)> 70):
                dev_arr[fileCount, sourceCount, 5] = -99
                dev_arr[fileCount, sourceCount, 6] = -99
                dev_arr[fileCount, sourceCount, 7] = -99
                dev_arr[fileCount, sourceCount, 8] = -99
                dev_arr[fileCount, sourceCount, 19 ] = -99
            else:
                dev_arr[fileCount, sourceCount, 5] = flux
                dev_arr[fileCount, sourceCount, 6] = sigxx
                dev_arr[fileCount, sourceCount, 7] = sigyy
                dev_arr[fileCount, sourceCount, 8] = sigxy
                dev_arr[fileCount, sourceCount, 19 ] = np.sqrt(  (guessmux - mux)**2 +  (guessmuy - muy)**2 )
                
                
                
            if(flux_f == None or flux_f <0 or psf_f == None or np.isnan(flux_f) or np.isnan(psf_f) or sigxx_f< 0 or np.abs(sigxx_f)> 70 or sigyy_f< 0 or np.abs(sigyy_f)> 70):
                dev_arr[fileCount, sourceCount, 9] = -99
                dev_arr[fileCount, sourceCount, 10] = -99
                dev_arr[fileCount, sourceCount, 11] = -99
                dev_arr[fileCount, sourceCount, 12] = -99
                dev_arr[fileCount, sourceCount, 20 ] = -99
            else:
                dev_arr[fileCount, sourceCount, 9] = flux_f
                dev_arr[fileCount, sourceCount, 10] = sigxx_f
                dev_arr[fileCount, sourceCount, 11] = sigyy_f
                dev_arr[fileCount, sourceCount, 12] = sigxy_f
                dev_arr[fileCount, sourceCount, 20 ] = np.sqrt(  (guessmux - mux_f)**2 +  (guessmuy - muy_f)**2 )
                
                
                
                
            if(flux_t == None or flux_t <0 or psf_t == None or np.isnan(flux_t) or np.isnan(psf_t) or sigxx_t< 0 or np.abs(sigxx_t)> 70 or sigyy_t< 0 or np.abs(sigyy_t)> 70):
                dev_arr[fileCount, sourceCount, 21] = -99
                dev_arr[fileCount, sourceCount, 22] = -99
                dev_arr[fileCount, sourceCount, 23] = -99
                dev_arr[fileCount, sourceCount, 24] = -99
                dev_arr[fileCount, sourceCount, 25 ] = -99
            else:
                dev_arr[fileCount, sourceCount, 21] = flux_t
                dev_arr[fileCount, sourceCount, 22] = sigxx_t
                dev_arr[fileCount, sourceCount, 23] = sigyy_t
                dev_arr[fileCount, sourceCount, 24] = sigxy_t
                dev_arr[fileCount, sourceCount, 25 ] = np.sqrt(  (guessmux - mux_t)**2 +  (guessmuy - muy_t)**2 )
                
            
            dev_arr[fileCount, sourceCount, 13] = flux_expected
            dev_arr[fileCount, sourceCount, 14] = guess_xx
            dev_arr[fileCount, sourceCount, 15] = guess_yy
            dev_arr[fileCount, sourceCount, 16] = guess_xy
            if(band == 'i'):
                dev_arr[fileCount, sourceCount, 17 ] == 1
            elif(band == 'r'):
                dev_arr[fileCount, sourceCount, 17 ] == 2
                
                
            
            
            
       
            

np.save('/scratch/bell/dutta26/test_data_i_t.npy', dev_arr)
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

