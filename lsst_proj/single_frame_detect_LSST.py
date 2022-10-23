#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 08:19:43 2022

@author: dutta26
"""


from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess

band = '2'
#Read the source catalog 
f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp_lsst.cat')
cat = f.readlines()
f.close()
raList =[]
decList=[]

for j in range(len(cat)):
    if((cat[j].split()[0]) == '#'):
     continue
    
    raList.append(float(cat[j].split()[5])) 
    decList.append(float(cat[j].split()[6]))
    
#Read IR coadd catalog 
band_coadd_df = pd.read_pickle('/scratch/halstead/d/dutta26/lsst/df1_'+band+'.pk1')
band_coadd_xx = np.array(band_coadd_df['xx'])
band_coadd_yy = np.array(band_coadd_df['yy'])
band_coadd_xy = np.array(band_coadd_df['xy'])
band_coadd_interpxx = np.array(band_coadd_df['interp_xx'])
band_coadd_interpyy = np.array(band_coadd_df['interp_yy'])
band_coadd_interpxy = np.array(band_coadd_df['interp_xy'])
band_coadd_flux = np.array(band_coadd_df['flux'])
star_bool=np.array(band_coadd_df['star_flag'])

band_coadd_df = np.array(band_coadd_df)



# =============================================================================
# star_arr = band_coadd_df[(np.where((band_coadd_df[:,2] == 1) & (band_coadd_df[:,3] > 2000) ))[0],  : ]
# #Now tuse k sigma clip to find usable stars. Just do for sigxx
# mean,median, std = sigma_clipped_stats(star_arr[:, 7])
# print (mean,median, std)
# 
# loc = np.where((band_coadd_df[:,2] == 1) & 
#                                       (band_coadd_df[:,7] >= mean-5*std) &
#                                       (band_coadd_df[:,7] <= mean+5*std) & 
#                                       (band_coadd_df[:,3] > 5000) )[0]
# for j in range(len(band_coadd_df)):
#     if(j in loc):
#         star_bool.append(1)
#     else:
#         star_bool.append(0)
#     
# star_bool = np.array(star_bool)
# =============================================================================

data_width = int(len(os.listdir('/scratch/halstead/d/dutta26/lsst/filter'+band +'/'))/2) +1
store = np.zeros((data_width, len(raList), 50), dtype = np.float32)
fileCount = -1
fileList = os.listdir('/scratch/halstead/d/dutta26/lsst/filter'+band +'/')


for file in fileList:
    
    if('10' not in file or 'flat' in file):
        continue
    fileCount += 1
    print (file, 'aa')
    
    f=fits.open('/scratch/halstead/d/dutta26/lsst/filter'+band +'/'+file+'/final_coadd.fits')
    data = np.array(f[0].data)
    f.close()
      
    if(len(data)<= 0):
        continue
    xList, yList = helper.convertToXY(raList, decList, '/scratch/halstead/d/dutta26/lsst/filter'+band +'/'+file+'/final_coadd.fits')
    
    ySize,xSize = np.shape(data)
    
    f=fits.open('/scratch/halstead/d/dutta26/lsst/filter'+band +'/'+file+'/sample.fits')
    temp_data = np.array(f[0].data)
    back_h = np.median(temp_data)
    mjd = float((f[0].header)['MJD-OBS'])
    airmass = float((f[0].header)['AIRMASS'])
    #Fix for moon phase is waxing gibbous
    
    mphase = float((f[0].header)['MOONPHS'])
    mAngle = float((f[0].header)['MOONDIS'])
    
    f.close()
    
    
    #Fist find sources and measure params . First pass
    cnt = cnt1 = cnt2 =cnt3 =0
    for j in range(len(xList)):
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(band_coadd_xx[j] <= 0 or band_coadd_yy[j]<=0 or np.isnan(band_coadd_xx[j]) or np.isnan(band_coadd_yy[j])):
            store[fileCount,j,0] = ra
            store[fileCount,j,1] = dec
            store[fileCount,j,2] = star_bool[j]
            store[fileCount,j,10] = x
            store[fileCount,j,11] = y
            continue
        
        size = np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
        if(size<4):
            size = 3.9
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        #Allow convergence only if (coadd flux >60)NO if S/N>4. Assume size ~ 4*seeing
        noise = np.sqrt( back_h * (size)) 
        SbyN = band_coadd_flux[j]/noise
        #if(band_coadd_flux[j]< 2):
        if(SbyN < 4):
            store[fileCount,j,0] = ra
            store[fileCount,j,1] = dec
            store[fileCount,j,2] = star_bool[j]
            store[fileCount,j,10] = x
            store[fileCount,j,11] = y
            continue
        
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
       
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        
       
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
        if(flux == None or flux<= 0 or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
            store[fileCount,j,0] = ra
            store[fileCount,j,1] = dec
            store[fileCount,j,2] = star_bool[j]
            store[fileCount,j,10] = x
            store[fileCount,j,11] = y
            continue
        cnt1 += 1
        #if(b_flag == 1 and cnt1 > 5000):
        #    sys.exit()
        store[fileCount,j,0:15] = ra, dec, star_bool[j], flux, mux, muy,  back_h, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag
    
# =============================================================================
#         if(j== 4):
#             sys.exit()
#         if(star_bool[j] == 1):
#             print (flux, band_coadd_flux[j])
#             hdu =fits.PrimaryHDU(cut)
#             hdu.writeto('/scratch/halstead/d/dutta26/lsst/lsst1/'+str(j)+'_'+str(flux)[0:4]+'.fits', overwrite = True)
#     sys.exit()
# =============================================================================
        
        
        
    #Find all good stars in the frame
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & (store[fileCount,:,12] == 99)  & (store[fileCount,:,14] == 0)))[0],  : ]

    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    print (sigma_clipped_stats(star_arr[:, 7]))
    threshold = 3 *median*3.14*back_h
    
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                          (store[fileCount,:,12] == 99)  & 
                                          (store[fileCount,:,3] > threshold)  & 
                                          (store[fileCount,:,14] == 0)))[0],  : ]
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    print (sigma_clipped_stats(star_arr[:, 7]))
    
    if(mean > 25 or std>3 or mean<1):
        store[fileCount, :,38] = -99
        store[fileCount, :,39] = -99
        store[fileCount, :,40] = -99
        continue
    
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                          (store[fileCount,:,7] >= mean-5*std) &
                                          (store[fileCount,:,7] <= mean+5*std) &
                                          (store[fileCount,:,12] == 99) & 
                                          (store[fileCount,:,3] > threshold) &
                                          (store[fileCount,:,14] == 0)))[0],  : ]

    
    
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr[:, 7]
    star_temp[:,1] = star_arr[:, 8]
    star_temp[:,2] = star_arr[:, 9]
    star_temp[:,3] = star_arr[:, 10]
    star_temp[:,4] = star_arr[:, 11]
    print (np.shape(star_arr))
    print (threshold, median, back_h)
    nStars = 10
    #Second pass. This pass we deconvolve with coadd PSF and recovolve with frame PSF of 10 nearest stars
    for j in range(len(xList)):
        
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        
        temp = np.copy(star_temp)
        temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
        temp = temp[temp[:,5].argsort()]
        
        #Check if same star. Then delete the entry
        if(temp[0,5]<5):
            temp = np.delete(temp, (0), axis = 0)
        
        #Checking for nans to avoid code from crashing
        if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
            continue
        if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
            continue
        if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
            continue
        
        #Conpute the average PSF values 
        avgSigxx_frame = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
        avgSigyy_frame = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
        avgSigxy_frame = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
        
        store[fileCount, j,38] = avgSigxx_frame
        store[fileCount, j,39] = avgSigyy_frame
        store[fileCount, j,40] = avgSigxy_frame
        
        
        
        
        
        
    
        #Use guess measure only when actual convergence fails 
        if(store[fileCount,j,12] == 99 ):
            continue
        #Find guess shapes
        guess_xx = band_coadd_xx[j]-band_coadd_interpxx[j] + avgSigxx_frame
        guess_yy = band_coadd_yy[j]-band_coadd_interpyy[j] + avgSigyy_frame
        guess_xy = band_coadd_xy[j]-band_coadd_interpxy[j] + avgSigxy_frame
        if(guess_xx<0 or guess_yy<0):
            continue

        #Make cutout
        ra = raList[j]
        dec = decList[j]
        x = int(round(xList[j]))
        y = int(round(yList[j]))
        mux_guess = xList[j] - x
        muy_guess = yList[j] - y
        if(band_coadd_xx[j] <= 0 or band_coadd_yy[j]<=0 or np.isnan(band_coadd_xx[j]) or np.isnan(band_coadd_yy[j])):
            continue
        size = np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
        if(size<4):
            size = 3.9
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        store[fileCount,j,2] =star_bool[j]
        store[fileCount,j,13] = v_flag
        store[fileCount,j,14] = b_flag
       
        
        #Measure cutout
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v3(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1, mux_guess, muy_guess)
        if(flux == None or e1== None or e2 == None or np.isnan(e1)):
            store[fileCount,j,12] = -99
            continue
        
        #print (j)    
        store[fileCount,j,31:38] = flux, mux, muy,  bkg, sigxx, sigyy, sigxy
        store[fileCount,j,12] = 1
        
        
        
        #sys.exit()

np.save('/scratch/halstead/d/dutta26/lsst/singleFrame_10star_' +str(band)+'_2.npy', store)        
        