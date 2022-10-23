#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 07:02:29 2022

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

band = 'i'
#band = str(sys.argv[1])
#Read the source catalog 
source_df = pd.read_pickle('/home/dutta26/codes/source_list.pk1')
#source_df = pd.read_pickle(str(sys.argv[2]))
star_bool = np.array(source_df['star_bool'])
raList = np.array(source_df['ra'])
decList = np.array(source_df['dec'])

#Read IR coadd catalog 
band_coadd_df = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
band_coadd_xx = np.array(band_coadd_df['xx'])
band_coadd_yy = np.array(band_coadd_df['yy'])
band_coadd_xy = np.array(band_coadd_df['xy'])
band_coadd_interpxx = np.array(band_coadd_df['interp_xx'])
band_coadd_interpyy = np.array(band_coadd_df['interp_yy'])
band_coadd_interpxy = np.array(band_coadd_df['interp_xy'])
band_coadd_flux = np.array(band_coadd_df['flux'])
band_coadd_df = np.array(band_coadd_df)



data_width = int(len(os.listdir('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'))/2.0) +1
store = np.zeros((data_width, len(raList), 50), dtype = np.float32)
fileCount = -1
for file in os.listdir('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'):
    
    if('.weight' in file):
        continue
    fileCount += 1
    
    
    #Run swarp to make a good image
    f= open('/home/dutta26/temp.ascii', 'w+')
    f.write('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'+file)
    f.close()
    
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
    #Read the swarp output
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
    data = np.array(f[0].data)
    f.close()        
    
    xList, yList = helper.convertToXY(raList, decList, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
    
    
    #Read the file data
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'+file)
    back_h = float((f[0].header)['SKYBG'])
    seeing = float((f[0].header)['SEEING'])
    zp = float((f[0].header)['MAGZERO'])
    fwhm = float((f[0].header)['FWHM_FLT'])
    mjd = float((f[0].header)['MJD-MID'])
    airmass = float((f[0].header)['AIRMASS'])
    #Fix for moon phase is waxing gibbous
    if(type((f[0].header)['MOONPHSE']) is str):
        mphase = -1.00000
    else:
        mphase = float((f[0].header)['MOONPHSE'])
    mAngle = float((f[0].header)['MOON_D'])
    expTime = float((f[0].header)['EXPTIME'])  
    focus = float((f[0].header)['TELFOCUS'])
    zp_n = float((f[0].header)['PHOTZP_N'])
    skymag = float((f[0].header)['SKYMAG']) 
    depth = float((f[0].header)['PHOTDPTH']) 
    mRa = float((f[0].header)['MOON_RA']) 
    mDec = float((f[0].header)['MOON_DEC']) 
    flux_scale = float((f[0].header)['FLXSCALE']) 
    f.close()
    
    
    ySize,xSize = np.shape(data)
    
    #Fist find sources and measure params . First pass
    cnt = cnt1 = cnt2 =cnt3 =0
    for j in range(len(xList)):
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(band_coadd_xx[j] <= 0 or band_coadd_yy[j]<=0 or np.isnan(band_coadd_xx[j]) or np.isnan(band_coadd_yy[j])):
            store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
            continue
         #Allow convergence only if coadd flux >60
        if(band_coadd_flux[j]< 30):
            store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
            continue
        
        size = np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
        if(size<4):
            size = 4
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        #If star then check if has weird lines through center or nans
        if(star_bool[j] == 1 ):
            v_flag = helper1.vert_stripe(cut)
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        
       
        t1, t2, t3, t4, t5, t6, t7 = measure_pythonV.measure_v2T(cut)
        if(t1 == None or t2<= 0 or t3 == None or t4==None or  np.isnan(t1) or np.isnan(t2) or np.isnan(t3)):
            store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
            continue
        
        store[fileCount,j,0:31] = ra, dec, star_bool[j], t1, t2, t3,  t4, t5, t6, t7, x, y, 99, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
    
    #Find all good stars in the frame
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & (store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) & (store[fileCount,:,14] == 0)))[0],  : ]

    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    if(mean > 25 or std>3 or mean<2):
        store[fileCount, :,38] = -99
        store[fileCount, :,39] = -99
        store[fileCount, :,40] = -99
        continue
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                          (store[fileCount,:,7] >= mean-5*std) &
                                          (store[fileCount,:,7] <= mean+5*std) &
                                          (store[fileCount,:,12] == 99) & 
                                          (store[fileCount,:,13] == 0) & 
                                          (store[fileCount,:,14] == 0)))[0],  : ]

    
    
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr[:, 7]
    star_temp[:,1] = star_arr[:, 8]
    star_temp[:,2] = star_arr[:, 9]
    star_temp[:,3] = star_arr[:, 10]
    star_temp[:,4] = star_arr[:, 11]
    
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
        if(np.isnan(np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            continue
        if(np.isnan(np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            continue
        if(np.isnan(np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5]))):
            continue
        
        #Conpute the average PSF values 
        avgSigxx_frame = np.sum( temp[0:10, 0] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigyy_frame = np.sum( temp[0:10, 1] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        avgSigxy_frame = np.sum( temp[0:10, 2] * 1/temp[0:10, 5])/np.sum(1/temp[0:10, 5])
        
        store[fileCount, j,38] = avgSigxx_frame
        store[fileCount, j,39] = avgSigyy_frame
        store[fileCount, j,40] = avgSigxy_frame
        

        #Use guess measure only when actual convergence fails 
        if(store[fileCount,j,12] == 99):
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
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(band_coadd_xx[j] <= 0 or band_coadd_yy[j]<=0 or np.isnan(band_coadd_xx[j]) or np.isnan(band_coadd_yy[j])):
            continue
        size = np.sqrt(band_coadd_xx[j] + band_coadd_yy[j])
        if(size<4):
            size = 4
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        #Measure cutout
        t1, t2, t3, t4, t5, t6, t7= measure_pythonV.measure_tver(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1)
        if(t1 == None or t2== None or t3 == None):
            store[fileCount,j,12] = -99
            continue
        
            
        store[fileCount,j,31:38] = t1, t2, t3,  t4, t5, t6, t7
        store[fileCount,j,12] = 1

            
        
# =============================================================================
# df_source = pd.DataFrame(store,  ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'e1', 'e2', 'bkg', 'size', 
#                                                               'xx', 'yy', 'xy','x', 'y','force_flag', 'vert_flag', 'bad_flag', 'back_sky', 'seeing', 'zp', 'fwhm', 'mjd' , 'airmass', 'mphase', 'mAngle' ,'expTime' ,'focus' ,'zp_n,skymag' ,'depth' ,'mRa' ,'mDec',
#                                                               'Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty', ])
# 
#     
# df_source.to_pickle('/home/dutta26/codes/frameSc_' + str(band)+'_.pk1')
# 
# =============================================================================
#np.save('/home/dutta26/codes/singleFrame_' +str(band)+'.npy', store)
np.save('/scratch/halstead/d/dutta26/abell_2390/test4t_' +str(band)+'.npy', store)
    