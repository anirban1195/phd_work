#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 17:41:30 2021

@author: dutta26
"""


from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
from scipy.ndimage import rotate
import pandas as pd
import os
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess
from scipy.signal import convolve2d as conv2
from skimage import color, data, restoration

#band = str(sys.argv[1])
#source_df_name = str(sys.argv[2])
#coadd_file = str(sys.argv[4])

band = 'ir'
source_df_name = '/home/dutta26/codes/source_list.pk1'
coadd_file = '/scratch/halstead/d/dutta26/abell_2390/abell_'+band+'_coadd_wted.fits'
ir_coadd_df_name = '/home/dutta26/codes/coaddSc_ir_.pk1'


#Read the catalog 
source_df = pd.read_pickle(source_df_name)
f=fits.open(coadd_file)
data = np.array(f[0].data)
ySize,xSize = np.shape(data)
f.close()      
  
star_bool = np.array(source_df['star_bool'])
raList = np.array(source_df['ra'])
decList = np.array(source_df['dec'])
xList, yList = helper.convertToXY(raList, decList, coadd_file)

store = np.zeros((len(xList), 50), dtype = np.float32)
cut = 0
#Fist find stars and measure params 
for j in range(len(xList)):
    #print (j)
    v_flag = b_flag= 0
    ra = raList[j]
    dec = decList[j]
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    size = np.sqrt(source_df['sex_xx'][j] + source_df['sex_yy'][j])
    if(size<4):
        size = 4
        
    if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
        continue
    
    del cut
    cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
    
    
    #If star then check if has weird lines through center or nans
    if(star_bool[j] == 1 ):
        v_flag = helper1.vert_stripe(cut)
        
    
    
    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
    mag = 0
    #if Measurement failed
    if(flux == None):
        if('ir_coadd' in coadd_file):
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(source_df['sex_xx'][j]) , np.sqrt(source_df['sex_yy'][j]) , source_df['sex_xy'][j], 1)
            if(flux == None or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
                continue
            else:
                print ('aa')
                store[j][0:15] = ra, dec,star_bool[j], flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 3, v_flag, b_flag

        else:
            ir_coadd_df = pd.read_pickle(ir_coadd_df_name)
            if(ir_coadd_df['xx'][j] == 0):
                continue
            else:
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*ir_coadd_df['xx'][j]) , np.sqrt(2*ir_coadd_df['yy'][j]) , 2*ir_coadd_df['xy'][j], 1)
                if(flux == None or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
                    continue
                else:
                    
                    store[j][0:15] = ra, dec, star_bool[j], flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 2, v_flag, b_flag
    #If measurement successful           
    else:
        #if(flux > 0):
            #mag = 25 - 2.5*np.log10(flux/60) - 4.445378
        store[j][0:15] = ra, dec, star_bool[j], flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 0, v_flag, b_flag,
        continue
    
#Now do a second pass to find interpolated values 
#Find all good stars in the frame
star_arr = store[(np.where((store[:,2] == 1) & (store[:,12] == 0) & (store[:,13] == 0) & (store[:,14] == 0)))[0],  : ]
#Now tuse k sigma clip to find usable stars. Just do for sigxx
mean,median, std = sigma_clipped_stats(star_arr[:, 7])
print (mean,median, std)
# =============================================================================
# if(mean > 25 or std>3 or mean<2):
#     store[fileCount, :,38] = -99
#     store[fileCount, :,39] = -99
#     store[fileCount, :,40] = -99
#     continue
# =============================================================================
star_arr = store[(np.where((store[:,2] == 1) & 
                                      (store[:,7] >= mean-5*std) &
                                      (store[:,7] <= mean+5*std) &
                                      (store[:,12] == 0) & 
                                      (store[:,13] == 0) & 
                                      (store[:,14] == 0)))[0],  : ]

q,r = np.shape(star_arr)
star_temp = np.zeros(( q , 6)   , dtype = np.float32)
star_temp[:,0] = star_arr[:, 7]
star_temp[:,1] = star_arr[:, 8]
star_temp[:,2] = star_arr[:, 9]
star_temp[:,3] = star_arr[:, 10]
star_temp[:,4] = star_arr[:, 11]       
nStars = 10
for j in range(len(xList)):
    #print (j)
    v_flag = b_flag = force_flag= 0
    ra = raList[j]
    dec = decList[j]
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    
    temp = np.copy(star_temp)
    temp[:,5] = (temp[:,3]-x)**2 + (temp[:,4]-y)**2
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
    
    avgSigxx = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigyy = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigxy = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    
    
# =============================================================================
#     avgSigxx = np.mean( temp[0:nStars, 0] )
#     avgSigyy = np.mean( temp[0:nStars, 1] )
#     avgSigxy = np.mean( temp[0:nStars, 2] )
# =============================================================================
    
    
    store[j,38] = avgSigxx
    store[j,39] = avgSigyy
    store[j,40] = avgSigxy
    store[j,41] = np.std(temp[0:nStars, 0])
    store[j,42] = np.std(temp[0:nStars, 1])
    store[j,43] = np.std(temp[0:nStars, 2])
    
# =============================================================================
#     fact = 0
#     flag = 0
#     while(((avgSigxx <store[j,7])   or (avgSigyy<store[j,8] )) and fact<3 and store[j,7]>0 and flag==0 and store[j,7]>0):
#         fact = fact+1
#         size = np.sqrt(source_df['sex_xx'][j] + source_df['sex_yy'][j])
#         if(size<4):
#             size = 4
#             
#         if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
#             continue
#         
#         del cut
#         cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
#         
#         guessAx = np.sqrt(2*(store[j,38]+ fact*store[j,41]) ) 
#         guessAy = np.sqrt(2*(store[j,39]+ fact*store[j,42]) )
#         guessAxy = 2*(store[j,39] + fact*store[j,43])
#         flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, guessAx, guessAy, guessAxy, 1)
#         
#         if((sigxx > avgSigxx) and (sigyy>avgSigyy)):
#             store[j,44] = sigxx
#             store[j,45] = sigyy
#             store[j,46] = sigxy
#             flag = 1
#             break
# =============================================================================
       
        
# =============================================================================
#     if(j == 17407):
#         sys.exit()
# =============================================================================
    del temp
    
 

df_source = pd.DataFrame(store,  columns = ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'bkg', 'xx', 'yy', 'xy','x', 'y', 'force_flag',  'vert_flag', 
 'bad_flag',  'back_sky',  'seeing',  'zp',  'fwhm',  'mjd' ,  'airmass',  'mphase',  'mAngle' , 'expTime' ,
 'focus' , 'zp_n' , 'skymag' , 'depth' , 'mRa' , 'mDec', 'scale_fact', 'flux_sf', 'mux_sf', 'muy_sf',
 'bkg_sf', 'xx_sf', 'yy_sf', 'xy_sf', 'interp_xx', 'interp_yy', 'interp_xy', 'Empty', 'Empty', 'Empty',
 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',  'Empty'])


df_source.to_pickle('/home/dutta26/codes/coaddSc1_' + str(band)+'_.pk1')
