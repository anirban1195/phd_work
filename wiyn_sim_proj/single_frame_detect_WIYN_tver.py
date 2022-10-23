#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 12:18:29 2022

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


folderLoc = "/scratch/halstead/d/dutta26/abell_2390/odi_img_realsee/"
swarpLoc = '/home/dutta26/apps/bin/bin/'

#Read the sextractor file to get ra and dec values
f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test_all.cat')
content = f.readlines()
f.close()
raList =[]
decList=[]
#Create a list of stars with ra and dec 
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue

    raList.append(float(content[j].split()[5])) 
    decList.append(float(content[j].split()[6])) 


coadd_df = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')
store = np.zeros((100,len(raList), 50), dtype = np.float32)
fileCount = -1
for files in os.listdir(folderLoc):
    print (files)
    fileCount += 1
    imgFile = folderLoc +files+'/final_coadd.fits'
    bkgFile = folderLoc +files+'/sample.fits'
    f=fits.open(bkgFile)
    data = np.array(f[0].data)
    f.close()
    bkg1 = sigma_clipped_stats(data)[1]
    
    f=fits.open(imgFile)
    data = np.array(f[0].data)
    f.close()
    
    
    
    ySize,xSize = np.shape(data)
    xList, yList = helper.convertToXY(raList, decList, imgFile)
    
    for j in range(len(xList)):
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(coadd_df[j,7] <= 0 or coadd_df[j,8]<=0 or np.isnan(coadd_df[j,7]) or np.isnan(coadd_df[j,8])):
            store[fileCount,j,0] = ra
            store[fileCount,j,1] = dec
            store[fileCount,j,2] = coadd_df[j,2]
            store[fileCount,j,10] = x
            store[fileCount,j,11] = y
            continue
        
        
        size = np.sqrt(coadd_df[j,7] + coadd_df[j,8])
        if(size<4):
            size = 3.9
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
    
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        
       
        t1,t2,t3,t4,t5,t6,t7 = measure_pythonV.measure_v2T(cut)
        mux = t3/t2
        muy = t4/t2
        flux = t2/t7
        sigxy = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy= (t6/t2) - (t4*t4) / (t2*t2)
        e1 = (sigxx - sigyy)/(sigxx + sigyy)
        e2 = 2*sigxy/(sigxx + sigyy)
        psf = np.sqrt(sigxx+sigyy)
        if(flux == None or flux<= 0 or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
            store[fileCount,j,0] = ra
            store[fileCount,j,1] = dec
            store[fileCount,j,2] = coadd_df[j,2]
            store[fileCount,j,10] = x
            store[fileCount,j,11] = y
            continue
        
        #if(b_flag == 1 and cnt1 > 5000):
        #    sys.exit()
        store[fileCount,j,0:19] = ra, dec, coadd_df[j,2], t1, t2, t3,  t4, t5, t6, t7, x, y, 99, 0, b_flag, sigxx, sigyy, sigxy,flux
    
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & (store[fileCount,:,12] == 99)  & (store[fileCount,:,14] == 0)))[0],  : ]
    mean,median, std = sigma_clipped_stats(star_arr[:, 15])
    print (sigma_clipped_stats(star_arr[:, 15]))
    threshold = 10000
    
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                          (store[fileCount,:,12] == 99)  & 
                                          (store[fileCount,:,18] > threshold)  & 
                                          (store[fileCount,:,14] == 0)))[0],  : ]
    mean,median, std = sigma_clipped_stats(star_arr[:, 15])
    print (sigma_clipped_stats(star_arr[:, 15]))
    
    if(mean > 25 or std>3 or mean<1):
        store[fileCount, :,38] = -99
        store[fileCount, :,39] = -99
        store[fileCount, :,40] = -99
        continue
    
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                          (store[fileCount,:,15] >= mean-5*std) &
                                          (store[fileCount,:,15] <= mean+5*std) &
                                          (store[fileCount,:,12] == 99) & 
                                          (store[fileCount,:,18] > threshold) &
                                          (store[fileCount,:,14] == 0)))[0],  : ]

    
    
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 6)   , dtype = np.float32)
    star_temp[:,0] = star_arr[:, 15]
    star_temp[:,1] = star_arr[:, 16]
    star_temp[:,2] = star_arr[:, 17]
    star_temp[:,3] = star_arr[:, 10]
    star_temp[:,4] = star_arr[:, 11]
    print (np.shape(star_arr))
    print (threshold, median)
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
        
        store[fileCount,j,41] = np.std(temp[0:nStars, 0])
        store[fileCount,j,42] = np.std(temp[0:nStars, 1])
        store[fileCount,j,43] = np.std(temp[0:nStars, 2])
        
        
        
        
    
        #Use guess measure only when actual convergence fails 
        if(store[fileCount,j,12] == 99 ):
            continue
        #Find guess shapes
        guess_xx = coadd_df[j,7]-coadd_df[j,38] + avgSigxx_frame
        guess_yy = coadd_df[j,8]-coadd_df[j,39] + avgSigyy_frame
        guess_xy = coadd_df[j,9]-coadd_df[j,40] + avgSigxy_frame
        if(guess_xx<0 or guess_yy<0):
            continue

        #Make cutout
        ra = raList[j]
        dec = decList[j]
        x = int(round(xList[j]))
        y = int(round(yList[j]))
        mux_guess = xList[j] - x
        muy_guess = yList[j] - y
        if(coadd_df[j,7] <= 0 or coadd_df[j,8]<=0 or np.isnan(coadd_df[j,7]) or np.isnan(coadd_df[j,8])):
            continue
        size = np.sqrt(guess_xx + guess_yy)
        if(size<4):
            size = 3.9
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        store[fileCount,j,2] =coadd_df[j,2]
        store[fileCount,j,13] = v_flag
        store[fileCount,j,14] = b_flag
       
        
        #Measure cutout
        t1,t2,t3,t4,t5,t6,t7 = measure_pythonV.measure_v2T(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1)
        mux = t3/t2
        muy = t4/t2
        flux = t2/t7
        sigxy = (t1/t2) - (t3*t4)/(t2*t2)
        sigxx = (t5/t2) - (t3*t3) / (t2*t2)
        sigyy= (t6/t2) - (t4*t4) / (t2*t2)
        e1 = (sigxx - sigyy)/(sigxx + sigyy)
        e2 = 2*sigxy/(sigxx + sigyy)
        psf = np.sqrt(sigxx+sigyy)
        
        #sys.exit()
        if(flux == None or e1== None or e2 == None or np.isnan(e1)):
            store[fileCount,j,12] = -99
            continue
        
        #print (j)    
        store[fileCount,j,31:38] = t1,t2,t3,t4,t5,t6,t7
        store[fileCount,j,12] = 1
    
    if(fileCount >= 100):
        break
    

np.save('/scratch/halstead/d/dutta26/abell_2390/test_10star_tver.npy', store)  