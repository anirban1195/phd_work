#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 17:10:22 2022

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
import wquantiles
import os
import helper

band = 'r'

ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
ir_coadd_data = np.array(ir_coadd_data)

band_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1'))

band_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test4t3_r.npy')
fileList = os.listdir('/scratch/bell/dutta26/abell_2390/'+band +'/temp/')
a,b,c = np.shape(band_sf_df)

for j in range(len(ir_coadd_data)):
    
    #Rejects stars 
    if(ir_coadd_data[j,2] == 1):
        continue
    if(ir_coadd_data[j,3] > 2 or ir_coadd_data[j,3]==0):
        continue
    
    
    loc = np.where((band_sf_df[:,j,13] == 0) & (band_sf_df[:,j,14] == 0) & (band_sf_df[:,j,46] == 0) & 
                   (band_sf_df[:,j,38] > 0) & (band_sf_df[:,j,39] > 0))[0]
    
    
    wtArr=[]
    loc1Arr=[]
    xArr=[]
    yArr=[]
    sizeArr=[]
    
    masterArr= np.ones((a, 15), dtype = np.float32)* 99999
    count = 0
    for k in loc:
        temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + band_sf_df[k,j,38] +1
        temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + band_sf_df[k,j,39] +1
        temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + band_sf_df[k,j,40]
        area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
        N = band_coadd_df [j,3] / band_sf_df[k,j,30]
        B = band_sf_df[k,j,34]
        error = np.sqrt(1/N + (4*area*B)/N**2 )
        wt = 1/error
        
        
        if(band_sf_df[k,j,10] < 15 or band_sf_df[k,j,11]<  15 or (temp_xx + temp_yy)<0):
            continue
        masterArr[k, 0 ] = error
        masterArr[k, 1 ] = band_sf_df[k,j,44]
        masterArr[k, 2 ] = band_sf_df[k,j,10]
        masterArr[k, 3 ] = band_sf_df[k,j,11]
        masterArr[k, 4 ] = np.sqrt(temp_xx + temp_yy)
        masterArr[k, 5 ] = band_sf_df[k,j,38]
        masterArr[k, 6 ] = band_sf_df[k,j,39]
        masterArr[k, 7 ] = band_sf_df[k,j,40]
        
        count += 1
        
    masterArr = masterArr[masterArr[:,0].argsort()]
    if(count <= 20):
        continue
    
    
    #Now find the files 
    while(k<len(masterArr)):
        
        #If cant form a set of 4 throw away
        if(masterArr[k,1] == 99999 or masterArr[k+1,1] == 99999 
           or masterArr[k+2,1] == 99999 or masterArr[k+3,1] == 99999):
            break
        #Make sure temp is clear 
        for files_temp in os.listdir('/scratch/bell/dutta26/temp/'):
                os.remove('/scratch/bell/dutta26/temp/'+files_temp)
        
        #Make the 4 cutouts
        for files in fileList:
            if(masterArr[k,1] in files):
                xPos = masterArr[k,2]
                yPos = masterArr[k,3]
                size = masterArr[k,4]
                helper.makeCutouts(band, files, 1, xPos, yPos, size)
                            
            
            if(masterArr[k+1,1] in files):  
                xPos = masterArr[k+1,2]
                yPos = masterArr[k+1,3]
                size = masterArr[k+1,4]
                helper.makeCutouts(band, files, 2, xPos, yPos, size)
                
            
            
            if(masterArr[k+2,1] in files):
                xPos = masterArr[k+2,2]
                yPos = masterArr[k+2,3]
                size = masterArr[k+2,4]
                helper.makeCutouts(band, files, 3, xPos, yPos, size)
                
                            
                                                
            if(masterArr[k+3,1] in files): 
                xPos = masterArr[k+3,2]
                yPos = masterArr[k+3,3]
                size = masterArr[k+3,4]
                helper.makeCutouts(band, files, 4, xPos, yPos, size)
                
        
        #Now make a file list and run swarp
        #Run swarp to make a good image
        f= open('/home/dutta26/temp.ascii', 'w+')
        f.write('/scratch/bell/dutta26/temp/cut1.fits \n')
        f.write('/scratch/bell/dutta26/temp/cut2.fits \n')
        f.write('/scratch/bell/dutta26/temp/cut3.fits \n')
        f.write('/scratch/bell/dutta26/temp/cut4.fits \n')
        f.close()
        
        swarpLoc = '/home/dutta26/apps/bin/bin/'
        
        swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/temp/coadd.fits  -WEIGHTOUT_NAME /scratch/bell/dutta26/temp/coadd.weight.fits -NTHREADS 1' 
        
        f=fits.open('/scratch/bell/dutta26/temp/coadd.fits')
        cut= np.array(f[0].data)
        f.close()
        
        #Now measure using forced 
        eff_psf_xx = 0.5 * (1/ ( ()**2 + ()**2 + ()**2 + ()**2  ))
        eff_psf_yy = 0.5 * (1/ ( ()**2 + ()**2 + ()**2 + ()**2  ))
        eff_psf_xy = 
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1)
        
        
        #If sigxx < eff psf then do monte carlo
        
        
        
        #Store the sigma values and weights
        
        
        
        
        k += 4
    sys.exit()
    
    
    
    
    
    
    
    
    
    
    
