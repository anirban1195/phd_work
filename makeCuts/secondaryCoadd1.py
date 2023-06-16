#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 19:52:13 2022

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
import subprocess
import measure_pythonV

#band = 'i'
band = str(sys.argv[1])
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
ir_coadd_data = np.array(ir_coadd_data)

band_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1'))

band_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test4t3_'+band+'.npy')
sf_loc = '/scratch/bell/dutta26/abell_2390/'+band +'/temp/'
fileList = os.listdir(sf_loc)

a,b,c = np.shape(band_sf_df)

sub_coadd_len_arr= [0, 9, 4, 0]
flux_bins = [0.3, 1.5, 2.5]
writeLoc = '/scratch/bell/dutta26/temp_'+band+'/'

for j in range(len(ir_coadd_data)):
#for j in [117]:
    print (j, '************************')
    #Rejects stars 
    if(ir_coadd_data[j,2] == 1 or ir_coadd_data[j,3] < 0.2 or  ir_coadd_data[j,3]> 2.5):
        continue
    if(ir_coadd_data[j,3] > 2 or ir_coadd_data[j,3]==0):
        continue
    
    flux_bin_index = np.digitize(ir_coadd_data[j,3], flux_bins)
    if(flux_bin_index<1 or flux_bin_index>2):
        continue
    sub_coadd_len = sub_coadd_len_arr[flux_bin_index]
    
    loc = np.where((band_sf_df[:,j,13] == 0) & (band_sf_df[:,j,14] == 0) & (band_sf_df[:,j,46] == 0) & 
                   (band_sf_df[:,j,38] > 0) & (band_sf_df[:,j,39] > 0) &  ( (band_sf_df[:,j,3] != 0) | (band_sf_df[:,j,31] != 0)))[0]
    
    
    masterArr= np.ones((a, 2), dtype = np.float32)* 99999
    count = 0
    for k in loc:
        if(band_sf_df[k,j,38] == -99 or ir_coadd_data [j,7]<0 or ir_coadd_data [j,8]<0 or ir_coadd_data [j,3]<0):
            continue
        
        #NOTE HERE THE fact we use psf for error estimation
        temp_xx = ir_coadd_data [j,7] - ir_coadd_data [j,38] + band_sf_df[k,j,38] 
        temp_yy = ir_coadd_data [j,8] - ir_coadd_data [j,39] + band_sf_df[k,j,39] 
        temp_xy = ir_coadd_data [j,9] - ir_coadd_data [j,40] + band_sf_df[k,j,40]
        
        #Use psf for error estimation 
        if(temp_xx< 0 or temp_yy< 0):
            temp_xx = band_sf_df[k,j,38]
            temp_yy = band_sf_df[k,j,39]
        area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
        N = band_coadd_df [j,3] / band_sf_df[k,j,30]
        B = band_sf_df[k,j,34]
        error = np.sqrt(1/N + (4*area*B)/N**2 )
        wt = 1/error
        
        
        if(np.isnan(error) or error == None or np.isinf(error) or error > 9999 or
           band_sf_df[k,j,10] < 15 or band_sf_df[k,j,11]<  15 or (temp_xx + temp_yy)<0):
            continue
        masterArr[k, 0 ] = error
        masterArr[k, 1 ] = k
        count += 1
        
    sortedIndices = masterArr[:,0].argsort()
    totSets = int(count / sub_coadd_len)
    leftOver = int(count % sub_coadd_len)
    
    #If too less viable cutouts then move on
    if(count < 20):
        continue
    
    current_set = 0
    while (current_set < totSets):
        startIndex = current_set*sub_coadd_len
        endIndex = (current_set+1)*sub_coadd_len - 1
        if((totSets - current_set ) == 1):
            endIndex = endIndex + leftOver
        
        #Make sure temp is clear 
        
        for files_temp in os.listdir(writeLoc):
                os.remove(writeLoc+files_temp)
        
        #Make cutouts
        
        seqNo = 0
        for k in np.arange(startIndex, endIndex+1):
            imageNo = sortedIndices[k]
            helper.make2Dsubcuts( band_sf_df[imageNo,j,:],  seqNo, writeLoc, sf_loc, band, j)
            seqNo += 1
        
        
        
        
        #Perform coadd 
        f= open('/home/dutta26/temp'+band+'.ascii', 'w+')
        for files in os.listdir(writeLoc):
            f.write(writeLoc + files + '\n')
        f.close()
        swarpLoc = '/home/dutta26/apps/bin/bin/'
        swarpCommand = './swarp @/home/dutta26/temp'+band+'.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/temp_'+band+'/coadd.fits  -WEIGHTOUT_NAME /scratch/bell/dutta26/temp_'+band+'/coadd.weight.fits -NTHREADS 1 -COMBINE_TYPE SUM'
        
        #Run swarp command 
        process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
        output, error = process.communicate()
        
       
        
        #Read the fits file 
        f=fits.open('/scratch/bell/dutta26/temp_'+band+'/coadd.fits')
        cut = np.array(f[0].data)
        f.close()
        
        
        #Find eff sigmaxx and sigmayy
        eff_psfsig_xx, eff_psfsig_yy, eff_psfsig_xy = helper.findEffSigmaxy(band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,38], 
                                                                   band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,39], 
                                                                   band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,40])
        err_psf_xx = np.sqrt(np.sum( band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,41] **2))
        err_psf_yy = np.sqrt(np.sum( band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,42] **2))
        err_psf_xy = np.sqrt(np.sum( band_sf_df[ sortedIndices[startIndex:endIndex+1] ,j,43] **2))
        
        
        if(ir_coadd_data[j,7]> 0 and ir_coadd_data[j,35] == 0):
            guess_xx = ir_coadd_data[j,7] -  ir_coadd_data[j,38] + eff_psfsig_xx
            guess_yy = ir_coadd_data[j,8] -  ir_coadd_data[j,39] + eff_psfsig_yy
            guess_xy = ir_coadd_data[j,9] -  ir_coadd_data[j,40] + eff_psfsig_xy
        elif(ir_coadd_data[j,7]> 0 and ir_coadd_data[j,35] >0):
            guess_xx = ir_coadd_data[j,35] -  ir_coadd_data[j,38] + eff_psfsig_xx
            guess_yy = ir_coadd_data[j,36] -  ir_coadd_data[j,39] + eff_psfsig_yy
            guess_xy = ir_coadd_data[j,37] -  ir_coadd_data[j,40] + eff_psfsig_xy
            
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1)
        
        
        #print (sigxx, sigyy ,sigxy)
        #Now run monte carlo if need be 
        corr_xx = sigxx - eff_psfsig_xx
        corr_yy = sigyy - eff_psfsig_yy
        corr_xy = sigxy - eff_psfsig_xy
        temp = corr_xx +corr_yy - 2*np.abs(corr_xy)
        if(corr_xx< 0 or corr_yy<0 or temp<0):
            area = 2*np.pi*np.sqrt(guess_xx*guess_yy - guess_xy**2)
            N = band_coadd_df[j,3]* np.sum(1/ band_sf_df[sortedIndices[startIndex:endIndex+1], j,30])
            B = np.sum(band_sf_df[sortedIndices[startIndex:endIndex+1], j,15])
            e1_measured = (guess_xx - guess_yy)/(guess_xx +guess_yy)
            e2_measured = 2*guess_xy/(guess_xx +guess_yy)
            
            s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * area
            s_xx = (1+e1_measured)*s
            s_yy = (1-e1_measured)*s
            s_xy = abs(e2_measured)*s
            
          
            
            corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                                 eff_psfsig_xx, eff_psfsig_yy, eff_psfsig_xy, 
                                                                 s_xx,s_yy,s_xy, 
                                                                 err_psf_xx, 
                                                                 err_psf_yy, 
                                                                 err_psf_xy)
            
            if(success<50):
                corr_xx = corr_yy = corr_xy = None
            
        #Now write the corrected sigmas to sf_df. Cannot add psf since we have eff psf here
        for index in sortedIndices[startIndex:endIndex+1]:
            if(corr_xx != None):
                band_sf_df[index, j , 50] = corr_xx 
                band_sf_df[index, j , 51] = corr_yy 
                band_sf_df[index, j , 52] = corr_xy 
            #print (corr_xx, corr_yy, corr_xy)
            
                
            
        
        #sys.exit()
        current_set += 1
    #sys.exit()
    
np.save('/scratch/bell/dutta26/abell_2390/test4t3_'+band+'_updated.npy', band_sf_df)    
    
    