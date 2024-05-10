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
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
#band = str(sys.argv[1])
#source_df_name = str(sys.argv[2])
#coadd_file = str(sys.argv[4])

#f1 = fits.open('/scratch/bell/dutta26/abell_2390/temp.fits', mode = 'update')

#txt_file = open('/home/dutta26/sig_list_ir.txt', 'w+')

def detect(source_df_name, coadd_file, ir_coadd_df_name, band, outFile, plotLoc, bandLoc ):
    a11 =[]
    a22=[]
    a33=[]
    
    if('ir_coadd' not in coadd_file):
        ir_coadd_df = np.load(ir_coadd_df_name)
    
    #Read the catalog 
    source_df = pd.read_pickle(source_df_name)
    f=fits.open(coadd_file)
    data = f[0].data
    ySize,xSize = np.shape(data)
    f.close()      
      
    star_bool = np.array(source_df['star_bool'])
    raList = np.array(source_df['ra'])
    decList = np.array(source_df['dec'])
    xList, yList = helper.convertToXY(raList, decList, coadd_file)
    
    
    arr=[]
    
    #Find total background and scale factors
    totBkg = 0
    totScale = 0
    fileList1 =[]
    if(band == 'ir'):
        f=open('/home/dutta26/codes/makeWeights/fileList_swarp_r.ascii')
        content_r = f.readlines()
        f.close()
        for j in range(len(content_r)):
            fileList1.append(content_r[j][:-1])
            
    
        f=open('/home/dutta26/codes/makeWeights/fileList_swarp_i.ascii')
        content_i = f.readlines()
        f.close()
        for j in range(len(content_i)):
            fileList1.append(content_i[j][:-1])
    
    else:
   
        f=open('/home/dutta26/codes/makeWeights/fileList_swarp_'+str(band)+ '.ascii')
        content_band = f.readlines()
        f.close()
        for j in range(len(content_band)):
            fileList1.append(content_band[j][:-1])
    
    for files in fileList1:
        if('temp' in files or 'weight' in files):
            continue
        f=fits.open(files)
        back_h = float((f[0].header)['SKYBG'])
        flux_scale = float((f[0].header)['FLXSCALE']) 
        f.close()
        totBkg += back_h
        print (flux_scale)
        totScale += 1/flux_scale
        arr.append(flux_scale)
        
    
    
    #print (totBkg, totScale)
    #return 
    cnt = 0
    store = np.zeros((len(xList), 100), dtype = np.float32)
    cut = 0
    cutArr=np.zeros((len(xList)), dtype = np.float32)
    #Fist find stars and measure params 
    for j in range(len(xList)):
        #print (j)
        v_flag = b_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        size = np.sqrt(source_df['sex_xx'][j] + source_df['sex_yy'][j])*1.25
        if(size<4):
            size = 4
        if(size > 12.5):
            size = 12.5
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            store[j,57] = 1
            continue
        
        del cut
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        
        #If star then check if has weird lines through center or nans
        if(star_bool[j] == 1 ):
            v_flag = helper1.vert_stripe(cut)
        
        print (j)
        flux = None
        cutx , cuty = np.shape(cut)
        while(cutx>=24 and cuty>= 24):
            
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
            if(flux == None or np.isnan(flux)):
                cut = cut[2:-2, 2:-2]
                cutx , cuty = np.shape(cut)
            else:
                break
        mag = 0
        
       
        
        mag = 0
        if(flux == None or np.isnan(flux)):  
            store[j,56] = 1
            cnt += 1
            
            #f1[0].data[y-20: y+20, x-20] = 1000*j
            #f1[0].data[y-20: y+20, x+20] = 1000*j
            #f1[0].data[y-20, x-20:x+20] = 1000*j
            #f1[0].data[y+20, x-20:x+20] = 1000*j
        
        #If measurement successful           
        else:
            #if(flux > 0):
            #    mag = 25 - 2.5*np.log10(flux/60) - 4.445378
            store[j][0:15] = ra, dec, star_bool[j], flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag
            cutx , cuty = np.shape(cut)
            cutArr[j] = cutx
            continue
    
    
    #f1.flush()
    print (cnt)
    
    
    
    
    if('ir' in band):
        boundArr, flagArr, minArr, angleArr = helper.findBoundRadius(store)
        bkgFlagArr = helper.bkgFlagCalculate(data, store)
        store[:,82] = bkgFlagArr
        store[:,78]= boundArr
        store[:,79]= flagArr
        store[:,80]= minArr
        store[:,81]= angleArr
    else:
        store[:,78]= ir_coadd_df[:,78]
        store[:,79]= ir_coadd_df[:,79]
        store[:,80]= ir_coadd_df[:,80]
        store[:,81]= ir_coadd_df[:,81]
        store[:,82]= ir_coadd_df[:,82]
        boundArr, flagArr, minArr, angleArr = store[:,78], store[:,79], store[:,80], store[:,81]
        bkgFlagArr = store[:,82]
    
    #bkgFlagArr = helper.bkgFlagCalculate(data, store)
    #store[:,82] = bkgFlagArr
    #DO a seconds pass with bkg fixed for sources .
    for j in range(len(xList)):
        if( cutArr[j]<=0):
            continue
        size = cutArr[j]
        x = int(round( xList[j]) + 3*np.cos(angleArr[j]+np.pi) )
        y = int(round(yList[j]) +  3*np.sin(angleArr[j]+np.pi))
        cut = data[y-int(size /2): y+int(size /2), x-int(size /2): x+int(size /2)]
        
        
        if(flagArr[j] == 0):
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, 100, - 3*np.cos(angleArr[j]+np.pi), - 3*np.cos(angleArr[j]+np.pi) , 3,3,0 ,100,0)
            #print (bkg, store[j,6],sigxx,  store[j,7])
            
            if(flux == None or np.isnan(flux)):  
                store[j,56] = 1
                cnt += 1
                
            else:
                
                store[j][3:10] = flux, mux, muy, bkg, sigxx, sigyy, sigxy
        else:
        
            angle1 = angleArr[j]+np.pi/2
            angle2 = angleArr[j]-np.pi/2
            angle3 = angleArr[j]+np.pi
            x1,y1 = int(x+(size/2)*np.cos(angle1)), int(y+(size/2)*np.sin(angle1))
            x2,y2 = int(x+(size/2)*np.cos(angle2)), int(y+(size/2)*np.sin(angle2))
            x3,y3 = int(x+(size/2)*np.cos(angle3)), int(y+(size/2)*np.sin(angle3))
            temp = np.vstack((data[y1-2:y1+2, x1-2:x1+2], data[y2-2:y2+2, x2-2:x2+2], data[y3-2:y3+2, x3-2:x3+2], np.zeros((4,4))))
            back_calc = np.median(temp)
            
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, 100, - 3*np.cos(angleArr[j]+np.pi), - 3*np.cos(angleArr[j]+np.pi) , 3,3,0 ,100,0, 1, back_calc)
            #print (bkg, store[j,6],sigxx,  store[j,7])
            if(flux == None or np.isnan(flux)):  
                store[j,56] = 1
                cnt += 1
                
            else:
                
                store[j][3:10] = flux, mux, muy, bkg, sigxx, sigyy, sigxy
    
    
    
    
    
    
    
    
    #Now do a second pass to find interpolated values 
    #Find all good stars in the frame
    star_arr = store[(np.where((store[:,2] == 1) & (store[:,12] == 99) & (store[:,13] == 0) & (store[:,14] == 0)))[0],  : ]
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    
    print (mean,median, std)
    
    
    star_arr = store[(np.where((store[:,2] == 1) & 
                                          (store[:,7] >= mean-3*std) &
                                          (store[:,7] <= mean+3*std) &
                                          (store[:,8] >= mean1-3*std1) &
                                          (store[:,8] <= mean1+3*std1) &
                                          (store[:,9] >= mean2-3*std2) &
                                          (store[:,9] <= mean2+3*std2) &
                                          (store[:,3] < 5000) &
                                          (store[:,12] == 99) & 
                                          (store[:,13] == 0) & 
                                          (store[:,14] == 0) &
                                          (store[:,79] == 0)))[0],  : ]
    
    print (np.shape(star_arr))
    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    median_size = np.sqrt(median+median1)
    
    
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 10)   , dtype = np.float32)
    star_temp[:,0] = star_arr[:, 7]
    star_temp[:,1] = star_arr[:, 8]
    star_temp[:,2] = star_arr[:, 9]
    star_temp[:,3] = star_arr[:, 10]
    star_temp[:,4] = star_arr[:, 11]     
    star_temp[:,5] = star_arr[:, 3]       
    
    #global_xx_err, global_yy_err, global_xy_err = helper.fitGaussian(star_arr[:, 7], star_arr[:, 8], star_arr[:, 9], star_arr[:, 3]*totScale, totBkg, plotLoc, band)
    #print (global_xx_err, global_yy_err, global_xy_err, 'Global errors')
    
    #Specific for ir coadd
# =============================================================================
#     global_xx_err =  global_xx_err*1.5
#     global_yy_err =  global_yy_err*1.5
# =============================================================================
    global_xx_err =  (np.percentile(star_arr[:, 7], 84) - np.percentile(star_arr[:, 7], 16))/2
    global_yy_err =  (np.percentile(star_arr[:, 8], 84) - np.percentile(star_arr[:, 8], 16))/2
    global_xy_err =  (np.percentile(star_arr[:, 9], 84) - np.percentile(star_arr[:, 9], 16))/2
    print (global_xx_err, global_yy_err, global_xy_err, 'Global errors')
    helper.plotStarGaussian(star_arr[:, 7], star_arr[:, 8], star_arr[:, 9], star_arr[:, 3]*totScale, totBkg, plotLoc, band, global_xx_err, global_yy_err, global_xy_err)
    #print (np.percentile(star_arr[:, 7], 84) - np.percentile(star_arr[:, 7], 16))
    #return
    nStars = 15
    cnt = 0
    
    for j in range(len(xList)):
        #print (j)
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        
        temp = np.copy(star_temp)
        temp[:,9] = (temp[:,3]-x)**2 + (temp[:,4]-y)**2
        temp = temp[temp[:,9].argsort()]
        
        #Check if same star. Then delete the entry
        if(temp[0,9]<5):
            temp = np.delete(temp, (0), axis = 0)
        
        #Checking for nans to avoid code from crashing
        if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9]))):
            store[j,58] = 1
            continue
        if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9]))):
            store[j,58] = 1
            continue
        if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9]))):
            store[j,58] = 1
            continue
        
        avgSigxx = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9])
        avgSigyy = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9])
        avgSigxy = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 9])/np.sum(1/temp[0:nStars, 9])
        
        
    
        
        store[j,38] = avgSigxx
        store[j,39] = avgSigyy
        store[j,40] = avgSigxy
        #Find average poisson term 
        area = 2*np.pi*np.sqrt(store[j,38]*store[j,39] - store[j,40]**2)
        median_psf_flux = np.median(temp[0:nStars, 5])*totScale
        B = totBkg
        #s_psf = np.sqrt( (area/(np.pi*median_psf_flux) + 4*area**2 * B/(np.pi * median_psf_flux**2)) ) * np.sqrt(2*store[j,38] + 2*store[j,39])
        e1_psf = (avgSigxx - avgSigyy)/(avgSigxx + avgSigyy)
        e2_psf = 2*avgSigxy/(avgSigxx + avgSigyy)
        
        s_psf = helper.getPSFErr(temp[0:nStars, 0], temp[0:nStars, 1], temp[0:nStars, 2], temp[0:nStars, 5]*totScale, 1/temp[0:nStars, 9], B)
        
        store[j,41] = np.sqrt(s_psf**2 + global_xx_err**2 )
        store[j,42] = np.sqrt(s_psf**2 + global_xx_err**2 )
        store[j,43] = np.sqrt((0.707*s_psf)**2 + global_xy_err**2)
        
        store[j,74] = np.sqrt(s_psf**2 )
        store[j,75] = np.sqrt(s_psf**2 )
        store[j,76] =  np.sqrt(0.5*s_psf**2 )
        
        #Run all sources through monte carlo
        #Do not use where source measurement was not successful
        if(store[j,7] == 0 or store[j,3]<= 0 or np.isnan(store[j,7]) or np.isnan(store[j,8]) or store[j,7] == None or store[j,8]== None):
            continue
        #corr_xx = store[j,7] - avgSigxx
        #corr_yy = store[j,8] - avgSigyy
        #corr_xy = store[j,9] - avgSigxy
        #temp = corr_xx +corr_yy - 2*np.abs(corr_xy)
        #print (j)
        e1 = (store[j,7] - store[j,8])/(store[j,7] + store[j,8])
        e2 = 2 * store[j,9]/(store[j,7] + store[j,8])
        area = 2* np.pi*np.sqrt(store[j,7]*store[j,8] - store[j,9]**2)
        size = np.sqrt(store[j,7]+store[j,8])
        if(area< 0 or np.isnan(area) or size<np.sqrt(store[j,38]+store[j,39])): #Use psf area if area is nan HOPEFULLY good appx at least
            area = 2*np.pi*np.sqrt(store[j,38]*store[j,39] - store[j,40]**2)
        N = store[j,3]*totScale
        B = totBkg
        #print (totScale, totBkg, size)
        #if(corr_xx< 0 or corr_yy<0 or temp<0):
        s = np.sqrt(size**4/N + 4*size**6 *np.pi* B/N**2) 
        
        s_xx = np.sqrt((s)**2 )
        s_yy = np.sqrt((s)**2 )
        s_xy = np.sqrt((s*0.707)**2)
        
    
        
        store[j,68] = s_xx
        store[j,69] = s_yy
        store[j,70] = s_xy
        
        
        store[j,71] = s
        store[j,72] = s
        store[j,73] = 0.707*s
        
        corr_xx, corr_yy, corr_xy , success = helper.correct(store[j,7], store[j,8], store[j,9], 
                                                             store[j,38], store[j,39], store[j,40], 
                                                             s_xx,s_yy,s_xy, 
                                                             store[j,41], 
                                                             store[j,42], 
                                                             store[j,43])
        
        store[j,77] = success
        if(success< 25 or np.isnan(corr_xx) or corr_xx<0 or corr_xx == None or np.isnan(corr_yy) or corr_yy<0 or corr_yy == None):
            #Skips mostly stars
            if(store[j,2] != 1 and store[j, 3]< 10000 and store[j, 3]>0):
                cnt += 1
                a11.append(store[j,7] -store[j,38])
                a22.append(success)
                a33.append(store[j, 3])
                
  
            store[j,59] = 1
            continue
        
        else:
            #print (j,store[j,3], corr_xx, corr_yy,corr_xy )
            #cnt += 1
            store[j,35] = corr_xx + store[j,38]
            store[j,36] = corr_yy + store[j,39]
            store[j,37] = corr_xy + store[j,40]
                
        
    
    
    
    
        del temp
        
     
   
    np.save(outFile, store)
    print (outFile, cnt )
    
    flux = store[:,3]
    flux = flux[flux>0]
    flux = flux[flux<25]
    #plt.figure(figsize=(13,10))
    n, bins, patches = plt.hist(x=flux, bins='auto',histtype=u'step')
    plt.xlabel('Flux')
    plt.savefig(plotLoc+'flux_hist_'+band+'.png')
    plt.close()
    
    size= np.sqrt(store[:,35] -store[:,38] + store[:,36] - store[:,39])
    print (size)
    size = size[size>0]
    print (size)
    size = size[size<10]
    print (size)
    #plt.figure(figsize=(13,10))
    n, bins, patches = plt.hist(x=size, bins='auto',histtype=u'step')
                            
    plt.xlabel('PSF Corrected Size')
    plt.ylabel('Frequency')
    plt.savefig(plotLoc+'corr_size_hist_'+band+'.png')
    plt.close()
    
    loc = np.where(store[:,2] == 1)[0]
    size= np.sqrt(store[loc,35] -store[loc,38] + store[loc,36] - store[loc,39])
    print (size)
    size = size[size>0]
    print (size)
    size = size[size<3]
    #plt.figure(figsize=(13,10))
    n, bins, patches = plt.hist(x=size, bins='auto', histtype=u'step')
                            
    plt.xlabel('Star PSF Corrected Size')
    plt.ylabel('Frequency')
    plt.savefig(plotLoc+'star_corr_size_hist_'+band+'.png')
    plt.close()
    
    
    loc = np.where(store[:,2] == 0)[0]
    size= np.sqrt(store[loc,35] -store[loc,38] + store[loc,36] - store[loc,39])
    print (size)
    size = size[size>0]
    print (size)
    size = size[size<10]
    #plt.figure(figsize=(13,10))
    n, bins, patches = plt.hist(x=size, bins='auto',histtype=u'step')
    plt.xlabel('Galaxy PSF Corrected Size')
    plt.ylabel('Frequency')
    plt.savefig(plotLoc+'galaxy_corr_size_hist_'+band+'.png')
    plt.close()
    
    xx = store[:,35] -store[:,38]
    yy = store[:,36] - store[:,39]
    xy = store[:,37] - store[:,40]
    e1 = (xx-yy)/(xx+yy)
    e2 = 2*xy/(xx+yy)
    ellip = np.sqrt(e1**2+e2**2)
    #plt.figure(figsize=(13,10))
    loc = np.where((store[:,2] == 0) &(ellip>0) &(store[:,12] == 99) & (store[:,13] == 0) & (store[:,14] == 0) &(store[:,79] == 0))[0]
    n, bins, patches = plt.hist(x=e1[loc], bins='auto',histtype=u'step', color='r', label='e1')
    n, bins, patches = plt.hist(x=e2[loc], bins='auto',histtype=u'step', color='b', label='e2')
    n, bins, patches = plt.hist(x=ellip[loc], bins='auto',histtype=u'step', color='k', label='Tot Ellip')
    plt.legend()
    plt.xlabel('Galaxy e1/e2/Tot ellip')
    plt.savefig(plotLoc+'galaxy_ellip_hist_'+band+'.png')
    plt.close()
    
    
    xx = store[:,35] 
    yy = store[:,36] 
    xy = store[:,37] 
    e1 = (xx-yy)/(xx+yy)
    e2 = 2*xy/(xx+yy)
    ellip = np.sqrt(e1**2+e2**2)
    #plt.figure(figsize=(13,10))
    loc = np.where((store[:,2] == 1)& (ellip>0) &(store[:,12] == 99) & (store[:,13] == 0) & (store[:,14] == 0) &(store[:,79] == 0) )[0]
    n, bins, patches = plt.hist(x=e1[loc], bins='auto',histtype=u'step', color='r', label='e1')
    n, bins, patches = plt.hist(x=e2[loc], bins='auto',histtype=u'step', color='b', label='e2')
    n, bins, patches = plt.hist(x=ellip[loc], bins='auto',histtype=u'step', color='k', label='Tot Ellip')
    plt.legend()
    plt.xlabel('Star e1/e2/Tot ellip')
    plt.savefig(plotLoc+'star_ellip_hist_'+band+'.png')
    plt.close()


# =============================================================================
# source_df_name = str(sys.argv[1])
# coadd_file = str(sys.argv[2]) 
# ir_coadd_df_name = str(sys.argv[3])
# band = str(sys.argv[4]) 
# outFile = str(sys.argv[5]) 
# plotLoc= str(sys.argv[6]) 
# bandLoc = str(sys.argv[7])
# =============================================================================
band = 'ir'
source_df_name = '/home/dutta26/codes/source_list.pk1'
coadd_file = '/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.fits'
ir_coadd_df_name =  '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy'
outFile = '/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.npy'
plotLoc= '/scratch/bell/dutta26/abell_2390/plot/'
bandLoc = '/scratch/bell/dutta26/abell_2390/'
detect(source_df_name, coadd_file, ir_coadd_df_name, band, outFile, plotLoc, bandLoc )