#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 11:07:35 2021

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
from scipy.signal import convolve2d as conv2
from skimage import color, data, restoration

band = 'r'
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

#Read band coadd catalog
ir_coadd_df = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
ir_coadd_xx = np.array(ir_coadd_df['xx'])
ir_coadd_yy = np.array(ir_coadd_df['yy'])
ir_coadd_xy = np.array(ir_coadd_df['xy'])


data_width = int(len(os.listdir('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'))/2.0) +1
store = np.zeros((data_width, len(raList), 50), dtype = np.float32)
fileCount = -1


low = [71, 1130, 3129, 3714, 4340, 4853, 4915, 6378, 6484, 9182, 11466, 
       11756, 12283, 13170, 13424, 13445, 13948, 15470, 16716, 17423, 
       18744, 18797, 21231, 21529, 21663, 23015, 23488, 23876, 24338,
       24962, 26647, 27633, 30153, 30508, 30699, 32003, 32224, 32355,
       35333, 35751, 36680, 37094, 38237, 39252, 39304, 39355, 40729, 
       40968, 41355, 41484, 42878, 43335, 43965, 44003, 45147, 45180,
       45296, 45736, 46020, 46652, 46653, 46735, 47147, 48130, 49555, 
       51223, 51866, 51893, 52763]

high = [816, 1114, 1715, 2532, 2788, 3270, 4693, 5023, 6112, 6713, 6883, 
        6941, 8422, 9534, 9544, 9724, 10149, 11820, 13548, 13753, 15143, 
        15626, 15665, 16062, 17092, 18289, 19246, 19495, 20774, 22774, 25733,
        26646, 26818, 26954, 27351, 27998, 29052, 29111, 30071, 31135, 
        31467, 32428, 32945, 33696, 34379, 34442, 35245, 35688, 36077, 
        37757, 38455, 38666, 39839, 39910, 40130, 40345, 40509, 42285, 
        42590, 42677, 42710, 42790, 43094, 43166, 43703, 44007, 44325, 
        45503, 45755, 45891, 46264, 46452, 46573, 46922, 47159, 47778, 
        47920, 48334, 48590, 49177, 49966, 50630, 51196, 51383]

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
    print (len(xList))
    
    #Read the file data
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'+file)
    back_h = float((f[0].header)['SKYBG'])
    seeing = float((f[0].header)['SEEING'])
    zp = float((f[0].header)['MAGZERO'])
    fwhm = float((f[0].header)['FWHM_FLT'])
    mjd = float((f[0].header)['MJD-MID'])
    airmass = float((f[0].header)['AIRMASS'])
    mphase = float((f[0].header)['MOONPHSE'])
    mAngle = float((f[0].header)['MOON_D'])
    expTime = float((f[0].header)['EXPTIME'])  
    focus = float((f[0].header)['TELFOCUS'])
    zp_n = float((f[0].header)['PHOTZP_N'])
    skymag = float((f[0].header)['SKYMAG']) 
    depth = float((f[0].header)['PHOTDPTH']) 
    mRa = float((f[0].header)['MOON_RA']) 
    mDec = float((f[0].header)['MOON_DEC']) 
    f.close()
    
    
    ySize,xSize = np.shape(data)
    sizeList =[]
    #Fist find sources and measure params 
    cnt = cnt1 = cnt2 =cnt3 =0
    for j in range(len(xList)):
        if(j not in high):
            continue
        
        
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        if(band_coadd_xx[j] == 0 or band_coadd_yy[j]==0 or np.isnan(band_coadd_xx[j]) or np.isnan(band_coadd_yy[j])):
            
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
        
        
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
        #Measurement failed. Try with guess values
        if(flux == None):
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*band_coadd_xx[j]), np.sqrt(2*band_coadd_yy[j]), 2*band_coadd_xy[j], 1)
            if(flux == None or e1 == None or e2==None):
                continue
            force_flag = 1
        
        store[fileCount,j,0:33] = ra, dec, star_bool[j], flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy, x, y, force_flag, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n,skymag ,depth ,mRa ,mDec 
        #hdu = fits.PrimaryHDU(cut)
        #print (fileCount)
        #name = str(file) + '_'+ str(star_bool[j])+'_'+ str(force_flag) +'_'+ str(v_flag) + '_'+ str(b_flag)+ '_'+str(flux)[0:4]+ '_'+str(j)+'_'+str(cat)+ 'we'
        #hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj2/'+name+'.fits', overwrite=True)
    #break
        
        
    arr_sig_xx=store[fileCount, high, 10]
    arr_sig_yy=store[fileCount, high, 11]
    arr_sig_xy=store[fileCount, high, 12]
    
    arr_sig_xx = arr_sig_xx[arr_sig_xx != 0]
    arr_sig_yy = arr_sig_yy[arr_sig_yy != 0]
    arr_sig_xy = arr_sig_xy[arr_sig_xy != 0]
    
    avg_sigxx = np.nanmedian(arr_sig_xx)
    avg_sigyy = np.nanmedian(arr_sig_yy)
    avg_sigxy = np.nanmedian(arr_sig_xy)
    print (avg_sigxx, avg_sigyy, avg_sigxy)
    
    for j in range(len(xList)):
        
        if(j not in low):
            continue
        
        
        v_flag = b_flag = force_flag= 0
        ra = raList[j]
        dec = decList[j]
        x = int(round( xList[j]))
        y = int(round(yList[j]))
        size = 6
        
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        #If star then check if has weird lines through center or nans
        if(star_bool[j] == 1 ):
            v_flag = helper1.vert_stripe(cut)
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        
        
        
        
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*avg_sigxx), np.sqrt(2*avg_sigyy), 2*avg_sigxy, 1)
        if(flux == None or e1 == None or e2==None):
            continue
        force_flag = 0
        
        store[fileCount,j,0:33] = ra, dec, star_bool[j], flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy, x, y, force_flag, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n,skymag ,depth ,mRa ,mDec 
           
    #break
# =============================================================================
# df_source = pd.DataFrame(store,  ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'e1', 'e2', 'bkg', 'size', 
#                                                               'xx', 'yy', 'xy','x', 'y','force_flag', 'vert_flag', 'bad_flag', 'back_sky', 'seeing', 'zp', 'fwhm', 'mjd' , 'airmass', 'mphase', 'mAngle' ,'expTime' ,'focus' ,'zp_n,skymag' ,'depth' ,'mRa' ,'mDec',
#                                                               'Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty', ])
# 
#     
# df_source.to_pickle('/home/dutta26/codes/frameSc_' + str(band)+'_.pk1')
# 
# =============================================================================
np.save('/scratch/halstead/d/dutta26/abell_2390/test1_' +str(band)+'.npy', store)
    