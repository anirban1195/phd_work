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
a11 =[]
a22=[]
a33=[]
band = 'ir'
source_df_name = '/home/dutta26/codes/source_list.pk1'
coadd_file = '/scratch/bell/dutta26/abell_2390/abell_'+band+'_coadd.fits'
ir_coadd_df_name = '/home/dutta26/codes/coaddSc_ir.pk1'
if('ir_coadd' not in coadd_file):
    ir_coadd_df = pd.read_pickle(ir_coadd_df_name)

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
# =============================================================================
#     fileList = os.listdir('/scratch/bell/dutta26/abell_2390/i/') 
#     for files in fileList:
#         fileList1.append('/scratch/bell/dutta26/abell_2390/i/' + files)
#     fileList = os.listdir('/scratch/bell/dutta26/abell_2390/r/') 
#     for files in fileList:
#         fileList1.append('/scratch/bell/dutta26/abell_2390/r/' + files)
# =============================================================================
else:
# =============================================================================
#     fileList = os.listdir( '/scratch/bell/dutta26/abell_2390/'+band+'/')
#     for files in fileList:
#         fileList1.append('/scratch/bell/dutta26/abell_2390/'+band+'/' + files)
# =============================================================================
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
    



cnt = 0
store = np.zeros((len(xList), 77), dtype = np.float32)
cut = 0
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
    while(cutx>=30 and cuty>= 30):
        
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
        store[j][0:15] = ra, dec, star_bool[j], flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag,
        continue

#f1.flush()
print (cnt)
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
                                      (store[:,3] < 2000) &
                                      (store[:,12] == 99) & 
                                      (store[:,13] == 0) & 
                                      (store[:,14] == 0)))[0],  : ]

print (np.shape(star_arr))
mean,median, std = sigma_clipped_stats(star_arr[:, 7])
mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
median_size = np.sqrt(median+median1)


q,r = np.shape(star_arr)
star_temp = np.zeros(( q , 6)   , dtype = np.float32)
star_temp[:,0] = star_arr[:, 7]
star_temp[:,1] = star_arr[:, 8]
star_temp[:,2] = star_arr[:, 9]
star_temp[:,3] = star_arr[:, 10]
star_temp[:,4] = star_arr[:, 11]     
star_temp[:,5] = star_arr[:, 3]       
  
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
    temp[:,5] = (temp[:,3]-x)**2 + (temp[:,4]-y)**2
    temp = temp[temp[:,5].argsort()]
    
    #Check if same star. Then delete the entry
    if(temp[0,5]<5):
        temp = np.delete(temp, (0), axis = 0)
    
    #Checking for nans to avoid code from crashing
    if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        store[j,58] = 1
        continue
    if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        store[j,58] = 1
        continue
    if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
        store[j,58] = 1
        continue
    
    avgSigxx = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigyy = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    avgSigxy = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
    
    

    
    store[j,38] = avgSigxx
    store[j,39] = avgSigyy
    store[j,40] = avgSigxy
    #Find average poisson term 
    area = 2*np.pi*np.sqrt(store[j,38]*store[j,39] - store[j,40]**2)
    median_psf_flux = np.median(temp[0:nStars, 5])*totScale
    B = totBkg
    s_psf = np.sqrt( (area/(np.pi*median_psf_flux) + 4*area**2 * B/(np.pi * median_psf_flux**2)) ) * np.sqrt(2*store[j,38] + 2*store[j,39])
    e1_psf = (avgSigxx - avgSigyy)/(avgSigxx + avgSigyy)
    e2_psf = 2*avgSigxy/(avgSigxx + avgSigyy)
    store[j,41] = np.sqrt(((median_size/30)*median_size*1.414)**2 + ((1+e1_psf)*s_psf)**2)
    store[j,42] = np.sqrt(((median_size/30)*median_size*1.414)**2 + ((1-e1_psf)*s_psf)**2)
    store[j,43] = np.sqrt(((median_size/30)*median_size*1.414 * 0.707)**2 + (0.707*(1+abs(e2_psf))*s_psf)**2)
    
    store[j,74] = ((1+e1_psf)*s_psf)
    store[j,75] = ((1-e1_psf)*s_psf)
    store[j,76] = (0.707*(1+abs(e2_psf))*s_psf)
    
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
    if(area< 0 or np.isnan(area)): #Use psf area if area is nan HOPEFULLY good appx at least
        area = 2*np.pi*np.sqrt(store[j,38]*store[j,39] - store[j,40]**2)
    N = store[j,3]*totScale
    B = totBkg
    
    #if(corr_xx< 0 or corr_yy<0 or temp<0):
    s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * np.sqrt(2*store[j,7] + 2*store[j,8])
    
    s_xx = np.sqrt(((1+e1)*s)**2 + ((median_size/30)*median_size*1.414)**2)
    s_yy = np.sqrt(((1-e1)*s)**2 + ((median_size/30)*median_size*1.414)**2)
    s_xy = np.sqrt((0.707*(1+abs(e2))*s)**2 + ((median_size/30)*median_size*1.414 * 0.707)**2)
    
# =============================================================================
#     s1 = np.sqrt( (1.3*area/(np.pi*N) + 2.2*4*area**2 * B/(np.pi * N**2)) ) * 2*store[j,7] 
#     s2 = np.sqrt( (1.3*area/(np.pi*N) + 2.2*4*area**2 * B/(np.pi * N**2)) ) * 2*store[j,8] 
#     s3 = np.sqrt( (1.3*area/(np.pi*N) + 2.2*4*area**2 * B/(np.pi * N**2))* 2*store[j,7] *2*store[j,8]  ) 
#     s=s1
#     s_xx = np.sqrt((s1)**2 + ((median_size/30)*median_size*1.414)**2)
#     s_yy = np.sqrt((s2)**2 + ((median_size/30)*median_size*1.414)**2)
#     s_xy = np.sqrt((s3)**2 + ((median_size/30)*median_size*1.414 * 0.707)**2)
# =============================================================================
    
    store[j,68] = s_xx
    store[j,69] = s_yy
    store[j,70] = s_xy
    
    
    store[j,71] = (1+e1)*s
    store[j,72] = (1-e1)*s
    store[j,73] = 0.707*(1+abs(e2))*s
    
    corr_xx, corr_yy, corr_xy , success = helper.correct(store[j,7], store[j,8], store[j,9], 
                                                         store[j,38], store[j,39], store[j,40], 
                                                         s_xx,s_yy,s_xy, 
                                                         store[j,41], 
                                                         store[j,42], 
                                                         store[j,43])
    
    if(success< 25 or np.isnan(corr_xx) or corr_xx<0 or corr_xx == None or np.isnan(corr_yy) or corr_yy<0 or corr_yy == None):
        #Skips mostly stars
        if(store[j,2] != 1 and store[j, 3]< 10000 and store[j, 3]>0):
            cnt += 1
            a11.append(store[j,7] -store[j,38])
            a22.append(success)
            a33.append(store[j, 3])
            
# =============================================================================
#             print (s,(store[j,7] -store[j,38]), (store[j,8] -store[j,39]),  (store[j,9] -store[j,40]) )
#             print (s,(store[j,7] - store[j,38]), (store[j,8] -store[j,39]))
#             txt_file.write('********************************  \n')
#             txt_file.write('s = '+str(s)[0:6]+ ' success = '+str(success) + '\n')
#             txt_file.write('Measured xx = '+str(store[j,7])[0:6] +  ' PSF xx' + str( store[j,38])[0:6]+ '\n')
#             txt_file.write('Measured yy = '+str(store[j,8])[0:6] +  ' PSF yy' + str( store[j,39])[0:6]+ '\n')
#             txt_file.write('Measured xy = '+str(store[j,9])[0:6] +  ' PSF yy' + str( store[j,40])[0:6]+ '\n')
#             txt_file.write('Err PSF xx = '+str(store[j,41])[0:6] +  ' Err PSF yy' + str( store[j,42])[0:6]+ ' Err PSF xy' + str( store[j,43])[0:6]+ '\n')
#             txt_file.write('Err xx = '+str(s_xx)[0:6] +  ' Err yy' + str( s_yy)[0:6]+ ' Err xy ' + str(s_xy)[0:6]+ '\n')
#             txt_file.write('Diff xy = '+str(store[j,9]- store[j,40])[0:6]+ '\n')
#             txt_file.write('x Loc = ' + str(store[j,10])+'y Loc = '+ str(store[j,11]) + '\n')
#             txt_file.write('j  = '+str(j) + '\n')
#             xLoca = int(store[j,10])
#             yLoca = int(store[j,11])
#             f1[0].data[yLoca-20: yLoca+20, xLoca-20] = 1000*j
#             f1[0].data[yLoca-20: yLoca+20, xLoca+20] = 1000*j
#             f1[0].data[yLoca-20, xLoca-20:xLoca+20] = 1000*j
#             f1[0].data[yLoca-20, xLoca-20:xLoca+20] = 1000*j
# =============================================================================
        store[j,59] = 1
        continue
    
    else:
        #print (j,store[j,3], corr_xx, corr_yy,corr_xy )
        #cnt += 1
        store[j,35] = corr_xx + store[j,38]
        store[j,36] = corr_yy + store[j,39]
        store[j,37] = corr_xy + store[j,40]
            
    




    del temp
    
 
#f1.flush()
#txt_file.close()
df_source = pd.DataFrame(store,  columns = ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'bkg', 'xx', 'yy', 'xy','x', 'y', 'force_flag',  'vert_flag', 
 'bad_flag',  'back_sky',  'seeing',  'zp',  'fwhm',  'mjd' ,  'airmass',  'mphase',  'mAngle' , 'expTime' ,
 'focus' , 'zp_n' , 'skymag' , 'depth' , 'mRa' , 'mDec', 'scale_fact', 'flux_sf', 'mux_sf', 'muy_sf',
 'bkg_sf', 'xx_sf', 'yy_sf', 'xy_sf', 'interp_xx', 'interp_yy', 'interp_xy', 'psf_xx_std', 'psf_yy_std', 'psf_xy_std',
 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',  'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',
 'Empty', 'Empty', 'Empty', 'Empty', 'Empty','Empty', 'Empty', 'Empty', 'Empty', 'Empty',
 'Empty', 'Empty', 'Empty', 'Empty', 'Empty', 'Empty','Empty'])


df_source.to_pickle('/home/dutta26/codes/coaddSc_' + str(band)+'.pk1')
