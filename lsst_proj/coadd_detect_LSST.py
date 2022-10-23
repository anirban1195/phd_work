#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 16:16:42 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np 
import os,sys
import subprocess
import gzip
import matplotlib.pyplot as plt
import shutil
from astropy.stats import sigma_clipped_stats
import helper
import measure_pythonV
import pandas as pd


swarpLoc = '/home/dutta26/apps/bin/bin/'
sexLoc = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/'

band ='ir'

if('ir' not in band):
    #First find flux to magnitude conversion i.e zp 
    raArr1 =np.arange(237.0083-(5*0.00555), 237.0083+(5*0.00555), 0.00555 )
    decArr1= np.arange(-15.2308 -(5*0.00555) ,  -15.2308 + (5*0.00555) , 0.005551) 
    raArr=[]
    decArr=[]
    magArr=[]
    mag = 16
    for j in range(len(raArr1)):
        mag = mag+1
        if(mag > 20):
            break
        for k in range(len(decArr1)):
            raArr.append(raArr1[j])
            decArr.append(decArr1[k])
            magArr.append(mag)
    
            
    coadd_file_flat = '/scratch/halstead/d/dutta26/lsst/filter'+band+ '/coadd_flat_all.fits'
    xList, yList = helper.convertToXY(raArr, decArr, coadd_file_flat)
    f=fits.open(coadd_file_flat)
    data=f[0].data
    f.close()
    fluxArr=[]
    for j in range(len(xList)):
        x = int(xList[j])
        y = int(yList[j])
        cut = data[y-20: y+20, x-20: x+20]
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
        fluxArr.append(flux)
        
    zpArr = magArr + 2.5*np.log10(fluxArr)

    zp = sigma_clipped_stats(zpArr)[1]
else:
    zp = 31



f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp_lsst.cat')
cat = f.readlines()
f.close()
raArr =[]
decArr=[]
fluxArr=[]
sex_xx=[]
sex_yy=[]
sex_xy=[]
for j in range(len(cat)):
    if((cat[j].split()[0]) == '#'):
     continue
    
    
    raArr.append(float(cat[j].split()[5])) 
    decArr.append(float(cat[j].split()[6]))
    fluxArr.append(float(cat[j].split()[1]))
    sex_xx.append(float(cat[j].split()[7]))
    sex_yy.append(float(cat[j].split()[8]))
    sex_xy.append(float(cat[j].split()[9]))
    
#zArr = np.load('/scratch/halstead/d/dutta26/lsst/lsst_zArr.npy')
#zArr=np.zeros(len(decArr), dtype=np.float32)

raArr = np.array(raArr)
decArr = np.array(decArr)
# =============================================================================
# min_dist_arr=[]
# for j in range(len(raArr)):
#     ra = raArr[j]
#     dec = decArr[j]
#     #print(j)
#     tem_dist = np.sqrt((raArr - ra) **2 + (decArr-dec)**2)
#     temp_dist = np.sort(tem_dist)
#     min_dist_arr.append(temp_dist[1]*3600/0.2)
# np.save('/scratch/halstead/d/dutta26/lsst/min_dist_arr.npy', min_dist_arr)
# sys.exit()
# =============================================================================
min_dist_arr = np.load('/scratch/halstead/d/dutta26/lsst/min_dist_arr.npy')


#Now measure sources 

if(band == 'ir'):
    coadd_file = '/scratch/halstead/d/dutta26/lsst/coadd1_ir.fits'
else:
    coadd_file = '/scratch/halstead/d/dutta26/lsst/filter'+band+ '/coadd_all.fits'
xList, yList = helper.convertToXY(raArr, decArr, coadd_file)
store = np.zeros((len(xList), 50), dtype = np.float32)
f=fits.open(coadd_file)
data= np.array(f[0].data)
f.close()
ySize,xSize = np.shape(data)
cut = 0
for j in range(len(xList)):
    print (j)
    v_flag = b_flag= 0
    ra = raArr[j]
    dec = decArr[j]
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    size = np.sqrt(sex_xx[j] + sex_yy[j])
    if(size<4):
        size = 3.9
    
    
    if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
        continue
    
    del cut
    cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
    
    flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut)
    
    
    
    
    if(flux == None):
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(sex_xx[j]) , np.sqrt(sex_yy[j]) , sex_xy[j], 1)
        if(flux == None or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
            store[j][0:2] = ra,dec
            store[j][10] = x
            store[j][11] = y
            continue
        else:
            print ('aa')
            store[j][0:15] = ra,dec, 0, flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 3, v_flag, b_flag
            
    else:
        e1 = (sigxx - sigyy)/(sigxx+sigyy)
        e2 = 2*sigxy/(sigxx+sigyy)
        ellip = np.sqrt(e1**2 +e2**2)
        star_bool = 1
        if(ellip < 0.1):
            star_bool = 1
        else:
            star_bool = 0
        store[j][0:15] = ra, dec, star_bool, flux, mux, muy, bkg, sigxx, sigyy, sigxy, x, y, 0, v_flag, b_flag
        
        continue

# =============================================================================
# df_source = pd.DataFrame(store,  columns = ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'bkg', 'xx', 'yy', 'xy','x', 'y', 'force_flag',  'vert_flag', 
#  'bad_flag',  'back_sky',  'seeing',  'zp',  'fwhm',  'mjd' ,  'airmass',  'mphase',  'mAngle' , 'expTime' ,
#  'focus' , 'zp_n' , 'skymag' , 'depth' , 'mRa' , 'mDec', 'scale_fact', 'flux_sf', 'mux_sf', 'muy_sf',
#  'bkg_sf', 'xx_sf', 'yy_sf', 'xy_sf', 'interp_xx', 'interp_yy', 'interp_xy', 'z', 'Empty', 'Empty',
#  'Empty', 'Empty', 'Empty', 'Empty', 'Empty',  'Empty'])
# 
# 
# df_source.to_pickle('/scratch/halstead/d/dutta26/lsst/df1_'+band+'.pk1')
# sys.exit()
# =============================================================================

star_arr = store[(np.where((store[:,2] == 1) & (store[:,3] > 10000) & (min_dist_arr >15) ))[0],  : ]
#Now tuse k sigma clip to find usable stars. Just do for sigxx
mean,median, std = sigma_clipped_stats(star_arr[:, 7])
print (mean,median, std)

star_arr = store[(np.where((store[:,2] == 1) & (min_dist_arr >15) &
                                      (store[:,7] >= mean-5*std) &
                                      (store[:,7] <= mean+5*std) &
                                      (store[:,3] > 20000) & (store[:,3] < 100000)))[0],  : ]

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
    ra = raArr[j]
    dec = decArr[j]
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
    
    

    
    
    store[j,38] = avgSigxx
    store[j,39] = avgSigyy
    store[j,40] = avgSigxy
    






df_source = pd.DataFrame(store,  columns = ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'bkg', 'xx', 'yy', 'xy','x', 'y', 'force_flag',  'vert_flag', 
 'bad_flag',  'back_sky',  'seeing',  'zp',  'fwhm',  'mjd' ,  'airmass',  'mphase',  'mAngle' , 'expTime' ,
 'focus' , 'zp_n' , 'skymag' , 'depth' , 'mRa' , 'mDec', 'scale_fact', 'flux_sf', 'mux_sf', 'muy_sf',
 'bkg_sf', 'xx_sf', 'yy_sf', 'xy_sf', 'interp_xx', 'interp_yy', 'interp_xy', 'z', 'Empty', 'Empty',
 'Empty', 'Empty', 'Empty', 'Empty', 'Empty',  'Empty'])


df_source.to_pickle('/scratch/halstead/d/dutta26/lsst/df1_'+band+'.pk1')

# =============================================================================
# f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp_lsst.cat')
# content = f.readlines()
# f.close()
# 
# for j in range(len(content)):
#     if(content[j].split()[0] == '#'):
#         continue
#     ra = float(content[j].split()[5])
#     dec = float(content[j].split()[6])
#     #if(ra > 237.275 or ra < 236.733):
#     #    continue
#     #if(dec > -15.4833 or dec < -15.9800):
#     #    continue
#     source_ra.append(ra)
#     source_dec.append(dec)
#     
# =============================================================================
# =============================================================================
# star_ra = np.array(star_ra)
# star_dec = np.array(star_dec)
# 
# 
# for folder in os.listdir(folderLoc):
#     
#     
#     bashCommand = './sex '+ folderLoc+folder 
#     process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
#     output, error = process.communicate()
#     
#     f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp.cat')
#     cat = f.readlines()
#     f.close()
#     raArr =[]
#     decArr=[]
#     fluxArr=[]
#     for j in range(len(cat)):
#         if((cat[j].split()[0]) == '#'):
#          continue
#         
#         if(float(cat[j].split()[1]) < 1000):
#             continue
#         raArr.append(float(cat[j].split()[5])) 
#         decArr.append(float(cat[j].split()[6]))
#         fluxArr.append(float(cat[j].split()[1]))
#     
#     fluxArr=np.array(fluxArr)
#     raArr = np.array(raArr)
#     decArr = np.array(decArr)    
#     if(len(fluxArr)> 80):
#         topIndices = np.argpartition(fluxArr, -80)[-80:]
#     elif(len(fluxArr) > 20):
#         topLen = len(fluxArr) - 2
#         topIndices = np.argpartition(fluxArr, -topLen)[-topLen:]
#     else:
#         continue
#         
#     raList = raArr[topIndices]
#     decList = decArr[topIndices]
#     #sys.exit()
#     
#     raShift=[]
#     decShift=[]
#     distShift=[]
#     #Now match indices 
#     for j in range(len(raList)):
#         ra = raList[j]
#         dec = decList[j]
#         delRa = (star_ra - ra)
#         delDec = (star_dec - dec)
#         dist = np.sqrt(delRa**2 + delDec**2)
#         loc = np.where(dist <0.001000 )[0]
#         print (len(loc), loc)
#         for l in range(len(loc)):
#             raShift.append(delRa[loc][l])
#             decShift.append(delDec[loc][l])
#             distShift.append(dist[loc][l])
#             
# # =============================================================================
# #         for l in range(len(delRa)):
# #             raShift.append(delRa[l])
# #             decShift.append(delDec[l])
# #             distShift.append(dist[l])
# # =============================================================================
#     print (folder)    
#     avgRaShift =  sigma_clipped_stats(raShift)[1] 
#     avgDecShift = sigma_clipped_stats(decShift)[1] 
#     #sys.exit()
#     print (sigma_clipped_stats(raShift))
#     print (sigma_clipped_stats(decShift))
#     shutil.copy(folderLoc+folder, destfolderLoc+folder)
#     f=fits.open(destfolderLoc+folder)
#     hdr = (f[0].header)
#     f.close()
#     
#     
#     print (hdr['CRVAL1'],avgRaShift, hdr['CRPIX1'])
#     fits.setval(destfolderLoc+folder, 'CRVAL1', value=hdr['CRVAL1']+avgRaShift)
#     fits.setval(destfolderLoc+folder, 'CRVAL2', value=hdr['CRVAL2']+avgDecShift)
#     #sys.exit()
# 
# 
# 
# =============================================================================
