#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 16:03:57 2022

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


swarpLoc = '/home/dutta26/apps/bin/bin/'
sexLoc = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/'
folderLoc ='/scratch/halstead/d/dutta26/abell_2390/lsst/'
fileList_name = '/home/dutta26/temp.ascii'
#Run swarp
fileArr =[]
for subfolder in os.listdir(folderLoc):
    

    
    
# =============================================================================
#     f=open(fileList_name, 'w+')
#     for files in os.listdir(folderLoc+subfolder):
#         if('C' in files or 'temp' in files):
#             continue
#         
#         if('gz' in files):
#             with gzip.open(folderLoc+subfolder+'/'+files, 'rb') as f_in:
#                 with open(folderLoc+subfolder+'/'+files[:-3], 'wb') as f_out:
#                     shutil.copyfileobj(f_in, f_out)
#             continue
#         f.write(folderLoc+subfolder+'/'+files + '\n')
#     f.close()
#     imageout = folderLoc+subfolder+'/temp.fits'
#     wtout=folderLoc+subfolder+'/temp.weight.fits'
#     bashCommand = './swarp @'+ fileList_name + ' -c /home/dutta26/default1.swarp -IMAGEOUT_NAME '+imageout+' -WEIGHTOUT_NAME '+wtout    
#     process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
#     output, error = process.communicate()
# =============================================================================
    fileArr.append(folderLoc+subfolder)


#Run sextractor
combined_ra =[]
combined_dec=[]
for folder in fileArr:
    
    bashCommand = './sex '+ folder+ '/temp.fits -WEIGHT_IMAGE '+folder+'/temp.weight.fits'  
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=sexLoc)
    output, error = process.communicate()
    
    f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/temp.cat')
    cat = f.readlines()
    f.close()
    raArr =[]
    decArr=[]
    fluxArr=[]
    for j in range(len(cat)):
        if((cat[j].split()[0]) == '#'):
         continue
        
        if(float(cat[j].split()[1]) < 100000):
            continue
        raArr.append(float(cat[j].split()[5])) 
        decArr.append(float(cat[j].split()[6]))
        fluxArr.append(float(cat[j].split()[1]))
    
    fluxArr=np.array(fluxArr)
    raArr = np.array(raArr)
    decArr = np.array(decArr)    
    topIndices = np.argpartition(fluxArr, -100)[-100:]
        
    combined_ra.append(raArr[topIndices])
    combined_dec.append(decArr[topIndices])
    #sys.exit()
    
combined_ra = np.array(combined_ra)
combined_dec= np.array(combined_dec)
np.save('/scratch/halstead/d/dutta26/abell_2390/ref_ra.npy', combined_ra[0])
np.save('/scratch/halstead/d/dutta26/abell_2390/ref_dec.npy', combined_dec[0])
np.save('/scratch/halstead/d/dutta26/abell_2390/comb_ra.npy', combined_ra)
np.save('/scratch/halstead/d/dutta26/abell_2390/comb_dec.npy', combined_dec)
#sys.exit()
#Perform match 
refCatRa = np.load('/scratch/halstead/d/dutta26/abell_2390/ref_ra.npy')
refCatDec = np.load('/scratch/halstead/d/dutta26/abell_2390/ref_dec.npy')
#combined_ra1 = np.load('/scratch/halstead/d/dutta26/abell_2390/comb_ra.npy')
#combined_dec = np.load('/scratch/halstead/d/dutta26/abell_2390/comb_dec.npy')
a = len(combined_dec)
finalRaShiftArr=[]
finalDecShiftArr=[]

tempArr_ra =[]
tempArr_dec =[]
for j in range(a):
    raShift = []
    decShift = []
    distShift =[]
    for k in range(len(combined_dec[j])):
        ra = combined_ra[j][k]
        dec = combined_dec[j][k]
        delRa = (refCatRa - ra)
        delDec = (refCatDec - dec)
        dist = np.sqrt(delRa**2 + delDec**2)
        loc = np.where(dist <0.001000 )[0]
        print (len(loc), loc)
        for l in range(len(loc)):
            raShift.append(delRa[loc][l])
            decShift.append(delDec[loc][l])
            distShift.append(dist[loc][l])
    
    raShift = np.array(raShift)
    decShift = np.array(decShift)
    distShift = np.array(distShift)
    htList = np.histogram(distShift, bins=10)[0]
    meanDist = np.histogram(distShift, bins=10)[1]
    maxHt = np.max(htList)
    index = np.where(htList== maxHt)[0][0]
    avgDist = (meanDist[index] +meanDist[index+1])*0.5
    
# =============================================================================
#     avgShiftRa = sigma_clipped_stats(raShift)[1]
#     avgShiftDec = sigma_clipped_stats(decShift)[1]
#     print (sigma_clipped_stats(raShift),  sigma_clipped_stats(decShift))
#     #sys.exit()
#     maxRa = max(1.2*avgShiftRa , 0.8*avgShiftRa)
#     minRa = min(1.2*avgShiftRa , 0.8*avgShiftRa)
#     
#     maxDec = max(1.2*avgShiftDec , 0.8*avgShiftDec)
#     minDec = min(1.2*avgShiftDec , 0.8*avgShiftDec)
#     loc = np.where( (raShift <= maxRa ) & (raShift >= minRa ) & (decShift <= maxDec) & (decShift >= minDec))
#     
# =============================================================================
    loc = np.where((distShift <= 1.05*avgDist) &(distShift >= 0.95*avgDist) )
    #if(j == 6):
    #    break
    if(len(loc[0]) ==0):
        finalRaShiftArr.append(0)
        finalDecShiftArr.append(0)
        continue
    finalRaShift = sigma_clipped_stats(raShift[loc])[1]
    finalDecShift = sigma_clipped_stats(decShift[loc])[1]  
# =============================================================================
#     n, bins, patches = plt.hist(x=raShift[loc], bins='auto', color='b',
#                             alpha=0.7, rwidth=0.85)
#     n, bins, patches = plt.hist(x=decShift[loc], bins='auto', color='r',
#                             alpha=0.7, rwidth=0.85)
# =============================================================================
    finalRaShiftArr.append(finalRaShift)
    finalDecShiftArr.append(finalDecShift)

#sys.exit()
#from astropy.wcs import WCS
#Correct wcs  
j = 0
for folder in fileArr:
    shutil.copy(folder+'/temp.fits','/scratch/halstead/d/dutta26/abell_2390/test'+str(j)+'.fits')
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/test'+str(j)+'.fits')
    hdr = (f[0].header)
    f.close()
    
    ra = finalRaShiftArr[j]
    dec = finalDecShiftArr[j]
    print (hdr['CRVAL1'],finalRaShiftArr[j], hdr['CRPIX1'])
    fits.setval('/scratch/halstead/d/dutta26/abell_2390/test'+str(j)+'.fits', 'CRVAL1', value=hdr['CRVAL1']+finalRaShiftArr[j])
    fits.setval('/scratch/halstead/d/dutta26/abell_2390/test'+str(j)+'.fits', 'CRVAL2', value=hdr['CRVAL2']+finalDecShiftArr[j])
    j += 1