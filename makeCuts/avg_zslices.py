#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:07:40 2024

@author: dutta26
"""

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

zminArr =[0.01, 0.35, 0.42, 0.51, 0.61, 0.70, 0.80]
zmaxArr =[0.35, 0.42, 0.51, 0.61, 0.70, 0.80, 8]
zminArr_subspaced = np.array([0.01, 0.3, 0.35, 0.38, 0.42, 0.47, 0.51, 0.55, 0.61, 0.65, 0.7, 0.75, 0.8, 0.85])

for j in range(len(zminArr)):
    if(j == 0):
        f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/temp/EMode_sf'+str(j+1)+'.fits')
        data = f[0].data
        f.close()
        hdu = fits.PrimaryHDU(data)  
        hdu.writeto('/scratch/bell/dutta26/abell_2390/zSlice/EMode_sf'+str(j+1)+'.fits', overwrite=True)
        continue
    loc = np.where(zminArr_subspaced == zminArr[j])[0][0]
    print (loc)
    imgArr = np.zeros((719, 622, 3), dtype = np.float32)
    f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/temp/EMode_sf'+str(j)+'.fits')
    imgArr[:,:,0] = f[0].data
    f.close()
    
    f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/temp/EMode_sf'+str(j+1)+'.fits')
    imgArr[:,:,1] = f[0].data
    f.close()
    
    f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/temp/EMode_sf'+str(j+2)+'.fits')
    imgArr[:,:,2] = f[0].data
    f.close()
    
    data = np.mean(imgArr[:,:,:], axis=2)
    hdu = fits.PrimaryHDU(data)  
    hdu.writeto('/scratch/bell/dutta26/abell_2390/zSlice/EMode_sf'+str(j+1)+'.fits', overwrite=True)
# =============================================================================
# imgArr = np.zeros((719, 622, 5))
# for j in np.arange(1,5):
#     print (j)
#     f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/temp/EMode_sf'+str(j)+'.fits')
#     imgArr[:,:,j] = f[0].data
#     f.close()
#     
# data = -imgArr[:,:,1] + np.mean(imgArr[:,:,2:5], axis=2)
# hdu = fits.PrimaryHDU(data)
# hdu.writeto('/scratch/bell/dutta26/abell_2390/zSlice/slices/Eslice_1.fits', overwrite=True)
# 
# =============================================================================
