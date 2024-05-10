#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:37:09 2024

@author: dutta26
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from scipy.ndimage import gaussian_filter
ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
master_frame=np.load('/scratch/bell/dutta26/abell_2390/master_arr_sf.npy')
zFile = '/home/dutta26/photz_eazy.zout'


f=fits.open('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits')
data = f[0].data
hdr =f[0].header
ySize,xSize = np.shape(data)
f.close() 
chopSize = 50
#Find the correct wcs
w = wcs.WCS(naxis=2)
w.wcs.crpix = [hdr['CRPIX1']/chopSize, hdr['CRPIX2']/chopSize]
w.wcs.cd = np.array([[hdr['CD1_1']*chopSize,hdr['CD1_2']], [hdr['CD2_1'], hdr['CD2_2']*chopSize]])
w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.cunit = [hdr['CUNIT1'], hdr['CUNIT2']]
#w.wcs.set_pv([(2, 1, 45.0)])
header = w.to_header()



if (zFile == None):
    redShiftArr = np.ones(len(ir_coadd_data))*9      
else:
    #Read Redshifts
    f=open(zFile)
    content = f.readlines()
    f.close()
    redShiftArr=[]
    for j in range(len(content)):
        
        if (content[j][0] == '#'):
            continue
        else:
            if('eazy' in zFile):
                if(float((content[j].split())[8]) >= 0.7 and float((content[j].split())[15])>=4):
                    redShiftArr.append(float((content[j].split())[7]))
                else:
                    redShiftArr.append(0)
            else:
                redShiftArr.append(float((content[j].split())[1]))
redShiftArr = np.array(redShiftArr)     





count = 1
zminArr =[0.01, 0.6]
zmaxArr =[0.6, 2.1]
n = len(zminArr)

matrix = np.zeros((n,n))
n_arr = np.zeros((n))
n = len(zminArr)
for j in range(n):
    zmin =zminArr[j]
    zmax = zmaxArr[j]
    print (zmin, zmax)
    cond = np.where((redShiftArr>zmin )& (redShiftArr<zmax) & (master_frame[:,2] !=0) &
                    (ir_coadd_data[:,2] == 0) &(master_frame[:,6] < 49) & 
                    (master_frame[:,7] < 49) & (ir_coadd_data[:,82] == 0) ) [0]
    print (len(cond))
    n_arr[j]=(len(cond))
    
for j in range(n):
    for k in range(n):
        if(k == (n-1)):
            matrix[j,k] = 1
        elif(k>= (n-1-j)):
            matrix[j,k] = (np.sum(n_arr[k+1:])) /np.sum(n_arr[k:])
        else:
            matrix[j,k] = 1
            
matrix = np.flip(matrix, axis =1 )          
inv_matrix = np.linalg.inv(matrix)

imgArr = np.zeros((719, 622, n))
for j in np.arange(n,0,-1):
    print (j)
    f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/EMode_sf'+str(j)+'.fits')
    if(j == 2):
        imgArr[:,:,n-j] = gaussian_filter(f[0].data, 0)
        f.close()
        continue
    elif(j == 1):
         imgArr[:,:,n-j] = gaussian_filter(f[0].data, 2)
         f.close() 
         continue
     
    elif(j == 3):
          imgArr[:,:,n-j] = gaussian_filter(f[0].data, 0)
          f.close() 
          continue
    imgArr[:,:,n-j] = f[0].data
    f.close()
    


for j in range(n):
    print (j)
    mulArr = inv_matrix[j,:]
    tot = np.zeros((719, 622), dtype = np.float32)
    print (mulArr)
    for k in range(n):
        loc = np.where(imgArr[:,:,k] < 0.014)
        imgArr[loc[0],loc[1],k] = 0
        tot += imgArr[:,:,k]*mulArr[k]
        
    hdu = fits.PrimaryHDU(tot,header=header)
    hdu.writeto('/scratch/bell/dutta26/abell_2390/zSlice/slices/Eslice_'+str(n-j)+'.fits', overwrite=True)
        

    


# =============================================================================
# count = 1
# zminArr =[0.01, 0.35, 0.42, 0.51, 0.61, 0.70, 0.80]
# zmaxArr =[0.35, 0.42, 0.51, 0.61, 0.70, 0.80, 8]
# matrix = np.zeros((7,7))
# n_arr = np.zeros((7))
# 
# for j in range(7):
#     zmin =zminArr[j]
#     zmax = zmaxArr[j]
#     print (zmin, zmax)
#     cond = np.where((redShiftArr>zmin )& (redShiftArr<zmax) & (master_frame[:,2] !=0) &
#                     (ir_coadd_data[:,2] == 0) &(master_frame[:,6] < 49) & 
#                     (master_frame[:,7] < 49)) [0]
#     print (len(cond))
#     n_arr[j]=(len(cond))
#     
# for j in range(7):
#     for k in range(7):
#         if(k == 6):
#             matrix[j,k] = 1
#         elif(k>= (6-j)):
#             matrix[j,k] = (np.sum(n_arr[k+1:])) /np.sum(n_arr[k:])
#         else:
#             matrix[j,k] = 1
#             
# matrix = np.flip(matrix, axis =1 )          
# inv_matrix = np.linalg.inv(matrix)
# 
# imgArr = np.zeros((719, 622, 7))
# for j in np.arange(7,0,-1):
#     f=fits.open('/scratch/bell/dutta26/abell_2390/zSlice/EMode_sf'+str(j)+'.fits')
#     imgArr[:,:,7-j] = f[0].data
#     f.close()
#     
# 
# 
# for j in range(7):
#     print (j)
#     mulArr = inv_matrix[j,:]
#     tot = np.zeros((719, 622), dtype = np.float32)
#     print (mulArr)
#     for k in range(7):
#         tot += imgArr[:,:,k]*mulArr[k]
#         
#     hdu = fits.PrimaryHDU(tot,header=header)
#     hdu.writeto('/scratch/bell/dutta26/abell_2390/zSlice/slices/Eslice_'+str(7-j)+'.fits', overwrite=True)
#         
# 
#     
# =============================================================================
