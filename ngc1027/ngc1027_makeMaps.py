#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 13:32:45 2020

@author: anirban
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import griddata

#This code creates residual PSFs and ellipticities
#First it boxes in a 1000x1000 sq the avg PSF, e1 and e2
#The uses linear bi interpolation to calculate PSF and e1 and e2 at any location 
#Read the data
catalog = '/home/dutta26/codes/ngc1027/1027_spatial_data'
with open(catalog) as f:
    content = f.readlines()

xInd = 0
yInd = 2
pInd = 4
e1Ind = 6
e2Ind = 8

xList=[]
yList =[]
pList= []
e1List=[]
e2List =[]

for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    xList.append(float(content[j].split()[xInd])) 
    yList.append( float(content[j].split()[yInd])) 
    pList.append(float(content[j].split()[pInd]))
    e1List.append(float(content[j].split()[e1Ind]))
    e2List.append(float(content[j].split()[e2Ind]))


#Open the original image and measure dimensions 
f=fits.open('/scratch/halstead/d/dutta26/ngc_1027/ngc_g_coadd.fits')
data = np.array(f[0].data)
f.close()

y, x=np.shape(data)
del data
x = round(x, -3)
y = round(y, -3)
xBins = int(x/1000)
yBins = int(y/1000)
binArr = np.zeros((100, yBins, xBins), dtype=np.float32)
binArrCounter  = np.zeros((yBins, xBins), dtype=np.int32)
#Put data in bins
for j in range(len(content)- 1):
    x = int(round( xList[j], -3)/1000)
    y = int(round( yList[j], -3)/1000)
    val = pList[j]
    #val = e1List[j]
    binArr[binArrCounter[y,x], y,x] += val
    binArrCounter[y,x] += 1
#Average each bin
img = np.zeros((yBins, xBins), dtype=np.float32)
for x in range(xBins):
    for y in range(yBins):
        if(len(np.where(np.abs(binArr[:, y,x])> 0 )[0]) >= 3):
            avg = sigma_clipped_stats(binArr[:, y,x], mask_value =0.0, sigma=3)[0]
            #print (avg)
            #avg = np.mean(binArr[:, y,x], mask_value =0.0, sigma=3)[0]
        else:
            avg=0
        img[y,x] = avg
#Write data
hdu = fits.PrimaryHDU(img) 
hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/test.fits', clobber=True)
    
    
#Open the original image and measure dimensions 
f=fits.open('/scratch/halstead/d/dutta26/ngc_1027/ngc_g_coadd.fits')
data = np.array(f[0].data)
f.close()

y, x=np.shape(data)
del data  
#Make the residual maps now
#First do linear interpolation to each pixel 
points=[]
values =[] 
grid_y, grid_x = np.meshgrid(np.linspace(0, y-1, y), np.linspace(0, x-1, x))

for x in range(xBins):
    for y in range(yBins):
        
        points.append([(1000*y)+ 500,(1000*x)+ 500])
        values.append(img[y,x])
        

grid_z1 = griddata(points, values, (grid_y, grid_x), method='cubic')    
#Write data
grid_z1 = np.array(grid_z1, dtype = np.float32).T
# =============================================================================
# hdu = fits.PrimaryHDU(grid_z1) 
# hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/residue.fits', clobber=True)   
# =============================================================================

#Now subtract the interpolated psfs from the actual psfs 
#Open the original image and measure dimensions 
f=fits.open('/scratch/halstead/d/dutta26/ngc_1027/ngc_g_coadd.fits')
data = np.array(f[0].data)
f.close()

y, x=np.shape(data)
del data  

points=[]
values =[] 
grid_y, grid_x = np.meshgrid(np.linspace(0, y-1, y), np.linspace(0, x-1, x))
for j in range(len(content)- 1):
    x = int(round(xList[j]))
    y = int(round(yList[j]))
    points.append([y,x])
    values.append( grid_z1[y,x] - e1List[j])
grid_z1 = griddata(points, values, (grid_y, grid_x), method='cubic')    
#Write data
grid_z1 = np.array(grid_z1, dtype = np.float32).T
hdu = fits.PrimaryHDU(grid_z1) 
hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/residue.fits', clobber=True)     



