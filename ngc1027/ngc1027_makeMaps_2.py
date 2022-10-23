#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 19:30:58 2021

@author: dutta26
"""


from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
from astropy.stats import sigma_clipped_stats
from scipy.interpolate import griddata
import helper1
#This code creates residual PSFs and ellipticities
#First it boxes in a 1000x1000 sq the avg PSF, e1 and e2
#The uses linear bi interpolation to calculate PSF and e1 and e2 at any location 
#Read the data
#catalog = '/scratch/halstead/d/dutta26/ngc_1027/1027_spatial_data'
#catalog = str(sys.argv[1])
#red_fact = int(sys.argv[3])

catalog = '/scratch/halstead/d/dutta26/ngc_1027/1027_spatial_data'
red_fact = 4
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
    if(np.isnan(float(content[j].split()[xInd])) or np.isnan(float(content[j].split()[yInd])) or np.isnan(float(content[j].split()[pInd])) or np.isnan(float(content[j].split()[e1Ind]))):
        continue
    xList.append(float(content[j].split()[xInd])) 
    yList.append( float(content[j].split()[yInd])) 
    pList.append(float(content[j].split()[pInd]))
    e1List.append(float(content[j].split()[e1Ind]))
    e2List.append(float(content[j].split()[e2Ind]))




#Open the original image and measure dimensions 
#originalFile = str(sys.argv[2])
originalFile = '/scratch/halstead/d/dutta26/m_38/coadd.fits'
f=fits.open(originalFile)
data = np.array(f[0].data)
f.close()

yShape, xShape =np.shape(data)
yShape = int(yShape / red_fact)+5
xShape = int(xShape / red_fact)+5
del data  

#Commented out to make histogram 
points=[]
values =[] 
grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
for j in range(len(xList)):
    x = int(round(xList[j]/red_fact))
    y = int(round(yList[j]/red_fact))
    
# =============================================================================
#     cacl_val = helper1.returnPolyVal(popt, x/50, y/50)
#     if(abs(cacl_val - e1List[j])) > 0.03:
#         continue
# =============================================================================
    if(pList[j] == np.nan or pList[j]>3 or pList[j]<2.5):
        continue
    points.append([y,x])
    
    values.append( pList[j])
grid_z1 = griddata(points, values, (grid_y, grid_x), method='linear')    
#Write data
grid_z1 = np.array(grid_z1, dtype = np.float32).T
hdu = fits.PrimaryHDU(grid_z1) 
hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/residue.fits', clobber=True)     


# =============================================================================
# import matplotlib.pyplot as plt
# 
# # An "interface" to matplotlib.axes.Axes.hist() method
# n, bins, patches = plt.hist(x=pList, bins='auto', color='#0504aa',
#                             alpha=0.7, rwidth=0.95)
# plt.grid(axis='y', alpha=0.75)
# plt.xlabel('PSF Value')
# plt.ylabel('Frequency')
# #plt.title('My Very Own Histogram')
# #plt.savefig('/home/dutta26/sine.png')
# =============================================================================


