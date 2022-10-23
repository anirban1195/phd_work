#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 11:42:23 2020

@author: dutta26
"""


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
import helper1
#This code creates residual PSFs and ellipticities
#First it boxes in a 1000x1000 sq the avg PSF, e1 and e2
#The uses linear bi interpolation to calculate PSF and e1 and e2 at any location 
#Read the data
#catalog = '/scratch/halstead/d/dutta26/ngc_1027/1027_spatial_data'
catalog = str(sys.argv[1])
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


#Reduce the data by factors of 100
xList1=[]
yList1=[]

for j in range(len(xList)):
    xList1.append(xList[j]/50.0)
    yList1.append(yList[j]/50.0)

#Get the rage of values of x and y
xMax = int(np.max(xList1)) + 5
yMax = int(np.max(yList1)) + 5
#End of getting range
#Make the arrays that will hold values for optimization
z = np.ones((yMax, xMax), dtype=np.float32)
sigma = np.ones((yMax, xMax), dtype=np.float32)
X = np.ones((yMax, xMax), dtype=np.int16)
Y = np.ones((yMax, xMax), dtype=np.int16)

#Where z is unmeasured it is mean vaue with std=mean value
z = z*np.mean(e1List)
sigma = sigma*100

# =============================================================================
# for j in range(len(xList1)):
#     x = int(round(xList1[j]))
#     y = int(round(yList1[j]))
#     z[y,x] = e1List[j]
#     sigma[y,x] = 0.0001
#     X[y,x]  = x
#     Y[y,x] = y
#     
#     
# #X,Y = np.meshgrid(xList1, yList1)
# popt, pcov = helper1.fitPoly2D(X,Y, z.ravel(), sigma.ravel())
# print (popt)
# 
# x = np.linspace(0, xMax, xMax-1)
# y = np.linspace(0, yMax, yMax-1)
# X, Y = np.meshgrid(x,y)
# img = helper1.returnPolyVal(popt, X, Y)
# #Write data
# hdu = fits.PrimaryHDU(img) 
# hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/test.fits', clobber=True)
# 
# 
# =============================================================================
#Open the original image and measure dimensions 
originalFile = str(sys.argv[2])
f=fits.open(originalFile)
data = np.array(f[0].data)
f.close()

yShape, xShape =np.shape(data)
del data  

# =============================================================================
# x = np.linspace(0, xShape, xShape-1)
# y = np.linspace(0, yShape, yShape-1)
# X, Y = np.meshgrid(x,y)
# residue = popt[0]*(X/50.0)**2 + popt[1]*(Y/50.0)**2 + popt[2]*(X/50.0) + popt[3]*(Y/50.0) + popt[4]*(X/50.0)*(Y/50.0) + popt[5]
# residue =np.array(residue, dtype=np.float32)
# hdu = fits.PrimaryHDU(residue) 
# hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/residue.fits', clobber=True)
# =============================================================================

points=[]
values =[] 
grid_y, grid_x = np.meshgrid(np.linspace(0, yShape-1, yShape), np.linspace(0, xShape-1, xShape))
for j in range(len(content)- 1):
    x = int(round(xList[j]))
    y = int(round(yList[j]))
    
# =============================================================================
#     cacl_val = helper1.returnPolyVal(popt, x/50, y/50)
#     if(abs(cacl_val - e1List[j])) > 0.03:
#         continue
# =============================================================================
    points.append([y,x])
    values.append( pList[j])
grid_z1 = griddata(points, values, (grid_y, grid_x), method='cubic')    
#Write data
grid_z1 = np.array(grid_z1, dtype = np.float32).T
hdu = fits.PrimaryHDU(grid_z1) 
hdu.writeto('/scratch/halstead/d/dutta26/ngc_1027/residue.fits', clobber=True)     





