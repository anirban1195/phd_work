#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 18:34:39 2024

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt

ir_sf_df=np.load( '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
zFile = '/home/dutta26/photz_eazy.zout'
redShiftArr =[]

f=open(zFile)
content = f.readlines()
f.close()
for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        if('eazy' in zFile):
            if(float((content[j].split())[8]) >= 0.8):
                redShiftArr.append(float((content[j].split())[7]))
            else:
                redShiftArr.append(0)
        else:
            redShiftArr.append(float((content[j].split())[1]))
def integrate(z):
    arr= np.arange(0, z, 0.01)
    #a= np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))
    #return a
    return np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))*(1/(1+z))

xArr=[]
yArr=[]
zArr=[]
for j in range(len(redShiftArr )):
    if(redShiftArr[j] <= 0.1 or redShiftArr[j]>0.5 or ir_sf_df[j,2]==1 or ir_sf_df[j,2]==2 or ir_sf_df[j,3]<= 1):
        continue
    dl =  4285.71*integrate(redShiftArr[j])
    xArr.append((ir_sf_df[j,10]-15000)*(0.11/3600)*(np.pi/180)*(1+integrate(redShiftArr[j]))**2*dl)
    yArr.append((ir_sf_df[j,11]-15000)*(0.11/3600)*(np.pi/180)*(1+integrate(redShiftArr[j]))**2*dl)
    zArr.append(dl)
    
#ax =plt.axes(projection='3d')
#ax.scatter3D(xArr, yArr,zArr, s=1)
xArr=np.array(xArr)
yArr=np.array(yArr)
zArr=np.array(zArr)

x_centers = np.linspace(np.min(xArr)-1, np.max(xArr)+1, 250)
y_centers = np.linspace(np.min(yArr)-1, np.max(yArr)+1, 250)
z_centers = np.linspace(np.min(zArr)-1, np.max(zArr)+1, 250)


densityArr = np.zeros((250, 250, 250), dtype=np.float32)
for x in range(250):
    print (x)
    for y in range(250):
        for z in range(250):
            loc = np.where((xArr>x_centers[x]-2) & (xArr<x_centers[x]+2) &
                           (yArr>y_centers[y]-2) & (yArr<y_centers[y]+2) &
                           (zArr>z_centers[z]-2) & (zArr<z_centers[z]+2) )[0]
            densityArr[x,y,z] = len(loc)
#np.save('/scratch/bell/dutta26/abell_2390/3d_density.npy', densityArr)


