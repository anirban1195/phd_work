#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 09:58:45 2020

@author: anirban
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
import wrapper
import numpy as np
from astropy.io import fits
import wrapper
from scipy.ndimage import rotate
import pandas as pd
from astropy.stats import sigma_clipped_stats
import math,sys


#Read the catalog 
catalog = '/home/anirban/test.cat'
with open(catalog) as f:
    content = f.readlines()
    
Xindex=3
Yindex=4

#Make an array containing x and y indices 
xList =[]
yList =[]
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    xList.append(float(content[j].split()[Xindex])) 
    yList.append( float(content[j].split()[Yindex])) 

#Discard any objects from bad regions (determined from visual)
xList1=[]
yList1=[]
for j in range(len(xList)):
    x = xList[j]
    y = yList[j]
    flag = 0
    if(x<2400 or x>22700):
        flag =1
    if(y<3500 or y>28500):
        flag = 1
    if (x<8600 and y>22500):
        flag = 1
    if(flag == 0):
        xList1.append(x)
        yList1.append(y)
        
#If there are two object closeby then consider NONE
xList2 =[]
yList2 =[]
for j in range(len(xList1)):
    x1 = xList1[j]
    y1 = yList1[j]
    flag =0
    minVal = j-2000
    maxVal = j+2000
    if(minVal < 0):
        minVal = 0
    if(maxVal > len(xList1)):
        maxVal = len(xList1)
    for k in range(minVal, maxVal):
        if(k == j):
            continue
        x2 = xList1[k]
        y2 = yList1[k]
        if(((x1-x2)**2 + (y1-y2)**2) < 1000):
            flag = 1
    #If flag=0 then it is a good object
    if(flag ==0 ):
        xList2.append(x1) 
        yList2.append(y1) 
    
del xList, yList, xList1, yList1
        
#Open the original image to create cutouts 
f=fits.open('/home/anirban/result/ngc_1027.odi_g.fits')
data = np.array(f[0].data)
f.close()
psfList =[]
e1List =[]
e2List =[]
filename = '/media/anirban/anirban/cpp_temp/1027_cut.fits'
txtf = open("/home/anirban/1027_spatial_data", "w")
txtf.write("# X Pos        Y Pos       PSF       e1       e2 \n")

#Create cutouts and run the measure algorithm 
for j in range(len(xList2)):
    x = int(round( xList2[j]))
    y = int(round(yList2[j]))
    cut = data[y-25: y+26, x-25: x+26]
    hdu = fits.PrimaryHDU(cut) 
    hdu.writeto(filename, clobber=True)
    #Run detection algorithm 
    filename_encode = filename.encode('utf-8')
    a=wrapper.measure([b"abc", filename_encode])
    #Measurement done. Now extract values
    f=open('/media/anirban/anirban/cpp_temp/test.txt')
    content = f.readlines()
    f.close()
    if(len(content) == 0):
        continue
    
    e1_m = float(content[len(content)-1].split()[5])
    e2_m = float(content[len(content)-1].split()[6])
    
    flux_measure =float(content[len(content)-1].split()[3])
    posx_measure = float(content[len(content)-1].split()[10])
    posy_measure = float(content[len(content)-1].split()[11])
    psf_measure = float(content[len(content)-1].split()[12])
    bkg_measure = float(content[len(content)-1].split()[13])
    #Definitelt merged object or something weird going on
    if(abs(e1_m)>0.1 or abs(e2_m)>0.1):
        continue
    line = str(x)+ ' , '+ str(y) + ' , ' + str(psf_measure) + ' , '+str(e1_m)+ ' , ' +str(e2_m)
    txtf.write(line + '\n')

del data
txtf.close()