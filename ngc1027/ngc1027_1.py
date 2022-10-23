#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 16:15:16 2020

@author: anirban
"""

from astropy.io import fits
from astropy import wcs
import numpy as np
import os
import sys
import numpy as np
from astropy.io import fits
#import wrapper
from scipy.ndimage import rotate
import pandas as pd
import helper
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys

filename = '/scratch/halstead/d/dutta26/m_38/coadd.fits'
#filename = str(sys.argv[1])
#Read the catalog 
catalog = '/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test1.cat'
#catalog = str(sys.argv[2])
with open(catalog) as f:
    content = f.readlines()
    
Raindex=5
Decindex=6

#Make an array containing x and y indices 
raList =[]
decList =[]
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue
    raList.append(float(content[j].split()[Raindex])) 
    decList.append( float(content[j].split()[Decindex])) 
    
indexList, xList, yList = helper.detNearby(raList, decList, 25 , filename)
xList,yList = helper.createNewList(xList, yList, indexList)
print (len(xList))
#Discard any objects from bad regions (determined from visual)
xList1=[]
yList1=[]
for j in range(len(xList)):
    x = xList[j]
    y = yList[j]
    flag = 0
    if(x<2200 or x>23500):
        flag =1
    if(y<3000 or y>29000):
        flag = 1
# =============================================================================
#     if (x<4390 and y>22300):
#         flag = 1
# =============================================================================
    if(flag == 0):
        xList1.append(x)
        yList1.append(y)
        

indexList = helper.detectBad_combined(filename, xList1, yList1)
xList,yList = helper.createNewList(xList1, yList1, indexList)
print (len(xList))
#Open the original image to create cutouts 
f=fits.open(filename)
data = np.array(f[0].data)
f.close()


#temp_filename = '/media/anirban/anirban/cpp_temp/1027_cut.fits'
#dataFile = str(sys.argv[3])
dataFile = '/scratch/halstead/d/dutta26/ngc_1027/1027_spatial_data'
txtf = open(dataFile, "w+")
txtf.write("# X Pos        Y Pos       PSF       e1       e2 \n")

#Create cutouts and run the measure algorithm 
for j in range(len(xList)):
    #print (j)
    x = int(round( xList[j]))
    y = int(round(yList[j]))
    cut = data[y-25: y+26, x-25: x+26]
    
    flux, mux, muy, e1, e2, bkg, psf = measure_pythonV.measure(cut) 
    if(flux == None or e1==None or e2 == None or mux == None or muy == None or psf== None):
        continue  
    if(flux == np.nan or e1==np.nan or e2 == np.nan or mux == np.nan or muy == np.nan or psf== np.nan):
        continue  
    if(math.isnan(flux) or math.isnan(e1) or math.isnan(e2) or math.isnan(mux) or math.isnan(muy) or math.isnan(psf)):
        continue  
    x=x+mux
    y= y +muy

# =============================================================================
#     hdu = fits.PrimaryHDU(cut) 
#     hdu.writeto(temp_filename, clobber=True)
# 
#     
# 
#     #Run detection algorithm 
#     filename_encode = temp_filename.encode('utf-8')
#     a=wrapper.measure([b"abc", filename_encode])
#     #Measurement done. Now extract values
#     f=open('/media/anirban/anirban/cpp_temp/test.txt')
#     content = f.readlines()
#     f.close()
#     if(len(content) == 0):
#         continue
#     
#     e1_m = float(content[len(content)-1].split()[5])
#     e2_m = float(content[len(content)-1].split()[6])
#     
#     flux_measure =float(content[len(content)-1].split()[3])
#     posx_measure = float(content[len(content)-1].split()[10])
#     posy_measure = float(content[len(content)-1].split()[11])
#     psf_measure = float(content[len(content)-1].split()[12])
#     bkg_measure = float(content[len(content)-1].split()[13])
# =============================================================================
    if(psf < 0):
        print (j, psf)
        break
# =============================================================================
#     if(np.sqrt(e1**2 + e2**2) > 0.2):
#         continue
# =============================================================================
    #if(flux < 10 or flux>500000):
    #    continue
    #Definitelt merged object or something weird going on
# =============================================================================
#     if(abs(e1)>0.15 or abs(e2)>0.15):
#         continue
# =============================================================================
    line = str(x)+ ' , '+ str(y) + ' , ' + str(psf) + ' , '+str(e1)+ ' , ' +str(e2)
    txtf.write(line + '\n')

del data
txtf.close()
