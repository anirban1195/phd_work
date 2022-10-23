#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 09:21:26 2022

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
folderLoc ='/scratch/halstead/d/dutta26/lsst/lsst1/'
destfolderLoc ='/scratch/halstead/d/dutta26/lsst/lsst2/'

#fileList_name = '/home/dutta26/temp.ascii'
catalogLoc = '/scratch/halstead/d/dutta26/lsst/lsst_catalog/0.txt_sheared'
#Read catalog 
f=open(catalogLoc)
content = f.readlines()
f.close()
source_ra =[]
source_dec = []
source_mag=[]
source_z=[]
count = 0

for j in range(len(content)):
    if(len(content[j].split()) == 15):
        if(float(content[j].split()[4]) > 25):
            continue
        ra = float(content[j].split()[2])
        dec = float(content[j].split()[3])
        if(ra > 237.275 or ra < 236.733):
            continue
        if(dec < -15.49833 or dec > -14.9700):
            continue
        source_ra.append(ra)
        source_dec.append(dec)
        source_mag.append(float(content[j].split()[4]))
        source_z.append(float(content[j].split()[6]))
    
    if(len(content[j].split()) == 30):
        
        if(float(content[j].split()[4]) > 25):
            continue
        if('.2' in (content[j].split()[1])):
            continue
        ra = float(content[j].split()[2])
        dec = float(content[j].split()[3]) 
        if(ra > 237.275 or ra < 236.733):
            continue
        if(dec < -15.4933 or dec > -14.9700):
            continue
        source_ra.append(ra)
        source_dec.append(dec)
        source_mag.append(float(content[j].split()[4]))
        source_z.append(float(content[j].split()[6]))
        
source_ra = np.array(source_ra)   
source_dec = np.array(source_dec) 
source_mag = np.array(source_mag)    
        
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
    

fluxArr=np.array(fluxArr)
raArr = np.array(raArr)
decArr = np.array(decArr)    
zArr=[]       
for j in range(len(raArr)):
    temp_ra = source_ra - raArr[j]
    temp_dec = source_dec - decArr[j]
    temp_dist = np.sqrt(temp_ra**2 + temp_dec**2)
    loc = np.where(temp_dist< 0.0002)[0]
    if(len(loc) <= 0):
        zArr.append(-99)
    elif(len(loc) == 1):
        zArr.append(source_z[loc[0]])
    else:
        #Find the brightest source 
        index = np.argmax(source_mag[loc])
        zArr.append(source_z[loc[index]])
        
zArr = np.array(zArr) 
np.save( '/scratch/halstead/d/dutta26/lsst/zArr.npy', zArr)