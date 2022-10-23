#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 16:40:23 2022

@author: dutta26
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os
import sys
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import math,sys
import subprocess

addCat = 0

catLoc = '/scratch/halstead/d/dutta26/lsst/lsst_catalog1/'
f=open(catLoc+ 'sheared.cat')
shear_cont = f.readlines()
f.close()

    
#Make sheared files 
for files in os.listdir(catLoc):
    if('cat' in files or 'sheared' in files or 'flat' in files):
        continue
    else:
        f=open(catLoc+files)
        cat_cont = f.readlines()
        f.close()
        
        f=open(catLoc+files+'_sheared', 'w+')
        for j in range(len(cat_cont)):
            if('stars' in cat_cont[j] or 'galaxies' in cat_cont[j]):
                continue
            f.write(cat_cont[j])
        
        for j in range(len(shear_cont)):
            f.write(shear_cont[j])
            
        f.close()
        
        
#Make the flat_txt files 
raArr =np.arange(237.0083-(5*0.00555), 237.0083+(5*0.00555), 0.00555 )
decArr = np.arange(-15.2308 -(5*0.00555) ,  -15.2308 + (5*0.00555) , 0.005551)    
for files in os.listdir(catLoc):
    if('cat' in files or 'sheared' in files or 'flat' in files):
        continue
    else:
        f=open(catLoc+files)
        cat_cont = f.readlines()
        f.close()
        f=open(catLoc+files+'_flat', 'w+')
        for j in range(len(cat_cont)):
            if('stars' in cat_cont[j] or 'galaxies' in cat_cont[j]):
                continue
            f.write(cat_cont[j])
            
        count = 1
        mag = 16
        for j in range(len(raArr)):
            mag = mag+1
            for k in range(len(decArr)):
                temp_str = str(count) + ' '+str(raArr[j])[0:9]+ ' '+ str(decArr[k])[0:9] + ' '+str(mag)
                f.write('object ' +temp_str + ' ../sky/sed_flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 star none none \n')
                count += 1
        f.close()
        

    
    