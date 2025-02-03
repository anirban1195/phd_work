#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 15:30:21 2024

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
import wquantiles
from astropy import wcs
import matplotlib.pyplot as plt

def integrate(z):
    arr= np.arange(0, z, 0.01)
    #a= np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))
    #return a
    return np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))*(1/(1+z))

#Read Redshifts
zMax = 2.0
zMin = 0.4
zFile = '/home/dutta26/photz_eazy.zout'
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
                if(float((content[j].split())[8]) >= 0.8 and float((content[j].split())[15])>=3):
                    redShiftArr.append(float((content[j].split())[7]))
                else:
                    redShiftArr.append(0)
            else:
                redShiftArr.append(float((content[j].split())[1]))
redShiftArr = np.array(redShiftArr)        

ir_coadd_data = np.load( '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
loc = np.where( (ir_coadd_data[:,2] == 0) & (redShiftArr<zMax) & (redShiftArr>zMin) )[0] 
tot = len(loc)
pArr =[]
dsArr = []
dlsArr = []
for j in range(10):
    print (zMin, zMin+0.2)
    loc = np.where( (ir_coadd_data[:,2] == 0) & (redShiftArr[:]<zMin+0.2) & (redShiftArr[:]>zMin)  )[0]
    pArr.append(len(loc)/tot)
    
    dsArr.append(integrate(zMin+.1))
    temp = (integrate(zMin+.1)*(zMin+.1+1) - integrate(0.23)*1.23)/(1+zMin+.1)
    dlsArr.append(temp)
    zMin += 0.2
    
pArr=np.array(pArr)
dsArr = np.array(dsArr)
dlsArr= np.array(dlsArr)
print (np.sum(dlsArr*pArr/dsArr))