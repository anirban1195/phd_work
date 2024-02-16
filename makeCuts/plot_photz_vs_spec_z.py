#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 10:22:47 2024

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
from astropy import wcs
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
df=pd.read_csv('/home/dutta26/table_objsearch_div.csv')
ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
redShiftArr=[]
zFile = '/home/dutta26/photz_eazy.zout'
if (zFile == None):
    redShiftArr = []
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
                if(float((content[j].split())[8]) >= 0.8):
                    redShiftArr.append(float((content[j].split())[7]))
                else:
                    redShiftArr.append(0)
            else:
                redShiftArr.append(float((content[j].split())[1]))

f=fits.open('/scratch/bell/dutta26/abell_2390/temp.fits')
data = f[0].data
f.close()
sizey, sizex = np.shape(data)
redShiftArr = np.array(redShiftArr)    
spec_z_arr =[]
phot_z_arr =[]
for j in range(len(df)):
    ra = df['RA'].iloc[j]
    dec = df['Dec'].iloc[j]
    spec_z = df['Redshift (z)'].iloc[j]
    #loc= np.where ()
    loc = np.where((ir_coadd_data[:,0] >= ra -0.0009) & (ir_coadd_data[:,0]<= ra +0.0009) &
                   (ir_coadd_data[:,1]>= dec -0.0009) & (ir_coadd_data[:,1]<= dec +0.0009) &
                   (ir_coadd_data[:,2]== 0))[0]
    #loc = np.where(ir_coadd_data[:,0] >= 3)[0]
    #print (len(loc), ra, dec)
    if(len(loc)>0):
        maxFlux = np.max(ir_coadd_data[loc, 3])
        index = np.where(ir_coadd_data[loc, 3] == maxFlux)[0]
        #print (index, len(loc))
        spec_z_arr.append(spec_z)
        phot_z_arr.append(redShiftArr[loc[index[0]]])
        continue
        
    xList, yList = helper.convertToXY([ra], [dec], '/scratch/bell/dutta26/abell_2390/temp.fits')
    #print (xList[0], yList[0])
    x=int(xList[0])
    y=int(yList[0])
    if(xList[0]< 0 or yList[0]<0 or xList[0]> sizex or yList[0]> sizey):
        continue
    else:
        print ('abc', x, y)
        data[y-20, x-20:x+20]=100
        data[y+20, x-20:x+20]=100
        data[y-20:y+20, x-20]=100
        data[y-20:y+20, x+20]=100
        
hdu = fits.PrimaryHDU(data)
hdu.writeto('/scratch/bell/dutta26/abell_2390/temp1.fits', overwrite=True)