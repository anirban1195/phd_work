#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:54:22 2024

@author: dutta26
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os
from astropy.io import fits
import numpy as np
import os
import subprocess
import gzip
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
from astropy.stats import sigma_clipped_stats

phosimLoc = '/home/dutta26/Downloads/phosim_release/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/bell/dutta26/abell_2390/'
#filt_list = ['g', 'r', 'u', 'z', 'i']
filt_list = ['g','r','z', 'i','u']

counter =0
zpArr =[]
flxScaleArr=[]
mjdArr =[]
for folder in filt_list:
    for files in os.listdir(dataSet+folder):
        #print (files)
        if('weight' in files or 'temp' in files):
            continue
        f=fits.open(dataSet+folder+'/'+files)
        #print (dataSet+folder+'/'+files)
        mjd = float((f[0].header)['MJD-MID'])
        ra_header = (f[0].header)['RA']
        dec_header = (f[0].header)['DEC']
        zp = (f[0].header)['MAGZERO']
        flx_scale = (f[0].header)['FLXSCALE']
        c=SkyCoord(ra_header, dec_header, unit=(u.hourangle, u.deg))
        ra,dec = c.ra.degree, c.dec.degree
        filt = (f[0].header)['FILTER']
        if(mjd < 8567):
            continue
        if(filt == 'odi_u'):
            filt_code =0
        elif(filt == 'odi_g'):
            filt_code = 1
        elif(filt == 'odi_r'):
            filt_code = 2
        elif(filt == 'odi_i'):
            filt_code = 3
        elif(filt == 'odi_z'):
            filt_code = 4
        else:
            filt_code = 5
        
        f.close()
        print (zp)
        zpArr.append(zp)
        flxScaleArr.append(flx_scale)
        mjdArr.append(mjd)

mjdArr= np.array(mjdArr)
flxScaleArr = np.array(flxScaleArr)
zpArr= np.array(zpArr)
f=open('/home/dutta26/comparison_abell2390_bkg1.txt')
content = f.readlines()
f.close()
colorArr= ['b', 'g', 'r', 'k', 'm']
magArr =[]
magArr1 =[]
for j in range(len(content)):
    
    if(j==0):
        continue
    mjd = float(content[j].split(',')[6][0:15])
    color = colorArr[int(content[j].split(',')[5])]
    if(int(content[j].split(',')[5]) != 4):
        continue
    
    sim_median =float(content[j].split(',')[1][0:12])
    sim_std = float(content[j].split(',')[2][0:12])
    
    
    actual_median = float(content[j].split(',')[3][0:12])
    actual_std = float(content[j].split(',')[4][0:12])
    
    #plt.errorbar(mjd, sim_median, yerr= sim_std, fmt = color+'.')
    #plt.errorbar(mjd, actual_median, yerr= actual_std, fmt = color+'^')
    
    loc = np.where(np.abs(mjdArr-mjd) <0.00001)[0]
    print (loc)
    zp = zpArr[loc]
    flxScale = flxScaleArr[loc]
    magArr1.append(zp[0] - 2.5*np.log10(actual_median*9.09*9.09/60))
    magArr.append(25 - 2.5*np.log10(actual_median*9.09*9.09*flxScale[0]))
    #mag_per_arcsec2 = zpArr[j-1] - 2.5*np.log10(actual_median*9.09*9.09/60)
    #plt.errorbar(mjd, mag_per_arcsec2, yerr= 0, fmt = color+'^')
    
    #mag_per_arcsec2 = 25 - 2.5*np.log10(actual_median*flxScaleArr[j-1]*9.09*9.09)
    #plt.errorbar(mjd, mag_per_arcsec2, yerr= 0, fmt = color+'.')
    
print (sigma_clipped_stats(magArr, cenfunc=np.median))
print (sigma_clipped_stats(magArr1, cenfunc=np.median))