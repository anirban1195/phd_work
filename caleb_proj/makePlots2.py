#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 09:39:27 2024

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

f=open('/home/dutta26/comparison_abell2390_bkg4.txt')
content = f.readlines()
f.close()
colorArr= ['b', 'g', 'r', 'k', 'm']
magArr_real =[[], [], [], [], []]
magArr_phosim =[[], [], [], [], []]
cts_real =[[], [], [], [], []]
cts_phosim =[[], [], [], [], []]
pan_starrs_mag =[21.9, 21, 20.3, 19.5 , 18]
pan_starrs_lambda =[ 480, 620, 750, 865, 960 ]

lambda_mean = [370, 500, 630, 750, 870 ]
for j in range(len(content)):
    
    if(j==0):
        continue
    idx = int(content[j].split(',')[5])
    
    
    sim_median =float(content[j].split(',')[1][0:12])
    sim_std = float(content[j].split(',')[2][0:12])
    
    
    actual_median = float(content[j].split(',')[3][0:12])
    actual_std = float(content[j].split(',')[4][0:12])
    
    star_flux = float(content[j].split(',')[8][0:12])
    flxScale = float(content[j].split(',')[7][0:8])
    print (flxScale)
    magArr_real[idx].append(25 - 2.5*np.log10(actual_median*9.09*9.09*flxScale) )
    cts_real[idx].append(actual_median)
    
    phosim_zp = 17 + 2.5*np.log10(star_flux) + 2.5*np.log10(500/lambda_mean[idx])
    
    magArr_phosim[idx].append(phosim_zp - 2.5*np.log10(sim_median*9.09*9.09))
    cts_phosim[idx].append(sim_median)
    

for j in range(5):
    plt.errorbar(lambda_mean[j], np.median(magArr_real[j]), yerr = np.std(magArr_real[j]), fmt = 'r.', markersize =10)
    plt.errorbar(lambda_mean[j], np.median(magArr_phosim[j]), yerr = np.std(magArr_phosim[j]), fmt = 'r*', markersize =10)
    #plt.errorbar(lambda_mean[j], np.median(cts_real[j]), yerr = np.std(cts_real[j]), fmt = 'r.', markersize =10)
    #plt.errorbar(lambda_mean[j], np.median(cts_phosim[j]), yerr = np.std(cts_phosim[j]), fmt = 'r*', markersize =10)
    #plt.errorbar(pan_starrs_lambda[j], pan_starrs_mag[j], fmt = 'b.', markersize =10)
    
    
plt.plot(0,0, 'r.', markersize =10, label ='Wiyn-ODI')
plt.plot(0,0, 'r*', markersize =10, label ='Phosim-WIYN')
#plt.errorbar(lambda_mean[j], np.median(cts_real[j]), yerr = np.std(cts_real[j]), fmt = colorArr[j]+'.', markersize =10)
#plt.errorbar(lambda_mean[j], np.median(cts_phosim[j]), yerr = np.std(cts_phosim[j]), fmt = colorArr[j]+'*', markersize =10)
#plt.plot(0,0, 'b.', markersize =10, label ='Pan-starrs')

    
plt.gca().invert_yaxis()
plt.legend()
plt.xlabel('Wavelength in nm')
plt.ylabel('Sky Brightness mag/arcsec^2')
#plt.ylabel('Counts')