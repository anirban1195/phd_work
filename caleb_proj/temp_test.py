#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 09:57:47 2024

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import os
import subprocess
import gzip
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
from astropy.stats import sigma_clipped_stats
import helper
import matplotlib.pyplot as plt
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')



# =============================================================================
# phosimLoc = '/home/dutta26/Downloads/phosim_release/'
# outLoc = phosimLoc + 'output/'
# 
# dataSet ='/scratch/bell/dutta26/abell_2390/'
# filt_list = ['g', 'r',  'z', 'i']
# #filt_list = [ 'u','z']
# colorArr = ['g', 'r',  'k', 'm']
# counter =0
# j =0
# temp_arr = np.zeros((100, 10))
# cloudArr =[]
# for folder in filt_list:
#     for files in os.listdir(dataSet+folder):
#         #print (files)
#         if('weight' in files or 'temp' in files ):
#             continue
#         
#         
#         f=fits.open(dataSet+folder+'/'+files)
#        # print (dataSet+folder+'/'+files)
#         mjd = float((f[0].header)['MJD-MID'])
#         airmass = float((f[0].header)['AIRMASS'])
#         zp = (f[0].header)['MAGZERO']
#         flx_scale = (f[0].header)['FLXSCALE']
#         ra_header = (f[0].header)['RA']
#         dec_header = (f[0].header)['DEC']
#         c=SkyCoord(ra_header, dec_header, unit=(u.hourangle, u.deg))
#         ra,dec = c.ra.degree, c.dec.degree
#         filt = (f[0].header)['FILTER']
#         skylevel = float((f[0].header)['SKYBG'])
#         skymag = float((f[0].header)['SKYMAG'])
#         print (skylevel)
# # =============================================================================
# #         if(mjd < 59187.036 or mjd > 59191.172):
# #             continue
# # =============================================================================
#         if(filt == 'odi_u'):
#             filt_code =0
#         elif(filt == 'odi_g'):
#             cloudArr.append( 26.4- float(zp))
#             filt_code = 1
#         elif(filt == 'odi_r'):
#             filt_code = 2
#             cloudArr.append(26.3 -float(zp))
#         elif(filt == 'odi_i'):
#             filt_code = 3
#             cloudArr.append(25.8 - float(zp))
#         elif(filt == 'odi_z'):
#             filt_code = 4
#             cloudArr.append(24.9 - float(zp))
#         else:
#             filt_code = 5
#         
#         f.close()
#         plt.plot(mjd, float(zp), colorArr[counter]+'.')
#         
#     counter += 1
#     
# 
# =============================================================================
        
        
master_arr  = np.load('/home/dutta26/2024/Apr/temp_5.npy')
colorArr = ['b','g', 'r','m', 'k']

a=[]
b=[]
for j in range(len(master_arr)):
    if(master_arr[j,0] == 0):
        continue
    if(int(master_arr[j,7]) != 0):
        continue
    color = colorArr[int(master_arr[j,7])]
    plt.plot(master_arr[j,6], master_arr[j,2], color+'+')
    plt.plot(master_arr[j,6], master_arr[j,3], color+'.')
    a.append(master_arr[j,2])
    b.append(master_arr[j,3])
    plt.xlabel('MJD')

print (np.median(a), np.median(b))