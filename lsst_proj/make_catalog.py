#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 16:24:26 2022

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
if(addCat == 1):
    f=open('/home/dutta26/apps/phosim-phosim_release-2479c7396c0c/output/catgen_1001.cat')
    content = f.readlines()
    f.close()
    
    
band = 'i'
ditherSeq = [[-0.0166]]
fileCount = 0
fileList = os.listdir('/scratch/halstead/d/dutta26/abell_2390/'+band +'/')
for file in fileList:

    if('.weight' in file or 'temp' in file):
        continue
    fileCount += 1
    print (file, fileCount)
    
    #if(fileCount > 1):
    #    break
    #Read the file data
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'+file)
    back_h = float((f[0].header)['SKYBG'])
    seeing = float((f[0].header)['SEEING'])
    zp = float((f[0].header)['MAGZERO'])
    fwhm = float((f[0].header)['FWHM_FLT'])
    mjd = float((f[0].header)['MJD-MID'])
    airmass = float((f[0].header)['AIRMASS'])
    #Fix for moon phase is waxing gibbous
    if(type((f[0].header)['MOONPHSE']) is str):
        mphase = 50
    else:
        mphase = float((f[0].header)['MOONPHSE'])*100
    mAngle = float((f[0].header)['MOON_D'])
    expTime = float((f[0].header)['EXPTIME'])  
    focus = float((f[0].header)['TELFOCUS'])
    zp_n = float((f[0].header)['PHOTZP_N'])
    skymag = float((f[0].header)['SKYMAG']) 
    depth = float((f[0].header)['PHOTDPTH']) 
    mRa = float((f[0].header)['MOON_RA']) 
    mDec = float((f[0].header)['MOON_DEC']) 
    flux_scale = float((f[0].header)['FLXSCALE']) 
    target_alt = float((f[0].header)['ELMAP'])*180/np.pi 
    target_az = float((f[0].header)['AZMAP'])*180/np.pi 
    
    moon_alt = float((f[0].header)['MOON_ALT'])
    moon_az = float((f[0].header)['MOON_AZ'])
    
    sun_alt = float((f[0].header)['SUN__ALT'])
    sun_az = float((f[0].header)['SUN__AZ'])
    
    f.close()
    
    f=open('/home/dutta26/Downloads/phosim_core/examples/odi_catalog/'+str(fileCount)+'_star.txt', 'w+')
    f.write('rightascension 328.4083 \n')    
    f.write('declination 17.6697 \n') 
    f.write('mjd '+ str(mjd)+' \n') 
    f.write('seed 100'+str(fileCount)+ '\n')
    f.write('vistime 60.0 \n')
    f.write('seeing '+str(seeing)+'\n')
    f.write('filter 4 \n')
    f.write('obshistid 100'+str(fileCount)+ '\n')
    f.write('stars 15.0 30.0 0.666 \n')
    #f.write('galaxies 15.0 30.0 0.666 \n')
    
# =============================================================================
#     f.write('altitude '+str(target_alt)+ '\n')
#     f.write('azimuth '+str(target_az)+ '\n')
#     f.write('moonalt '+str(moon_alt)+ '\n')
#     f.write('moonaz '+str(moon_az)+ '\n')
#     f.write('moonphase '+str(mphase)+ '\n')
#     f.write('sunalt '+str(sun_alt)+ '\n')
#     f.write('sunaz '+str(sun_az)+ '\n')
#     f.write('seed 100'+str(fileCount)+ '\n')
#     f.write('vistime 60.0 \n')
#     f.write('filter 4 \n')
#     f.write('seeing '+str(seeing)+'\n')
#     f.write('obshistid 100'+str(fileCount)+ '\n')
#     if(addCat == 1):
#         for j in range(len(content)):
#             mag = float((content[j].split())[4])
#             if(mag > 20):
#                 f.write(content[j])
#     else:
#         f.write('stars 15.0 30.0 0.666 \n')
#         f.write('galaxies 15.0 30.0 0.666 \n')
# =============================================================================
    f.close()
    
    
    