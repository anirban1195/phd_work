#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 16:46:50 2022

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

band = 'z'
fileList = os.listdir('/scratch/bell/dutta26/abell_2390/'+band +'/')
f= open('/home/dutta26/codes/Merged_A2390_proj/make_sf_' +band +'.sub', 'w+')
f.write('#!/bin/sh -l \n' )
f.write('# FILENAME: make_sf_' +band +'.sub \n')
f.write('source ~/.bashrc \n')
for file in fileList:
    
    if('.weight' in file or 'temp' in file):
        continue
    print (file)
    
    
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    #swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    swarpCommand = 'swarp /scratch/bell/dutta26/abell_2390/'+band +'/'+file +' -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file)
    f.write(swarpCommand + '\n')
    
f.close()    