#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:32:28 2023

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
#Read Redshifts
f=open('/home/dutta26/zphot_2390.out')
content = f.readlines()
f.close()
redShiftArr=[]
for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        redShiftArr.append(float((content[j].split())[1]))
        
coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
redShiftArr = np.array(redShiftArr)
# =============================================================================
# loc = np.where((coadd_data[:,10]> 10274) & (coadd_data[:,10]< 13283) &(coadd_data[:,11]> 17164) & (coadd_data[:,11]< 19073) &
#                (coadd_data[:,2]== 0)& (redShiftArr>0))[0]
# 
# =============================================================================

loc = np.where((coadd_data[:,10]> 13890) & (coadd_data[:,10]< 15240) &(coadd_data[:,11]> 17165) & (coadd_data[:,11]< 19000) &
               (coadd_data[:,2]== 0)& (redShiftArr>0))[0]

n, bins, patches = plt.hist(x=redShiftArr[loc], bins=40,histtype=u'step', color='b', label='e1')