#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 22:50:36 2021

@author: dutta26
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt
from shutil import copy
import sys

bandList = ['g', 'i', 'r']
#folder =str(sys.argv[1])
folder = '/scratch/bell/dutta26/abell_2390/'

segment_stat = []
for j in range(1,31):
    for y in range(8):
        for x in range(8):
            seg_id = "%02d" % j+str(y) +str(x)
            segment_stat.append([int(seg_id), 0, 0])

segment_stat = np.array(segment_stat)


for band in bandList:
    loc = folder+band
    for file in os.listdir(folder+band+'/'):
        if('.weight' in file or '2023' in file or 'temp' in file):
            continue
        f=fits.open(folder+band+'/'+file)
        for j in range(1,31):
            
            img = np.array(f[j].data)
            #img[np.isnan(img)] = 0
            for y in range(8):
                for x in range(8):
                    seg_id = "%02d" % j+str(y) +str(x)
                    xStart = 5+(x*508)
                    xEnd = xStart + 470
                    yStart = 5+ (y*505)
                    yEnd = yStart + 484
                    data = img[yStart:yEnd, xStart:xEnd]
                    index = np.where(segment_stat[:,0] == int(seg_id))[0][0]
                    
                    mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
                    if(mean < 0 or median<0 or stddev < 0 ):
                        segment_stat[index, 2] += 1
                        continue
                    
                    if(stddev**2 < (median+300) and stddev**2 > median*0.9): #Changed here
                        segment_stat[index, 1] += 1
                    else:
                        segment_stat[index, 2] += 1
                    
        
        
        
        f.close()

f = open ('/home/dutta26/codes/makeWeights/segment_stat_new1.txt', 'w+')        
for j in range(len(segment_stat)):
    line = str(segment_stat[j,0]) + ', ' + str(segment_stat[j,1]) + ', ' + str(segment_stat[j,2]) + '\n'
    f.write(line)
    
f.close()
