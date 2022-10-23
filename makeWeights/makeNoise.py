#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 08:51:10 2020

@author: dutta26
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt
from shutil import copy
import sys

#folder =str(sys.argv[1])
folder ='/scratch/halstead/d/dutta26/abell_2390/r_noise/'
segDataLoc ='/home/dutta26/codes/makeWeights/segment_stat_new.txt'


#Read segement stats from file
f=open(segDataLoc)
content = f.readlines()
f.close()
segment_stat=[]
for j in range(len(content)):
    temp = content[j].split(',')
    segment_stat.append([int(temp[0]), int(temp[1]), int(temp[2])])
segment_stat = np.array(segment_stat)       
meanArr=[]

        
for files in os.listdir(folder)[100:212]:
    print (files)
    if( '.weight' in files):
         continue
    f=fits.open(folder+files, mode='update')
    back = float((f[0].header)['SKYBG'])
    flag = 0
    for j in range(1,31):
         img = np.array(f[j].data)
         a = []
         d=[]
         for y in range(8):
             for x in range(8):
                 xStart = 5+(x*508)
                 xEnd = xStart + 470
                 yStart = 5+ (y*505)
                 yEnd = yStart + 484
                 data = img[yStart:yEnd, xStart:xEnd]
                 mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.median)
                                  
                 seg_id =int( "%02d" % j+str(y) +str(x))
                 index = np.where(segment_stat[:,0] == seg_id)[0][0]
                 #If half the data is bad or half of total (80) data is accounted for make wt 0
                 if(segment_stat[index,2] > 0.5*(segment_stat[index,1]+segment_stat[index,2]) or (segment_stat[index,2]+segment_stat[index,1])< 100 or flag==1):
                     f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
                     continue
                 else:
                     if( stddev > 0 ):
                         if(mean<0.5*back):
                             f[j].data[yStart:yEnd, xStart:xEnd] =0.0
                             
                         else:
                             f[j].data[yStart:yEnd, xStart:xEnd] =np.random.normal(mean, stddev, (484,470))
                             
                     else:
                         f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
                             

f.flush()
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         #weight = np.zeros((4096,4096), dtype =np.float16)
# =============================================================================
#          gain = float((f[j].header)['GAINX'])
#          for y in range(8):
#              for x in range(8):
#                  xStart = 5+(x*508)
#                  xEnd = xStart + 470
#                  yStart = 5+ (y*505)
#                  yEnd = yStart + 484
#                  data = img[yStart:yEnd, xStart:xEnd]
#                  data = data * gain
#                  mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
#                  if(mean<0):
#                      mean=0.1
#                  stddev = np.sqrt(mean)
#                  #Convert back to ADU
#                  mean = mean /gain
#                  stddev = stddev/ gain
#                  seg_id =int( "%02d" % j+str(y) +str(x))
#                  index = np.where(segment_stat[:,0] == seg_id)[0][0]
#                  
#                  if(segment_stat[index,2] > 0.5*(segment_stat[index,1]+segment_stat[index,2]) or (segment_stat[index,2]+segment_stat[index,1])<40):
#                      f[j].data[yStart:yEnd, xStart:xEnd] = np.nan
#                      continue
#                  else:
#                      if( stddev > 0 ):
#                          if(mean<0.5*back):
#                              f[j].data[yStart:yEnd, xStart:xEnd] =np.nan
#                              
#                          else:
#                              f[j].data[yStart:yEnd, xStart:xEnd] =np.random.normal(mean, stddev, (484,470))
#                      else:
#                          f[j].data[yStart:yEnd, xStart:xEnd] = np.nan
#                              
# 
# f.flush()
#         
# =============================================================================
        