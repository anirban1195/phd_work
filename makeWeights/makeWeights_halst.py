#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 08:04:07 2019

@author: anirban
"""
from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt
from shutil import copy
import sys

folder ='/home/anirban/test/'
segDataLoc ='/home/anirban/test/segment_stat.txt'

#Check if the weight exist. If so then delete them and create new
for files in os.listdir(folder):
    if('.weight' in files):
        os.remove(folder+files)

#Create new weights
for files in os.listdir(folder):
    if('.fits' in files):
        new_filename = files[:-5] + '.weight.fits'
        copy(folder+files , folder+new_filename)
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
        
for files in os.listdir(folder):
    if( '.weight' not in files):
         continue
    f=fits.open(folder+files, mode='update')
    for j in range(1,31):
         img = np.array(f[j].data)
         a = []
         d=[]
         #weight = np.zeros((4096,4096), dtype =np.float16)
         for y in range(8):
             for x in range(8):
                 xStart = 5+(x*508)
                 xEnd = xStart + 470
                 yStart = 5+ (y*505)
                 yEnd = yStart + 484
                 data = img[yStart:yEnd, xStart:xEnd]
                 mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
                 seg_id =int( "%02d" % j+str(y) +str(x))
                 index = np.where(segment_stat[:,0] == seg_id)[0][0]
                 
                 if(segment_stat[index,2] > 0.5*(segment_stat[index,1]+segment_stat[index,2]) or (segment_stat[index,2]+segment_stat[index,1])<40):
                     f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
                     continue
                 else:
                     if( stddev > 0 ):
                         if(len(meanArr) < 10):
                             f[j].data[yStart:yEnd, xStart:xEnd] =1/stddev
                             meanArr.append(mean)
                         else:
                             if((mean> 1.1*np.mean(meanArr)) or (mean< 0.9*np.mean(meanArr))):
                                 f[j].data[yStart:yEnd, xStart:xEnd] =0.0
                             else:
                                 f[j].data[yStart:yEnd, xStart:xEnd] =1/stddev
                                 meanArr.append(mean)
                     else:
                         f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
                             

f.flush()
                         
#                    arr[~ np.isnan(data)] = 1/(stddev**2)
# #                  f[j].data[yStart:yEnd, xStart:xEnd] = arr




# =============================================================================
# c=[]
# d_final=[]
# a_final =[]
# for files in os.listdir(folder):
#     if( 'abc' in files):
#         continue
#     f=fits.open(folder+files)
#     #print (arr)
#     for j in range(1,31):
#         img = np.array(f[j].data)
#         a = []
#         d=[]
#         #weight = np.zeros((4096,4096), dtype =np.float16)
#         for y in range(8):
#             for x in range(8):
#                 xStart = 5+(x*508)
#                 xEnd = xStart + 470
#                 yStart = 5+ (y*505)
#                 yEnd = yStart + 484
#                 data = img[yStart:yEnd, xStart:xEnd]
#                 mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
#                 #arr= np.zeros(np.shape(data))
#                 wt_inv = np.sqrt(abs(stddev**2 - mean))
#                 if(wt_inv< 2):
#                     wt_inv = 2
#                 a.append(mean)
#                 d.append(1/stddev)
#                 c.append(stddev**2)
#                 #a[x,y] = 1000/(stddev**3)
#         a=np.array(a)
#         d=np.array(d)
#         d_final.append(d)
#         a_final.append(a)
#     f.close()
# # =============================================================================
# #                 if( stddev > 0):
# #                     arr[~ np.isnan(data)] = 1/(stddev**2)
# #                 f[j].data[yStart:yEnd, xStart:xEnd] = arr
# # =============================================================================
# c=np.array(c)
# a_final = np.array(a_final)
# a_final = a_final.flatten()
# d_final = np.array(d_final)
# d_final = d_final.flatten()
# 
# plt.subplot(111)        
# plt.plot(a_final,c, 'b.', markersize =1 )
# plt.plot(np.arange(0,300,1), np.arange(0,300,1), 'r-')
# plt.plot(np.arange(0,300,1), np.arange(0,300,1)+70, 'g-')
# #plt.plot(d_final,c, 'r.', markersize =1 )
# plt.xlabel('Mean')
# plt.ylabel('Total Standadrd Deviation')
# plt.show()
# =============================================================================
