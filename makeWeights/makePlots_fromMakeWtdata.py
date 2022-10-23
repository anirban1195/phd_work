#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 10:14:07 2020

@author: anirban
"""

from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# =============================================================================
# g= 1/1.2308620689
# bandList = ['u','g','r','z']
# #bandList = ['u']
# content=[]
# for band in bandList:
#     f=open('/home/anirban/testbash/odi_'+band+'_data.txt')
#     content.append( f.readlines())
#     f.close()    
# mean=[]
# mean1=[]
# variance=[]
# variance1=[]
# mjdArr=[]
# bad =0
# good = 0 
# const= 100
# 
# segment_stat = []
# for j in range(1,31):
#     for y in range(8):
#             for x in range(8):
#                 seg_id = "%02d" % j+str(y) +str(x)
#                 segment_stat.append([int(seg_id), 0, 0])
# 
# 
# segment_stat = np.array(segment_stat)
# for k in range(len(content)):
#     for j in range(len(content[k])):
#         mjd=float(content[k][j].split(',')[2])
#         mjdArr.append(mjd)
#         temp_mean = float(content[k][j].split(',')[0])
#         temp_var = float(content[k][j].split(',')[1])
# # =============================================================================
# #         if(temp_var > 10000 or temp_var<0 or temp_mean< 20):
# #             continue
# # =============================================================================
#         mean.append(float(content[k][j].split(',')[0]))
#         variance.append(float(content[k][j].split(',')[1]))
#         
#         seg_id= int((content[k][j].split(',')[3])[1:5])
#         index = np.where(segment_stat[:,0] == seg_id)[0][0]
#         
#         if((g*temp_mean+ const)> temp_var and (g*temp_mean)< temp_var):
#             good += 1
#             segment_stat[index,1] += 1
#         else:
#             bad += 1
#             segment_stat[index,2] += 1
#             
# 
# mean=np.array(mean) 
# variance=np.array(variance) 
# 
# hist_data=[]
# for j in range(len(segment_stat)):
#     tot = segment_stat[j, 1] +segment_stat[j, 2]
#     if(tot< -40):
#         continue
#     else:
#         goodPercent = segment_stat[j, 1]/tot
#         #print (tot)
#         hist_data.append(goodPercent)
# 
# n, bins, patches=plt.hist(hist_data, 20, facecolor='blue', alpha=0.8)
# plt.xlabel('Fraction of cases when segement is good')
# plt.ylabel('No of segments')
# 
# =============================================================================
# =============================================================================
# plt.subplot(111)       
# plt.plot(mean,variance, 'b.', markersize =1 )
# plt.plot(np.arange(0,3000,1), g*np.arange(0,3000,1), 'r-')
# plt.plot(np.arange(0,3000,1), g*np.arange(0,3000,1)+const, 'r-')
# plt.xlabel('Mean')
# plt.ylabel('Standadrd Deviation^2')
# plt.show()
# percent = good/(good+bad)
# =============================================================================
#print ('Percent good data is ',percent)


# =============================================================================
# from astropy.io import fits 
# f=fits.open('/home/anirban/testbash/test1.fits', mode="update")
# f1 = fits.open('/home/anirban/testbash/20191006T191934.3_Abell2390_odi_g.8445.fits')
# good=0
# bad = 0
# meanArr=[]
# 
# for j in range(1,31):
#         img = np.array(f1[j].data)
#         for y in range(8):
#             for x in range(8):
#                 xStart = 5+(x*508)
#                 xEnd = xStart + 470
#                 yStart = 5+ (y*505)
#                 yEnd = yStart + 484
#                 data = img[yStart:yEnd, xStart:xEnd]
#                 mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
#                 if(mean>0):
#                     meanArr.append(mean)
# 
# global_mean = np.median(meanArr)
# 
# 
# 
# 
# for j in range(1,31):
#         img = np.array(f1[j].data)
#         for y in range(8):
#             for x in range(8):
#                 xStart = 5+(x*508)
#                 xEnd = xStart + 470
#                 yStart = 5+ (y*505)
#                 yEnd = yStart + 484
#                 data = img[yStart:yEnd, xStart:xEnd]
#                 
#                 seg_id =int( "%02d" % j+str(y) +str(x))
#                 
#                 index = np.where(segment_stat[:,0] == seg_id)[0][0]
#                 mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
#                 if(segment_stat[index,2] > 0.5*(segment_stat[index,1]+segment_stat[index,2]) or (segment_stat[index,2]+segment_stat[index,1])<10
#                    ):
#                     f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
#                     bad += 1
#                 else:
#                     
#                     if(stddev**2 >mean):
#                         f[j].data[yStart:yEnd, xStart:xEnd] = 1/np.sqrt(stddev**2-mean)
#                     else:
#                         f[j].data[yStart:yEnd, xStart:xEnd] = 0.0
#                     good += 1
# f.flush()
# f1.close()              
# =============================================================================
# =============================================================================
# f=open('/home/anirban/segment_stat.txt', 'w+')
# for j in range(len(segment_stat)):  
#     seg_id = segment_stat[j][0]
#     good = segment_stat[j][1]
#     bad = segment_stat[j][2]
#     f.write(str(seg_id) +', '+str(good) +', '+str(bad) + ' \n')
# f.close()
# =============================================================================
f1 = fits.open('/home/anirban/testbash/20191006T191934.3_Abell2390_odi_g.8445.fits')
img  =np.array(f1[1].data)
x=7
y=6
xStart = 5+(x*508)
xEnd = xStart + 470
yStart = 5+ (y*505)
yEnd = yStart + 484
data1= img[yStart:yEnd, xStart:xEnd]
x=4
y=1
xStart = 5+(x*508)
xEnd = xStart + 470
yStart = 5+ (y*505)
yEnd = yStart + 484
data2= img[yStart:yEnd, xStart:xEnd]

print(sigma_clipped_stats(data1, cenfunc=np.mean))
print(sigma_clipped_stats(data2, cenfunc=np.mean))

data1 = data1.flatten()
data2 = data2.flatten()
data1 = np.clip(data1, -250, 3500)
data2 = np.clip(data2, -250, 3500)
num_bins = 100
n, bins, patches=plt.hist(data1, num_bins, facecolor='blue', alpha=0.3)
n, bins, patches=plt.hist(data2, num_bins,facecolor='red', alpha=0.3)
plt.show()
f1.close()
