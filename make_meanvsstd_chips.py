#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 21:59:13 2021

@author: dutta26
"""


# =============================================================================
# from astropy.io import fits 
# import numpy as np 
# import matplotlib.pyplot as plt
# f=fits.open('/scratch/halstead/d/dutta26/abell_2390/i/20171108T183209.4_abell_2390_odi_i.7190.fits')
# data = np.array(f[1].data)
# f.close()
# 
# bad = data[2100:2400 , 1100:1400]
# good = data[2100:2400 , 2100:2400]
# 
# bad = bad.flatten()
# good = good.flatten()
# 
# weightsb = np.ones_like(bad)/float(len(bad))
# weightsg = np.ones_like(good)/float(len(good))
# 
# #plt.hist(myarray, weights=weights)
# 
# n, bins, patches = plt.hist(x=bad, bins=100, color='r',histtype= 'step',
#                             density= True, label = 'Stripe')
# n, bins, patches = plt.hist(x=good, bins=100, color='b',histtype= 'step',
#                             density= True, label = 'Normal')
# #plt.grid(axis='y', alpha=0.75)
# plt.xlabel('Value')
# plt.ylabel('Frequency')
# plt.legend()
# 
# =============================================================================


from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt
from shutil import copy
import sys

bandList = ['r']
#folder =str(sys.argv[1])
folder = '/scratch/halstead/d/dutta26/abell_2390/'

meanArr = []
stdArr=[]
yearArr =[]
for band in bandList:
    loc = folder+band
    for file in os.listdir(folder+band+'/'):
        if('.weight' in file):
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
                    
                    
                    mean, median, stddev = sigma_clipped_stats(data, cenfunc=np.mean)
                    #if(mean > 1500 and (mean-stddev)< 300  and (mean-stddev)>100):
                    #    sys.exit()
                    meanArr.append(median)
                    stdArr.append(stddev)
                    yearArr.append(int(file[0:4]))
        
        
        f.close()

sys.exit()        
meanArr =np.array(meanArr)
stdArr = np.array(stdArr)
a = np.arange(1,2000, 1)
for j in range(len(meanArr)):
    if(yearArr[j] == 2017):
        plt.plot(meanArr, stdArr**2, 'b.', markersize = 1)
    elif(yearArr[j] == 2019):
        plt.plot(meanArr, stdArr**2, 'b+', markersize = 1)
    else:
        plt.plot(meanArr, stdArr**2, 'b*', markersize = 1)
plt.plot(a, a+100, 'r--', label = 'y= x+100')
#plt.plot(a, a, 'g--')
plt.xlim(-1000, 2000)
plt.ylim(-1000, 2000)
plt.ylabel('Variance')
plt.xlabel('Mean')
plt.legend()