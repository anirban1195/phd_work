#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  3 08:38:48 2024

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
#import wquantiles
from astropy import wcs
import matplotlib.pyplot as plt

f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_sf.fits')
data= f[0].data 
f.close()

f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_error_sf.fits')
data_err= f[0].data 
f.close()
# =============================================================================
# data=data[350-100:350+100, 320-100:320+100]
# 
# print (np.std(data.flatten()), np.mean(data.flatten()), np.median(data.flatten()))
# mean, med, std = sigma_clipped_stats(data.flatten())
# print (mean,med,std)
# =============================================================================
a=[]
b=[]
a_err=[]
x_mid = 232
y_mid = 357
for j in range(3,75, 2):
    a.append(np.sum(data[y_mid-j:y_mid+j+1, x_mid-j:x_mid+j+1])*4.23*5.59e-4)
    print (np.sum(data[y_mid-j:y_mid+j+1, x_mid-j:x_mid+j+1]))
    a_err.append(np.sum(data_err[y_mid-j:y_mid+j+1, x_mid-j:x_mid+j+1])*4.23*5.59e-4 )
    b.append(j*50*0.11*3.7)
a=np.array(a)    
plt.errorbar(b, a, yerr=a_err, fmt='r.', label='This work')


f=open('/home/dutta26/2390_kaiser.txt')
content = f.readlines()
f.close()
ksb_mass =[]
ksb_mass_up =[]
ksb_mass_down =[]
ksb_radius =[]
for j in range(11):
    
    
    #print (float(temp[0][0:6]), float(temp[1][0:6]))
    temp=content[11+j].split()
    ksb_mass_up.append((float(temp[1][0:6]))/0.7)
    temp=content[22+j].split()
    ksb_mass_down.append((float(temp[1][0:6]))/0.7)
    temp=content[j+1].split()
    temp=content[j].split()
    ksb_mass.append(float(temp[1][0:6])/0.7)
    ksb_radius.append(float(temp[0][0:6])*3.7)

downErr = np.array(ksb_mass) - np.array(ksb_mass_down)
upErr = np.array(ksb_mass_up) - np.array(ksb_mass)
plt.errorbar(ksb_radius, ksb_mass,  fmt='.k', label='SKF', yerr =[upErr, downErr])
plt.ylabel(r'Mass in units of $10^{15} $M$_{\odot}$')
plt.xlabel('Radius in kpc')
plt.legend()