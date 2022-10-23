#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 10:47:23 2022

@author: dutta26
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import os,shutil
from astropy.stats import sigma_clipped_stats

main = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_data.npy')

cut12 = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_datasz13.npy')

cut15 = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_datasz16.npy')

cut20 = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_datasz20.npy')

cut45 = np.load('/scratch/halstead/d/dutta26/abell_2390/wiyn_datasz45.npy')


f=open('/home/dutta26/Downloads/sextractor-2.19.5/share/sextractor/test_all.cat')
content = f.readlines()
f.close()
raList =[]
decList=[]
sex_xx=[]
sex_yy=[]
sex_xy =[]
elonList=[]
Raindex=5
Decindex=6
elonIndex = 12
xxIndex = 7
yyIndex = 8
xyIndex = 9
#Create a list of stars with ra and dec 
for j in range(len(content)):
    #Skip the comment lines
    if((content[j].split()[0]) == '#'):
        continue

    raList.append(float(content[j].split()[Raindex])) 
    decList.append(float(content[j].split()[Decindex])) 
    sex_xx.append(float(content[j].split()[xxIndex]))
    sex_yy.append(float(content[j].split()[yyIndex]))
    sex_xy.append(float(content[j].split()[xyIndex]))
    elonList.append(float(content[j].split()[elonIndex]))

size12 = np.sqrt(cut12[:,7] + cut12[:,8])
size15 = np.sqrt(cut15[:,7] + cut15[:,8])
size20 = np.sqrt(cut20[:,7] + cut20[:,8])
size45 = np.sqrt(cut45[:,7] + cut45[:,8])

size = np.sqrt(main[:,7] + main[:,8])
size_sex =  np.sqrt(sex_xx +sex_yy)

#sizeRatio = sexSize/size
loc = np.where((main[:,3]>20000) & (main[:,3]<500000) & (main[:,2] == 1 ))[0]
print (len(loc))
temp1=[]
temp2=[]
mean,median,std = sigma_clipped_stats(size[loc])
for j in loc:
    if(size15[j]> 8 or size15[j]<1 or size[j]<(median-2*std) or size[j]>(median+2*std)):
        continue
    else:
        temp1.append(size[j])
        temp2.append(size15[j])
        
a=np.polyfit(temp1, temp2, deg=1)
print (len(temp1))
b=np.arange(1,10,1)
plt.plot(size[loc], size15[loc], 'r+', markersize = 1)
plt.plot(b, a[0]*b+ a[1], 'r--',label='30x30 cut')
print (median,std)

temp1=[]
temp2=[]
mean,median,std = sigma_clipped_stats(size[loc])
for j in loc:
    if(size20[j]> 8 or size20[j]<1 or size[j]<(median-2*std) or size[j]>(median+2*std)):
        continue
    else:
        temp1.append(size[j])
        temp2.append(size20[j])
        
a=np.polyfit(temp1, temp2, deg=1)
print (len(temp1))
b=np.arange(1,10,1)
plt.plot(size[loc], size20[loc], 'b+', markersize = 1)
plt.plot(b, a[0]*b+ a[1], 'b--', label='40x40 cut')
print (median,std)



temp1=[]
temp2=[]
mean,median,std = sigma_clipped_stats(size[loc])
for j in loc:
    if(size12[j]> 8 or size12[j]<1 or size[j]<(median-2*std) or size[j]>(median+2*std)):
        continue
    else:
        temp1.append(size[j])
        temp2.append(size12[j])
        
a=np.polyfit(temp1, temp2, deg=1)
print (len(temp1))
b=np.arange(1,10,1)
plt.plot(size[loc], size12[loc], 'g+', markersize = 1)
plt.plot(b, a[0]*b+ a[1], 'g--', label='24x24 cut')
print (median,std)




temp1=[]
temp2=[]
mean,median,std = sigma_clipped_stats(size[loc])
for j in loc:
    if(size45[j]> 8 or size45[j]<1 or size[j]<(median-2*std) or size[j]>(median+2*std)):
        continue
    else:
        temp1.append(size[j])
        temp2.append(size45[j])
        
a=np.polyfit(temp1, temp2, deg=1)
print (len(temp1))
b=np.arange(1,10,1)
plt.plot(size[loc], size45[loc], 'k+', markersize = 1)
plt.plot(b, a[0]*b+ a[1], 'k--', label='90x90 cut')
print (median,std)

#plt.plot(size, size, 'k.', markersize = 1)
plt.legend()
plt.xlabel('Size with best estimates')
plt.ylabel('Size with fixed cutout size')
plt.title('Stars with flux 20k - 500k ')
