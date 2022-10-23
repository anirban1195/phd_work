#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 17:32:19 2022

@author: dutta26
"""
import numpy as np
raArr =np.arange(237.0083-(5*0.00555), 237.0083+(5*0.00555), 0.00555 )
decArr = np.arange(-15.2308 -(5*0.00555) ,  -15.2308 + (5*0.00555) , 0.005551)    
for j in range(30):
    fileName = '/scratch/halstead/d/dutta26/lsst/wiyn_cat1/'+str(j)+'.txt1_flat'
    f=open(fileName, 'w+')
    f.write('rightascension 237.0083  \n')    
    f.write('declination -15.2308  \n') 
    f.write('seed 100'+str(j)+ '\n')
    f.write('vistime 100.0 \n')
    f.write('filter 4 \n')
    f.write('obshistid 100'+str(j)+ '\n')
    seeing = np.random.uniform(0.6, 0.7)
    f.write('seeing '+str(seeing)[0:5]+'\n')    
    count = 1
    mag = 15
    for j in range(len(raArr)):
        mag = mag+1.2
        for k in range(len(decArr)):
            temp_str = str(count) + ' '+str(raArr[j])[0:9]+ ' '+ str(decArr[k])[0:9] + ' '+str(mag)[0:5]
            f.write('object ' +temp_str + ' ../sky/sed_flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 star none none \n')
            count += 1
    f.close()