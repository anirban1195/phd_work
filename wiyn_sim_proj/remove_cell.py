#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 09:06:35 2022

@author: dutta26
"""

import os
loc = '/scratch/halstead/d/dutta26/abell_2390/odi_img_star1/'
#loc = '/scratch/halstead/d/dutta26/lsst/lsst2/'
swarpFile = '/home/dutta26/fileList_wiyn_star.ascii'
f=open(swarpFile, 'w+')
for files in os.listdir(loc):
    f.write(loc+files+'\n')
    
f.close()


loc = '/scratch/halstead/d/dutta26/abell_2390/odi_img1/'
swarpFile = '/home/dutta26/fileList_wiyn.ascii'
f=open(swarpFile, 'w+')
for files in os.listdir(loc):
    f.write(loc+files+'\n')
    
f.close()

loc = '/scratch/halstead/d/dutta26/abell_2390/odi_img_long1/'
swarpFile = '/home/dutta26/fileList_wiyn_long.ascii'
f=open(swarpFile, 'w+')
for files in os.listdir(loc):
    f.write(loc+files+'\n')
    
f.close()
    
    