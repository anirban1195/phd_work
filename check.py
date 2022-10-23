#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 23:42:29 2020

@author: anirban
"""

import os
from astropy.io import fits
import numpy as np
path ='/scratch/halstead/d/dutta26/'

txtPath =[]
fitsPath =[]
npPath =[]
for r, d, f in os.walk(path):
    for file in f:
        if('.gz' in file):
            continue
        if '.txt' in file:
            txtPath.append(os.path.join(r, file))
        elif '.fits' in file:
            fitsPath.append(os.path.join(r, file))
        elif '.npy' in file:
            npPath.append(os.path.join(r, file))
        
            
print (len(fitsPath))
print (len(txtPath))
for files in fitsPath:
    f=fits.open(files)
    data= np.array(f[0].data)
    #mean=np.nanmean(data)
    f.close()

for files in txtPath:
    f=open(files)
    content = f.readlines()
    f.close()

for files in npPath:
    frame_data = np.load(files)
    del frame_data
