#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 23:42:29 2020

@author: anirban
"""

import os
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
path ='/scratch/bell/dutta26/wiyn_sim (copy 1)/'

txtPath =[]
fitsPath =[]
npPath =[]
pkPath =[]
catPath =[]
pngPath=[]
for r, d, f in os.walk(path):
    for file in f:
        if('.gz' in file):
            continue
        if ( '.txt' in file or '.cat' in file or '.swarp' in file or '.ascii' in file):
            txtPath.append(os.path.join(r, file))
        elif '.fits' in file:
            fitsPath.append(os.path.join(r, file))
        elif '.npy' in file:
            npPath.append(os.path.join(r, file))
        elif '.pk' in file:
            pkPath.append(os.path.join(r, file))
        elif '.png' in file:
            pngPath.append(os.path.join(r, file))
       
        
            
print (len(fitsPath), len(txtPath), len(npPath), len(pkPath))

for files in fitsPath:
    f=fits.open(files)
    data= f[0].data
    #mean=np.nanmean(data)
    f.close()

for files in txtPath:
    f=open(files)
    content = f.readlines()
    f.close()

for files in npPath:
    frame_data = np.load(files)
    del frame_data

for files in pkPath:
    source_df = pd.read_pickle(files)
    
    
for files in pngPath:
    img= mpimg.imread(files)