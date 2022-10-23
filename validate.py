#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 22:46:05 2021

@author: dutta26
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
import subprocess
from astropy.io import fits
band = 'i'
frame_data = np.load('/home/dutta26/codes/singleFrame_'+band+'.npy')
coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
coadd_data = np.array(coadd_data)
coadd_flux = coadd_data[:,3]
a,b,c = np.shape(frame_data)
arr = []
mjd = []
airmass =[]
bkg = []
cnt = 0
for file in os.listdir('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'):
    
    if('.weight' in file):
        continue
    
    print (file)
    #Run swarp to make a good image
    f= open('/home/dutta26/temp.ascii', 'w+')
    f.write('/scratch/halstead/d/dutta26/abell_2390/'+band +'/'+file)
    f.close()
    
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
    #Read the swarp output
    f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
    data = np.array(f[0].data)
    f.close()  
    
    flux_dev =  (frame_data[cnt,:,3] ) /coadd_flux  
    valid = np.where((flux_dev>0.5)& (flux_dev<1.5) & (frame_data[cnt,:,2] == 1) & 
                     (frame_data[cnt,:,15] == 0) &(frame_data[cnt,:,16] == 0) & 
                     (frame_data[cnt,:,17] == 0) & (coadd_data[:,3] > 80) & (coadd_data[:,3] < 100))
    
    
    for index in valid[0]:
        x = frame_data[cnt,index,13]
        y = frame_data[cnt,index,14]
        cut = data[int(y)-30:int(y)+30, int(x)-30: int(x)+30]
        print (frame_data[cnt,index,9], frame_data[cnt,index,8])
        hdu = fits.PrimaryHDU(cut)
        hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/'+str(index)+'_'+str(coadd_data[index,3]) + '.fits', overwrite= True)
    
    
    break
    
    
    