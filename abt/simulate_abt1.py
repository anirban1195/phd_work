#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 08:41:48 2021

@author: anirban
"""

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
import subprocess

sample_img_loc  = '/home/dutta26/codes/abt/ztf.fits'
sample_exposure = 0
sample_aper = 0

f=fits.open(sample_img_loc)
data = np.array(f[0].data)
w = WCS(f[0].header)
ra = (f[0].header)['CRVAL1']
dec = (f[0].header)['CRVAL2']
exp = (f[0].header)['TOTEXPT']
zp = (f[0].header)['MAGZP']
f.close()
(y,x) =np.shape(data)
central_wcs = w.wcs_pix2world([[y/2, x/2]], 0)
#Conver nans to 0
data[np.isnan(data)] = 0
data[data<0] = 0

mean, median, sigma = sigma_clipped_stats(data, cenfunc= np.median)
bkg = mean
gain = (sigma**2)/mean

counts = np.sum(data)
mag=-2.5*np.log10(counts/gain/exp) + zp
print (mag)

data=data- bkg
data[data<0] = 0

counts = np.sum(data)
mag=-2.5*np.log10(counts/gain/exp) + zp
print (mag)


hdu = fits.PrimaryHDU(data)  
hdu.writeto('/home/dutta26/apps/johnrpeterson-phosim_core-7ba62a798e55/data/images/temp.fits', overwrite=True)


#Create run file 
f=open('/home/dutta26/codes/abt/run_temp.txt', 'w+')
f.write('rightascension '+str(central_wcs[0,0]) +'\n')
f.write('declination '+str(central_wcs[0,1]) +'\n')
f.write('rotskypos 0.0 '+'\n')
f.write('rottelpos 0.0 '+'\n')
f.write('vistime 300.0 '+'\n')
f.write('nsnap 1 '+'\n')
f.write('moonalt -90.0 '+'\n')
f.write('obshistid 1001 '+'\n')
f.write('object 0 193.85918 52.2645 13.00 ../sky/sed_flat.txt 0 0 0 0 0 0 temp.fits 1.0 0.0 none none '+'\n')
f.write('object 1 193.859 52.2654 16.8 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none '+'\n')

f.close()


phosimLoc = '/home/dutta26/apps/johnrpeterson-phosim_core-7ba62a798e55/'
outLoc = '/home/dutta26/codes/abt/'
bashCommand = './phosim /home/dutta26/codes/abt/run_temp.txt -i abt -c '+phosimLoc+'examples/quickbackground --thread=15 -o'+outLoc
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
output, error = process.communicate()








