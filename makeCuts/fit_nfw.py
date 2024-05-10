#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 16:39:33 2024

@author: dutta26
"""

from astropy.io import fits
import numpy as np
from astropy import wcs 
import helper 
import os,shutil
from scipy.optimize import curve_fit

f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_sf.fits')
data = f[0].data
w = wcs.WCS(f[0].header)
f.close()


x_cent = 234
y_cent = 359

size_y = 12
size_x = 10
cut = data[y_cent-size_y:y_cent+size_y, x_cent-size_x:x_cent+size_x]

y = np.arange(y_cent-size_y,y_cent+size_y)
x = np.arange(x_cent-size_x,x_cent+size_x)
X,Y = np.meshgrid(x,y)

temp = np.array(w.all_pix2world(X,Y, 0))
raList = temp[0,:,:].flatten()
decList = temp[1,:,:].flatten()

y_cut = cut.flatten()

#helper.get_gamma_kappa1([raList,decList], temp[0,int(size/2),int(size/2)], temp[1,int(size/2),int(size/2)], 0.5, 3)

popt,pcov = curve_fit(helper.get_kappa_fit, [raList,decList], cut.flatten(),
                      [328.40034536, 17.69826912, 0.5, 3])

print (popt)
#Find the residual

x_cent = 230
y_cent = 358

size_y = 24
size_x = 24
cut = data[y_cent-size_y:y_cent+size_y, x_cent-size_x:x_cent+size_x]

y = np.arange(y_cent-size_y,y_cent+size_y)
x = np.arange(x_cent-size_x,x_cent+size_x)
X,Y = np.meshgrid(x,y)

temp = np.array(w.all_pix2world(X,Y, 0))
raList = temp[0,:,:].flatten()
decList = temp[1,:,:].flatten()
residual = cut - helper.get_kappa_fit([raList,decList], popt[0], popt[1], popt[2], popt[3]).reshape(48,48)

shutil.copyfile('/scratch/bell/dutta26/abell_2390/kappa_sf.fits', '/scratch/bell/dutta26/abell_2390/kappa_sf_nfwResidual.fits')

f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_sf_nfwResidual.fits', mode='update')
f[0].data[y_cent-size_y:y_cent+size_y, x_cent-size_x:x_cent+size_x] = residual
f.flush()




# =============================================================================
# hdu = fits.PrimaryHDU(helper.get_kappa_fit([raList,decList], popt[0], popt[1], popt[2], popt[3]).reshape(56,48))
# hdu.writeto('/scratch/bell/dutta26/abell_2390/test.fits', overwrite=True)
# 
# =============================================================================


# =============================================================================
# #For the E mode
# 
# f=fits.open('/scratch/bell/dutta26/abell_2390/EMode_sf2.fits')
# data = f[0].data
# w = wcs.WCS(f[0].header)
# f.close()
# 
# 
# x_cent = 233
# y_cent = 358
# 
# size = 12
# cut = data[y_cent-size:y_cent+size, x_cent-size:x_cent+size]*4
# 
# y = np.arange(y_cent-size,y_cent+size)
# x = np.arange(x_cent-size,x_cent+size)
# X,Y = np.meshgrid(x,y)
# 
# temp = np.array(w.all_pix2world(X,Y, 0))
# raList = temp[0,:,:].flatten()
# decList = temp[1,:,:].flatten()
# 
# y_cut = cut.flatten()
# 
# #helper.get_gamma_kappa1([raList,decList], temp[0,int(size/2),int(size/2)], temp[1,int(size/2),int(size/2)], 0.5, 3)
# 
# popt,pcov = curve_fit(helper.get_gamma_fit, [raList,decList], cut.flatten(),
#                       [temp[0,int(size/2),int(size/2)], temp[1,int(size/2),int(size/2)], 0.5, 3])
# 
# print (popt)
# =============================================================================



