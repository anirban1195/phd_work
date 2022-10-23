#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 07:51:05 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os, sys
from astropy.nddata import Cutout2D
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS

#file ='/scratch/halstead/d/dutta26/abell_2390/abell_ir_noise.fits'
file = '/scratch/halstead/d/dutta26/abell_2390/odi_sim_long.fits'
new_file = '/scratch/halstead/d/dutta26/abell_2390/odi_sim_long1.fits'
#f=fits.open(file, mode='update')
hdu = fits.open(file)[0]
wcs = WCS(hdu.header)
centralra = 328.4083
centraldec = 17.6697
temp = wcs.wcs_world2pix([[centralra, centraldec]], 0)
cent_x = temp[:,0]
cent_y = temp[:,1]
print (temp)
# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=(cent_x,cent_y), size=(10000,10000), wcs=wcs)

# Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = new_file
hdu.writeto(cutout_filename, overwrite=True)
