#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 11:54:13 2020

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
file = str(sys.argv[1])
new_file = str(sys.argv[2])
#f=fits.open(file, mode='update')
hdu = fits.open(file)[0]
wcs = WCS(hdu.header)

# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=(16000,17000), size=(23500,21500), wcs=wcs)

# Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = new_file
hdu.writeto(cutout_filename, overwrite=True)
