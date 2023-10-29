#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 07:17:54 2023

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

file ='/scratch/bell/dutta26/wiyn_sim/median_coadds/ir_coadd_med.fits'
new_file = '/scratch/bell/dutta26/wiyn_sim/median_coadds/ir_coadd_med_crop.fits'
#file = str(sys.argv[1])
#new_file = str(sys.argv[2])
#f=fits.open(file, mode='update')
hdu = fits.open(file)[0]
wcs = WCS(hdu.header)

# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=(10500,10500), size=(10000,10000), wcs=wcs)

# Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = new_file
hdu.writeto(cutout_filename, overwrite=True)
