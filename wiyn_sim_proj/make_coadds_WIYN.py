#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 08:06:14 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import os, sys,shutil, gzip
from astropy.nddata import Cutout2D
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS
import subprocess


folderLoc = "/scratch/halstead/d/dutta26/abell_2390/odi_img_realsee/"
swarpLoc = '/home/dutta26/apps/bin/bin/'
tempLoc = "/scratch/halstead/d/dutta26/abell_2390/odi_img1/"
for files in os.listdir(folderLoc):
    for files1 in os.listdir(folderLoc+"/"+files):
        if('CELL' in files1 or "coadd" in files1):
            continue
        else:
            with gzip.open(folderLoc+"/"+files+'/'+files1, 'rb') as f_in:
                with open(tempLoc+files1[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
    f=open("/home/dutta26/fileList_wiyn.ascii", 'w+')
    for files1 in os.listdir(tempLoc):
        f.write(tempLoc+'/'+files1 +'\n')
    f.close()
    
    #bashCommand = './swarp @/home/dutta26/fileList_wiyn.ascii -c /home/dutta26/default1.swarp -IMAGEOUT_NAME '+folderLoc+"/"+files+'/final_coadd.fits -WEIGHTOUT_NAME '+folderLoc+"/"+files+'/final_coadd.weight.fits -RESAMPLE_DIR /scratch/halstead/d/dutta26/lsst/' 
    #process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    #output, error = process.communicate()
    
    a=os.listdir(tempLoc)
    shutil.copy(tempLoc+'/'+a[0], folderLoc+"/"+files+'/sample.fits')
                          
    for files1 in os.listdir(tempLoc):
        os.remove(tempLoc+'/'+files1)
        
    
        
    