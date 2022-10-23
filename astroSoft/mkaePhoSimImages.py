#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 09:49:27 2020

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import os
import subprocess
import gzip
from astropy import units as u
from astropy.coordinates import SkyCoord
import time
from astropy.stats import sigma_clipped_stats


phosimLoc = '/home/dutta26/apps/phosim-phosim_release_5_3_10/'
outLoc = phosimLoc + 'output/'
newOutLoc = '/scratch/halstead/d/dutta26/images/'
visList = [10]
bkgList =[0]
for bkg in bkgList:
    for vis in visList:
        #Make test file 
        f=open(phosimLoc+'examples/star1', 'w+')
        f.write('rightascension 0'+' \n')
        f.write('declination 0'+' \n')
        f.write('filter	2'+'\n')
        f.write('nsnap	1' + '\n')
        f.write('vistime ' +str(vis) +  '\n')
        f.write('seeing 0.0'+' \n')
        f.write('obshistid 999'+' \n')
        #seed=np.random.randint(0,1000)
        seed=100
        f.write('seed '+str(seed)+' \n')
        f.write('object  0 0.0 0.0 17 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none' + '\n')
        f.close()
        
        for j in range(20):
            #Clean output loc
            for files in os.listdir(outLoc):
                os.remove(outLoc+files)
            center = 512
            #Run phosim 
            bashCommand = './phosim '+phosimLoc+'examples/star1 -c '+phosimLoc+'examples/nobackground -e 0 -o'+outLoc
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
            output, error = process.communicate()
            #Read and cut images
            files = os.listdir(outLoc)
            fileHandle = gzip.open(outLoc+files[0])
            f=fits.open(fileHandle)
            data = np.array(f[0].data)
            f.close()
            fileHandle.close()
            z = data[center-25:center+25, center-25:center+25]
            print (np.sum(data))
            #Rewrite to new loc
            totNoise =int(bkg*50*50)
            x= np.random.randint(0, 50, size=totNoise)
            y= np.random.randint(0, 50, size=totNoise)
            np.add.at(z,(x,y),1)
            hdu = fits.PrimaryHDU(z)  
            hdu.writeto(newOutLoc + 'vis_'+str(vis)+ '_bkg_'+str(bkg)+'_no_'+ str(j)+'.fits', overwrite=True)
            del z,data,x,y
        
    
        
        
        