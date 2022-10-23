#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 13:10:53 2020

@author: anirban
"""

from astropy.io import fits
import numpy as np
import os
import subprocess
import gzip
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
import time


phosimLoc = '/home/dutta26/apps/phosim-phosim_release-60db50aae970/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/halstead/d/dutta26/abell_2390/'
folder = [ 'g']
outFile= open('/home/dutta26/comparison_starflux.txt', 'w+')
f1=open('/home/dutta26/codes/caleb_proj/starList_radec.txt')
content = f1.readlines()
f1.close()



for sets in os.listdir(dataSet):
    if(sets not in folder):
        continue
    for files in os.listdir(dataSet+sets):
        if('20171108T190519.3_abell_2390_odi_g.7190.fits' not in files):
            continue
        print (files)
        if('weight' in files):
            continue
        f=fits.open(dataSet+sets+'/'+files)
        mjd = float((f[0].header)['MJD-MID'])
        print (dataSet+sets+'/'+files)
        outFile.write(dataSet+sets+'/'+files+'\n')
        for j in range(1,31):
            wcs_val = wcs.WCS(f[j].header)
            for k in range(len(content)):
                flag =0 
                ra = float(content[k].split(',')[0])
                dec = float(content[k].split(',')[1])
                [[x,y]] = wcs_val.wcs_world2pix([[ra,dec]], 0)
                x = int(round(x))
                y = int(round(y))
                tot =0 
                quality = 9999
                if(x>30 and x<4066 and y>30 and y<4066):
                    flag =1 
                    data = np.array(f[j].data)
                    cut = data[x-25:x+25, y-25:y+25]
                    #print (cut[25,25], x,y,k)
                    tot = np.nansum(cut) - (50*50*np.nanmean(data))
                    quality = np.nansum(np.isnan(cut))
                if(flag == 1):
                    print (cut[25,25], x,y,k, quality, tot, np.nanmean(data))
                    outFile.write(str(tot)+ ', ' +str(quality)+ ', ' +str(ra) + ', ' +str(dec) +', '+str(mjd) + ', ' +str(k) +'\n')
                
                
outFile.close()
                

                    
                    
            