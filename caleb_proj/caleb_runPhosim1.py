#!/usr/bin/env python30
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 09:57:31 2020

@author: anirban
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


phosimLoc = '/home/dutta26/apps/phosim-phosim_release_5_3/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/halstead/d/dutta26/abell_115_1/'
folder = ['images']
outFile= open('/home/dutta26/comparison_abell115_1_z0a0e0.txt', 'w+')
outFile.write(' Zodaical at 70% and 23 and 23.5 ' + '\n')
counter =0


for sets in os.listdir(dataSet):
    if(sets not in folder):
        continue
    for files in os.listdir(dataSet+sets):
        print (files)
        if('weight' in files):
            continue
        f=fits.open(dataSet+sets+'/'+files)
        print (dataSet+sets+'/'+files)
        mjd = float((f[0].header)['MJD-MID'])
        ra_header = (f[0].header)['RA']
        dec_header = (f[0].header)['DEC']
        c=SkyCoord(ra_header, dec_header, unit=(u.hourangle, u.deg))
        ra,dec = c.ra.degree, c.dec.degree
        filt = (f[0].header)['FILTER']
        #if(mjd > 58067):
        #    continue
        if(filt == 'odi_u'):
            filt_code =0
        elif(filt == 'odi_g'):
            filt_code = 1
        elif(filt == 'odi_r'):
            filt_code = 2
        elif(filt == 'odi_i'):
            filt_code = 3
        elif(filt == 'odi_z'):
            filt_code = 4
        else:
            filt_code = 5
        medianArr=[]
        stdArr=[]
        for j in range(1,31):
            data = np.array(f[j].data)
            mean, median, std = sigma_clipped_stats(data)
            medianArr.append(median)
            stdArr.append(std)
        #skyBkg =(f[0].header)['SKYBG']
        #skystd = (f[0].header)['SKYBGSTD']
        #print (medianArr)
        #print (stdArr)
        skyBkg = np.nanmean(medianArr)
        skystd = np.nanmean(stdArr)
        f.close()
        
        #Make test file 
        f=open(phosimLoc+'examples/test2', 'w+')
        f.write('rightascension  '+str(ra) + '\n')
        f.write('declination  '+str(dec) + '\n')
        f.write('mjd  '+str(mjd) + '\n')
        f.write('filter	'+str(filt_code) +'\n')
        f.write('nsnap	1' + '\n')
        f.write('vistime	60' + '\n')
        f.write('object  0  '+str(ra)+'  '+str(dec)+'  ' +' 35 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none' + '\n')
        f.close()
        
        #os.chmod(phosimLoc+'examples/test1.txt', 0o777)
        #First make sure the outputLoc is empty
        for files in os.listdir(outLoc):
            os.remove(outLoc+files)
        
        #Make the bash command 
        files = os.listdir(outLoc)
        attempts = 1
        while(len(files) == 0 and attempts< 3):
            bashCommand = './phosim '+phosimLoc+'examples/test2 -i wiyn_odi -c '+phosimLoc+'examples/quickbackground -e 0 -o'+outLoc
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
            output, error = process.communicate()
            time.sleep(10)
            files = os.listdir(outLoc)
            attempts = attempts + 1
        if(len(files) == 0 and attempts>= 2):
            continue
        fileHandle = gzip.open(outLoc+files[0])
        f=fits.open(fileHandle)
        data = np.array(f[0].data)
        mean, median, std = np.nanmean(data), np.nanmedian(data), np.nanstd(data)
        f.close()
        fileHandle.close()
        outLine = str(mean)+ ' , ' +str(median)+ ' , ' +str(std)+ ' , ' +str(skyBkg)+ ' , ' + str(skystd)+ ' , ' +str(filt_code) + ' , '+str(mjd)
        outFile.write(outLine + '\n')
        
        counter += 1
        #if(counter > 1):
        #    break
outFile.close()
        
