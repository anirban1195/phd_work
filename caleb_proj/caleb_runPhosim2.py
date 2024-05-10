#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:48:05 2024

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
import helper
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')



phosimLoc = '/home/dutta26/Downloads/phosim_release/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/bell/dutta26/abell_2390/'
#filt_list = ['g', 'r', 'u', 'z', 'i']
filt_list = [ 'i']
outFile= open('/home/dutta26/comparison_abell2390_bkg5.txt', 'w+')
outFile.write(' Zodaical at 70% and 21.7 and 22.2 ' + '\n')
counter =0


for folder in filt_list:
    for files in os.listdir(dataSet+folder):
        print (files)
        if('weight' in files or 'temp' in files ):
            continue
        
        
        f=fits.open(dataSet+folder+'/'+files)
        print (dataSet+folder+'/'+files)
        mjd = float((f[0].header)['MJD-MID'])
        airmass = float((f[0].header)['AIRMASS'])
        zp = (f[0].header)['MAGZERO']
        flx_scale = (f[0].header)['FLXSCALE']
        ra_header = (f[0].header)['RA']
        dec_header = (f[0].header)['DEC']
        c=SkyCoord(ra_header, dec_header, unit=(u.hourangle, u.deg))
        ra,dec = c.ra.degree, c.dec.degree
        filt = (f[0].header)['FILTER']
        if(mjd < 8567):
            continue
        
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
        for j in range(1,30):
            data = np.array(f[j].data)
            mean, median, std = sigma_clipped_stats(data)
            medianArr.append(median)
            stdArr.append(std)
        
        skyBkg = np.nanmedian(medianArr)
        skystd = np.nanmedian(stdArr)
        f.close()
        
        #Make test file 
        f=open(phosimLoc+'examples/test2', 'w+')
        f.write('rightascension  '+str(ra) + '\n')
        f.write('declination  '+str(dec) + '\n')
        f.write('mjd  '+str(mjd) + '\n')
        f.write('filter	'+str(filt_code) +'\n')
        f.write('nsnap	1' + '\n')
        f.write('seeing	2' + '\n')
        f.write('vistime 10' + '\n')
        f.write('object  0  '+str(ra)+'  '+str(dec)+'  ' +' 16 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none' + '\n')
        f.close()
        
        #os.chmod(phosimLoc+'examples/test1.txt', 0o777)
        #First make sure the outputLoc is empty
        for files in os.listdir(outLoc):
            os.remove(outLoc+files)
        
        #Make the bash command 
        files = os.listdir(outLoc)
        attempts = 1
        while(len(files) == 0 and attempts< 2):
            bashCommand = './phosim '+phosimLoc+'examples/test2 -i wiyn_odi -c '+phosimLoc+'examples/quickbackground -e 0 -o '+outLoc
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
        
        print (mean, median, std)
        loc = np.where(data>  (median+4.5*std))
        y = int(np.median(loc[0]))
        x = int(np.median(loc[1]))
        cut = data[y-50:y+50, x-50:x+50]
        flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
        if(flux == None or np.isnan(flux)):
            flux = 0.0
        
        print (flux)
        f.close()
        fileHandle.close()
        
        outLine = str(mean)+ ' , ' +str(median)+ ' , ' +str(std)+ ' , ' +str(skyBkg)+ ' , ' + str(skystd)+ ' , ' +str(filt_code) + ' , '+str(mjd)+' , '+str(flx_scale)+ ' , '+str(flux) + ' , '+str(airmass)
        outFile.write(outLine + ' \n')
        print (outLine)
        counter += 1
        if(counter >= 2):
            break
        
outFile.close()
        
