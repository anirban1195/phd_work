#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 12:31:11 2024

@author: dutta26
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:48:05 2024

@author: dutta26
"""



from astropy.io import fits
import numpy as np
import os,shutil
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



phosimLoc = '/home/dutta26/Downloads/phosim_core/'
outLoc = phosimLoc + 'output/'

dataSet ='/scratch/bell/dutta26/abell_2390/'
filt_list = ['g', 'z']

counter =0
lambdaArr = [350, 500, 600, 750, 875]
temp_arr = np.zeros((1000, 11))
j = 0
moon_alt_arr=np.zeros((1000))
for folder in filt_list:
    for file_name in os.listdir(dataSet+folder):
        
        if('weight' in file_name or 'temp' in file_name ):
            continue
        
        
        f=fits.open(dataSet+folder+'/'+file_name)
        
        mjd = float((f[0].header)['MJD-MID'])
        airmass = float((f[0].header)['AIRMASS'])
        zp = float((f[0].header)['MAGZERO'])
        flx_scale = float((f[0].header)['FLXSCALE'])
        ra_header = (f[0].header)['RA']
        dec_header = (f[0].header)['DEC']
        c=SkyCoord(ra_header, dec_header, unit=(u.hourangle, u.deg))
        ra,dec = c.ra.degree, c.dec.degree
        filt = (f[0].header)['FILTER']
        moon_alt = float((f[0].header)['MOON_ALT'])
# =============================================================================
#         if(mjd < 59189.036 or mjd > 59189.172):
#             continue
# =============================================================================
# =============================================================================
#         if(mjd < 59527.055 or mjd > 59527.255):
#             continue
# =============================================================================
        
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
        skylevel = float((f[0].header)['SKYBG'])
        skymag = float((f[0].header)['SKYMAG'])
        f.close()
        if(moon_alt>0):
            continue
        
        #Make test file 
        f=open(phosimLoc+'examples/test2', 'w+')
        f.write('rightascension  '+str(ra) + '\n')
        f.write('declination  '+str(dec) + '\n')
        f.write('mjd  '+str(mjd) + '\n')
        f.write('filter	'+str(filt_code) +'\n')
        f.write('nsnap	1' + '\n')
        f.write('seeing	1.6' + '\n')
        f.write('vistime 30' + '\n')
        f.write('object  0  '+str(ra)+'  '+str(dec)+'  ' +' 16 ../sky/sed_flat.txt 0 0 0 0 0 0 star none none' + '\n')
        #f.write('cloudcover 1' + '\n')
        #f.write('clouddepth 0.0' + '\n')
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
        
        phosim_zp = 16 + 2.5*np.log10(flux/30) - 2.5*np.log10(lambdaArr[filt_code]/500)
        phosim_skymag = phosim_zp - 2.5*np.log10(median*9.09*9.09/30)
        phosim_counts = median * 9.09*9.09/30
        print (flux)
        print (phosim_zp, zp, phosim_skymag, skymag, phosim_counts, skylevel*9.09*9.09/60, mjd, filt_code, airmass)
        f.close()
        fileHandle.close()
        
        #shutil.copy(outLoc+files[0], '/scratch/bell/dutta26/phosim_055045/neutral_0504_'+str(file_name))
        
        #outLine = str(mean)+ ' , ' +str(median)+ ' \, ' +str(std)+ ' , ' +str(skyBkg)+ ' , ' + str(skystd)+ ' , ' +str(filt_code) + ' , '+str(mjd)+' , '+str(flx_scale)+ ' , '+str(flux) + ' , '+str(airmass)
        #outLine = str(phosim_zp)+ ' , ' +str(phosim_skymag)+ ' , ' +str(phosim_counts)+ ' , ' + str(mjd)+ ' , ' +str(filt_code) + ' , '+str(mjd)+' , '+str(skylevel*9.09*9.09/60)+ ' , '+str(skymag) + ' , '+str(airmass) + ' , '+str(zp)  
        temp_arr[j,0:10] = phosim_zp, zp, phosim_skymag, skymag, phosim_counts, skylevel*9.09*9.09/60, mjd, filt_code, airmass, moon_alt
        moon_alt_arr[j] = moon_alt
        j += 1
        
        
np.save('/home/dutta26/2024/May/temp_neutral.npy', temp_arr)
#np.save('/home/dutta26/2024/Apr/moon_alt_arr.npy', moon_alt_arr)        
