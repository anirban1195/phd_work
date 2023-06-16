#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 09:41:02 2023

@author: anirban
"""
import os
import gzip
import higer_test
import numpy as np
from astropy.modeling.models import Sersic2D, Moffat2D, Gaussian2D
import matplotlib.pyplot as plt
import measure_pythonV_dist
from scipy.signal import convolve2d
from astropy.io import fits
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
import shutil
import subprocess
phosimLoc = '/home/dutta26/Downloads/phosim_core/'
id_n = '6'
sersicList= [0.5, 2, 4, 6]
eDiffArr =[]
eErrArr =[]
kArr =[]
kErrArr=[]
for sersicIndex in sersicList:
    tempeArr= []
    tempkArr =[]
    for j in range(20):
        fileList = os.listdir('/home/dutta26/Downloads/phosim_core/output/')
        if('/home/dutta26/Downloads/phosim_core/output/data'+id_n+'.fits' in fileList ):
            print ('aa')
            os.remove('/home/dutta26/Downloads/phosim_core/output/data'+id_n+'.fits')
            os.remove('/home/dutta26/Downloads/phosim_core/output/generic_e_'+id_n+'_f3_chip_E000.fits.gz')
            os.remove('/home/dutta26/Downloads/phosim_core/output/generic_a_'+id_n+'_f3_chip_amplifier_E000.fits.gz')

        #Make catalog
        f=open(phosimLoc+'examples/small_catalog'+id_n, 'w+')
        f.write('rightascension 0 \n')
        f.write('declination 0 \n')
        f.write('vistime 30 \n')
        f.write('moonalt -90.0 \n')
        f.write('obshistid '+id_n+' \n')
        f.write('filter 3 \n')
        f.write('object 0 0.0 0.0 14 ../sky/sed_flat.txt 0 0 0 0 0 0 sersic2d 12.0 8.0 0.0 '+str(sersicIndex)+' none none \n')
        f.write('object 1 0.03 0.03 16 ../sky/sed_flat.txt 0 0 0 0 0 0 point none none \n')
        f.close()
        
        
        #Run Phosim
        bashCommand = './phosim '+phosimLoc+'examples/small_catalog'+id_n+' -c '+phosimLoc+'examples/nobackground -t4'
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, cwd=phosimLoc)
        output, error = process.communicate()
        

        #Extract files
        with gzip.open('/home/dutta26/Downloads/phosim_core/output/generic_e_'+id_n+'_f3_chip_E000.fits.gz', 'rb') as f_in:
            with open('/home/dutta26/Downloads/phosim_core/output/data'+id_n+'.fits', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        
        
        print ('Success')
        f =fits.open('/home/dutta26/Downloads/phosim_core/output/data'+id_n+'.fits')
        img = f[0].data
        f.close()
        
        
        
        psf = img[302-30:302+30, 722-30: 722+30]
        if(np.sum(psf) <10000):
            psf = img[722-30: 722+30, 302-30:302+30]
        img =img[512-50:512+50, 512-50:512+50]
        
        
        psf = psf + np.random.normal(300, np.sqrt(300), np.shape(psf))
        img = img + np.random.normal(300, np.sqrt(300), np.shape(img))
        flux_measure, mux_measure, muy_measure, e1_measure, e2_measure, bkg_measure, size_img, sigxx_img, sigyy_img, sigxy_img = measure_pythonV_dist.measure(img, lut1, lut2, 100, 0 ,0, 5, 5, 0, 100, 0)
        flux_measure, mux_measure, muy_measure, e1_measure, e2_measure, bkg_measure, size_psf, sigxx_psf, sigyy_psf, sigxy_psf = measure_pythonV_dist.measure(psf, lut1, lut2, 100, 0 ,0, 5, 5, 0, 100, 0)
        k1,k2 = higer_test.measure(img , lut1, lut2, 100, 0 ,0, 5, 5, 0, 100, 0)
        
        e1 =  ((sigxx_img- sigxx_psf) -((sigyy_img- sigyy_psf)))/((sigxx_img- sigxx_psf) +((sigyy_img- sigyy_psf)))
        tempeArr.append(e1+ 0.3846)
        tempkArr.append(k2 - 2)
        
    eDiffArr.append(np.median(tempeArr))
    eErrArr.append(np.std(tempeArr))
    kArr.append(np.median(tempkArr))
    kErrArr.append(np.std(tempkArr))

a = np.zeros ((4, len(eDiffArr)))
a[0,:] = eDiffArr
a[1,:] = eErrArr
a[2,:] = kArr
a[3,:] = kErrArr

np.save('/home/dutta26/codes/sersic_test/data'+id_n+'.npy' , a)
#plt.errorbar(kArr, eDiffArr, yerr= eErrArr, xerr =kErrArr ,fmt = 'b.')
#plt.xlabel('Excess Kurtosis')
#plt.ylabel('Delta Elliptcity')
#plt.title('True Elliptcity =0.385')