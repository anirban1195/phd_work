#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 14:52:57 2023

@author: dutta26
"""



from astropy.io import fits
import numpy as np

sizex = sizey = 80
sigxArr= np.array([ 2,3,4,5])
sigyArr= np.array([ 2,3,4,5])
fluxArr = np.hstack((np.arange(1, 10, 1), np.arange(10, 100, 10), np.arange(100, 1000, 100), 
                     np.arange(1000, 10000, 1000), np.arange(10000, 110000, 10000) ))
bkgArr= [50, 100, 300, 500, 1000, 2000]
# =============================================================================
# sizex = sizey = 80
# sigxArr= np.array([2,3,4])
# sigyArr= np.array([2,3,4])
# fluxArr = np.hstack((np.arange(1, 10, 1), np.arange(10, 100, 20), np.arange(100, 1000, 200)
#                      ))
# bkgArr= [100, 300]
# =============================================================================
import measure_pythonV_dist
print (len(fluxArr)*len(bkgArr)* len(sigxArr)*len(sigyArr))
finalArr = []
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
for sigx in sigxArr :
    for sigy in sigyArr:
        if(sigx>5 or sigy>5):
            sizex = sizey = 100
        else:
            sizex = sizey = 80

        print (sigy)
        val = sigx*sigy*0.15
        sigxyArr = np.linspace(-val, val, 5)
        for sigxy in sigxyArr:
            for flux in fluxArr:
                
                for bkg in bkgArr:
                    temp_fluxArr =[]
                    temp_sigxxArr =[]
                    temp_sigxyArr =[]
                    temp_muxArr =[]
                    count = 0
                    area = np.pi * np.sqrt(sigx**2 * sigy**2 - sigxy**2)
                    snr = flux/np.sqrt(4*area*bkg + flux)
                    if(snr<1 or snr>8):
                        continue
                    
                    
                    
                    finalImg = np.zeros((sizex, sizey))
                    for j in range(200):
                        randx = randy = randxy = 0
                        while((randx*randy) <= randxy**2):
                            randxy = sigxy+np.random.normal(0, abs(0.15*sigxy))
                            randx = (sigx**2)+ np.random.normal(0, 0.15*sigx**2)
                            randy = (sigy**2)+ np.random.normal(0, 0.15*sigy**2)
                         
                        guessx = guessy = guessxy = 0
                        while((guessx*guessy) <= abs(guessxy) ):
                            guessx = np.sqrt((sigx**2)+ np.random.normal(0.0*sigx**2, 0.15*sigx**2))
                            guessy = np.sqrt((sigy**2)+ np.random.normal(0.0*sigy**2, 0.15*sigy**2))
                            guessxy = sigxy + np.random.normal(0.0*sigxy, abs(0.15*sigxy))
                        
                            
                        muArr= [sizex/2.0-0.5+ np.random.normal(0, 1), sizey/2.0-0.5+ np.random.normal(0, 1)]
                        cov = [[randx,randxy], [randxy, randy]]
                        const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
                        if(const<1):
                            const = 1
                        x, y = np.random.multivariate_normal(muArr, cov, const).T
                        x = np.int32(np.round(x))
                        y = np.int32(np.round(y))
                        obj = np.zeros((sizex,sizey))
                        np.add.at(obj, (y,x), 1)
                        
                        
                        noise = np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
                        tot = np.add(obj,noise)
                        finalImg = finalImg + tot 
                        
                        flux_measure, mux_measure, muy_measure, e1_measure, e2, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure_pythonV_dist.measure(tot, lut1, lut2, flux+np.random.normal(0.0*flux, 0.0*flux), 0 ,0, guessx, guessy, guessxy, 1, 1)
                        if(flux_measure == None or np.isnan(flux_measure) or psf_measure == None or np.isnan(psf_measure)):
                            count += 1
                            continue
                        if(np.isinf(sigxx_measure) or np.isinf(sigyy_measure) or np.isinf(sigxy_measure) ):
                            count += 1
                            continue
                        if(sigxx_measure == None or np.isnan(sigxx_measure) or sigyy_measure == None or np.isnan(sigyy_measure)):
                            count += 1
                            continue
                        if(sigxx_measure<0 or sigxx_measure>70 or sigyy_measure<0 or sigyy_measure>70 ):
                            count += 1
                            continue
                        temp_fluxArr.append(flux_measure)
                        temp_sigxxArr.append(sigxx_measure)
                        temp_sigxyArr.append(sigxy_measure)
                        temp_muxArr.append(mux_measure)
                    percent_flux_err = (np.median(temp_fluxArr) - flux) / flux
                    percent_sigx_err = (np.median(temp_sigxxArr) - (sigx**2/2)) / (sigx**2/2)
                    percent_sigxy_err = (np.median(temp_sigxyArr) - (sigxy/2)) / (sigxy/2)
                    #percent_sigxy_err = 0
                    failRate = count / 2
                    area = np.pi * np.sqrt(sigx**2 * sigy**2 - sigxy**2)
                    e1 = (sigx**2 - sigy**2)/(sigx**2 + sigy**2)
                    e2 = 2*sigxy/ (sigx**2 + sigy**2)
                    snr = flux/np.sqrt(4*area*bkg + flux)
                    
                    finalImg = finalImg/200
                    flux_measure, mux_measure, muy_measure, e1_measure, e2_measure, bkg_measure, psf_measure, sigxx_measure,sigyy_measure, sigxy_measure = measure_pythonV_dist.measure(finalImg, lut1, lut2, flux, 0 ,0, 3, 3,0, 100, 0)
                    
                    percent_flux_err1 = (flux_measure - flux) / flux
                    percent_sigx_err1 = (sigxx_measure - (sigx**2/2)) / (sigx**2/2)
                    percent_sigxy_err1 = (sigxy_measure - (sigxy/2)) / (sigxy/2)
                    #percent_sigxy_err1 = 0
                    
                    
                    
                    finalArr.append([percent_flux_err, percent_sigx_err, percent_sigxy_err, percent_flux_err1, percent_sigx_err1,
                                     percent_sigxy_err1, failRate, flux, area, e1, e2, bkg, sigx, sigy, sigxy, np.median(temp_muxArr), mux_measure])
                    
                    
finalArr= np.array(finalArr)                    
np.save('/scratch/bell/dutta26/abell_2390/paper_plot_2ndlast_mux.npy', finalArr)    
                    
                    
                    
                    
                    