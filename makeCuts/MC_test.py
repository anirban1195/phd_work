#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:06:11 2023

@author: dutta26
"""

import numpy as np
import helper 
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
import sys
import matplotlib.pyplot as plt

coadd_e1=[]
sf_e1= []
true_e1_arr =[]
for samples in range(400):
    print (samples)
    sig_xx = np.random.uniform(3, 6) #4
    sig_yy = np.random.uniform(3, 6) #5
    sig_xy = 0
    true_e1 = (sig_xx-sig_yy)/(sig_xx+sig_yy)
    true_e1_arr.append(true_e1)
    
    coadd_xx = []
    coadd_yy = []
    coadd_xy = []
    
    
    sf_xx =[]
    sf_yy =[]
    sf_xy =[]
    
    
    bkg = 200
    flux = 250
    sizex =sizey = 100
    sigxx_psf = np.random.uniform(6, 8)#6
    sigyy_psf = np.random.uniform(6, 8)#6
    sigxy_psf  = 0
    delmu = 1
    delSig = 1
    for j in range(1):
        coadd = np.zeros((sizex, sizey))
        temp_xx =[]
        temp_yy= []
        temp_xy =[] 
        img_cube = np.zeros((200, 100, 100))
        for k in range(200):
            muRandX = np.random.normal(0, delmu)
            muRandY = np.random.normal(0, delmu)
            muArr= [sizex/2.0-0.5+ muRandX, sizey/2.0-0.5+muRandY]
            
            cov = [[(sig_xx+sigxx_psf)+np.random.normal(0,delSig), sig_xy], [sig_xy, (sig_yy+sigyy_psf)+np.random.normal(0,delSig)]]
            const = int(round(np.random.normal(flux, np.sqrt(flux)))) 
            x, y = np.random.multivariate_normal(muArr, cov, const).T
            x = np.int32(np.round(x))
            y = np.int32(np.round(y))
            obj = np.zeros((sizex,sizey))
            np.add.at(obj, (y,x), 1)
            obj = obj + np.random.normal(bkg, np.sqrt(bkg), (sizex,sizey))
            img_cube[k,:,:] = obj
            coadd = coadd + obj
            
            #sys.exit()
            
        flux_c, mux_c, muy_c, e1, e2, bkg_c, psf_c, sigxx_c, sigyy_c, sigxy_c = helper.measure_new(coadd, lut1, lut2)
        e1_guess = (sigxx_c - sigyy_c)/(sigxx_c + sigyy_c)
        e2_guess = 2*sigxy_c/(sigxx_c+sigyy_c)
        area = 2* np.pi*np.sqrt(sigxx_c*sigyy_c - sigxy_c**2)
        N =flux_c
        B=bkg_c
        s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * np.sqrt(sigxx_c+ sigyy_c) * 1.414
        s_xx = np.sqrt(((1+e1_guess)*s)**2 )
        s_yy = np.sqrt(((1-e1_guess)*s)**2 )
        s_xy = np.sqrt((0.5*(1+abs(e2_guess))*s)**2 )
        corr_xx_coadd, corr_yy_coadd, corr_xy_coadd , success = helper.correct(sigxx_c, sigyy_c, sigxy_c, 
                                                             sigxx_psf/2, sigyy_psf/2, sigxy_psf/2, 
                                                             s_xx,s_yy,s_xy, 
                                                             0,  0, 0)
        #print (s_xx, corr_xx_coadd, sigxx_c, sigxx_psf)
        coadd_xx.append(corr_xx_coadd)
        
        coadd_e1.append((corr_xx_coadd - corr_yy_coadd)/(corr_xx_coadd+corr_yy_coadd))
        
        #sys.exit()
        
        
        for k in range(200):
            guess_xx = corr_xx_coadd + sigxx_psf/2
            guess_yy = corr_yy_coadd + sigyy_psf/2
            guess_xy = corr_xy_coadd
            
            #guess_xx = (sig_xx + sigxx_psf)/2
            #guess_yy = (sig_yy + sigyy_psf)/2
            #guess_xy = 0
            guessmux = 0.5
            guessmuy = 0.5
            obj = img_cube[k,:,:]
            flux_c, mux_c, muy_c, e1, e2, bkg_c, psf_c, sigxx_c, sigyy_c, sigxy_c = helper.measure_new(obj, lut1, lut2, flux, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,1)
            
            #Run MC ALWAYS
            e1_guess = (guess_xx - guess_yy)/(guess_xx + guess_yy)
            e2_guess = 2*guess_xy/(guess_xx+guess_yy)
            area = 2* np.pi*np.sqrt(guess_xx*guess_yy - guess_xy**2)
            N =flux_c
            B=bkg_c
            s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * np.sqrt(guess_xx + guess_yy) * 1.414
            s_xx = np.sqrt(((1+e1_guess)*s)**2 )
            s_yy = np.sqrt(((1-e1_guess)*s)**2 )
            s_xy = np.sqrt((0.5*(1+abs(e2_guess))*s)**2 )
            
            
            corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx_c, sigyy_c, sigxy_c, 
                                                                 sigxx_psf/2, sigyy_psf/2, sigxy_psf/2, 
                                                                 0,0, 0, 
                                                                 0,  0, 0)
            temp_xx.append(corr_xx)
            temp_yy.append(corr_yy)
            temp_xy.append(corr_xy)
            
            #print (s_xx, corr_xx, sigxx_c, sigxx_psf)
        
        
        sf_xx.append(np.median(temp_xx))
        sf_e1.append(  (np.median(temp_xx) - np.median(temp_yy)) / ((np.median(temp_xx) + np.median(temp_yy))))
        
    #print (sf_xx, coadd_xx)
    #print (true_e1, coadd_e1, sf_e1)
    
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(sf_e1, true_e1_arr, 'b.')
ax1.plot(np.arange(-0.3, 0.3,0.01), np.arange(-0.3, 0.3,0.01), 'k--')
ax1.set(xlabel='Single Frames e1', ylabel='True e1')

ax2.plot(coadd_e1, true_e1_arr, 'r.')
ax2.plot(np.arange(-0.3, 0.3,0.01), np.arange(-0.3, 0.3,0.01), 'k--')
ax2.set(xlabel='Coadd Frames e1')



