#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 08:01:10 2023

@author: dutta26
"""

from astropy.io import fits
import numpy as np
import helper
import pandas as pd
import getgamma 
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats



def getWt(ellip_err):
    
    if(ellip_err>0.25):
        return 0
    if(ellip_err<0.001):
        return 1000
    else:
        return (1/ellip_err)


def get_avg_obs(datafile, img_file, cent_ra, cent_dec, radius, zLens, ir_coadd_file):
   
    #Load the master file      
    master_frame = np.load(datafile)
    ir_coadd_data = np.load(ir_coadd_file)
    
    #ir_coadd_data = np.array(ir_coadd_data)
    distArr = np.sqrt( (cent_ra-ir_coadd_data[:,0])**2 + (cent_dec-ir_coadd_data[:,1])**2)
    #print (distArr[10000])
    size = np.sqrt(ir_coadd_data[:,7] + ir_coadd_data[:,8])
    loc = np.where( (distArr <= radius) & (master_frame[:,4]> 1) & (ir_coadd_data[:,2]==0) & (ir_coadd_data[:,3]>10)
                   & (ir_coadd_data[:,82]==0))[0]
    print (len(loc), '***************', np.median(master_frame[loc, 15]))
    med_err = np.median(master_frame[loc, 15])
    totShear =0
    ellipArr=[]
    wtArr=[]
    errArr=[]
    [x_mid], [y_mid] = helper.convertToXY(cent_dec, cent_ra, img_file)
    
    cnt =0
    arr=[]
    th_arr =[]
    z_arr=[]
    signalArr =[]
    #print (ir_coadd_data[0,1],ir_coadd_data[0,1])
    for ind in loc:
        dx = ir_coadd_data[ind,10] -x_mid
        dy = ir_coadd_data[ind,11] -y_mid
        #dx = ir_coadd_data[ind,0] -cent_dec
        #dy = ir_coadd_data[ind,1] - cent_ra   
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        e1 = master_frame[ind,2]
        #e1 = (ir_coadd_data[ind,7]-ir_coadd_data[ind,8])/(ir_coadd_data[ind,7]+ir_coadd_data[ind,8])
        e2 = master_frame[ind,3]
        #e2 = 2*ir_coadd_data[ind,9]/(ir_coadd_data[ind,7]+ir_coadd_data[ind,8])
        snr = ir_coadd_data[ind,3]*800/np.sqrt(4*np.pi*(master_frame[ind,6] + master_frame[ind,7])*18000)
        size = np.sqrt(master_frame[ind,6] + master_frame[ind,7])
        size_coadd =  np.sqrt(ir_coadd_data[ind,7]+ir_coadd_data[ind,8])
        #ellip_err = 1.414*np.sqrt( 1/(ir_coadd_data[ind,3]*800) + 4*np.pi*18000*size_coadd**2/(ir_coadd_data[ind,3]*800)**2  )
        ellip_err = master_frame[ind, 15]
        tot_ellip = np.sqrt(e1**2 + e2**2)
        
        
        wt=(1/ellip_err)**2
        if(wt>1/(med_err/3)**2):
            wt=1/(med_err/3)**2
        #print (wt)
        if(abs(e1)>5 or abs(e2)>5 or tot_ellip==0 or np.isnan(size) or size_coadd<0.2 or size_coadd>10 ):
            wt = 0
            
        
        epar= - (e1*cos2phi + e2*sin2phi)
        if(abs(e1)>5 or abs(e2)>5 or tot_ellip==0 or np.isnan(size) or size<0.2 or size>10 or wt<=0 ):
            continue
        else:
            #print (master_frame[ind,6] , ir_coadd_data[ind,7], ir_coadd_data[ind,38])
            #print (e1)
            
            totShear += e2*wt#(e1*wt)
            cnt += (1*wt)
            signalArr.append(e2)
            wtArr.append(wt)
            errArr.append(master_frame[ind, 15])
            ellipArr.append(tot_ellip)
            
            #gamma1, gamma2, kappa = getgamma.getGamma(328.3944, 17.6695, 0.3, ir_coadd_data[ind,0], ir_coadd_data[ind,1], master_frame[ind, 4], 0.25)
            #th_arr.append(- (gamma1*cos2phi + gamma2*sin2phi))
            #z_arr.append(master_frame[ind, 4])
    ellipArr =np.array(ellipArr)
    signalArr = np.array(signalArr)
    wtArr=np.array(wtArr)
    errArr=np.array(errArr)
    wtArr=wtArr/np.sum(wtArr)
    ellipArr =np.array(ellipArr)
    mean,med,mode = sigma_clipped_stats(ellipArr**2, cenfunc= np.median)
    print (len(ellipArr), np.median(ellipArr**2), cnt, totShear, np.sum(ellipArr**2 * wtArr**0.5)/np.sum(wtArr**0.5), np.sum(ellipArr* wtArr),np.sum(ellipArr**2 * wtArr)) 
    #n, bins, patches = plt.hist(x=ellipArr, bins='auto',histtype=u'step', density=True)
    #plt.xlabel('Galaxy Corrected Size')-np.median(ellipArr**2)
    #print (np.median(errArr))
    wt_tild = wtArr**0.5/np.sum(wtArr**0.5)
    shear = totShear /(2*cnt*(1-np.sum(ellipArr**2 * wt_tild)/np.sum(wt_tild) ))
    error = shear / np.sqrt(len(ellipArr))
    
    wtArr1 = wtArr /(2*cnt*(1-np.sum(ellipArr**2 * wt_tild)/np.sum(wt_tild) ))
    wtArr1 = wtArr1/np.sum(wtArr1)
    err1 = np.sqrt(np.sum(wtArr1**2 * errArr**2))
    #print( np.sqrt(1.414*np.sum(errArr**2))/len(errArr), np.sqrt(np.sum(wtArr**2 * errArr**2))/(2*(1-np.mean(ellipArr**2))))
    #print( np.sqrt(np.sum(2*errArr**2))/len(errArr), np.sqrt(np.sum(wtArr**2 * errArr**2))/(2*(1-np.mean(ellipArr**2))))
    #return shear, np.std(ellipArr)/np.sqrt(len(ellipArr))
    #print (np.median(th_arr))
    #print (np.mean(th_arr), np.mean(z_arr), np.median(z_arr), np.sum(z_arr*wtArr))
    print (np.sum(wtArr**2 * errArr**2)/np.sum(wtArr * signalArr)**2)
    print (2*np.sum(wt_tild**2* errArr**2)/(1-np.sum(ellipArr**2 * wt_tild)/np.sum(wt_tild) )**2)
    a = (np.sum(wtArr**2 * errArr**2)/np.sum(wtArr * signalArr)**2) +(2*np.sum(wt_tild**2* errArr**2)/(1-np.sum(ellipArr**2 * wt_tild)/np.sum(wt_tild) )**2)
    print (np.sqrt(a)*shear, shear, np.sqrt(a))
    #a = np.sqrt()
    #return shear, np.sqrt(np.sum(wtArr**2 * errArr**2) + (np.sqrt(2*np.sum(errArr**2))/len(errArr))**2) , np.median(th_arr)
    return shear, np.sqrt(a)*shear , np.median(th_arr), ellipArr, wtArr
    #return shear, err1 , np.median(th_arr), ellipArr, wtArr
    


def get_avg_th(cent_ra, cent_dec, radius, zLens):
    
    bins =np.arange(50/3600, radius, 1/3600)
    totGamma =0
    cnt =0
    for radius in bins:
        gamma1, gamma2, kappa = getgamma.getGamma(cent_ra, cent_dec, 0.3, cent_ra+radius, cent_dec, 0.7, 0.25)
        totGamma += (np.sqrt(gamma1**2 + gamma2**2) * radius * (1/3600))
        cnt += (radius * (1/3600))
    return totGamma/cnt


    
 
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
# =============================================================================
# f=fits.open('/scratch/bell/dutta26/wiyn_sim/EMode_sf.fits')
# data= f[0].data
# f.close()
# 
# cut = data[212-15:212+15, 212-15:212+15]
# flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
# cent_dec, cent_ra = helper.convertToRaDec([211.5+mux], [211.5+muy], '/scratch/bell/dutta26/wiyn_sim/EMode_sf.fits')
# print (mux, muy,cent_dec, cent_ra)
# =============================================================================

# =============================================================================
# cut = data[212-4:212+3, 212-4:212+3]
# num_x=0
# num_y=0
# for j in range(7):
#     for k in range(7):
#         num_x += cut[j,k]*(k-3)
#         num_y += cut[j,k]*(j-3)
#         
# mux, muy = num_x/np.sum(cut),num_y/np.sum(cut)
# 
# dec, ra = helper.convertToRaDec([212+mux], [212+muy], '/scratch/bell/dutta26/wiyn_sim/EMode_sf.fits')
# print (mux, muy,dec, ra)
# 
# =============================================================================
# =============================================================================
# finalDataSet_loc='/scratch/bell/dutta26/backup/master_arr_coaddMC_sfMC_z.npy'
# imgLoc = '/scratch/bell/dutta26/backup/wted_coadds/ir_coadd_wt.fits'
# ir_coadd_file = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'
# =============================================================================

finalDataSet_loc='/scratch/bell/dutta26/wiyn_sim1/master_arr_coaddMC_sfMC.npy'
imgLoc = '/scratch/bell/dutta26/wiyn_sim1/wted_coadds/ir_coadd_wt.fits'
ir_coadd_file = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy'


#Now go in radial bins and determine shear
bins = np.arange(60/3600, 600/3600, 30/3600 ) 
th_avg_shear = []
obs_avg_shear =[]
obs_shear_err=[]
for radius in bins:
    #th_avg_shear.append(get_avg_th(328.3941 , 17.6697, radius, 0.3))
    th_avg_shear.append(-0.05)
    shear_obs, error, th_med, ellipArr, wtArr = get_avg_obs(finalDataSet_loc, imgLoc, [328.3944], [16.6695], radius, 0.3, ir_coadd_file )
    obs_avg_shear.append(shear_obs)
    obs_shear_err.append(error)
    #break
    #th_avg_shear.append(-th_med)
    
#yerr = np.linspace(0.5, 1, 24).T
plt.plot(bins*3600, th_avg_shear, 'k--')
plt.errorbar(bins*3600, obs_avg_shear,  yerr= obs_shear_err, fmt='b.', label = 'Single Frame',capsize=3)

print (np.mean(obs_avg_shear[-10:]))




finalDataSet_loc='/scratch/bell/dutta26/wiyn_sim1/master_arr_coaddMC.npy'
imgLoc = '/scratch/bell/dutta26/wiyn_sim1/wted_coadds/ir_coadd_wt.fits'
ir_coadd_file = '/home/dutta26/codes/wiyn_wl_sim1/coaddSc_ir.npy'


#Now go in radial bins and determine shear
bins = np.arange(60/3600, 600/3600, 30/3600 ) 
th_avg_shear = []
obs_avg_shear =[]
obs_shear_err=[]
for radius in bins:
    #th_avg_shear.append(get_avg_th(328.3941 , 17.6697, radius, 0.3))
    th_avg_shear.append(-0.05)
    shear_obs, error, th_med, ellipArr, wtArr = get_avg_obs(finalDataSet_loc, imgLoc, [328.3944], [16.6695], radius, 0.3, ir_coadd_file )
    obs_avg_shear.append(shear_obs)
    obs_shear_err.append(error)
    #break
    #th_avg_shear.append(-th_med)
    
#yerr = np.linspace(0.5, 1, 24).T
plt.plot(bins*3600, th_avg_shear, 'k--', label = r'Input $\gamma_2$')
plt.errorbar(bins*3600, obs_avg_shear,  yerr= obs_shear_err, fmt='r.', label = 'Coadd', capsize=3)

print (np.mean(obs_avg_shear[-10:]))
plt.xlabel('Radius (in arcseconds)')
plt.ylabel(r'Average  $\gamma_2$')
plt.legend()
plt.savefig('/scratch/bell/dutta26/wiyn_sim1/comparison.png')
#plt.close()




