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
    print (len(loc), '***************')
    totShear =0
    ellipArr=[]
    wtArr=[]
    errArr=[]
    [x_mid], [y_mid] = helper.convertToXY(cent_dec, cent_ra, img_file)
    print (x_mid, y_mid)
    cnt =0
    arr=[]
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
        
       
        wt=(size/ellip_err)
        if(wt>100):
            wt=100
        #print (wt)
        
        epar= - (e1*cos2phi + e2*sin2phi)
        if(abs(e1)>1 or abs(e2)>1 or tot_ellip==0 or np.isnan(size) or size<0.05 or size>8 or wt<=0):
            continue
        else:
            #print (master_frame[ind,6] , ir_coadd_data[ind,7], ir_coadd_data[ind,38])
            #print (e1)
            
            totShear += epar*wt#(e1*wt)
            cnt += (1*wt)
            wtArr.append(wt)
            errArr.append(master_frame[ind, 15])
            ellipArr.append(tot_ellip)
    print (len(ellipArr), np.mean(ellipArr), cnt, totShear) 
    wtArr=np.array(wtArr)
    errArr=np.array(errArr)
    wtArr=wtArr/np.sum(wtArr)
    ellipArr =np.array(ellipArr)
    #n, bins, patches = plt.hist(x=ellipArr, bins='auto',histtype=u'step', density=True)
    #plt.xlabel('Galaxy Corrected Size')
    print (np.median(errArr))
    shear = totShear /(2*cnt*(1-np.mean(ellipArr**2)))
    error = shear / np.sqrt(len(ellipArr))
    print( np.sqrt(1.414*np.sum(errArr**2))/len(errArr), np.sqrt(np.sum(wtArr**2 * errArr**2))/(2*(1-np.mean(ellipArr**2))))
    print( np.sqrt(np.sum(2*errArr**2))/len(errArr), np.sqrt(np.sum(wtArr**2 * errArr**2))/(2*(1-np.mean(ellipArr**2))))
    #return shear, np.std(ellipArr)/np.sqrt(len(ellipArr))
    return shear, np.sqrt(np.sum(wtArr**2 * errArr**2) + (np.sqrt(2*np.sum(errArr**2))/len(errArr))**2)
    
def get_avg_th(cent_ra, cent_dec, radius, zLens):
    
    bins =np.arange(50/3600, radius, 1/3600)
    totGamma =0
    cnt =0
    for radius in bins:
        gamma1, gamma2, kappa = getgamma.getGamma(cent_ra, cent_dec, 0.3, cent_ra+radius, cent_dec, 0.8, 0.25)
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
finalDataSet_loc='/scratch/bell/dutta26/backup/master_arr_coaddMC_z.npy'
imgLoc = '/scratch/bell/dutta26/backup/wted_coadds/ir_coadd_wt.fits'
ir_coadd_file = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'

# =============================================================================
# finalDataSet_loc='/scratch/bell/dutta26/wiyn_sim/master_arr_coaddNoMC_sfNoMC.npy'
# imgLoc = '/scratch/bell/dutta26/wiyn_sim/wted_coadds/ir_coadd_wt.fits'
# ir_coadd_file = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy'
# =============================================================================


#Now go in radial bins and determine shear
bins = np.arange(60/3600, 600/3600, 30/3600 ) 
th_avg_shear = []
obs_avg_shear =[]
obs_shear_err=[]
for radius in bins:
    th_avg_shear.append(get_avg_th(328.3941 , 17.6697, radius, 0.3))
    #th_avg_shear.append(-0.1)
    shear_obs, error = get_avg_obs(finalDataSet_loc, imgLoc, [328.3941], [17.6697], radius, 0.3, ir_coadd_file )
    obs_avg_shear.append(shear_obs)
    obs_shear_err.append(error)
    
    
#yerr = np.linspace(0.5, 1, 24).T
plt.plot(bins*3600, th_avg_shear, 'r--')
plt.errorbar(bins*3600, obs_avg_shear,  yerr= obs_shear_err, fmt='r.', label = 'Single Frame No MC')


plt.xlabel('In arcesconds.')
plt.ylabel('Average shear')







