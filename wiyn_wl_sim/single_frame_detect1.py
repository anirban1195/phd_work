#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 15:07:31 2023

@author: dutta26
"""


from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os,sys
import pandas as pd
import helper, helper1
from astropy.stats import sigma_clipped_stats
import subprocess

lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')
def run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc , plotloc):
   
    #Read the source catalog 
    
    source_df = pd.read_pickle(source_df_name)
    #print (band, str(sys.argv[2]))
    star_bool = np.array(source_df['star_bool'])
    raList = np.array(source_df['ra'])
    decList = np.array(source_df['dec'])
    
    #Read IR coadd catalog 
    band_coadd_df = np.load(band_coadd_df_name)
    band_coadd_xx = band_coadd_df[:,35]
    band_coadd_yy = band_coadd_df[:,36]
    band_coadd_xy = band_coadd_df[:,37]
    band_coadd_interpxx = band_coadd_df[:,38]
    band_coadd_interpyy = band_coadd_df[:,39]
    band_coadd_interpxy = band_coadd_df[:,40]
    band_coadd_flux = band_coadd_df[:,3]
    
    
    ir_coadd_df = np.load(ir_coadd_df_name)
    
    
    
    data_width = int(len(os.listdir(bandLoc +band+'/'))/2.0) +1
    data_width = 40
    store = np.zeros((data_width, len(raList), 77), dtype = np.float32)
    fileCount = -1
    fileList = os.listdir(bandLoc +band+'/')
    for file in fileList:
        
        if('star' in file):
            continue
        
        
       
        
        fileCount += 1
        print (file)
        #temp_arr= [10,11,13,16,19,68,70,104,105,106,122,123,124, 128]
        #if(fileCount not in temp_arr):
        #    continue
        fileArr = os.listdir(bandLoc +band+'/'+file)
        print (fileArr)
        for index in range(len(fileArr)):
            if('output' in fileArr[index]):
                continue
            else:
                suppFile= fileArr[index]
                
        #Read the swarp output
        #f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
        f=fits.open(bandLoc +band+'/'+file+'/output.fits')
        data = f[0].data
        f.close()   
        
        print (np.shape(data), '*************')     
        
        #xList, yList = helper.convertToXY(raList, decList, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
        
        xList, yList = helper.convertToXY(decList, raList, bandLoc +band+'/'+file+'/output.fits')
        
        #Read the file data
        f=fits.open(bandLoc +band+'/'+file+'/output.fits')
        back_h = float((f[0].header)['BKG'])
        seeing = float((f[0].header)['SIZE'])
        zp = float((f[0].header)['ZP'])
        mjd = float((f[0].header)['MJD-OBS'])
        expTime = 60
        flux_scale = float((f[0].header)['FLXSCALE']) 
        f.close()
        #Read supplementary file 
        f=fits.open(bandLoc +band+'/'+file+'/'+suppFile)
        airmass = float((f[0].header)['AIRMASS'])
        f.close()
        mphase=mAngle =focus =zp_n= skymag =depth =mRa =mDec= fwhm=0 
        #Store chip positions in store
        #store[fileCount,:, 47],store[fileCount,:, 48], store[fileCount,:, 49] = helper.getLoc(raList, decList, '/scratch/bell/dutta26/abell_2390/'+band +'/'+file)
        
        store[fileCount,:, 44] = float(file)
        
        ySize,xSize = np.shape(data)
        
        #Fist measure stars
        cnt = cnt1 = cnt2 =cnt3 =0
        for j in range(len(xList)):
            #print (j)
            v_flag = b_flag = force_flag= 0
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            
            #If not star skip
            if(star_bool[j] == 0):
                store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
                store[fileCount,j,0] = ra
                store[fileCount,j,1] = dec
                store[fileCount,j,2] = star_bool[j]
                store[fileCount,j,10] = x
                store[fileCount,j,11] = y
                continue
            
            cut = data[y-25: y+25, x-25: x+25]
            #print (x,y)
            
            
           
            #If star then check if the image can pass measure
            b_flag = helper1.detectBad(cut)
            
            
            
           
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
            if(flux == None or flux<= 0 or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
                store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
                store[fileCount,j,0] = ra
                store[fileCount,j,1] = dec
                store[fileCount,j,2] = star_bool[j]
                store[fileCount,j,10] = x
                store[fileCount,j,11] = y
                continue
            #print (flux)
            store[fileCount,j,0:31] = ra, dec, star_bool[j], flux, mux, muy,  bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
        
        
        import pathlib
        path = pathlib.Path(plotloc+'/'+band+'/'+file)
        path.mkdir(parents=True, exist_ok=True)
        
        #Find all good stars in the frame
        #Now tuse k sigma clip to find usable stars. Just do for sigxx
        star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & (store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) 
                                              & (store[fileCount,:,14] == 0) & (ir_coadd_df[:,79] == 0)))[0],  : ]
    
        mean,median, std = sigma_clipped_stats(star_arr[:, 7])
        mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
        mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
        median_size = np.sqrt(median+median1)
        
        print ('*************************')
        print (mean, median, std)
        
        if(mean > 25 or std>3 or mean<1):
            store[fileCount, :,38] = -99
            store[fileCount, :,39] = -99
            store[fileCount, :,40] = -99
            continue
        
        threshold = np.sqrt(4* np.pi*back_h* (median+ median1 )) *100
        loc = np.where((store[fileCount,:,2] == 1) & 
                                              (store[fileCount,:,7] >= mean-3*std) &
                                              (store[fileCount,:,7] <= mean+3*std) &
                                              (store[fileCount,:,8] >= mean1-3*std1) &
                                              (store[fileCount,:,8] <= mean1+3*std1) &
                                              (store[fileCount,:,9] >= mean2-3*std2) &
                                              (store[fileCount,:,9] <= mean2+3*std2) &
                                              (store[fileCount,:,12] == 99) & 
                                              (store[fileCount,:,13] == 0) & 
                                              (store[fileCount,:,14] == 0) & 
                                              (store[fileCount,:,3] > threshold)&
                                              (ir_coadd_df[:,79] == 0))[0]
    
        print (len(loc))
        if(len(loc) == 0):
            store[fileCount, :,38] = -99
            store[fileCount, :,39] = -99
            store[fileCount, :,40] = -99
            continue
        
        
        
        star_arr = store[fileCount, loc,:]
        q,r = np.shape(star_arr)
        star_temp = np.zeros(( q , 10)   , dtype = np.float32)
        star_temp[:,0] = star_arr[:, 7]
        star_temp[:,1] = star_arr[:, 8]
        star_temp[:,2] = star_arr[:, 9]
        star_temp[:,3] = star_arr[:, 10]
        star_temp[:,4] = star_arr[:, 11]
        star_temp[:,6] = (star_arr[:, 3])/band_coadd_flux[loc]
        star_temp[:,7] = (star_arr[:, 3])
        print (star_arr[:, 3]/band_coadd_flux[loc])
        print (np.shape(star_arr))
        #print (threshold, median)
        nStars = 10
        global_xx_err, global_yy_err, global_xy_err = helper.fitGaussian(star_arr[:, 7], star_arr[:, 8], star_arr[:, 9], star_arr[:, 3], back_h, plotloc+'/'+band+'/'+file+'/', band)
        print (global_xx_err, global_yy_err, global_xy_err, 'Global errors')
        
        
        #Second pass. This pass we deconvolve with coadd PSF and recovolve with frame PSF of 10 nearest stars
        for j in range(len(xList)):
            
# =============================================================================
#             if(j != 5067):
#                 continue
# =============================================================================
            
            
            v_flag = b_flag = force_flag= 0
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            
            temp = np.copy(star_temp)
            temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
            temp = temp[temp[:,5].argsort()]
            
            #Check if same star. Then delete the entry
            if(temp[0,5]<5):
                temp = np.delete(temp, (0), axis = 0)
            
            #Checking for nans to avoid code from crashing
            if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            
            #Conpute the average PSF values 
            avgSigxx_frame = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            avgSigyy_frame = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            avgSigxy_frame = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            
            store[fileCount, j,38] = avgSigxx_frame
            store[fileCount, j,39] = avgSigyy_frame
            store[fileCount, j,40] = avgSigxy_frame
            #Find the average Poisson noise for stars used for PSF
            area = 2*np.pi*np.sqrt(store[fileCount,j,38]*store[fileCount, j,39] - store[fileCount, j,40]**2)
            median_psf_flux = np.median(temp[0:nStars, 7])
            size_psf = np.sqrt(store[fileCount,j,38]+ store[fileCount,j,39])
            B = back_h
            
            
            #s_psf =np.sqrt( (size_psf**4/median_psf_flux + 4*size_psf**6*np.pi * B/(median_psf_flux**2)) )
            s_psf=helper.getPSFErr(temp[0:nStars, 0], temp[0:nStars, 1], temp[0:nStars, 2], temp[0:nStars, 7], 1/temp[0:nStars, 5], B)
            
            
            
            store[fileCount,j,41] =  np.sqrt(s_psf**2 + global_xx_err**2) #np.std(temp[0:nStars, 0])
            store[fileCount,j,42] =  np.sqrt(s_psf**2 + global_yy_err**2) #np.std(temp[0:nStars, 1])
            store[fileCount,j,43] =  np.sqrt((s_psf* 0.707)**2 + global_xy_err**2)  #np.std(temp[0:nStars, 2])
            
            store[fileCount,j,74] = s_psf
            store[fileCount,j,75] = s_psf
            store[fileCount,j,76] = s_psf
            
            #Find flux ratio of these 10 nearest avg 
            #ratio = np.nanmean(temp[0:nStars, 6])
            ratio = np.mean(np.ma.masked_invalid(temp[0:nStars, 6]))
            
            flux_expected = band_coadd_flux[j]*ratio
            if(np.isnan(flux_expected) or np.isinf(flux_expected) or flux_expected<=0):
                store[fileCount, j,67] = 1
                if(band_coadd_df[j,3]>0):
                    sys.exit()###########################
                continue
            store[fileCount,j,60] = flux_expected
            #Use guess measure only when actual convergence fails 
            if(store[fileCount,j,12] == 99 ):
                continue
            #Find guess shapes. Use monte carlo shapes if available
            
            if(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] > 0):
                guess_xx = ir_coadd_df[j,35]  - ir_coadd_df[j,38] + avgSigxx_frame
                guess_yy = ir_coadd_df[j,36]  - ir_coadd_df[j,39] + avgSigyy_frame
                guess_xy = ir_coadd_df[j,37]  - ir_coadd_df[j,40] + avgSigxy_frame
            elif(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] == 0):
                guess_xx = ir_coadd_df[j,7]  - ir_coadd_df[j,38] + avgSigxx_frame
                guess_yy = ir_coadd_df[j,8]  - ir_coadd_df[j,39] + avgSigyy_frame
                guess_xy = ir_coadd_df[j,9]  - ir_coadd_df[j,40] + avgSigxy_frame
            else:
                store[fileCount, j,62] = 1
                
            
            if(guess_xx<0 or guess_yy<0):
                store[fileCount, j,63] = 1
                continue
    
            #Make cutout
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            guessmux = xList[j] - x + 0.5  #CHECK THIS(Checked and changed -0.5 to +0.5 on Aprl 19 2023)
            guessmuy = yList[j] - y + 0.5
            if(ir_coadd_df[j,7] <= 0 or ir_coadd_df[j,8]<=0 or np.isnan(ir_coadd_df[j,7]) or np.isnan(ir_coadd_df[j,8])):
                store[fileCount, j,64] = 1
                continue
            #size = np.sqrt(ir_coadd_df[j,7] + ir_coadd_df[j,8])
            size = np.sqrt(guess_xx + guess_yy)
            if(size<4):
                size = 3.9
            if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
                continue
            
            cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
            
            #If star then check if the image can pass measure
            b_flag = helper1.detectBad(cut)
            store[fileCount,j,2] =star_bool[j]
            store[fileCount,j,13] = 0
            store[fileCount,j,14] = b_flag
           
            snr = flux_expected / np.sqrt(flux_expected + 4*3.14* (guess_xx+guess_yy)*back_h )
            forced_measure_do = 0
            if(snr > 15): #Try convergence
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
                if(flux == None or flux<= 0 or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
                    forced_measure_do = 1
                else:
                    store[fileCount,j,0:31] = ra, dec, star_bool[j], flux, mux, muy,  bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
                    #continue
            else:
                forced_measure_do = 1
    
    
            if(forced_measure_do ==1 ):
                #Measure cutout
                
                #if(( 2*guess_xx * 2*guess_yy - 4*guess_xy**2)<=0 or  ( 2*guess_xx * 2*guess_yy - 4*guess_xy**2)> 6096):
                #    cnt1 += 1
                #    print (guess_xx, guess_yy, guess_xy, j)
                #flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = measure_pythonV.measure_v2(cut, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1)
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,1)
        
                
                if(flux == None or e1== None or e2 == None):
                    store[fileCount,j,12] = -99
                    store[fileCount, j,65] = 1
                    continue
                
                    
                store[fileCount,j,31:38] = flux, mux, muy,  bkg, sigxx, sigyy, sigxy
                store[fileCount,j,12] = 1
            
            
        
        #Run MC ALWAYS
        for j in range(len(xList)):
# =============================================================================
#             if(j != 5067):
#                 continue
#             
# =============================================================================
            if(store[fileCount, j,12] == 1):
                bkg = store[fileCount, j,34]
                
                sigxx, sigyy, sigxy = store[fileCount, j,35], store[fileCount, j,36], store[fileCount, j,37]
                if(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] > 0):
                    guess_xx = ir_coadd_df[j,35]  - ir_coadd_df[j,38] + store[fileCount, j,38]
                    guess_yy = ir_coadd_df[j,36]  - ir_coadd_df[j,39] + store[fileCount, j,39]
                    guess_xy = ir_coadd_df[j,37]  - ir_coadd_df[j,40] + store[fileCount, j,40]
                elif(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] == 0):
                    guess_xx = ir_coadd_df[j,7]  - ir_coadd_df[j,38] + store[fileCount, j,38]
                    guess_yy = ir_coadd_df[j,8]  - ir_coadd_df[j,39] + store[fileCount, j,39]
                    guess_xy = ir_coadd_df[j,9]  - ir_coadd_df[j,40] + store[fileCount, j,40]
                else:
                    store[fileCount, j,62] = 1
            elif(store[fileCount, j,12] == 99):
                sigxx, sigyy, sigxy = store[fileCount, j,7], store[fileCount, j,8], store[fileCount, j,9]
                bkg = store[fileCount, j,6]
                guess_xx = store[fileCount, j,7]
                guess_yy = store[fileCount, j,8]
                guess_xy = store[fileCount, j,9]
            else:
                continue
            
            N = store[fileCount, j,60]
            if(N<=0):
                continue
            #print (guess_xx, guess_yy, guess_xy, '****', ir_coadd_df[j,7], ir_coadd_df[j,8], '****',band_coadd_df[j,7], band_coadd_df[j,8])
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            #print (x,y, "x and y coords")
            #area = 2* np.pi*np.sqrt(guess_xx*guess_yy - guess_xy**2)
            size = np.sqrt(guess_xx + guess_yy)
            
            if(size< 0 or np.isnan(size)): #Use psf area if area is nan HOPEFULLY good appx at least
                size = np.sqrt(store[fileCount,j,38]+store[fileCount, j,39])
            
            cut = data[y-int(5*size): y+int(5*size), x-int(5*size): x+int(5*size)]    
            
            
            if(np.isnan(bkg) or bkg == None or bkg<=0 ):
                B = np.nanmedian(cut)
                if(np.isnan(B) or B == None or B<=0):
                    B = back_h
            else:
                B= bkg
            correction_fact = 0
            sqrtba = np.sqrt(B)*3.14*size**2
            #print (B, size, N, j)
            fbysqrtba = np.log10(N/sqrtba)
            index = int((fbysqrtba + 2.8)/0.1)
            if(index<=0):
                index = 1
            if(index>= 48):
                index = 47
            correction_fact = np.interp(fbysqrtba,lut_forcedDist[index-1:index+2,0], lut_forcedDist[index-1:index+2,2] )
            if(correction_fact >1):
                correction_fact = 1
            #print (correction_fact, N, sqrtba, store[fileCount, j,12])
            if(store[fileCount, j,12] == 1):
                s = np.sqrt( size**4/N + 4*B*np.pi*size**6/N**2)*correction_fact*0.5  
            elif(store[fileCount, j,12] == 99):
                s = np.sqrt( size**4/N + 4*B*np.pi*size**6/N**2)
            s_xx = np.sqrt((s)**2 )
            s_yy = np.sqrt((s)**2 )
            s_xy = np.sqrt((0.707*s)**2 )
            #print (N,B,size, s, s_xx)
            store[fileCount,j,68] = s_xx
            store[fileCount,j,69] = s_yy
            store[fileCount,j,70] = s_xy
            
            store[fileCount,j,71] = s_xx
            store[fileCount,j,72] = s_yy
            store[fileCount,j,73] = s_xy
            
            corr_xx = corr_yy = corr_xy = 0
            
            corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                                 store[fileCount,j,38], store[fileCount,j,39], store[fileCount,j,40], 
                                                                 s_xx,s_yy,s_xy, 
                                                                 store[fileCount,j,41], 
                                                                 store[fileCount,j,42], 
                                                                 store[fileCount,j,43], 30000)
            
# =============================================================================
#             corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
#                                                                  store[fileCount,j,38], store[fileCount,j,39], store[fileCount,j,40], 
#                                                                  0,0,0, 
#                                                                  0, 
#                                                                  0, 
#                                                                  0,100)
# =============================================================================
            
            if(success< 50 or np.isnan(corr_xx) or corr_xx<0 or corr_xx == None):
                store[fileCount, j,66] = 1
                continue
            
            else:
                #print (j,store[j,3], corr_xx, corr_yy,corr_xy )
                #cnt += 1
                store[fileCount,j,51] = corr_xx +store[fileCount,j,38]
                store[fileCount,j,52] = corr_yy + store[fileCount,j,39]
                store[fileCount,j,53] = corr_xy +store[fileCount,j,40]
                
        #Make validation plots for every frame 
        
        
        loc = np.where((store[fileCount,:,3]>0) & (store[fileCount,:,3]<100000) )[0]
        
        flux = store[fileCount,loc,3]
        n, bins, patches = plt.hist(x=flux, bins='auto',histtype=u'step')
        plt.xlabel('Flux')
        plt.savefig(plotloc+'/'+band+'/'+file+'/flux_hist_'+band+'.png')
        plt.close()
        
        
        size= np.sqrt(store[fileCount,:,51] -store[fileCount,:,38] + store[fileCount,:,52] - store[fileCount,:,39])
        loc= np.where( (store[fileCount,:,3] >0))[0]
        size=size[loc]
        
        size = size[size>0]
        
        size = size[size<15]
        print (size)
        n, bins, patches = plt.hist(x=size, bins='auto',histtype=u'step')                       
        plt.xlabel('Corrected Size')
        plt.savefig(plotloc+'/'+band+'/'+file+'/corr_size_hist_'+band+'.png')
        plt.close()
        
        
        loc= np.where((store[fileCount,:,2] == 1) &(store[fileCount,:,3] >0))[0]
        size= np.sqrt(store[fileCount,loc,51] -store[fileCount,loc,38] + store[fileCount,loc,52]  -store[fileCount,loc,39])
        
        
        n, bins, patches = plt.hist(x=size, bins='auto', histtype=u'step')                       
        plt.xlabel('Star Corrected Size')
        plt.savefig(plotloc+'/'+band+'/'+file+'/star_corr_size_hist_'+band+'.png')
        plt.close()
        
        
        
        loc= np.where((store[fileCount,:,2] == 0) &(store[fileCount,:,3] >0))[0]
        size= np.sqrt(store[fileCount,loc,51] - store[fileCount,loc,38] + store[fileCount,loc,52] - store[fileCount,loc,39])
        
        size = size[size>0]
        
        size = size[size<10]
        n, bins, patches = plt.hist(x=size, bins='auto', histtype=u'step')                    
        plt.xlabel('Galaxy Corrected Size')
        plt.savefig(plotloc+'/'+band+'/'+file+'/galaxy_corr_size_hist_'+band+'.png')
        plt.close()
        
        
        
        xx = store[fileCount,:,51] -store[fileCount,:,38]
        yy = store[fileCount,:,52] -store[fileCount,:,39]
        xy = store[fileCount,:,53] -store[fileCount,:,40]
        e1 = (xx-yy)/(xx+yy)
        e2 = 2*xy/(xx+yy)
        ellip = np.sqrt(e1**2+e2**2)
        loc = np.where((store[fileCount,:,2] == 0) &(ellip>0) &(store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) & (store[fileCount,:,14] == 0) &(ir_coadd_df[:,79] == 0))[0]
        n, bins, patches = plt.hist(x=e1[loc], bins='auto',histtype=u'step', color='r', label='e1')
        n, bins, patches = plt.hist(x=e2[loc], bins='auto',histtype=u'step', color='b', label='e2')
        n, bins, patches = plt.hist(x=ellip[loc], bins='auto',histtype=u'step', color='k', label='Tot Ellip')
        plt.legend()
        plt.xlabel('Galaxy e1/e2/Tot ellip')
        plt.savefig(plotloc+'/'+band+'/'+file+'/galaxy_ellip_hist_'+band+'.png')
        plt.close()
        
        loc = np.where((store[fileCount,:,2] == 1) &(ellip>0) &(store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) & (store[fileCount,:,14] == 0) &(ir_coadd_df[:,79] == 0))[0]
        n, bins, patches = plt.hist(x=e1[loc], bins='auto',histtype=u'step', color='r', label='e1')
        n, bins, patches = plt.hist(x=e2[loc], bins='auto',histtype=u'step', color='b', label='e2')
        n, bins, patches = plt.hist(x=ellip[loc], bins='auto',histtype=u'step', color='k', label='Tot Ellip')
        plt.legend()
        plt.xlabel('Star e1/e2/Tot ellip')
        plt.savefig(plotloc+'/'+band+'/'+file+'/star_ellip_hist_'+band+'.png')
        plt.close()
        
        
        
        
        
        
        
    
    outFile = '/home/dutta26/codes/wiyn_wl_sim/' +str(band)+'_withoutMC.npy'
    #outFile = '/scratch/bell/dutta26/backup/' +str(band)+'_withMC.npy'
    #np.save(outFile, store)


band = 'i'
bandLoc= '/scratch/bell/dutta26/wiyn_sim/'
band_coadd_df_name = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_'+ str(band)+'.npy'
ir_coadd_df_name = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy'
outLoc= '/home/dutta26/codes/wiyn_wl_sim/sf_' + str(band)+'.npy'
source_df_name = '/home/dutta26/codes/wiyn_wl_sim/source_list.pk1'
plotLoc = '/scratch/bell/dutta26/wiyn_sim/plot/phosim/'
run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc, plotLoc )


band = 'r'
bandLoc= '/scratch/bell/dutta26/wiyn_sim/'
band_coadd_df_name = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_'+ str(band)+'.npy'
ir_coadd_df_name = '/home/dutta26/codes/wiyn_wl_sim/coaddSc_ir.npy'
outLoc= '/home/dutta26/codes/wiyn_wl_sim/sf_' + str(band)+'.npy'
source_df_name = '/home/dutta26/codes/wiyn_wl_sim/source_list.pk1'
plotLoc = '/scratch/bell/dutta26/wiyn_sim/plot/phosim/'
run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc, plotLoc )

# =============================================================================
# band = 'i'
# bandLoc= '/scratch/bell/dutta26/backup/'
# band_coadd_df_name = '/scratch/bell/dutta26/backup/coaddSc_'+ str(band)+'.npy'
# ir_coadd_df_name = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'
# outLoc= '/scratch/bell/dutta26/backup/sf_' + str(band)+'.npy'
# source_df_name = '/scratch/bell/dutta26/backup/source_list.pk1'
# plotLoc = '/scratch/bell/dutta26/backup/plot/phosim/'
# run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc, plotLoc )
# 
# 
# band = 'r'
# bandLoc= '/scratch/bell/dutta26/backup/'
# band_coadd_df_name = '/scratch/bell/dutta26/backup/coaddSc_'+ str(band)+'.npy'
# ir_coadd_df_name = '/scratch/bell/dutta26/backup/coaddSc_ir.npy'
# outLoc= '/scratch/bell/dutta26/backup/sf_' + str(band)+'.npy'
# source_df_name = '/scratch/bell/dutta26/backup/source_list.pk1'
# plotLoc = '/scratch/bell/dutta26/backup/plot/phosim/'
# run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc, plotLoc )
# 
# =============================================================================
