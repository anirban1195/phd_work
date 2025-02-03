#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 19:21:19 2023

@author: dutta26
"""

#Import libraries 
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os,sys
import pandas as pd
import helper, helper1
from astropy.stats import sigma_clipped_stats
import subprocess


#Read lookup tables for forced measurement and adjustment factors 
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')



def run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc , plotloc):
   
    #Read the source catalog 
    source_df = pd.read_pickle(source_df_name)
    star_bool = np.array(source_df['star_bool']) #A boolean array containing 1 for pt sources else 0
    raList = np.array(source_df['ra'])
    decList = np.array(source_df['dec'])
    
    #Read filter coadd catalog 
    band_coadd_df = np.load(band_coadd_df_name)
    band_coadd_xx = band_coadd_df[:,35] 
    band_coadd_yy = band_coadd_df[:,36]
    band_coadd_xy = band_coadd_df[:,37]
    band_coadd_interpxx = band_coadd_df[:,38]
    band_coadd_interpyy = band_coadd_df[:,39]
    band_coadd_interpxy = band_coadd_df[:,40]
    band_coadd_flux = band_coadd_df[:,3] 
    
    #Read ir_coadd_catalog
    ir_coadd_df = np.load(ir_coadd_df_name)
    
    
    #Load the bad chip locations 
    #Format is (chip_serial_no, pos_x, pos_x)
    segDataLoc ='/home/dutta26/codes/makeWeights_A2261/segment_stat_new.txt'
    f=open(segDataLoc)
    content = f.readlines()
    f.close()
    segment_stat=[]
    for j in range(len(content)):
        temp = content[j].split(',')
        segment_stat.append([int(temp[0]), int(temp[1]), int(temp[2])])
    segment_stat = np.array(segment_stat)   
    #Custom defined bad pixel cells in the same format    
    badLoc = np.load('/home/dutta26/codes/stripeLoc.npy')

    
    #store variable is for storing all the single frame data. 3d array
    data_width = int(len(os.listdir(bandLoc +band+'/'))/2.0) + 5
    store = np.zeros((data_width, len(raList), 80), dtype = np.float32)
    fileCount = -1
    fileList = os.listdir(bandLoc +band+'/')
    
    #Iterate over all single frame images
    for file in fileList:
        
        #Reject stars, weights and temp files 
        if('star' in file or '.weight' in file or 'temp' in file):
            continue
        
        fileCount += 1

        
        print (file)
        
        #Run swarp
        f= open('/home/dutta26/temp.ascii', 'w+')
        f.write('/scratch/bell/dutta26/abell_2261/'+band +'/'+file)
        f.close()
        
        #COMMENT OUT this section IF FILES ALREADY PRESENT
        swarpLoc = '/home/dutta26/apps/bin/bin/'
        swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/abell_2261/'+band+'/temp/temp_'+str(file)+' -WEIGHTOUT_NAME /scratch/bell/dutta26/abell_2261/'+band+'/temp/temp_'+str(file)[0:-5]+'.weight.fits'
        process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
        output, error = process.communicate()
                
        #Read the swarp output
        f=fits.open('/scratch/bell/dutta26/abell_2261/'+band+'/temp/temp_'+str(file))
        data = f[0].data
        f.close() 
        
        print (np.shape(data), '*************')     
        
        #Convert the ra dec list to x,y coords to create cutouts
        xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2261/'+band+'/temp/temp_'+str(file))
        
        #Read the file specific data (bkg, mjd, zp, seeing etc) from headers.
        f=fits.open('/scratch/bell/dutta26/abell_2261/'+band +'/'+file)
        back_h = float((f[0].header)['SKYBG'])
        seeing = float((f[0].header)['SEEING'])
        zp = float((f[0].header)['MAGZERO'])
        fwhm = float((f[0].header)['FWHM_FLT'])
        mjd = float((f[0].header)['MJD-MID'])
        airmass = float((f[0].header)['AIRMASS'])
        #Fix for moon phase is waxing gibbous
        if(type((f[0].header)['MOONPHSE']) is str):
            mphase = -1.00000
        else:
            mphase = float((f[0].header)['MOONPHSE'])
        mAngle = float((f[0].header)['MOON_D'])
        expTime = float((f[0].header)['EXPTIME'])  
        focus = float((f[0].header)['TELFOCUS'])
        zp_n = float((f[0].header)['PHOTZP_N'])
        skymag = float((f[0].header)['SKYMAG']) 
        depth = float((f[0].header)['PHOTDPTH']) 
        mRa = float((f[0].header)['MOON_RA']) 
        mDec = float((f[0].header)['MOON_DEC']) 
        flux_scale = float((f[0].header)['FLXSCALE']) 
        f.close()
        
        #Store chip positions
        #getLoc return the location of the source in which (CCD, pixelCell_x, pixelCell_y) format
        store[fileCount,:, 47],store[fileCount,:, 48], store[fileCount,:, 49], bkg_variance_arr = helper.getLoc(raList, decList, '/scratch/bell/dutta26/abell_2261/'+band +'/'+file)
        store[fileCount,:, 44] = float(file[10:17])
        ySize,xSize = np.shape(data)
        
        
        
        #Fist measure stars
        cnt = cnt1 = cnt2 =cnt3 =0  #Use these variables to check for issues 
        for j in range(len(xList)):
            
            #Check if in excluded region
            ssindex = (store[fileCount,j, 47]- 1)*64 + store[fileCount,j, 48]*8 + store[fileCount,j, 49]
            ssindex = int(ssindex)
            if(segment_stat[ssindex,2] > 0.5*(segment_stat[ssindex,1]+segment_stat[ssindex,2]) or (segment_stat[ssindex,2]+segment_stat[ssindex,1])< 100 ):
                store[fileCount,j, 45] = 1
                
            #Check if in bad regions
            ccdNo = int(store[fileCount,j, 47])
            chip_x = int(store[fileCount,j, 48])
            chip_y = int(store[fileCount,j, 49])
            if(badLoc[ccdNo, chip_x, chip_y] == 0):
                store[fileCount,j, 46] = 1
                
            
            v_flag = b_flag = force_flag= 0 #Vflag = vertical stripe, bFlag = Bad region , force_flag = foced measurement
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            
            #If not star skip
            if(star_bool[j] != 1):
                store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
                store[fileCount,j,0] = ra
                store[fileCount,j,1] = dec
                store[fileCount,j,2] = star_bool[j]
                store[fileCount,j,10] = x
                store[fileCount,j,11] = y
                continue
            
            #A 80x80 cutout to measue star
            cut = data[y-40: y+40, x-40: x+40]
            
            
            #If star then check if has weird lines through center or nans
            if(star_bool[j] == 1 ):
                v_flag = helper1.vert_stripe(cut)
                
            #If star then check if the image can pass measure
            b_flag = helper1.detectBad(cut, 8)
            
            
            
            #Measure using AM alogrithm
            flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
            #Check if measurement failed i.e got nan
            if(flux == None or flux<= 0 or e1 == None or e2==None or  np.isnan(e1) or np.isnan(e2) or np.isnan(psf)):
                store[fileCount,j,15:31]=back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
                store[fileCount,j,0] = ra
                store[fileCount,j,1] = dec
                store[fileCount,j,2] = star_bool[j]
                store[fileCount,j,10] = x
                store[fileCount,j,11] = y
                continue
            
            #Save images to visually inspect stars 
# =============================================================================
#             hdu= fits.PrimaryHDU(cut)
#             hdu.writeto('/scratch/bell/dutta26/abell_2390/obj1/'+str(j)+'.fits')
# =============================================================================
            store[fileCount,j,0:31] = ra, dec, star_bool[j], flux, mux, muy,  bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
        
        
        #Make plot folder. If exists dont make
        import pathlib
        path = pathlib.Path(plotloc+'/'+band+'/'+file)
        path.mkdir(parents=True, exist_ok=True)
        
        #Find all good stars in the frame
        #Now use k sigma clip to find usable stars. Just do for sigxx
        #Conditions are >threshold flux, not in bad regions, must be identified as pt source and
        #not too close to other sources
        threshold = np.sqrt(4* np.pi*back_h* 5*5) *100
        star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & 
                                              (store[fileCount,:,12] == 99) & 
                                              (store[fileCount,:,13] == 0) & 
                                              (store[fileCount,:,14] == 0) &
                                              (store[fileCount,:,3] >threshold) & 
                                              (ir_coadd_df[:,79] == 0) ))[0],  : ]
        
        mean,median, std = sigma_clipped_stats(star_arr[:, 7])
        mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
        mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
        median_size = np.sqrt(median+median1)
        
        print (np.shape(star_arr))
        print ('*************************')
        print (mean, median, std, 'Mean, Med Std of xx')
        print (mean1, median1, std1, 'Mean, Med Std of yy')
        print (threshold)
        print (star_arr[:,7])
        print (star_arr[:,8])
        
        #The condition below usually indicates some sort of issue with the image
        #Reject those 
        if(mean > 30 or std>3 or mean<1 or mean1 > 30 or std1 >3 or mean1<1):
            store[fileCount, :,38] = -99
            store[fileCount, :,39] = -99
            store[fileCount, :,40] = -99
            continue
        #Recacluate the thereshold based on the median star size and S/N >= 100
        #Reject all stars that are 3sigma away from median second moment values 
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
        
        #If no stars are found in the previous condition,
        #Indicats weird issues with the image
        print (len(loc))
        if(len(loc) == 0):
            store[fileCount, :,38] = -99
            store[fileCount, :,39] = -99
            store[fileCount, :,40] = -99
            continue
        
        
        #Array star_arr and star_temp contain the stars selcted using the previous cuts.
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
        print (np.shape(star_arr), 'printing shape of star arr')
        
        nStars = 10 #Number of stars used for PSF interpolation 
        
        #Find the average variation of PSF across the frame 
        #Defined as half of 16th-84th percentile . These are global errors
        totScale= 1
        totBkg = back_h
        global_xx_err =  (np.percentile(star_arr[:, 7], 84) - np.percentile(star_arr[:, 7], 16))/2
        global_yy_err =  (np.percentile(star_arr[:, 8], 84) - np.percentile(star_arr[:, 8], 16))/2
        global_xy_err =  (np.percentile(star_arr[:, 9], 84) - np.percentile(star_arr[:, 9], 16))/2
        print (global_xx_err, global_yy_err, global_xy_err, 'Global errors')
        helper.plotStarGaussian(star_arr[:, 7], star_arr[:, 8], star_arr[:, 9], star_arr[:, 3]*totScale, totBkg, plotloc+'/'+band+'/'+file+'/', band, global_xx_err, global_yy_err, global_xy_err)
        
        
        #Second pass. This pass we deconvolve with coadd PSF and recovolve with frame PSF of 10 nearest stars
        #This is needed to porduce guess shape for forced measuremnent 
        
        for j in range(len(xList)):
            
            
            
            v_flag = b_flag = force_flag= 0
            
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round( yList[j]))
            
            #I used copy because sort was messing up the star_temp
            temp = np.copy(star_temp)
            temp[:,5] = np.sqrt((temp[:,3]-x)**2 + (temp[:,4]-y)**2)
            temp = temp[temp[:,5].argsort()]
            
            #Check if same star. Then delete the entry
            if(temp[0,5]<5):
                temp = np.delete(temp, (0), axis = 0)
            
            #Checking for nans to avoid code from crashing
            #If nan's found raise flag and skip this source 
            if(np.isnan(np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            if(np.isnan(np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            if(np.isnan(np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5]))):
                store[fileCount, j,61] = 1
                continue
            
            #Conpute the average PSF values . Inverse distance weighted
            avgSigxx_frame = np.sum( temp[0:nStars, 0] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            avgSigyy_frame = np.sum( temp[0:nStars, 1] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            avgSigxy_frame = np.sum( temp[0:nStars, 2] * 1/temp[0:nStars, 5])/np.sum(1/temp[0:nStars, 5])
            
            store[fileCount, j,38] = avgSigxx_frame
            store[fileCount, j,39] = avgSigyy_frame
            store[fileCount, j,40] = avgSigxy_frame
            
            
            
            #Find the average Poisson noise for stars used for PSF. helper.getPSFErr does this
            area = 2*np.pi*np.sqrt(store[fileCount,j,38]*store[fileCount, j,39] - store[fileCount, j,40]**2)
            median_psf_flux = np.median(temp[0:nStars, 7])
            size_psf = np.sqrt(store[fileCount,j,38]+ store[fileCount,j,39])
            B = back_h 
            #s_psf =np.sqrt( (size_psf**4/median_psf_flux + 4*size_psf**6*np.pi * B/(median_psf_flux**2)) )
            s_psf=helper.getPSFErr(temp[0:nStars, 0], temp[0:nStars, 1], temp[0:nStars, 2], temp[0:nStars, 7], 1/temp[0:nStars, 5], B)
            
            
            #Combine the Poisson error with the global error in quadrature. And store them in output array
            store[fileCount,j,41] =  np.sqrt(s_psf**2 + global_xx_err**2) #np.std(temp[0:nStars, 0])
            store[fileCount,j,42] =  np.sqrt(s_psf**2 + global_yy_err**2) #np.std(temp[0:nStars, 1])
            store[fileCount,j,43] =  np.sqrt((s_psf* 0.707)**2 + global_xy_err**2)  #np.std(temp[0:nStars, 2])
            
            store[fileCount,j,74] = s_psf
            store[fileCount,j,75] = s_psf
            store[fileCount,j,76] = s_psf
            
            #Find flux ratio of these 10 nearest avg. This is used to calculate guess flux for forced measuremnt
            #ratio = np.nanmean(temp[0:nStars, 6])
            ratio = np.mean(np.ma.masked_invalid(temp[0:nStars, 6]))
            
            flux_expected = band_coadd_flux[j]*ratio  #This is guess flux
            #If guess flux is nan, inf or negative, something is wrong with this source. Skip it
            if(np.isnan(flux_expected) or np.isinf(flux_expected) or flux_expected<=0):
                store[fileCount, j,67] = 1
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
                
            #If the guess shapes are nagtive, then skip the source
            if(guess_xx<0 or guess_yy<0):
                store[fileCount, j,63] = 1
                continue
    
            #Make cutout
            ra = raList[j]
            dec = decList[j]
            x = int(round( xList[j]))
            y = int(round( yList[j]))
            guessmux = xList[j] - x + 0.5  #CHECK THIS(Checked and changed -0.5 to +0.5 on Aprl 19 2023)
            guessmuy = yList[j] - y + 0.5
            store[fileCount,j,77] = guessmux
            store[fileCount,j,78] = guessmuy
            
            #If the second moments of the sources in i+r coadd is nagative or nan, then skip the source
            if(ir_coadd_df[j,7] <= 0 or ir_coadd_df[j,8]<=0 or np.isnan(ir_coadd_df[j,7]) or np.isnan(ir_coadd_df[j,8])):
                store[fileCount, j,64] = 1
                continue
            
            #Cutout size should be optimal. Hence the miinimum is 32x32 
            size = np.sqrt(guess_xx + guess_yy)
            if(size<4):
                size = 3.9
            #Check if cutout size outside the edge of the image
            if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
                continue
            
            cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
            
            #If star then check if has weird lines through center or nans
            if(star_bool[j] == 1 ):
                v_flag = helper1.vert_stripe(cut)
            #Check if the image can pass measure i.e if bad conditions are triggered. 
            b_flag = helper1.detectBad(cut, size)
            store[fileCount,j,2] =star_bool[j]
            store[fileCount,j,13] = v_flag
            store[fileCount,j,14] = b_flag
           
            
            #Caluclate the expected SNR based on guess flux and guess size. If SNR > 15 attempt normal measuremnt
            #Else if SNR<15 or convergence fails attempt forced measurement. 
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
    
            #Forced measurement section. Forced measurement done if S/N<15 or converge fails
            if(forced_measure_do ==1 ):
                # Force Measure cutout
                flux, mux, muy, e1, e2, bkg, psf, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2, flux_expected, guessmux, guessmuy, np.sqrt(2*guess_xx), np.sqrt(2*guess_yy), 2*guess_xy, 1,1)
                #If the measurement produces nan's then skip the sources. Indicative of issues with the source
                if(flux == None or e1== None or e2 == None or np.isnan(flux) or flux == 0):
                    store[fileCount,j,12] = -99
                    store[fileCount, j,65] = 1
                    continue
                
                    
                store[fileCount,j,31:38] = flux, mux, muy,  bkg, sigxx, sigyy, sigxy
                store[fileCount,j,12] = 1
            
            
        #Run MC PSF correction on all sources
        for j in range(len(xList)):
            #MC PSC correction depends on if the source has convered using AM algorithm or not store in store[fileCount, j,12]
            #Because error made by the algorithm is differnnt in two cases 
            if(store[fileCount, j,12] == 1):  #Forced measurement case
                bkg = store[fileCount, j,34]
                
                sigxx, sigyy, sigxy = store[fileCount, j,35], store[fileCount, j,36], store[fileCount, j,37]
                if(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] > 0):  #If uncorrected i+r coadd measuremtn was used
                    guess_xx = ir_coadd_df[j,35]  - ir_coadd_df[j,38] + store[fileCount, j,38]
                    guess_yy = ir_coadd_df[j,36]  - ir_coadd_df[j,39] + store[fileCount, j,39]
                    guess_xy = ir_coadd_df[j,37]  - ir_coadd_df[j,40] + store[fileCount, j,40]
                elif(ir_coadd_df[j,7] > 0 and ir_coadd_df[j,35] == 0):#If MC corrected i+r coadd measuremtn was used
                    guess_xx = ir_coadd_df[j,7]  - ir_coadd_df[j,38] + store[fileCount, j,38]
                    guess_yy = ir_coadd_df[j,8]  - ir_coadd_df[j,39] + store[fileCount, j,39]
                    guess_xy = ir_coadd_df[j,9]  - ir_coadd_df[j,40] + store[fileCount, j,40]
                else:
                    store[fileCount, j,62] = 1
                    
            elif(store[fileCount, j,12] == 99): #When AM converged 
                sigxx, sigyy, sigxy = store[fileCount, j,7], store[fileCount, j,8], store[fileCount, j,9]
                bkg = store[fileCount, j,6]
                guess_xx = store[fileCount, j,7]
                guess_yy = store[fileCount, j,8]
                guess_xy = store[fileCount, j,9]
            else:
                continue
            
            N = store[fileCount, j,60] #guess flux
            if(N<=0):
                continue
            #print (guess_xx, guess_yy, guess_xy, '****', ir_coadd_df[j,7], ir_coadd_df[j,8], '****',band_coadd_df[j,7], band_coadd_df[j,8])
            x = int(round( xList[j]))
            y = int(round(yList[j]))
            
            size = np.sqrt(guess_xx + guess_yy)
            
            if(size< 0 or np.isnan(size)): #Use psf area if area is nan HOPEFULLY good appx at least
                size = np.sqrt(store[fileCount,j,38]+store[fileCount, j,39])
            
            cut = data[y-int(5*size): y+int(5*size), x-int(5*size): x+int(5*size)]    
            #print (j, x, y, store[fileCount,j, 47],store[fileCount,j, 48], store[fileCount,j, 49])
            
            #Often the bkg variance is more important and bkg since correlated noise leads to higher variance 
            #than Poisson statistics. Henc sqrt(variance) is better than median bkg. This this is needed 
            #This calculated the bkg variance of the particulal pixel cell. 
            bkg_var = bkg_variance_arr[int(store[fileCount,j, 47]),int(store[fileCount,j, 48]), int(store[fileCount,j, 49])]
            
            if(np.isnan(bkg) or bkg == None or bkg<=0 ):
                B = np.nanmedian(cut)
                if(np.isnan(B) or B == None or B<=0):
                    B = back_h
            else:
                B= bkg
                
            if(np.isnan(bkg_var) or bkg_var == None or np.isinf(bkg_var) or bkg_var <=0 ):
               bkg_var = B**2
                
                
            correction_fact = 0
            sqrtba = np.sqrt(B)*3.14*size**2
            #print (B, size, N, j)
            #Calculate the lookup table index where the correction factor is found
            #Its stores as  a function of N/sqrt(BA). See Dutta 2024 for details. 
            fbysqrtba = np.log10(N/sqrtba)
            index = int((fbysqrtba + 2.8)/0.1)
            if(index<=0):
                index = 1
            if(index>= 48):
                index = 47
                
            # Calculate the correciton factor. p factor lookup tabel from Dutta et al 2024
            correction_fact = np.interp(fbysqrtba,lut_forcedDist[index-1:index+2,0], lut_forcedDist[index-1:index+2,2] )
            if(correction_fact >1):
                correction_fact = 1
            #print (correction_fact, N, sqrtba, store[fileCount, j,12])
            
            #This condition calculates the error in size made by the AM algorithm as s. It depends on if converged or not 
            if(store[fileCount, j,12] == 1):
                s = np.sqrt( size**4/N + 4*np.sqrt(bkg_var)*np.pi*size**6/N**2)*correction_fact*0.5  
            elif(store[fileCount, j,12] == 99):
                s = np.sqrt( size**4/N + 4*np.sqrt(bkg_var)*np.pi*size**6/N**2)
            s_xx = np.sqrt((s)**2 ) #error in xx
            s_yy = np.sqrt((s)**2 )  #error in yy
            s_xy = np.sqrt((0.707*s)**2 ) #error in xy
            #print (N,B,size, s, s_xx)
            
            #Not sure why I stored same variabe in two different places ??
            store[fileCount,j,68] = s_xx
            store[fileCount,j,69] = s_yy
            store[fileCount,j,70] = s_xy
            
            store[fileCount,j,71] = s_xx
            store[fileCount,j,72] = s_yy
            store[fileCount,j,73] = s_xy
            
            #Do the MC psf correction by calling helper.correct()
            #It needs the measured size , error in measured size, PSF and error in PSF
            corr_xx = corr_yy = corr_xy = 0
            corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                                 store[fileCount,j,38], store[fileCount,j,39], store[fileCount,j,40], 
                                                                 s_xx,s_yy,s_xy, 
                                                                 store[fileCount,j,41], 
                                                                 store[fileCount,j,42], 
                                                                 store[fileCount,j,43])
            
            
            #If nans are returned or only a few success then MC was not successful
            if(success< 50 or np.isnan(corr_xx) or corr_xx<0 or corr_xx == None):
                store[fileCount, j,66] = 1
                continue
            
            else: #If MC correction successful then store them in array
                #print (j,store[j,3], corr_xx, corr_yy,corr_xy )
                #cnt += 1
                store[fileCount,j,51] = corr_xx +store[fileCount,j,38]
                store[fileCount,j,52] = corr_yy + store[fileCount,j,39]
                store[fileCount,j,53] = corr_xy +store[fileCount,j,40]
                
        #Make validation/diagnostic plots for every frame 
        #Plot histogram of fluxes
        loc = np.where((store[fileCount,:,3]>0) & (store[fileCount,:,3]<100000) )[0]
        flux = store[fileCount,loc,3]
        n, bins, patches = plt.hist(x=flux, bins='auto',histtype=u'step')
        plt.xlabel('Flux')
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/flux_hist_'+band+'.png')
        plt.close()
        
        #Plot histogram of sizes. Cutoff at 15pix
        size= np.sqrt(store[fileCount,:,51] -store[fileCount,:,38] + store[fileCount,:,52] - store[fileCount,:,39])
        loc= np.where( (store[fileCount,:,3] >0))[0]
        size=size[loc] 
        size = size[size>0]
        size = size[size<15]
        print (size)
        n, bins, patches = plt.hist(x=size, bins='auto',histtype=u'step')                       
        plt.xlabel('Corrected Size')
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/corr_size_hist_'+band+'.png')
        plt.close()
        
        
        #Plot hisgotgram of star sizes after PSF correction
        loc= np.where((store[fileCount,:,2] == 1) &(store[fileCount,:,3] >0))[0]
        size= np.sqrt(store[fileCount,loc,51] -store[fileCount,loc,38] + store[fileCount,loc,52]  -store[fileCount,loc,39])
        n, bins, patches = plt.hist(x=size, bins='auto', histtype=u'step')                       
        plt.xlabel('Star Corrected Size')
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/star_corr_size_hist_'+band+'.png')
        plt.close()
        
        
        #Plot histogram of galaxy sizes after PSF correction. Cutoff at 10 pixels 
        loc= np.where((store[fileCount,:,2] == 0) &(store[fileCount,:,3] >0))[0]
        size= np.sqrt(store[fileCount,loc,51] - store[fileCount,loc,38] + store[fileCount,loc,52] - store[fileCount,loc,39])
        size = size[size>0]
        size = size[size<10]
        n, bins, patches = plt.hist(x=size, bins='auto', histtype=u'step')                    
        plt.xlabel('Galaxy Corrected Size')
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/galaxy_corr_size_hist_'+band+'.png')
        plt.close()
        
        
        #Plot ellipticty histogram of galaxies
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
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/galaxy_ellip_hist_'+band+'.png')
        plt.close()
        
        
        
        #Plot elliptcity histogram of stars. e1, e2 and totalellip
        xx = store[fileCount,:,51] 
        yy = store[fileCount,:,52] 
        xy = store[fileCount,:,53] 
        e1 = (xx-yy)/(xx+yy)
        e2 = 2*xy/(xx+yy)
        ellip = np.sqrt(e1**2+e2**2)
        loc = np.where((store[fileCount,:,2] == 1) &(ellip>0) &(store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) & (store[fileCount,:,14] == 0) &(ir_coadd_df[:,79] == 0))[0]
        n, bins, patches = plt.hist(x=e1[loc], bins='auto',histtype=u'step', color='r', label='e1')
        n, bins, patches = plt.hist(x=e2[loc], bins='auto',histtype=u'step', color='b', label='e2')
        n, bins, patches = plt.hist(x=ellip[loc], bins='auto',histtype=u'step', color='k', label='Tot Ellip')
        plt.legend()
        plt.xlabel('Star e1/e2/Tot ellip')
        plt.ylabel('Frequency')
        plt.savefig(plotloc+'/'+band+'/'+file+'/star_ellip_hist_'+band+'.png')
        plt.close()
        
        
        
        
        
        
        
    
    outFile = '/scratch/bell/dutta26/abell_2261/' +str(band)+'_sf.npy'
    #outFile = '/scratch/bell/dutta26/backup/' +str(band)+'_withMC.npy'
    np.save(outFile, store)


band = str(sys.argv[1])
bandLoc= '/scratch/bell/dutta26/abell_2261/'
band_coadd_df_name = '/scratch/bell/dutta26/abell_2261/abell_'+ str(band)+'_coadd.npy'
ir_coadd_df_name = '/scratch/bell/dutta26/abell_2261/abell_ir_coadd.npy'
outLoc= '/scratch/bell/dutta26/abell_2261/abell_'+band+'_coadd.npy'
source_df_name = str(sys.argv[2])
plotLoc = '/scratch/bell/dutta26/abell_2261/plot/'
star_arr=run_sf_detect(band, bandLoc, band_coadd_df_name, ir_coadd_df_name, source_df_name, outLoc, plotLoc )


