# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 17:14:22 2022

@author: dutta26
"""


from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import wcs
import numpy as np
import os,sys
import pandas as pd
import helper, helper1
import measure_pythonV
from astropy.stats import sigma_clipped_stats
import subprocess

lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
segDataLoc ='/home/dutta26/codes/makeWeights/segment_stat_new.txt'

f=open(segDataLoc)
content = f.readlines()
f.close()
segment_stat=[]
for j in range(len(content)):
    temp = content[j].split(',')
    segment_stat.append([int(temp[0]), int(temp[1]), int(temp[2])])
segment_stat = np.array(segment_stat)       



#band = 'r'
band = str(sys.argv[1])
#Read the source catalog 
#source_df = pd.read_pickle('/home/dutta26/codes/source_list.pk1')
source_df = pd.read_pickle(str(sys.argv[2]))
#print (band, str(sys.argv[2]))
star_bool = np.array(source_df['star_bool'])
raList = np.array(source_df['ra'])
decList = np.array(source_df['dec'])

#Read IR coadd catalog 
band_coadd_df = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'.pk1')
band_coadd_xx = np.array(band_coadd_df['xx'])
band_coadd_yy = np.array(band_coadd_df['yy'])
band_coadd_xy = np.array(band_coadd_df['xy'])
band_coadd_interpxx = np.array(band_coadd_df['interp_xx'])
band_coadd_interpyy = np.array(band_coadd_df['interp_yy'])
band_coadd_interpxy = np.array(band_coadd_df['interp_xy'])
band_coadd_flux = np.array(band_coadd_df['flux'])
band_coadd_df = np.array(band_coadd_df)


ir_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1'))


badLoc = np.load('/home/dutta26/codes/stripeLoc.npy')



data_width = int(len(os.listdir('/scratch/bell/dutta26/abell_2390/'+band +'/'))/2.0) +1
data_width = 150
store = np.zeros((data_width, len(raList), 77), dtype = np.float32)
fileCount = -1
fileList = os.listdir('/scratch/bell/dutta26/abell_2390/'+band +'/')
for file in fileList:
    
    if('.weight' in file or 'temp' in file):
        continue
    
    
   
    
    fileCount += 1
    print (file)
    #temp_arr= [10,11,13,16,19,68,70,104,105,106,122,123,124, 128]
    #if(fileCount not in temp_arr):
    #    continue
    #Run swarp to make a good image
    f= open('/home/dutta26/temp.ascii', 'w+')
    f.write('/scratch/bell/dutta26/abell_2390/'+band +'/'+file)
    f.close()
    
    #swarpLoc = '/home/dutta26/apps/bin/bin/'
    #swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits'
    #swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_coadd_'+str(fileCount)+'.fits  -WEIGHTOUT_NAME /scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_coadd_'+str(fileCount)+'.weight.fits' 
    #swarpCommand = './swarp @/home/dutta26/temp.ascii -c /home/dutta26/default_1.swarp -IMAGEOUT_NAME /scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file)+' -WEIGHTOUT_NAME /scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file)[0:-5]+'.weight.fits'

    #process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    #output, error = process.communicate()
    
    #Read the swarp output
    #f=fits.open('/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
    f=fits.open('/scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file))
    data = f[0].data
    f.close()   
    
    print (np.shape(data))     
    
    #xList, yList = helper.convertToXY(raList, decList, '/scratch/halstead/d/dutta26/abell_2390/temp_coadd.fits')
    xList, yList = helper.convertToXY(raList, decList, '/scratch/bell/dutta26/abell_2390/'+band+'/temp/temp_'+str(file))
    
    #Read the file data
    f=fits.open('/scratch/bell/dutta26/abell_2390/'+band +'/'+file)
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
    #Store chip positions in store
    store[fileCount,:, 47],store[fileCount,:, 48], store[fileCount,:, 49] = helper.getLoc(raList, decList, '/scratch/bell/dutta26/abell_2390/'+band +'/'+file)
    
    store[fileCount,:, 44] = float(file[10:17])
    
    ySize,xSize = np.shape(data)
    
    #Fist measure stars
    cnt = cnt1 = cnt2 =cnt3 =0
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
        
        
        
        #If star then check if has weird lines through center or nans
        if(star_bool[j] == 1 ):
            v_flag = helper1.vert_stripe(cut)
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
        
        store[fileCount,j,0:31] = ra, dec, star_bool[j], flux, mux, muy,  bkg, sigxx, sigyy, sigxy, x, y, 99, v_flag, b_flag, back_h ,seeing ,zp ,fwhm, mjd , airmass, mphase,mAngle ,expTime ,focus ,zp_n, skymag ,depth ,mRa ,mDec, flux_scale
    
    
    #Find all good stars in the frame
    #Now tuse k sigma clip to find usable stars. Just do for sigxx
    star_arr = store[fileCount, (np.where((store[fileCount,:,2] == 1) & (store[fileCount,:,12] == 99) & (store[fileCount,:,13] == 0) & (store[fileCount,:,14] == 0)
                                          & (store[fileCount,:,14] == 0)))[0],  : ]

    mean,median, std = sigma_clipped_stats(star_arr[:, 7])
    mean1,median1, std1 = sigma_clipped_stats(star_arr[:, 8])
    mean2,median2, std2 = sigma_clipped_stats(star_arr[:, 9])
    median_size = np.sqrt(median+median1)
    
    print ('*************************')
    print (mean, median, std)
    
    if(mean > 25 or std>3 or mean<2):
        store[fileCount, :,38] = -99
        store[fileCount, :,39] = -99
        store[fileCount, :,40] = -99
        store[fileCount, :,60] = 1
        continue
    
    threshold = np.sqrt(back_h* 4* 3.14 * 2* np.sqrt(median* median )) *100
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
                                          (store[fileCount,:,3] > threshold))[0]

    print (len(loc))
    if(len(loc) == 0):
        store[fileCount, :,38] = -99
        store[fileCount, :,39] = -99
        store[fileCount, :,40] = -99
        store[fileCount, :,60] = 1
        continue
    
    
    
    star_arr = store[fileCount, loc,:]
    q,r = np.shape(star_arr)
    star_temp = np.zeros(( q , 8)   , dtype = np.float32)
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
    
    #Second pass. This pass we deconvolve with coadd PSF and recovolve with frame PSF of 10 nearest stars
    for j in range(len(xList)):
        
        
        
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
        B = back_h
        s_psf = np.sqrt( (area/(np.pi*median_psf_flux) + 4*area**2 * B/(np.pi * median_psf_flux**2)) ) * np.sqrt(store[fileCount,j,38] + store[fileCount, j,39]) * 1.414
        e1_psf = (avgSigxx_frame - avgSigyy_frame)/(avgSigxx_frame + avgSigyy_frame)
        e2_psf = 2*avgSigxy_frame/(avgSigxx_frame + avgSigyy_frame)
        
        
        store[fileCount,j,41] = np.sqrt(((median_size/30)*median_size*1.414)**2 +  ((1+e1_psf)*s_psf)**2 ) #np.std(temp[0:nStars, 0])
        store[fileCount,j,42] = np.sqrt(((median_size/30)*median_size*1.414)**2 +  ((1-e1_psf)*s_psf)**2 )#np.std(temp[0:nStars, 1])
        store[fileCount,j,43] = np.sqrt(((median_size/30)*median_size*1.414 * 0.707)**2 +  (0.707*(1+abs(e2_psf))*s_psf)**2 )  #np.std(temp[0:nStars, 2])
        
        store[fileCount,j,74] = ((1+e1_psf)*s_psf)
        store[fileCount,j,75] = ((1-e1_psf)*s_psf)
        store[fileCount,j,76] = (0.707*(1+abs(e2_psf))*s_psf)
        
        #Find flux ratio of these 10 nearest avg 
        #ratio = np.nanmean(temp[0:nStars, 6])
        ratio = np.mean(np.ma.masked_invalid(temp[0:nStars, 6]))
        
        flux_expected = band_coadd_flux[j]*ratio
        if(np.isnan(flux_expected) or np.isinf(flux_expected) or flux_expected<=0):
            store[fileCount, j,67] = 1
            if(band_coadd_df[j,3]>0):
                sys.exit()###########################
            continue
        
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
        size = np.sqrt(ir_coadd_df[j,7] + ir_coadd_df[j,8])
        if(size<4):
            size = 3.9
        if(y-int(4*size) < 0 or x-int(4*size)<0 or y+int(4*size)> ySize or x+int(4*size)> xSize):
            continue
        
        cut = data[y-int(4*size): y+int(4*size), x-int(4*size): x+int(4*size)]
        
        #If star then check if has weird lines through center or nans
        if(star_bool[j] == 1 ):
            v_flag = helper1.vert_stripe(cut)
        #If star then check if the image can pass measure
        b_flag = helper1.detectBad(cut)
        store[fileCount,j,2] =star_bool[j]
        store[fileCount,j,13] = v_flag
        store[fileCount,j,14] = b_flag
       
        snr = flux_expected / np.sqrt(flux_expected + 4*3.14* np.sqrt(guess_xx*guess_yy)*back_h )
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
        e1_guess = (guess_xx - guess_yy)/(guess_xx + guess_yy)
        e2_guess = 2*guess_xy/(guess_xx+guess_yy)
        area = 2* np.pi*np.sqrt(guess_xx*guess_yy - guess_xy**2)
        if(area< 0 or np.isnan(area)): #Use psf area if area is nan HOPEFULLY good appx at least
            area = 2*np.pi*np.sqrt(store[fileCount,j,38]*store[fileCount, j,39] - store[fileCount, j,40]**2)
        N = flux_expected
        if(np.isnan(bkg) or bkg == None):
            B = np.nanmedian(cut)
            if(np.isnan(B) or B == None):
                B = back_h
        else:
            B= bkg
        s = np.sqrt( (area/(np.pi*N) + 4*area**2 * B/(np.pi * N**2)) ) * np.sqrt(guess_xx + guess_yy) * 1.414
        s_xx = np.sqrt(((1+e1_guess)*s)**2 + ((median_size/30)*median_size*1.414)**2)
        s_yy = np.sqrt(((1-e1_guess)*s)**2 + ((median_size/30)*median_size*1.414)**2)
        s_xy = np.sqrt((0.707*(1+abs(e2_guess))*s)**2 + ((median_size/30)*median_size*1.414*0.707)**2)
        
        store[fileCount,j,68] = s_xx
        store[fileCount,j,69] = s_yy
        store[fileCount,j,70] = s_xy
        
        store[fileCount,j,71] = ((1+e1_guess)*s)
        store[fileCount,j,72] = ((1-e1_guess)*s)
        store[fileCount,j,73] = (0.707*(1+abs(e2_guess))*s)
        
        corr_xx = corr_yy = corr_xy = 0
        corr_xx, corr_yy, corr_xy , success = helper.correct(sigxx, sigyy, sigxy, 
                                                             store[fileCount,j,38], store[fileCount,j,39], store[fileCount,j,40], 
                                                             s_xx,s_yy,s_xy, 
                                                             store[fileCount,j,41], 
                                                             store[fileCount,j,42], 
                                                             store[fileCount,j,43])
        if(success< 50 or np.isnan(corr_xx) or corr_xx<0 or corr_xx == None):
            store[fileCount, j,66] = 1
            continue
        
        else:
            #print (j,store[j,3], corr_xx, corr_yy,corr_xy )
            #cnt += 1
            store[fileCount,j,51] = corr_xx +store[fileCount,j,38]
            store[fileCount,j,52] = corr_yy + store[fileCount,j,39]
            store[fileCount,j,53] = corr_xy +store[fileCount,j,40]
    #sys.exit()
                
# =============================================================================
#         if(j == 21745):
#             sys.exit()
# =============================================================================
    
        
# =============================================================================
# df_source = pd.DataFrame(store,  ['ra', 'dec', 'star_flag','flux', 'mux', 'muy', 'e1', 'e2', 'bkg', 'size', 
#                                                               'xx', 'yy', 'xy','x', 'y','force_flag', 'vert_flag', 'bad_flag', 'back_sky', 'seeing', 'zp', 'fwhm', 'mjd' , 'airmass', 'mphase', 'mAngle' ,'expTime' ,'focus' ,'zp_n,skymag' ,'depth' ,'mRa' ,'mDec',
#                                                               'Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty','Empty', ])
# 
#     
# df_source.to_pickle('/home/dutta26/codes/frameSc_' + str(band)+'_.pk1')
# 
# =============================================================================
#np.save('/home/dutta26/codes/singleFrame_' +str(band)+'.npy', store)
np.save('/scratch/bell/dutta26/abell_2390/test2_' +str(band)+'.npy', store)
    