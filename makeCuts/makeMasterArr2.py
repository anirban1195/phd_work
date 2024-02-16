#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 11:50:50 2023

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
import wquantiles
from astropy import wcs


#f1 = fits.open('/scratch/bell/dutta26/abell_2390/temp.fits', mode = 'update')

def weighted_quantiles_interpolate(values, weights, quantiles=0.5):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, np.array(quantiles) * c[-1])]]

badLocs = np.load('/home/dutta26/codes/stripeLoc.npy')

#Read Redshifts
f=open('/home/dutta26/zphot.out')
content = f.readlines()
f.close()
redShiftArr=[]
for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        redShiftArr.append(float((content[j].split())[1]))
redShiftArr = np.ones(54141)*9      
        
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir.pk1')
ir_coadd_data = np.array(ir_coadd_data)
redShiftArr = np.array(redShiftArr)        
z_min = 0.3
z_max = 10
bandList =['g', 'r', 'i']

#g_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_g.npy')
r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_i.npy')
#z_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_z.npy')

g_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_g.pk1'))
r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r.pk1'))
i_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_i.pk1'))
#z_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1'))

bkgIndices = np.where( (redShiftArr>z_min )& (redShiftArr<z_max)& (ir_coadd_data[:,3] > 0.1) & (ir_coadd_data[:,3] < 100000)
                      & (ir_coadd_data[:,2] == 0))[0]

temp_coadd = np.zeros((40,40))

a1,b1 = np.shape(ir_coadd_data)
master_frame = np.zeros((a1,20), dtype = np.float32)
cnt1 = cnt2=cnt3= cnt4=cnt5=cnt6=cnt7=cnt8=cnt9=0



#Find total background 
totBkg = 0
totScale = 0
for value in r_sf_df[:,0,15]:
    if(value > 0):
        totBkg += value
for value in i_sf_df[:,0,15]:
    if(value > 0):
        totBkg += value

#Find eff scale factor
for value in r_sf_df[:,0,30]:
    if(value > 0):
        totScale += 1/value
for value in i_sf_df[:,0,30]:
    if(value > 0):
        totScale += 1/value

failedArr=[]
print ('aa')


for j in range(a1):
#for j in [169]:   
    print (j)
    
    if(j not in bkgIndices or ir_coadd_data [j,3] == 0):
        continue
    
    sigxx_arr=[]
    sigyy_arr=[]
    sigxy_arr=[]
    
    sigxx_err_arr=[]
    sigyy_err_arr=[]
    sigxy_err_arr=[]
    
    psfsigxx_arr=[]
    psfsigyy_arr=[]
    psfsigxy_arr=[]
    
    psfsigxx_uncertain_arr=[]
    psfsigyy_uncertain_arr=[]
    psfsigxy_uncertain_arr=[]
    
    size_sf_arr =[]
    
    wtArr= []
    area_arr=[]
    B_arr=[]
    N_arr=[]
    indexArr=[]
    coadd_measurement_flag = 0.0
    
    #If super faint use coadd measurements
    if(ir_coadd_data [j,3] < 0.25):
        
        sigxx_arr.append(ir_coadd_data [j,35])
        sigyy_arr.append(ir_coadd_data [j,36])
        sigxy_arr.append(ir_coadd_data [j,37])
        
        psfsigxx_arr.append(ir_coadd_data [j,38])
        psfsigyy_arr.append(ir_coadd_data [j,39])
        psfsigxy_arr.append(ir_coadd_data [j,40])
        
        psfsigxx_uncertain_arr.append(ir_coadd_data [j,41])
        psfsigyy_uncertain_arr.append(ir_coadd_data [j,42])
        psfsigxy_uncertain_arr.append(ir_coadd_data [j,43])
        
        area = 2*np.pi*np.sqrt(ir_coadd_data [j,7]* ir_coadd_data [j,8] - ir_coadd_data [j,9]**2)
        B = totBkg
        N = ir_coadd_data [j,3] *totScale
        e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
        e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
        e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
        error = np.sqrt((1/N + (4*area*B)/N**2)*( 1/(1-e_temp**2)) )
        if(error <= 0 or error== None or np.isnan(error)):
            error = 100000000
        wtArr.append(1/error)
        
        area_arr.append(area)
        N_arr.append(N)
        B_arr.append(B)
        indexArr.append(j)
        coadd_measurement_flag = 1
        
        sigxx_err_arr.append( np.sqrt(ir_coadd_data [j,71]**2 + ir_coadd_data [j,74]**2 ))
        sigyy_err_arr.append( np.sqrt(ir_coadd_data [j,72]**2 + ir_coadd_data [j,75]**2 ))
        sigxy_err_arr.append( np.sqrt(ir_coadd_data [j,73]**2 + ir_coadd_data [j,76]**2 ))
        
        size_sf_arr.append(np.sqrt(ir_coadd_data [j,38] + ir_coadd_data [j,39]))
        
    else:
        
        #First do r band 
        a,b,c = np.shape(r_sf_df)
        for k in range(a):
            
            
            if(r_sf_df[k,j,12] <= 0 or r_sf_df[k,j,13]> 0 or r_sf_df[k,j,14]> 0 or r_sf_df[k,j,38]==-99):
                continue
            if(np.sum(r_sf_df[k,j,60:67])>=1):
                continue
            if(r_sf_df[k,j,12] == 99):
                if(r_sf_df[k,j,7] < 0 or r_sf_df[k,j,7]>70 or r_sf_df[k,j,8] < 0 or r_sf_df[k,j,8]>70):
                    continue
                if(r_sf_df[k,j,7] == None or np.isnan(r_sf_df[k,j,7]) or r_sf_df[k,j,8] == None or np.isnan(r_sf_df[k,j,8])):
                    continue
                if(r_sf_df[k,j,3]<= 0 or r_sf_df[k,j,3]== None or np.isnan(r_sf_df[k,j,3]) ):
                    continue
                    
                
            if(r_sf_df[k,j,12] == 1):
                if(r_sf_df[k,j,35] < 0 or r_sf_df[k,j,35]>70 or r_sf_df[k,j,36] < 0 or r_sf_df[k,j,36]>70):
                    continue
                if(r_sf_df[k,j,35] == None or np.isnan(r_sf_df[k,j,35]) or r_sf_df[k,j,36] == None or np.isnan(r_sf_df[k,j,36])):
                    continue
                if(r_sf_df[k,j,31]<= 0 or r_sf_df[k,j,31]== None or np.isnan(r_sf_df[k,j,31]) ):
                    continue
            
            
            if(r_sf_df[k,j,12] == 99):
                sigxx_arr.append(r_sf_df[k,j,7])
                sigyy_arr.append(r_sf_df[k,j,8])
                sigxy_arr.append(r_sf_df[k,j,9])
                area = 2*np.pi*np.sqrt(r_sf_df[k,j,7]* r_sf_df[k,j,8] - r_sf_df[k,j,9]**2)
                #N =r_sf_df[k,j,3]
                N = r_coadd_df [j,3] / r_sf_df[k,j,30]
                B = r_sf_df[k,j,6]
                area_arr.append(area)
                if(N< 0 or N>1e7 or area > 600):
                    N = 0.00001
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            elif(r_sf_df[k,j,12] == 1):
                sigxx_arr.append(r_sf_df[k,j,35])
                sigyy_arr.append(r_sf_df[k,j,36])
                sigxy_arr.append(r_sf_df[k,j,37])
                temp_xx = r_sf_df[k,j,35]
                temp_yy = r_sf_df[k,j,36]
                temp_xy = r_sf_df[k,j,37]
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                N =r_sf_df[k,j,31]
                if(N< 0 or N>1e7 or area > 600):
                    N = 0.00001
                B = r_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            
            psfsigxx_arr.append(r_sf_df[k,j,38])
            psfsigyy_arr.append(r_sf_df[k,j,39])
            psfsigxy_arr.append(r_sf_df[k,j,40])
            
            psfsigxx_uncertain_arr.append(r_sf_df[k,j,41])
            psfsigyy_uncertain_arr.append(r_sf_df[k,j,42])
            psfsigxy_uncertain_arr.append(r_sf_df[k,j,43])
            
            sigxx_err_arr.append( np.sqrt(r_sf_df[k,j,71]**2 + r_sf_df[k,j,74]**2 ))
            sigyy_err_arr.append( np.sqrt(r_sf_df[k,j,72]**2 + r_sf_df[k,j,75]**2 ))
            sigxy_err_arr.append( np.sqrt(r_sf_df[k,j,73]**2 + r_sf_df[k,j,76]**2 ))
            size_sf_arr.append(np.sqrt(r_sf_df[k,j,38] + r_sf_df[k,j,39]))
            
            
            e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
            error = np.sqrt((1/N + (4*area*B)/N**2) * (1/ (1-e_temp**2)) )
            if(error == 0 or error == None or np.isnan(error)):
                error = 10000000000
            if(badLocs[int(r_sf_df[k,j,47]), int(r_sf_df[k,j,48]), int(r_sf_df[k,j,49])] == 0):
                wtArr.append(0)
                #print ('aa')
            else:
                
                wtArr.append(1/error)
                
         
            
        #Now do I band
        a,b,c = np.shape(i_sf_df)
        for k in range(a):
            
            
            if(i_sf_df[k,j,12] <= 0 or i_sf_df[k,j,13]> 0 or i_sf_df[k,j,14]> 0) or i_sf_df[k,j,38]==-99:
                continue
            if(np.sum(i_sf_df[k,j,60:67])>=1):
                continue
            if(i_sf_df[k,j,12] == 99):
                if(i_sf_df[k,j,7] < 0 or i_sf_df[k,j,7]>70 or i_sf_df[k,j,8] < 0 or i_sf_df[k,j,8]>70):
                    continue
                if(i_sf_df[k,j,7] == None or np.isnan(i_sf_df[k,j,7]) or i_sf_df[k,j,8] == None or np.isnan(i_sf_df[k,j,8])):
                    continue
                if(i_sf_df[k,j,3]<= 0 or i_sf_df[k,j,3]== None or np.isnan(i_sf_df[k,j,3]) ):
                    continue
                    
                
            if(i_sf_df[k,j,12] == 1):
                if(i_sf_df[k,j,35] < 0 or i_sf_df[k,j,35]>70 or i_sf_df[k,j,36] < 0 or i_sf_df[k,j,36]>70):
                    continue
                if(i_sf_df[k,j,35] == None or np.isnan(i_sf_df[k,j,35]) or i_sf_df[k,j,36] == None or np.isnan(i_sf_df[k,j,36])):
                    continue
                if(i_sf_df[k,j,31]<= 0 or i_sf_df[k,j,31]== None or np.isnan(i_sf_df[k,j,31]) ):
                    continue
            
            
            if(i_sf_df[k,j,12] == 99):
                sigxx_arr.append(i_sf_df[k,j,7])
                sigyy_arr.append(i_sf_df[k,j,8])
                sigxy_arr.append(i_sf_df[k,j,9])
                area = 2*np.pi*np.sqrt(i_sf_df[k,j,7]* i_sf_df[k,j,8] - i_sf_df[k,j,9]**2)
                #N =i_sf_df[k,j,3]
                N = i_coadd_df [j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,6]
                area_arr.append(area)
                if(N< 0 or N>1e7 or area > 600):
                    N = 0.00001
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
                
            elif(i_sf_df[k,j,12] == 1):
                sigxx_arr.append(i_sf_df[k,j,35])
                sigyy_arr.append(i_sf_df[k,j,36])
                sigxy_arr.append(i_sf_df[k,j,37])
                temp_xx = i_sf_df[k,j,35] 
                temp_yy = i_sf_df[k,j,36] 
                temp_xy = i_sf_df[k,j,37] 
                area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
                N = i_sf_df[k,j,31]
                if(N< 0 or N>1e7 or area > 600):
                    N = 0.00001
                #N = i_coadd_df[j,3] / i_sf_df[k,j,30]
                B = i_sf_df[k,j,34]
                area_arr.append(area)
                N_arr.append(N)
                B_arr.append(B)
                indexArr.append(j)
            
            psfsigxx_arr.append(i_sf_df[k,j,38])
            psfsigyy_arr.append(i_sf_df[k,j,39])
            psfsigxy_arr.append(i_sf_df[k,j,40])
            
            psfsigxx_uncertain_arr.append(i_sf_df[k,j,41])
            psfsigyy_uncertain_arr.append(i_sf_df[k,j,42])
            psfsigxy_uncertain_arr.append(i_sf_df[k,j,43])
            
            sigxx_err_arr.append( np.sqrt(i_sf_df[k,j,71]**2 + i_sf_df[k,j,74]**2 ))
            sigyy_err_arr.append( np.sqrt(i_sf_df[k,j,72]**2 + i_sf_df[k,j,75]**2 ))
            sigxy_err_arr.append( np.sqrt(i_sf_df[k,j,73]**2 + i_sf_df[k,j,76]**2 ))
            size_sf_arr.append(np.sqrt(i_sf_df[k,j,38] + i_sf_df[k,j,39]))
            
            e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
            e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
            error = np.sqrt((1/N + (4*area*B)/N**2) * (1/ (1-e_temp**2)) )
            
            if(error == 0 or error == None or np.isnan(error)):
                error = 10000000000
            if(badLocs[int(i_sf_df[k,j,47]), int(i_sf_df[k,j,48]), int(i_sf_df[k,j,49])] == 0):
                wtArr.append(0)
                #print ('aa')
            else:
                wtArr.append(1/error)
                
        
        

                
    #print (j,len(wtArr), len(indexArr))
    if(len(wtArr) <= 2):
        continue
    #Now correct for PSF and use monte carlo if needed
    corr_sigxx_arr=[]
    corr_sigyy_arr=[]
    corr_sigxy_arr=[]
    wtArr = np.array(wtArr)
    wtArr = wtArr**2 #(Inverse variance weight)
    wtArr = wtArr/np.sum(wtArr)
    
    
    median_size= np.nanmedian(size_sf_arr)
    sigxx_err_arr = np.array(sigxx_err_arr)
    sigyy_err_arr = np.array(sigyy_err_arr)
    sigxy_err_arr = np.array(sigxy_err_arr)
    
    
    
    validIndices = np.where(wtArr >0)[0]
    for k in range(len(sigxx_arr)):
        corr_xx = sigxx_arr[k] - psfsigxx_arr[k]
        corr_yy = sigyy_arr[k] - psfsigyy_arr[k]
        corr_xy = sigxy_arr[k] - psfsigxy_arr[k]
        corr_sigxx_arr.append(corr_xx)
        corr_sigyy_arr.append(corr_yy)
        corr_sigxy_arr.append(corr_xy)
    corr_sigxx = np.sum(corr_sigxx_arr*wtArr)/np.sum(wtArr)
    corr_sigyy = np.sum(corr_sigyy_arr*wtArr)/np.sum(wtArr)
    corr_sigxy = np.sum(corr_sigxy_arr*wtArr)/np.sum(wtArr)
    
    #corr_sigxx = weighted_quantiles_interpolate(corr_sigxx_arr, wtArr)
    #corr_sigyy = weighted_quantiles_interpolate(corr_sigyy_arr, wtArr)
    #corr_sigxy = weighted_quantiles_interpolate(corr_sigxy_arr, wtArr)
    #If any of them 0 or e2>1 then do monte carlo 
    temp = corr_sigxx +corr_sigyy - 2*np.abs(corr_sigxy)
    
    err_xx = np.sqrt(np.sum(wtArr**2 * sigxx_err_arr**2 ))
    err_yy = np.sqrt(np.sum(wtArr**2 * sigyy_err_arr**2 ))
    err_xy = np.sqrt(np.sum(wtArr**2 * sigxy_err_arr**2 ))
    
    err_xx = np.sqrt(err_xx**2 + ((median_size/30)*median_size*1.414)**2)
    err_yy = np.sqrt(err_yy**2 + ((median_size/30)*median_size*1.414)**2)
    err_xy = np.sqrt(err_xy**2 + ((median_size/30)*median_size*1.414*0.707)**2)
    
    #if(corr_sigxx<0 or corr_sigyy<0 or temp<0):
    if(True):
        cnt1 += 1
        
        a1 ,b1, c1, success = helper.correct(corr_sigxx, corr_sigyy, corr_sigxy, 0,0,0,
                                                            err_xx, err_yy, err_xy, 0,0,0)
        
        if(success < 50):
            #print (corr_sigxx,  err_xx, np.sqrt(ir_coadd_data[j,68]**2 + ir_coadd_data[j,41]**2))
            #print (corr_xx, corr_yy, err_xx, err_yy)
            failedArr.append(j)
            corr_sigxx =corr_sigyy=1
            corr_sigxy=0
            cnt2 += 1
            
            x = int(ir_coadd_data [j,10])
            y = int(ir_coadd_data [j,11])
            #f1[0].data[y-20: y+20, x-20] = 1000*j
            #f1[0].data[y-20: y+20, x+20] = 1000*j
            #f1[0].data[y-20, x-20:x+20] = 1000*j
            #f1[0].data[y+20, x-20:x+20] = 1000*j
        else:
            corr_sigxx = a1
            corr_sigyy = b1
            corr_sigxy = c1
            
        
    
   
    #Combine optimally 
    e1 = (corr_sigxx - corr_sigyy)/ (corr_sigxx+corr_sigyy)
    e2 = 2*corr_sigxy/(corr_sigxx+corr_sigyy)
    if(np.isnan(e1) or np.isnan(e2)):
        #sys.exit()
        e1 = e2 = 0.0
    
    
    
    
    master_frame[j, 0] = ir_coadd_data[j,10]
    master_frame[j, 1] = ir_coadd_data[j,11]
    master_frame[j, 2] = e1
    master_frame[j, 3] = e2
    master_frame[j, 4] = redShiftArr[j]
    master_frame[j,5] = len(validIndices)
    master_frame[j,6] = corr_sigxx
    master_frame[j,7] = corr_sigyy
    master_frame[j,8] = corr_sigxy
    master_frame[j,9] = ir_coadd_data [j,3]
    size = np.sqrt(corr_sigxx + corr_sigyy)
    #Condition for stars 
    if(np.log10(ir_coadd_data [j,3])> 1 and  size<1.5):
        master_frame[j, 2] = 0
        master_frame[j, 3] = 0
# =============================================================================
#     if(len(wtArr) > 50 and j>5000 and ir_coadd_data[j,3]> 400 ):
#         print ('aa')
#         sys.exit()
# =============================================================================
#sys.exit()
np.save('/scratch/bell/dutta26/abell_2390/master_arr_sf.npy', master_frame)
#f1.flush()

sys.exit()




















master_frame = np.load('/scratch/bell/dutta26/abell_2390/master_arr_sf.npy')
#Make the make wrt to ir coadd 
coadd_file = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits'
f=fits.open(coadd_file)
data = f[0].data
hdr =f[0].header
ySize,xSize = np.shape(data)
f.close()  
del data 



    
chopSize = 50
alphax = 1000



#Find the correct wcs
w = wcs.WCS(naxis=2)
w.wcs.crpix = [hdr['CRPIX1']/chopSize, hdr['CRPIX2']/chopSize]
w.wcs.cd = np.array([[hdr['CD1_1']*chopSize,hdr['CD1_2']], [hdr['CD2_1'], hdr['CD2_2']*chopSize]])
w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
w.wcs.cunit = [hdr['CUNIT1'], hdr['CUNIT2']]
#w.wcs.set_pv([(2, 1, 45.0)])
header = w.to_header()

mod_ySize = int(round(ySize/chopSize)) + 1
mod_xSize = int(round(xSize/chopSize)) + 1
imgE = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
imgB = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
count_img = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)

for j in range(mod_ySize):
    print (int(j*chopSize + chopSize/2))
    for k in range(mod_xSize):
        
        x_mid = int(k*chopSize + chopSize/2)
        y_mid = int(j*chopSize + chopSize/2)
        cond = np.where((master_frame[:,0] > x_mid-3000) & (master_frame[:,0] < x_mid+3000) & 
                        (master_frame[:,1] > y_mid-3000) & (master_frame[:,1] < y_mid+3000) 
                        & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,2] !=0)  ) 
        

        temp = np.copy(master_frame[cond])
        dx = temp[:,0] -x_mid
        dy = temp[:,1] -y_mid
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        
        e1 = temp[:,2]
        e2 = temp[:,3]
        tot_ellip = np.sqrt(e1**2 + e2**2)
        epar= - (e1*cos2phi+e2*sin2phi)
        eper= (e2*cos2phi-e1*sin2phi)
        #goodEllipIndices = np.where((np.abs(e1)<0.8) & (np.abs(e2)<0.8) & (r2<(3*alphax)**2) & (r2>100))
        goodEllipIndices = np.where((tot_ellip < 0.9)  & (r2<(3*alphax)**2) & (r2>100) & (tot_ellip > 0.0) )
        wt = np.exp(-(r2/(2*alphax**2) ))
        
        mean,median,std = sigma_clipped_stats(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2)
        e1sum = np.sum(epar[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        e2sum = np.sum(eper[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        count = np.sum(wt[goodEllipIndices[0]]**2 *(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2) )
        
        if(count == 0):
            e1sum = 0
            e2sum = 0
            #sys.exit()
            #continue
        count = np.sqrt(0.5*count)
        e1sum = e1sum/count
        e2sum = e2sum/count
        #print (count, e1sum, e2sum)
        
        
        imgE[j, k] = e1sum
        imgB[j ,k] = e2sum
        count_img[j ,k] = len(epar[goodEllipIndices[0]])
        del temp
        
hdu = fits.PrimaryHDU(imgE,header=header)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/EMode_sf.fits', overwrite=True)

hdu = fits.PrimaryHDU(imgB,header=header)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/BMode_sf.fits', overwrite=True)

hdu = fits.PrimaryHDU(count_img, header=header)  
hdu.writeto('/scratch/bell/dutta26/abell_2390/count_img_sf.fits', overwrite=True)
            

            





