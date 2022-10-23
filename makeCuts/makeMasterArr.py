#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:23:53 2022

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import wquantiles,sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime

f=open('/home/dutta26/zphot.out')
content = f.readlines()
f.close()

redShiftArr=[]

for j in range(len(content)):
    if (content[j][0] == '#'):
        continue
    else:
        redShiftArr.append(float((content[j].split())[1]))
ir_coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_ir_.pk1')
ir_coadd_data = np.array(ir_coadd_data)
redShiftArr = np.array(redShiftArr)        
z_min = 0.4
z_max = 10
bandList =['r','i']
coadd_data = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r_.pk1'))
bkgIndices = np.where( (redShiftArr>z_min )& (redShiftArr<z_max)& (coadd_data[:,3] > 0.01) & (coadd_data[:,3] < 1000))

temp_coadd = np.zeros((40,40))

a1,b1 = np.shape(ir_coadd_data)
master_frame = np.zeros((5,a1,10), dtype = np.float32)
master_frame1 = np.zeros((a1,10), dtype = np.float32)
upLt =1.3
for band in bandList:
    if(band =='i'):
        bandNdx = 0
    if(band == 'r'):
        bandNdx = 1
    
    
    frame_data = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_'+band+ '.npy')
    coadd_data = pd.read_pickle('/home/dutta26/codes/coaddSc_'+band+'_.pk1')
    coadd_data = np.array(coadd_data)
    a,b,c = np.shape(frame_data)
    a1,b1 = np.shape(coadd_data)
    
    zp = []
    seeing = []
    zpAvg = 25
    
    for j in range(a):
        useableIndices = np.where(frame_data[j,:,16]> 0)
        seeing.append( np.mean(frame_data[j, useableIndices[0], 16]) )
        zp.append( np.mean(frame_data[j, useableIndices[0], 17]) )
    
    
    for j in bkgIndices[0]:
        #Dont consider stars
        if(coadd_data[j, 2] == 1):
            continue
        sigxx_arr=[]
        sigyy_arr=[]
        sigxy_arr=[]
        wtArr= []
        
        
        
        lt = 1.0
        while ((lt*coadd_data[j,7] -coadd_data[j, 38])<0.5 or  (lt*coadd_data[j,8] -coadd_data[j, 38])<0.5):
            lt = lt +0.05
            if(lt>upLt):
                break
        if(lt>upLt):
            continue
        sigxx_arr.append( lt*coadd_data[j,7] -coadd_data[j, 38] )
        sigyy_arr.append( lt*coadd_data[j,8] -coadd_data[j, 39])
        sigxy_arr.append( lt*coadd_data[j,9] -coadd_data[j, 40])
        wtArr.append(1)
        
        for k in range(a):
            #Check if measurement was successful 
            if(frame_data[k,j,7] > 0 and frame_data[k,j,7] < 40 
               and frame_data[k,j,8] > 0 and frame_data[k,j,9] < 40 and frame_data[k,j,12] == 99 
               and frame_data[k,j,38] != -99 and frame_data[k,j,14] == 0):
                lt = 1.0
                while ((lt*frame_data[k,j,7] -frame_data[k,j,38])<0.5 or  (lt*frame_data[k,j,8] -frame_data[k,j,39])<0.5):
                    lt = lt +0.05
                    if(lt>upLt):
                        break
                if(lt>upLt):
                    continue
                sigxx_arr.append(lt*frame_data[k,j,7] - frame_data[k,j,38])
                sigyy_arr.append(lt*frame_data[k,j,8] - frame_data[k,j,39])
                sigxy_arr.append(lt*frame_data[k,j,9] - frame_data[k,j,40] )
                
                wtArr.append(0.5)
            #Check if FORCED measurement was successful     
            elif(frame_data[k,j,35] > 0 and frame_data[k,j,35] < 40 
               and frame_data[k,j,36] > 0 and frame_data[k,j,36] < 40 and frame_data[k,j,12] == 1
               and frame_data[k,j,38] != -99 and frame_data[k,j,14] == 0):
                
                lt = 1.0
# =============================================================================
#                 while ((lt*frame_data[k,j,35] -frame_data[k,j,38])<0.5 or  (lt*frame_data[k,j,36] -frame_data[k,j,39])<0.5):
#                     lt = lt +0.05
#                     if(lt>upLt):
#                         break
#                 if(lt>upLt):
#                     continue
# =============================================================================
                
                sigxx_arr.append(lt*frame_data[k,j,35] - frame_data[k,j,38]+1)
                sigyy_arr.append(lt*frame_data[k,j,36] - frame_data[k,j,39]+1)
                sigxy_arr.append(lt*frame_data[k,j,37] - frame_data[k,j,40] )
                #wtArr.append(0.1)
                
                wtArr.append(0.1*np.power(10,(zp[k]-zpAvg)/2.5)/(frame_data[k,j,15]* frame_data[k,j,30]*seeing[k]**2))
# =============================================================================
#                 if(j == bkgIndices[0][563]):
#                     #print(datetime.now(), 1)
#                     f=fits.open('/scratch/halstead/d/dutta26/abell_2390/'+band+'/temp/temp_coadd_'+str(k)+'.fits')
#                     data = np.array(f[0].data)
#                     f.close()
#                     #print(datetime.now(), 2)
#                     xPos = int(frame_data[k,j,10])
#                     yPos = int(frame_data[k,j,11])
#                     cut = data[yPos-20:yPos+20, xPos-20:xPos+20 ]
#                     hdu = fits.PrimaryHDU(cut) 
#                     #print (cut)
#                     hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/obj1/' + str(k)+'.fits', clobber=True)     
#                     zeoLoc = np.where(cut <= 0)
#                     if(len(zeoLoc[0]) == 0):
#                         temp_coadd = temp_coadd + cut*frame_data[k,j,30]
#                     print (k, 0.1*np.power(10,(zp[k]-zpAvg)/2.5)/(frame_data[k,j,15]* frame_data[k,j,30]*seeing[k]**2), frame_data[k,j,35] - frame_data[k,j,38]+1 )
#                     #print(datetime.now(), 3)
# =============================================================================
        
        sigxx_arr = np.array(sigxx_arr)   
        sigyy_arr = np.array(sigyy_arr)   
        sigxy_arr = np.array(sigxy_arr)   
        wtArr = np.array(wtArr)   
        master_frame[bandNdx, j, 0] = wquantiles.median(sigxx_arr ,wtArr)
        master_frame[bandNdx, j, 1] = wquantiles.median(sigyy_arr ,wtArr)
        master_frame[bandNdx, j, 2] = wquantiles.median(sigxy_arr ,wtArr)
        
# =============================================================================
#         if(j == bkgIndices[0][563]):
#             print (wtArr, sigxx_arr, wquantiles.median(sigxx_arr ,wtArr))
#             a1= wtArr
#             b1= sigxx_arr
# =============================================================================
        
        
# =============================================================================
#         master_frame[bandNdx, j, 0] = sigma_clipped_stats(sigxx_arr)[1]
#         master_frame[bandNdx, j, 1] = sigma_clipped_stats(sigyy_arr)[1]
#         master_frame[bandNdx, j, 2] = sigma_clipped_stats(sigxy_arr)[1]
# =============================================================================
        master_frame[bandNdx, j, 3] = ir_coadd_data[j, 10]
        master_frame[bandNdx, j, 4] = ir_coadd_data[j, 11]

#sys.exit()
#Converge the measurements of 5 frames to 1 single frame 
for j in bkgIndices[0]:
    success_indices = np.where((master_frame[:,j,0] != 0) & (master_frame[:,j,1] != 0 ))
    master_frame1[j,0] = np.median(master_frame[success_indices[0], j ,0])
    master_frame1[j,1] = np.median(master_frame[success_indices[0], j ,1])
    master_frame1[j,2] = np.median(master_frame[success_indices[0], j ,2])
# =============================================================================
#     master_frame1[j,0] = ir_coadd_data[j, 7] - ir_coadd_data[j, 38 ] +1
#     master_frame1[j,1] = ir_coadd_data[j, 8] - ir_coadd_data[j, 39] +1
#     master_frame1[j,2] = ir_coadd_data[j, 9] - ir_coadd_data[j, 40] 
# =============================================================================
    master_frame1[j,3] = ir_coadd_data[j, 10]
    master_frame1[j,4] = ir_coadd_data[j, 11]

       
np.save('/scratch/halstead/d/dutta26/abell_2390/masterArr.npy', master_frame1)

#Redefine master_frame1 as master frame
master_frame = np.load('/scratch/halstead/d/dutta26/abell_2390/masterArr.npy')        
#Make the make wrt to ir coadd 
coadd_file = '/scratch/halstead/d/dutta26/abell_2390/abell_ir_coadd_wted.fits'
f=fits.open(coadd_file)
data = np.array(f[0].data)
ySize,xSize = np.shape(data)
f.close()  
del data 
    
chopSize = 50
alphax = 1000

mod_ySize = int(round(ySize/chopSize)) + 5
mod_xSize = int(round(xSize/chopSize)) + 5
imgE = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
imgB = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)


#master_i = master_frame[0,:,:]
#master_r = master_frame[1,:,:]
        
        
for j in range(mod_ySize):
    print (int(j*chopSize + chopSize/2))
    for k in range(mod_xSize):
        
        x_mid = int(k*chopSize + chopSize/2)
        y_mid = int(j*chopSize + chopSize/2)
        cond = np.where((ir_coadd_data[:,10] > x_mid-3000) & (ir_coadd_data[:,10] < x_mid+3000) & 
                        (ir_coadd_data[:,11] > y_mid-3000) & (ir_coadd_data[:,11] < y_mid+3000) 
                        & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,0]> 0)
                        & (master_frame[:,1]> 0) ) 
        
        temp1 = np.copy(ir_coadd_data[cond])
        temp = np.copy(master_frame[cond])
        dx = temp1[:,10] -x_mid
        dy = temp1[:,11] -y_mid
        r2 = dx**2 + dy**2
        cos2phi = (dx**2-dy**2)/r2
        sin2phi = 2.0*dx*dy/r2
        
        e1 = (temp[:,0]-temp[:,1])/(temp[:,0]+temp[:,1])
        e2 = 2*temp[:,2]/(temp[:,0]+temp[:,1])
        
        epar= - (e1*cos2phi+e2*sin2phi)
        eper= (e2*cos2phi-e1*sin2phi)
        goodEllipIndices = np.where((np.abs(e1)<0.8) & (np.abs(e2)<0.8) & (r2<(3*alphax)**2) & (r2>100))
        wt = np.exp(-(r2/(2*alphax**2) ))
        
        mean,median,std = sigma_clipped_stats(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2)
        e1sum = np.sum(epar[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        e2sum = np.sum(eper[goodEllipIndices[0]]*wt[goodEllipIndices[0]])
        count = np.sum(wt[goodEllipIndices[0]]**2 *(e1[goodEllipIndices[0]]**2+e2[goodEllipIndices[0]]**2) )
        #print (k, len(goodEllipIndices[0]))
        #print (np.max(cos2phi), np.max(sin2phi), np.sqrt(np.max(r2)))
# =============================================================================
#         eTan = eRad =0
#         print (k, len(cond[0]))
#         for index in cond[0]:
#             success_indices = np.where((master_frame[:,index,0]>0) & (master_frame[:,index,1]>0 ) )
#             sigxx = np.nanmedian(master_frame[success_indices,index,0])
#             sigyy = np.nanmedian(master_frame[success_indices,index,1])
#             sigxy = np.nanmedian(master_frame[success_indices,index,2])
#             e1 = (sigxx - sigyy)/(sigxx+sigyy)
#             e2 = 2*sigxy/(sigxx+sigyy)
#             if(np.abs(e1)>0.8 or np.abs(e2)>0.8):
#                 continue
#             #print ('aa')
#             epar= - (e1*cos2phi[index]+e2*sin2phi[index])
#             eper= (e2*cos2phi[index]-e1*sin2phi[index])
#             
#             
#             wt = np.exp(-(r2[index]/(2*alphax**2) ))
#             if(r2[index]> 3000**2):
#                 wt = 0
#                 
#             e1sum=e1sum+epar*wt
#             e2sum=e2sum+eper*wt
#             count += (wt**2 * (e1**2+e2**2) )
# =============================================================================
            
        if(count == 0):
            e1sum = 0
            e2sum = 0
            sys.exit()
            continue
        count = np.sqrt(0.5*count)
        e1sum = e1sum/count
        e2sum = e2sum/count
        #print (count, e1sum)
        
        
        imgE[j, k] = e1sum
        imgB[j ,k] = e2sum
        del temp1, temp
        
hdu = fits.PrimaryHDU(imgE)  
hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/EMode5.fits', overwrite=True)

hdu = fits.PrimaryHDU(imgB)  
hdu.writeto('/scratch/halstead/d/dutta26/abell_2390/BMode5.fits', overwrite=True)
            

            
            
            
        
    
    
    



    
    
    
    