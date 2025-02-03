#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 18:13:11 2024

@author: dutta26
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 20:28:33 2023

@author: dutta26
"""


from astropy.io import fits
import numpy as np
import sys
import pandas as pd
from astropy.stats import sigma_clipped_stats
from datetime import datetime
import helper
from astropy import wcs
import matplotlib.pyplot as plt
lut_forcedDist = np.load('/home/dutta26/codes/forced_truncated_calib.npy')


#f1 = fits.open('/scratch/bell/dutta26/abell_2390/temp.fits', mode = 'update')

def weighted_quantiles_interpolate(values, weights, quantiles=0.5):
    i = np.argsort(values)
    c = np.cumsum(weights[i])
    return values[i[np.searchsorted(c, np.array(quantiles) * c[-1])]]


def EandB(ir_coadd_data_name, r_coadd_npy_name, i_coadd_npy_name, zFile, outFile, coadd_img, r_sf_name, i_sf_name):
    
    flag_count_arr = np.zeros(10)
    
    #redShiftArr = np.ones(15106)*9      
    r_sf_npy = np.load(r_sf_name)
    i_sf_npy = np.load(i_sf_name)
    r_coadd_npy = np.load(r_coadd_npy_name)
    i_coadd_npy = np.load(i_coadd_npy_name)
    #ir_coadd_data = pd.read_pickle(ir_coadd_data)
    #ir_coadd_data = np.array(ir_coadd_data)
    ir_coadd_data = np.load(ir_coadd_data_name)
    
    
    z_min = 0.35 - 0.25
    z_max = 0.75 - 0.25
    bandList =['g', 'r', 'i']
    if (zFile == None):
        redShiftArr = np.ones(len(ir_coadd_data))*9      
    else:
        #Read Redshifts
        f=open(zFile)
        content = f.readlines()
        f.close()
        redShiftArr=[]
        for j in range(len(content)):
            
            if (content[j][0] == '#'):
                continue
            else:
                if('eazy' in zFile):
                    if(float((content[j].split())[8]) >= 0.8 and float((content[j].split())[15])>=3):
                        redShiftArr.append(float((content[j].split())[7]))
                    else:
                        redShiftArr.append(0)
                else:
                    redShiftArr.append(float((content[j].split())[1]))
    redShiftArr = np.array(redShiftArr)     
    #g_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_g.npy')
    #r_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_r.npy')
    #i_sf_df = np.load('/scratch/bell/dutta26/abell_2390/test2_i.npy')
    #z_sf_df = np.load('/scratch/halstead/d/dutta26/abell_2390/test4t3_z.npy')
    
    #g_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_g.pk1'))
    #r_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_r.pk1'))
    #i_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_i.pk1'))
    #z_coadd_df = np.array(pd.read_pickle('/home/dutta26/codes/coaddSc_z_.pk1'))
    
    bkgIndices = np.where(  (ir_coadd_data[:,3] > 0.1) & (ir_coadd_data[:,3] < 1e10) & (ir_coadd_data[:,2] == 0))[0]
    print (len(bkgIndices))
    #sys.exit()
    
    a1,b1 = np.shape(ir_coadd_data)
    master_frame = np.zeros((a1,20), dtype = np.float32)
    cnt1 = cnt2=cnt3= cnt4=cnt5=cnt6=cnt7=cnt8=cnt9=0
    
    
    
    #Find total background 
    totBkg = 53293.7880
    totScale = 24083.1857
    
    
    failedArr=[]
    print ('aa')
    arr1=[]
    arr2=[]
    arr3=[]
    arr_sf_wt =[]
    arr_sf_wt1 =[]
# =============================================================================
#     for j in range(a1):
#         
#     #for j in [17420]:   
#         #print (j)
#         
#         if(j not in bkgIndices or ir_coadd_data [j,3] == 0  or ir_coadd_data [j,13] == 1 or ir_coadd_data [j,14] == 1 or  ir_coadd_data [j,59] == 1):
#             continue
#         
#         sigxx_arr=[]
#         sigyy_arr=[]
#         sigxy_arr=[]
#         
#         sigxx_err_arr=[]
#         sigyy_err_arr=[]
#         sigxy_err_arr=[]
#         
#         psfsigxx_arr=[]
#         psfsigyy_arr=[]
#         psfsigxy_arr=[]
#         
#         psfsigxx_uncertain_arr=[]
#         psfsigyy_uncertain_arr=[]
#         psfsigxy_uncertain_arr=[]
#         
#         size_sf_arr =[]
#         
#         wtArr= []
#         area_arr=[]
#         B_arr=[]
#         N_arr=[]
#         N_corr_arr=[]
#         indexArr=[]
#         coadd_measurement_flag = 0.0
#     
#         #If super faint use coadd measurements
#         if(ir_coadd_data [j,3] < 1.0 or r_coadd_npy[j,3]<=0 or i_coadd_npy[j,3]<=0 ):
#         
#             sigxx_arr.append(ir_coadd_data [j,35])
#             sigyy_arr.append(ir_coadd_data [j,36])
#             sigxy_arr.append(ir_coadd_data [j,37])
#             
#             psfsigxx_arr.append(ir_coadd_data [j,38])
#             psfsigyy_arr.append(ir_coadd_data [j,39])
#             psfsigxy_arr.append(ir_coadd_data [j,40])
#             
#             psfsigxx_uncertain_arr.append(ir_coadd_data [j,41])
#             psfsigyy_uncertain_arr.append(ir_coadd_data [j,42])
#             psfsigxy_uncertain_arr.append(ir_coadd_data [j,43])
#             
#             area = 2*np.pi*np.sqrt(ir_coadd_data [j,7]* ir_coadd_data [j,8] - ir_coadd_data [j,9]**2)
#             size = np.sqrt(ir_coadd_data [j,7]+ ir_coadd_data [j,8])
#             B = totBkg
#             N = ir_coadd_data [j,3] *totScale
#             e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#             e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#             e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
#             error = ir_coadd_data[j,68] #np.sqrt((size**4/N + (4*np.pi*size**6 *B)/N**2)) 
#             if(error <= 0 or error== None or np.isnan(error)):
#                 error = 100000000
#                 
#             
#             wtArr.append(1/error)
#             
#             area_arr.append(area)
#             N_arr.append(N)
#             B_arr.append(B)
#             indexArr.append(j)
#             coadd_measurement_flag = 1
#             
#             sigxx_err_arr.append( np.sqrt(ir_coadd_data [j,71]**2 + ir_coadd_data [j,41]**2 ))
#             sigyy_err_arr.append( np.sqrt(ir_coadd_data [j,72]**2 + ir_coadd_data [j,42]**2 ))
#             sigxy_err_arr.append( np.sqrt(ir_coadd_data [j,73]**2 + ir_coadd_data [j,43]**2 ))
#             
#             size_sf_arr.append(np.sqrt(ir_coadd_data [j,38] + ir_coadd_data [j,39]))
#     
#         else:
#             
#             #First do r band 
#             a,b,c = np.shape(r_sf_npy)
#             for k in range(a):
#                 if( r_coadd_npy[j,7]<=0 or r_coadd_npy[j,8]<=0):
#                     flag_count_arr[7]+= 1
#                     continue
#                 if( abs((ir_coadd_data[j,7] - r_coadd_npy[j,7])/ir_coadd_data[j,7] ) >2  ):
#                     flag_count_arr[6]+= 1
#                     continue
#                 if( abs((ir_coadd_data[j,8] - r_coadd_npy[j,8])/ir_coadd_data[j,8] ) >2  ):
#                     flag_count_arr[5]+= 1
#                     continue
#                 
#                 
#                 if(r_sf_npy[k,j,12] <= 0  or r_sf_npy[k,j,14]> 0 or r_sf_npy[k,j,38]==-99):
#                     flag_count_arr[3]+= 1
#                     continue
#                 if(np.sum(r_sf_npy[k,j,61:67])>=1 or r_sf_npy[k,j,38]==-99):
#                     flag_count_arr[4]+= 1
#                     continue
#                 if(r_sf_npy[k,j,12] == 99):
#                     if(r_sf_npy[k,j,7] < 0 or r_sf_npy[k,j,7]>70 or r_sf_npy[k,j,8] < 0 or r_sf_npy[k,j,8]>70):
#                         flag_count_arr[0]+= 1
#                         continue
#                     if(r_sf_npy[k,j,7] == None or np.isnan(r_sf_npy[k,j,7]) or r_sf_npy[k,j,8] == None or np.isnan(r_sf_npy[k,j,8])):
#                         flag_count_arr[1]+= 1
#                         continue
#                     if(r_sf_npy[k,j,3]<= 0 or r_sf_npy[k,j,3]== None or np.isnan(r_sf_npy[k,j,3]) ):
#                         flag_count_arr[2]+= 1
#                         continue
#                         
#                     
#                 if(r_sf_npy[k,j,12] == 1):
#                     if(r_sf_npy[k,j,35] < 0 or r_sf_npy[k,j,35]>70 or r_sf_npy[k,j,36] < 0 or r_sf_npy[k,j,36]>70):
#                         flag_count_arr[0]+= 1
#                         continue
#                     if(r_sf_npy[k,j,35] == None or np.isnan(r_sf_npy[k,j,35]) or r_sf_npy[k,j,36] == None or np.isnan(r_sf_npy[k,j,36])):
#                         flag_count_arr[1]+= 1
#                         continue
#                     if(r_sf_npy[k,j,31]<= 0 or r_sf_npy[k,j,31]== None or np.isnan(r_sf_npy[k,j,31]) ):
#                         flag_count_arr[2]+= 1
#                         continue
#                 
#                 
#                 if(r_sf_npy[k,j,12] == 99):
#                     sigxx_arr.append(r_sf_npy[k,j,51])
#                     sigyy_arr.append(r_sf_npy[k,j,52])
#                     sigxy_arr.append(r_sf_npy[k,j,53])
#                     area = 2*np.pi*np.sqrt(r_sf_npy[k,j,7]* r_sf_npy[k,j,8] - r_sf_npy[k,j,9]**2)
#                     size = np.sqrt(r_sf_npy[k,j,7]+r_sf_npy[k,j,8])
#                     N =r_sf_npy[k,j,3]
#                     #N = r_coadd_df [j,3] / r_sf_npy[k,j,30]
#                     B = r_sf_npy[k,j,6]
#                     area_arr.append(area)
#                     if(N< 0 or N>1e5 ):
#                         N = 0.00001
#                     N_arr.append(N)
#                     N_corr_arr.append(N*r_sf_npy[k,j,30])
#                     B_arr.append(B)
#                     indexArr.append(j)
#                    
#                     
#                 elif(r_sf_npy[k,j,12] == 1):
#                     sigxx_arr.append(r_sf_npy[k,j,51])
#                     sigyy_arr.append(r_sf_npy[k,j,52])
#                     sigxy_arr.append(r_sf_npy[k,j,53])
#                     temp_xx = r_sf_npy[k,j,51]
#                     temp_yy = r_sf_npy[k,j,52]
#                     temp_xy = r_sf_npy[k,j,53]
#                     area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
#                     size = np.sqrt(temp_xx+temp_yy)
#                     N =r_sf_npy[k,j,31]
#                     if(N< 0 or N>1e5):
#                         N = 0.00001
#                     B = r_sf_npy[k,j,34]
#                     area_arr.append(area)
#                     N_arr.append(N)
#                     N_corr_arr.append(N*r_sf_npy[k,j,30])
#                     B_arr.append(B)
#                     indexArr.append(j)
#                     
#                 
#                 psfsigxx_arr.append(r_sf_npy[k,j,38])
#                 psfsigyy_arr.append(r_sf_npy[k,j,39])
#                 psfsigxy_arr.append(r_sf_npy[k,j,40])
#                 
#                 psfsigxx_uncertain_arr.append(r_sf_npy[k,j,41])
#                 psfsigyy_uncertain_arr.append(r_sf_npy[k,j,42])
#                 psfsigxy_uncertain_arr.append(r_sf_npy[k,j,43])
#                 
#                 sigxx_err_arr.append( np.sqrt(r_sf_npy[k,j,71]**2 + r_sf_npy[k,j,41]**2 ))
#                 sigyy_err_arr.append( np.sqrt(r_sf_npy[k,j,72]**2 + r_sf_npy[k,j,42]**2 ))
#                 sigxy_err_arr.append( np.sqrt(r_sf_npy[k,j,73]**2 + r_sf_npy[k,j,43]**2 ))
#                 size_sf_arr.append(np.sqrt(r_sf_npy[k,j,38] + r_sf_npy[k,j,39]))
#                 
#                 
#                 e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#                 e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#                 e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
#                 
#                 turb_xx = np.sqrt(r_sf_npy[k,j,41]**2 )
#                 turb_yy = np.sqrt(r_sf_npy[k,j,42]**2 )
#                 turb_avg = 0.5*(turb_xx+turb_yy)
#                 
#                
#                 
#                 
#                 error =np.sqrt(size**4/N + (4*np.pi*size**6*B)/N**2 +turb_avg)
#                 if(error == 0 or error == None or np.isnan(error)):
#                     error = 10000000000
#                     wtArr.append(1/error)
#                 else:
#                     
#                     wtArr.append(1/error)
#                     
#              
#             print (len(wtArr))   
#             #Now do I band
#             a,b,c = np.shape(i_sf_npy)
#             for k in range(a):
#                 if( i_coadd_npy[j,7]<=0 or i_coadd_npy[j,8]<=0):
#                     continue
#                 
#                 if( abs((ir_coadd_data[j,7] - i_coadd_npy[j,7])/ir_coadd_data[j,7] ) >2  ):
#                     continue
#                 if( abs((ir_coadd_data[j,8] - i_coadd_npy[j,8])/ir_coadd_data[j,8] ) >2  ):
#                     continue
#                
#                 if(i_sf_npy[k,j,12] <= 0  or i_sf_npy[k,j,14]> 0) or i_sf_npy[k,j,38]==-99:
#                     continue
#                 if(np.sum(i_sf_npy[k,j,61:67])>=1 or i_sf_npy[k,j,38]==-99):
#                     continue
#                 if(i_sf_npy[k,j,12] == 99):
#                     if(i_sf_npy[k,j,7] < 0 or i_sf_npy[k,j,7]>70 or i_sf_npy[k,j,8] < 0 or i_sf_npy[k,j,8]>70):
#                         continue
#                     if(i_sf_npy[k,j,7] == None or np.isnan(i_sf_npy[k,j,7]) or i_sf_npy[k,j,8] == None or np.isnan(i_sf_npy[k,j,8])):
#                         continue
#                     if(i_sf_npy[k,j,3]<= 0 or i_sf_npy[k,j,3]== None or np.isnan(i_sf_npy[k,j,3]) ):
#                         continue
#                         
#                     
#                 if(i_sf_npy[k,j,12] == 1):
#                     if(i_sf_npy[k,j,35] < 0 or i_sf_npy[k,j,35]>70 or i_sf_npy[k,j,36] < 0 or i_sf_npy[k,j,36]>70):
#                         continue
#                     if(i_sf_npy[k,j,35] == None or np.isnan(i_sf_npy[k,j,35]) or i_sf_npy[k,j,36] == None or np.isnan(i_sf_npy[k,j,36])):
#                         continue
#                     if(i_sf_npy[k,j,31]<= 0 or i_sf_npy[k,j,31]== None or np.isnan(i_sf_npy[k,j,31]) ):
#                         continue
#                 
#                 
#                 if(i_sf_npy[k,j,12] == 99):
#                     sigxx_arr.append(i_sf_npy[k,j,51])
#                     sigyy_arr.append(i_sf_npy[k,j,52])
#                     sigxy_arr.append(i_sf_npy[k,j,53])
#                     area = 2*np.pi*np.sqrt(i_sf_npy[k,j,7]* i_sf_npy[k,j,8] - i_sf_npy[k,j,9]**2)
#                     size = np.sqrt(i_sf_npy[k,j,7]+ i_sf_npy[k,j,8])
#                     N =i_sf_npy[k,j,3]
#                     #N = i_coadd_df [j,3] / i_sf_npy[k,j,30]
#                     B = i_sf_npy[k,j,6]
#                     area_arr.append(area)
#                     if(N< 0 or N>1e5 ):
#                         N = 0.00001
#                     N_arr.append(N)
#                     N_corr_arr.append(N*i_sf_npy[k,j,30])
#                     B_arr.append(B)
#                     indexArr.append(j)
#                     
#                    
#                     
#                 elif(i_sf_npy[k,j,12] == 1):
#                     sigxx_arr.append(i_sf_npy[k,j,51])
#                     sigyy_arr.append(i_sf_npy[k,j,52])
#                     sigxy_arr.append(i_sf_npy[k,j,53])
#                     temp_xx = i_sf_npy[k,j,51] 
#                     temp_yy = i_sf_npy[k,j,52] 
#                     temp_xy = i_sf_npy[k,j,53] 
#                     area = 2*np.pi*np.sqrt(temp_xx*temp_yy - temp_xy**2)
#                     size = np.sqrt(temp_xx+temp_yy)
#                     N = i_sf_npy[k,j,31]
#                     if(N< 0 or N>1e5 ):
#                         N = 0.00001
#                     #N = i_coadd_df[j,3] / i_sf_npy[k,j,30]
#                     B = i_sf_npy[k,j,34]
#                     area_arr.append(area)
#                     N_arr.append(N)
#                     N_corr_arr.append(N*i_sf_npy[k,j,30])
#                     B_arr.append(B)
#                     indexArr.append(j)
#                     
#                     
#                 
#                 psfsigxx_arr.append(i_sf_npy[k,j,38])
#                 psfsigyy_arr.append(i_sf_npy[k,j,39])
#                 psfsigxy_arr.append(i_sf_npy[k,j,40])
#                 
#                 psfsigxx_uncertain_arr.append(i_sf_npy[k,j,41])
#                 psfsigyy_uncertain_arr.append(i_sf_npy[k,j,42])
#                 psfsigxy_uncertain_arr.append(i_sf_npy[k,j,43])
#                 
#                 sigxx_err_arr.append( np.sqrt(i_sf_npy[k,j,71]**2 + i_sf_npy[k,j,41]**2 ))
#                 sigyy_err_arr.append( np.sqrt(i_sf_npy[k,j,72]**2 + i_sf_npy[k,j,42]**2 ))
#                 sigxy_err_arr.append( np.sqrt(i_sf_npy[k,j,73]**2 + i_sf_npy[k,j,43]**2 ))
#                 size_sf_arr.append(np.sqrt(i_sf_npy[k,j,38] + i_sf_npy[k,j,39]))
#                 
#                 e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#                 e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
#                 e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
#                 
#                 turb_xx = np.sqrt(i_sf_npy[k,j,41]**2 )
#                 turb_yy = np.sqrt(i_sf_npy[k,j,42]**2 )
#                 turb_avg = 0.5*(turb_xx+turb_yy)
#                 
#                 
#                 
#                 error = np.sqrt(size**4/N + (4*np.pi*size**6*B)/N**2 +turb_avg)
#                 
#                 if(error == 0 or error == None or np.isnan(error)):
#                     error = 10000000000
#                     wtArr.append(1/error)
#                 else:
#                     wtArr.append(1/error)
#                     
#             
#             
#     
#                     
#         print (j,len(wtArr), len(indexArr))
#         
# # =============================================================================
# #         if(len(wtArr) <= 0):
# #             
# #             sigxx_arr=[]
# #             sigyy_arr=[]
# #             sigxy_arr=[]
# #             
# #             sigxx_err_arr=[]
# #             sigyy_err_arr=[]
# #             sigxy_err_arr=[]
# #             
# #             psfsigxx_arr=[]
# #             psfsigyy_arr=[]
# #             psfsigxy_arr=[]
# #             
# #             psfsigxx_uncertain_arr=[]
# #             psfsigyy_uncertain_arr=[]
# #             psfsigxy_uncertain_arr=[]
# #             
# #             size_sf_arr =[]
# #             
# #             wtArr= []
# #             area_arr=[]
# #             B_arr=[]
# #             N_arr=[]
# #             N_corr_arr=[]
# #             indexArr=[]
# #             sigxx_arr.append(ir_coadd_data [j,35])
# #             sigyy_arr.append(ir_coadd_data [j,36])
# #             sigxy_arr.append(ir_coadd_data [j,37])
# #             
# #             psfsigxx_arr.append(ir_coadd_data [j,38])
# #             psfsigyy_arr.append(ir_coadd_data [j,39])
# #             psfsigxy_arr.append(ir_coadd_data [j,40])
# #             
# #             psfsigxx_uncertain_arr.append(ir_coadd_data [j,41])
# #             psfsigyy_uncertain_arr.append(ir_coadd_data [j,42])
# #             psfsigxy_uncertain_arr.append(ir_coadd_data [j,43])
# #             
# #             area = 2*np.pi*np.sqrt(ir_coadd_data [j,7]* ir_coadd_data [j,8] - ir_coadd_data [j,9]**2)
# #             size = np.sqrt(ir_coadd_data [j,7]+ ir_coadd_data [j,8])
# #             B = totBkg
# #             N = ir_coadd_data [j,3] *totScale
# #             e1_temp = (sigxx_arr[len(sigxx_arr)-1] - sigyy_arr[len(sigyy_arr)-1])/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
# #             e2_temp = 2* sigxy_arr[len(sigxy_arr)-1]/(sigxx_arr[len(sigxx_arr)-1] + sigyy_arr[len(sigyy_arr)-1])
# #             e_temp = np.sqrt(e1_temp**2 + e2_temp**2)
# #             turb_xx = np.sqrt(ir_coadd_data[j,41]**2 )
# #             turb_yy = np.sqrt(ir_coadd_data[j,42]**2 )
# #             turb_avg = 0.5*(turb_xx+turb_yy)
# #             
# #             
# #             error = np.sqrt(ir_coadd_data[j,68] + turb_avg) #np.sqrt((size**4/N + (4*np.pi*size**6 *B)/N**2)) 
# #             if(error <= 0 or error== None or np.isnan(error)):
# #                 error = 100000000
# #                 
# #             
# #             wtArr.append(1/error)
# #             
# #             area_arr.append(area)
# #             N_arr.append(N)
# #             B_arr.append(B)
# #             indexArr.append(j)
# #             coadd_measurement_flag = 1
# #             
# #             sigxx_err_arr.append( np.sqrt(ir_coadd_data [j,71]**2 + ir_coadd_data [j,41]**2 ))
# #             sigyy_err_arr.append( np.sqrt(ir_coadd_data [j,72]**2 + ir_coadd_data [j,42]**2 ))
# #             sigxy_err_arr.append( np.sqrt(ir_coadd_data [j,73]**2 + ir_coadd_data [j,43]**2 ))
# #             
# #             size_sf_arr.append(np.sqrt(ir_coadd_data [j,38] + ir_coadd_data [j,39]))
# # =============================================================================
#             
# # =============================================================================
# #         if(len(wtArr)<=50):
# #             continue
# # =============================================================================
#         
#         if(len(wtArr)== 1):
#             cnt3+= 1
#             
#         N_corr_arr = np.array(N_corr_arr)
#         mean_N,med_N, std_N = sigma_clipped_stats(N_corr_arr, cenfunc=np.median)
#         
#              
#         #Now correct for PSF and use monte carlo if needed
#         corr_sigxx_arr=[]
#         corr_sigyy_arr=[]
#         corr_sigxy_arr=[]
#         wtArr = np.array(wtArr)
#         wtArr = wtArr/np.sum(wtArr)
#         #print (wtArr)
#         wtArr = wtArr**2#(Inverse variance weight)
#         wtArr = wtArr/np.sum(wtArr)
#         #print (list(wtArr))
#         
#         
#         median_size= np.nanmedian(size_sf_arr)
#         sigxx_err_arr = np.array(sigxx_err_arr)
#         sigyy_err_arr = np.array(sigyy_err_arr)
#         sigxy_err_arr = np.array(sigxy_err_arr)
#         sigxx_arr = np.array(sigxx_arr)
#         sigyy_arr = np.array(sigyy_arr)
#         sigxy_arr = np.array(sigxy_arr)
#         sigxx_err_arr = np.array(sigxx_err_arr)
#         sigyy_err_arr = np.array(sigyy_err_arr)
#         sigxy_err_arr = np.array(sigxy_err_arr)
#         psfsigxx_arr = np.array(psfsigxx_arr)
#         psfsigyy_arr = np.array(psfsigyy_arr)
#         psfsigxy_arr = np.array(psfsigxy_arr)
#         
#         mean_xx, med_xx, std_xx = sigma_clipped_stats(sigxx_arr - psfsigxx_arr, cenfunc=np.median)
#         mean_yy, med_yy, std_yy = sigma_clipped_stats(sigyy_arr - psfsigyy_arr, cenfunc=np.median)
#         mean_xy, med_xy, std_xy = sigma_clipped_stats(sigxy_arr - psfsigxy_arr, cenfunc=np.median)
#         
#         validIndices = np.where((wtArr >0) & (N_corr_arr< (med_N+3*std_N)) & (N_corr_arr> (med_N-3*std_N)) 
#                                 & (N_corr_arr> 0 ) )[0]
#         
#         for k in range(len(sigxx_arr)):
#             corr_xx = sigxx_arr[k] - psfsigxx_arr[k]
#             corr_yy = sigyy_arr[k] - psfsigyy_arr[k]
#             corr_xy = sigxy_arr[k] - psfsigxy_arr[k]
#             corr_sigxx_arr.append(corr_xx)
#             corr_sigyy_arr.append(corr_yy)
#             corr_sigxy_arr.append(corr_xy)
#             if(k not in validIndices and len(wtArr)>10):
#                 wtArr[k] = 0
#             
#         print (wtArr)
#         wtArr = wtArr/np.sum(wtArr)
#         #corr_sigxx = np.sum(corr_sigxx_arr*1/sigxx_err_arr**4)/np.sum(1/sigxx_err_arr**4)
#         #corr_sigyy = np.sum(corr_sigyy_arr*1/sigxx_err_arr**4)/np.sum(1/sigxx_err_arr**4)
#         #corr_sigxy = np.sum(corr_sigxy_arr*1/sigxx_err_arr**4)/np.sum(1/sigxx_err_arr**4)
#         corr_sigxx = np.sum(corr_sigxx_arr*wtArr)/np.sum(wtArr)
#         corr_sigyy = np.sum(corr_sigyy_arr*wtArr)/np.sum(wtArr)
#         corr_sigxy = np.sum(corr_sigxy_arr*wtArr)/np.sum(wtArr)
#         #print (corr_sigxx, ir_coadd_data[j,35]-ir_coadd_data[j,38])
#         #corr_sigxx = weighted_quantiles_interpolate(corr_sigxx_arr, wtArr)
#         #corr_sigyy = weighted_quantiles_interpolate(corr_sigyy_arr, wtArr)
#         #corr_sigxy = weighted_quantiles_interpolate(corr_sigxy_arr, wtArr)
#         #print(corr_sigxx_arr)
#         #corr_sigxx = np.median(corr_sigxx_arr)
#         #corr_sigxy = np.median(corr_sigxy_arr)
#         #If any of them 0 or e2>1 then do monte carlo 
#         temp = corr_sigxx +corr_sigyy - 2*np.abs(corr_sigxy)
#         
#         #print (N_arr, N_corr_arr)
#         
#         err_xx = np.sqrt(np.sum(wtArr**2 * sigxx_err_arr**2 ))
#         err_yy = np.sqrt(np.sum(wtArr**2 * sigyy_err_arr**2 ))
#         err_xy = np.sqrt(np.sum(wtArr**2 * sigxy_err_arr**2 ))
#         
#         err_xx = np.sqrt(err_xx**2 )
#         err_yy = np.sqrt(err_yy**2 )
#         err_xy = np.sqrt(err_xy**2 )
# # =============================================================================
# #         if(corr_sigxx<0 or corr_sigyy<0 or temp<0):
# #             #print ('wtf',corr_sigxx, ir_coadd_data [j,68], ir_coadd_data [j,41],ir_coadd_data [j,7], ir_coadd_data [j,3],corr_sigxx_arr, j)
# #         #if(True):
# #             #print (corr_sigxx,  corr_sigyy, corr_sigxy)
# #             cnt1 += 1
# #             
# #             a1 ,b1, c1, success = helper.correct(corr_sigxx, corr_sigyy, corr_sigxy, 0,0,0,
# #                                                                 err_xx, err_yy, err_xy, 0,0,0)
# #             
# #             if(success < 500000):
# #                 #print (corr_sigxx,  err_xx, np.sqrt(ir_coadd_data[j,68]**2 + ir_coadd_data[j,41]**2))
# #                 #print (corr_xx, corr_yy, err_xx, err_yy)
# #                 failedArr.append(j)
# #                 corr_sigxx =corr_sigyy=0.01
# #                 corr_sigxy=0
# #                 cnt2 += 1
# #                 
# #                 x = int(ir_coadd_data [j,10])
# #                 y = int(ir_coadd_data [j,11])
# #                 #f1[0].data[y-20: y+20, x-20] = 1000*j
# #                 #f1[0].data[y-20: y+20, x+20] = 1000*j
# #                 #f1[0].data[y-20, x-20:x+20] = 1000*j
# #                 #f1[0].data[y+20, x-20:x+20] = 1000*j
# #             else:
# #                 corr_sigxx = a1
# #                 corr_sigyy = b1
# #                 corr_sigxy = c1
# # =============================================================================
#                 
#             
#         
#         
#         #Combine optimally 
#         e1 = (corr_sigxx - corr_sigyy)/ (corr_sigxx+corr_sigyy)
#         e2 = 2*corr_sigxy/(corr_sigxx+corr_sigyy)
#         if(np.isnan(e1) or np.isnan(e2)):
#             #sys.exit()
#             e1 = e2 = 0.0
#         
#         xx_coadd = ir_coadd_data [j,35]-ir_coadd_data [j,38]
#         yy_coadd = ir_coadd_data [j,36]-ir_coadd_data [j,39]
#         xy_coadd = ir_coadd_data [j,37]+ir_coadd_data [j,40]
#         e1_coadd = (xx_coadd-yy_coadd)/(xx_coadd+yy_coadd)
#         e2_coadd = 2*xy_coadd/(xx_coadd+yy_coadd)
#         
#         #print (e1, e1_coadd, j)
#         
#         #print (corr_sigxx-ir_coadd_data [j,35]+ir_coadd_data [j,38], corr_sigyy-ir_coadd_data [j,36]+ir_coadd_data [j,39])
#         arr1.append(corr_sigxx - ir_coadd_data [j,35]+ir_coadd_data [j,38] )
#         arr2.append(corr_sigyy - ir_coadd_data [j,36]+ir_coadd_data [j,39] )
#         arr3.append(corr_sigxy-  ir_coadd_data [j,37]+ir_coadd_data [j,40] )
#         master_frame[j, 0] = ir_coadd_data[j,10]
#         master_frame[j, 1] = ir_coadd_data[j,11]
#         master_frame[j, 2] = e1
#         master_frame[j, 3] = e2
#         master_frame[j, 4] = redShiftArr[j]
#         master_frame[j,5] = len(validIndices)
#         master_frame[j,6] = corr_sigxx
#         master_frame[j,7] = corr_sigyy
#         master_frame[j,8] = corr_sigxy
#         master_frame[j,9] = ir_coadd_data [j,3]
#         master_frame[j,10] = err_xx
#         master_frame[j,11] = err_yy
#         master_frame[j,12] = err_xy
#         master_frame[j,13] = np.sqrt( err_xx**2 * 4* (corr_sigxx**2+corr_sigyy**2)/(corr_sigxx+corr_sigyy)**4)
#         master_frame[j,14] = np.sqrt( (err_xx**2 * 8* (corr_sigxy**2)/(corr_sigxx+corr_sigyy)**4) + (err_xy**2 * 4/(corr_sigxx+corr_sigyy)**2))
#         #master_frame[j,15] = np.sqrt((e1_coadd**2 *master_frame[j,13]**2 / (e1_coadd**2 + e2_coadd**2)) + (e2_coadd**2 *master_frame[j,14]**2 / (e1_coadd**2 + e2_coadd**2)) )
#         master_frame[j,15] = np.sqrt((e1**2 *master_frame[j,13]**2 / (e1**2 + e2**2)) + (e2**2 *master_frame[j,14]**2 / (e1**2 + e2**2)) )
# 
#         size = np.sqrt(corr_sigxx + corr_sigyy)
# # =============================================================================
# #         if(j==1387):
# #             print (xx_coadd, yy_coadd,corr_sigxx, master_frame[j,6], corr_sigyy_arr, corr_sigxx_arr, N_corr_arr)
# #             print (np.polyfit(corr_sigyy_arr, corr_sigxx_arr, deg=1, w=wtArr))
# #             return
# # =============================================================================
#         #Condition for stars 
#         if(np.log10(ir_coadd_data [j,3])> 1 and  size<=0):
#             print ('abc', ir_coadd_data [j,2], size)
#             
#     np.save(outFile, master_frame)
#     np.save('/scratch/bell/dutta26/abell_2390/fail_stat.npy',flag_count_arr)
#     #f1.flush()
#     
#     #sys.exit()
#     arr1=np.array(arr1)
#     arr2=np.array(arr2)
#     arr3=np.array(arr3)
#     arr1=arr1[arr1<3]
#     arr1=arr1[arr1>-3]
#     arr2=arr2[arr2<3]
#     arr2=arr2[arr2>-3]
#     arr3=arr3[arr3<3]
#     arr3=arr3[arr3>-3]
#     n, bins, patches = plt.hist(x=arr1, bins=60, histtype=u'step', color='r', label='xx') 
#     n, bins, patches = plt.hist(x=arr2, bins=60, histtype=u'step', color='b', label='yy')  
#     n, bins, patches = plt.hist(x=arr3, bins=60, histtype=u'step', color='k', label='xy')                      
#     plt.xlabel('Deviation of sigmas from I+R coadd sigmas')
#     plt.legend()
#     plt.savefig('/scratch/bell/dutta26/abell_2390/sig_dev.png')
#     plt.close()
#     print (arr1, arr2, arr3)
#     print (cnt1, cnt3, len(bkgIndices))
#     
# =============================================================================

    
    
    #return
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    master_frame = np.load(outFile)
    #Make the make wrt to ir coadd 
    
    f=fits.open(coadd_img)
    data = f[0].data
    hdr =f[0].header
    ySize,xSize = np.shape(data)
    f.close()  
    del data 
    
    
    
        
    chopSize = 50
    alphax = 400
    
    
    
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
    err_img = np.zeros((mod_ySize, mod_xSize), dtype = np.float32)
    
    loc =np.where((master_frame[:,15] >0) & (master_frame[:,2] !=0) &
                   (ir_coadd_data[:,2] == 0))[0] 
    print(len(loc), "median determinin length")
    median_ellip_err = np.median(master_frame[loc,15])
    for j in range(80, 600):
        print (int(j*chopSize + chopSize/2))
        for k in range(80, 600):
            
            x_mid = int(k*chopSize + chopSize/2)
            y_mid = int(j*chopSize + chopSize/2)
            cond = np.where((master_frame[:,0] > x_mid-2000) & (master_frame[:,0] < x_mid+2000) & 
                            (master_frame[:,1] > y_mid-2000) & (master_frame[:,1] < y_mid+2000) 
                            & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,2] !=0) &
                            (ir_coadd_data[:,2] == 0) & (master_frame[:,6] < 49) & 
                            (master_frame[:,7] < 49)) [0]
            
    
            temp = np.copy(master_frame[cond,:])
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
            goodEllipIndices = np.where( (r2<(3*alphax)**2) & (r2>400) & (tot_ellip > 0.0) )[0]
            wt = np.exp(-(r2/(2*alphax**2) )) * ir_coadd_data[cond, 3] * 1/(2*np.pi*(alphax*0.001833)**2)
            #wt = np.sqrt(1/r2)*(1- (1+ r2/(2*200**2))*np.exp(-(r2/(2*200**2) )))
            #wt = (1+ r2/(2*alphax**2))*np.exp(-(r2/(2*alphax**2) ))
            
            loc = np.where(temp[:,15]< median_ellip_err/3)
            temp[loc,15] = median_ellip_err/3
            wt_ellip_err = 1/temp[:,15]**2
            
            
            wt_tild = wt_ellip_err[goodEllipIndices]**0.5/np.sum(wt_ellip_err[goodEllipIndices]**0.5)
            fudge_fact = 1 /(2*(1-np.sum(tot_ellip[goodEllipIndices]**2 * wt_tild)/np.sum(wt_tild) ))
            
            e1sum = np.sum(epar[goodEllipIndices]*wt[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
            e2sum = np.sum(eper[goodEllipIndices]*wt[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
            count = np.sum(wt[goodEllipIndices]**2 * wt_ellip_err[goodEllipIndices]**1)
            #count = np.sum(wt[goodEllipIndices]**2 * (e1[goodEllipIndices]**2+e2[goodEllipIndices]**2))
# =============================================================================
#             e1sum = np.sum(epar[goodEllipIndices]*wt[goodEllipIndices])
#             e2sum = np.sum(eper[goodEllipIndices]*wt[goodEllipIndices])
#             count = np.sum(wt[goodEllipIndices]**2 *(e1[goodEllipIndices]**2+e2[goodEllipIndices]**2) )
# =============================================================================
            
            if(count == 0):
                e1sum = 0
                e2sum = 0
                #sys.exit()
                #continue
            count = np.sqrt(0.5*count)
            if(len(epar[goodEllipIndices])>0):
                e1sum = e1sum*fudge_fact/count
                e2sum = e2sum*fudge_fact/count
            else:
                e1sum, e2sum = 0,0
            
            
            
            imgE[j, k] = np.sum(wt)
            imgB[j ,k] = e2sum
            count_img[j ,k] = len(epar[goodEllipIndices])
            err_img[j ,k] = count
            del temp
            
    hdu = fits.PrimaryHDU(imgE,header=header)  
    hdu.writeto('/scratch/bell/dutta26/abell_2390/proj_density_low.fits', overwrite=True)
    
    
                
    
                
img = '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.fits'
ebModeFile = '/scratch/bell/dutta26/abell_2390/master_arr_sf.npy'
#zFile = '/home/dutta26/zphot_2390.out'
zFile = '/home/dutta26/photz_eazy.zout'

EandB( '/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy',  '/scratch/bell/dutta26/abell_2390/abell_r_coadd.npy', '/scratch/bell/dutta26/abell_2390/abell_i_coadd.npy', zFile, ebModeFile, img, '/scratch/bell/dutta26/abell_2390/r_sf.npy', '/scratch/bell/dutta26/abell_2390/i_sf.npy')
   
    
    
