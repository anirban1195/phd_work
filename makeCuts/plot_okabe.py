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

f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_sf.fits')
kappa_map= f[0].data
f.close()


zArr=[]
z_min = 0.4
z_max = 2
zFile = '/home/dutta26/photz_eazy.zout'
if (zFile == None):
    redShiftArr = 9
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
                if(float((content[j].split())[8]) >= 0.8):
                    redShiftArr.append(float((content[j].split())[7]))
                else:
                    redShiftArr.append(0)
            else:
                redShiftArr.append(float((content[j].split())[1]))
redShiftArr = np.array(redShiftArr)      
master_frame = np.load('/scratch/bell/dutta26/abell_2390/master_arr_sf.npy')
ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
loc =np.where((master_frame[:,15] >0) & (master_frame[:,2] !=0) &
               (ir_coadd_data[:,2] == 0))[0] 
print(len(loc), "median determinin length")
median_ellip_err = np.median(master_frame[loc,15])

x_mid = 11600
y_mid = 17850
gamma_t_arr =[]
gamma_err_arr=[]
distArr=[]
dist_min = 0
dist_max = 200
x_mid_kappa = int(round(x_mid/50))
y_mid_kappa = int(round(y_mid/50))
kappa_dist_arr =np.zeros(np.shape(kappa_map))       
y,x = np.shape(kappa_map)
for j in range(x):
    for k in range(y):
        dist = np.sqrt( (j-x_mid_kappa)**2 + (k-y_mid_kappa)**2)
        kappa_dist_arr[k,j] = dist

r_min =dist_max= 50
r_max= 6000
bins = 10
delta = (r_max**2 - r_min**2)/10
for j in range(10):
    dist_min = dist_max
    dist_max = np.sqrt(dist_max**2 + delta)
    print (dist_min, dist_max)
    width = 7000
    cond = np.where((master_frame[:,0] > x_mid-width) & (master_frame[:,0] < x_mid+width) & 
                    (master_frame[:,1] > y_mid-width) & (master_frame[:,1] < y_mid+width) 
                    & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,2] !=0) &
                    (ir_coadd_data[:,2] == 0) & (master_frame[:,6] < 64) & 
                     (master_frame[:,7] < 64)) [0]
    z_temp = redShiftArr[cond]
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
    goodEllipIndices = np.where( (r2<(dist_max)**2) & (r2>(dist_min**2)) & (tot_ellip > 0.0) )[0]
    
    
    loc = np.where(temp[:,15]< median_ellip_err/3)
    temp[loc,15] = median_ellip_err/3
    wt_ellip_err = 1/temp[:,15]**2
    
    
    wt_tild = wt_ellip_err[goodEllipIndices]**0.5/np.sum(wt_ellip_err[goodEllipIndices]**0.5)
    fudge_fact = 1 /(2*(1-np.sum(tot_ellip[goodEllipIndices]**2 * wt_tild)/np.sum(wt_tild) ))
    
    #Process kappa
    radius = dist_max/50
    loc = np.where(kappa_dist_arr<radius)
    avg_kappa = np.mean(kappa_map[loc[0], loc[1]])
    #print (avg_kappa, 'Avg kappa')
    
    
    e1sum = np.sum(epar[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
    e2sum = np.sum(eper[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
    count = np.sum(wt_ellip_err[goodEllipIndices])
    print(len(goodEllipIndices))
    gamma_t_arr.append(e1sum*fudge_fact/(count*(1-avg_kappa)))
    distArr.append(dist_max*0.11/60)
    gamma_err_arr.append(  np.sqrt((np.sqrt(count)*fudge_fact/(count*(1-avg_kappa)))**2 + 0.005**2))

gamma_err_arr = np.array(gamma_err_arr)
gamma_t_arr = np.array(gamma_t_arr)
plt.loglog(distArr, gamma_t_arr, 'k.', label = 'Our work z>0.5')  

#plt.loglog(distArr, gamma_t_arr+gamma_err_arr, 'r.', label = 'Our work z>0.5')  
#plt.loglog(distArr, gamma_t_arr-gamma_err_arr, 'b.', label = 'Our work z>0.5')  



#Plot okabe 
f=open('/home/dutta26/okabe_2390.txt')
content = f.readlines()
f.close()

okabe_shear =[]
okabe_radius =[]
for j in range(len(content)):
    temp = content[j].split()
    okabe_shear.append((float(temp[1]) ))
    okabe_radius.append((float(temp[0])))
    
    
plt.loglog(okabe_radius, okabe_shear, 'r+', label='Okabe Et al')  
plt.ylabel('Avg Tangential shear')
plt.xlabel('Radius in arcmin')















# =============================================================================
# def get_angularDistance(z):
#     arr= np.arange(0, z, 0.01)
#     #a= np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))
#     #return a
#     return np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))*(1/(1+z))*4285.71
# 
# def get_comvDistance(z):
#     arr= np.arange(0, z, 0.01)
#     #a= np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))
#     #return a
#     return np.sum( 0.01/np.sqrt(0.3*(1+arr)**3 + 0.7))*4285.71
# 
# 
# def get_angualrDistDiff(D1, D2, z1, z2):
#     #temp = (D2 - D1)/(1+z2)
#     return (D2 - D1)/(1+z2)
# 
# 
# a=[]
# b=[]
# c=[]
# for k in np.arange(0.4, 2, 0.1):
#     width = 7000
#     cond = np.where((master_frame[:,0] > x_mid-width) & (master_frame[:,0] < x_mid+width) & 
#                     (master_frame[:,1] > y_mid-width) & (master_frame[:,1] < y_mid+width) 
#                     & (redShiftArr>k )& (redShiftArr<k+0.1) & (master_frame[:,2] !=0) &
#                     (ir_coadd_data[:,2] == 0) & (master_frame[:,6] < 64) & 
#                      (master_frame[:,7] < 64)) [0]
#     print (len(cond))
#     
#     
#     temp = np.copy(master_frame[cond,:])
#     loc = np.where(temp[:,15]< median_ellip_err/3)
#     temp[loc,15] = median_ellip_err/3
#     wt_ellip_err = 1/temp[:,15]**2
#     
#     
#     
#     zl = 0.23
#     zs = (2*k+0.1)/2
#     dl = get_angularDistance(zl)
#     ds = get_angularDistance(zs)
#     d_comv_l = get_comvDistance(zl)
#     d_comv_s = get_comvDistance(zs)
#     dls = get_angualrDistDiff(d_comv_l, d_comv_s, zl, zs)
#     #print (dl, ds, dls)
#     a.append(dls/ds)
#     b.append(np.sum(wt_ellip_err))
#     c.append(zs)
#     #print (len(wt_ellip_err))
#     
# a=np.array(a)
# b=np.array(b)
# b = b/np.sum(b)
# print (np.sum(a*b))
# 
# 
# 
# #Fit curve 
# x = np.log10(distArr)
# y = np.log10(gamma_t_arr)
# y_err = np.log10(gamma_err_arr)
# param = np.polyfit(x,y,2, w = 1/y_err)
# 
# x_range = np.arange(-2, 2, 0.01)
# y_pred = 10** (param[0]*x_range**2 + param[1]*x_range+ param[2])
# plt.plot(10**x_range, y_pred, 'k-')
# 
# 
# term_a = 0
# term_b = 0
# theta_m_arr = np.array([1,2,3,4,5,6,7])
# theta_1 = 20
# theta_2 = 25
# for theta_m in theta_m_arr:
#     d_slice1 = np.arange(theta_m, theta_1, 0.01)
#     d_slice1_log = np.log10(d_slice1)
#     d_slice2 = np.arange(theta_1, theta_2, 0.01)
#     d_slice2_log = np.log10(d_slice2)
#     term1 = 2*np.sum( (0.01/d_slice1)*  10** (param[0]*d_slice1_log**2 + param[1]*d_slice1_log+ param[2]))
#     term2 = 2*np.sum( (0.01/d_slice2)*  10** (param[0]*d_slice2_log**2 + param[1]*d_slice2_log+ param[2]))/ (1- theta_1**2/theta_2**2)
#     print (term1, term2, (term1-term2))
# 
# 
# =============================================================================

# =============================================================================
# f=fits.open('/scratch/bell/dutta26/abell_2390/kappa_sf.fits')
# kappa_map= f[0].data
# f.close()
# 
# 
# 
# z_min = 0.7
# z_max = 3
# zFile = '/home/dutta26/photz_eazy.zout'
# if (zFile == None):
#     redShiftArr = 9
# else:
#     #Read Redshifts
#     f=open(zFile)
#     content = f.readlines()
#     f.close()
#     redShiftArr=[]
#     for j in range(len(content)):
#         
#         if (content[j][0] == '#'):
#             continue
#         else:
#             if('eazy' in zFile):
#                 if(float((content[j].split())[8]) >= 0.8):
#                     redShiftArr.append(float((content[j].split())[7]))
#                 else:
#                     redShiftArr.append(0)
#             else:
#                 redShiftArr.append(float((content[j].split())[1]))
# redShiftArr = np.array(redShiftArr)      
# master_frame = np.load('/scratch/bell/dutta26/abell_2390/master_arr_sf.npy')
# ir_coadd_data = np.load('/scratch/bell/dutta26/abell_2390/abell_ir_coadd.npy')
# loc =np.where((master_frame[:,15] >0) & (master_frame[:,2] !=0) &
#                (ir_coadd_data[:,2] == 0))[0] 
# print(len(loc), "median determinin length")
# median_ellip_err = np.median(master_frame[loc,15])
# 
# x_mid = 11600
# y_mid = 17850
# gamma_t_arr =[]
# distArr=[]
# dist_min = 0
# dist_max = 200
# x_mid_kappa = int(round(x_mid/50))
# y_mid_kappa = int(round(y_mid/50))
# kappa_dist_arr =np.zeros(np.shape(kappa_map))       
# y,x = np.shape(kappa_map)
# for j in range(x):
#     for k in range(y):
#         dist = np.sqrt( (j-x_mid_kappa)**2 + (k-y_mid_kappa)**2)
#         kappa_dist_arr[k,j] = dist
# for j in range(25):
#     #dist_min = dist_min+ 200
#     dist_max = dist_max + 200
#     width = 8000
#     cond = np.where((master_frame[:,0] > x_mid-width) & (master_frame[:,0] < x_mid+width) & 
#                     (master_frame[:,1] > y_mid-width) & (master_frame[:,1] < y_mid+width) 
#                     & (redShiftArr>z_min )& (redShiftArr<z_max) & (master_frame[:,2] !=0) &
#                     (ir_coadd_data[:,2] == 0) & (master_frame[:,6] < 64) & 
#                      (master_frame[:,7] < 64)) [0]
#     temp = np.copy(master_frame[cond,:])
#     dx = temp[:,0] -x_mid
#     dy = temp[:,1] -y_mid
#     r2 = dx**2 + dy**2
#     cos2phi = (dx**2-dy**2)/r2
#     sin2phi = 2.0*dx*dy/r2
#     
#     
#     e1 = temp[:,2]
#     e2 = temp[:,3]
#     tot_ellip = np.sqrt(e1**2 + e2**2)
#     epar= - (e1*cos2phi+e2*sin2phi)
#     eper= (e2*cos2phi-e1*sin2phi)
#     #goodEllipIndices = np.where((np.abs(e1)<0.8) & (np.abs(e2)<0.8) & (r2<(3*alphax)**2) & (r2>100))
#     goodEllipIndices = np.where( (r2<(dist_max)**2) & (r2>(dist_min**2)) & (tot_ellip > 0.0) )[0]
#     
#     
#     loc = np.where(temp[:,15]< median_ellip_err/3)
#     temp[loc,15] = median_ellip_err/3
#     wt_ellip_err = 1/temp[:,15]**2
#     
#     
#     wt_tild = wt_ellip_err[goodEllipIndices]**0.5/np.sum(wt_ellip_err[goodEllipIndices]**0.5)
#     fudge_fact = 1 /(2*(1-np.sum(tot_ellip[goodEllipIndices]**2 * wt_tild)/np.sum(wt_tild) ))
#     
#     #Process kappa
#     radius = dist_max/50
#     loc = np.where(kappa_dist_arr<radius)
#     avg_kappa = np.mean(kappa_map[loc[0], loc[1]])
#     print (avg_kappa, 'Avg kappa')
#     
#     
#     e1sum = np.sum(epar[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
#     e2sum = np.sum(eper[goodEllipIndices]* wt_ellip_err[goodEllipIndices])
#     count = np.sum(wt_ellip_err[goodEllipIndices])
#     print(len(goodEllipIndices))
#     gamma_t_arr.append(e1sum*fudge_fact/(count*(1-avg_kappa)))
#     distArr.append(dist_max*0.11/60)
# plt.loglog(distArr, gamma_t_arr, 'b.', label ='Our Work z>0.4')  
# plt.legend()
# 
# =============================================================================
