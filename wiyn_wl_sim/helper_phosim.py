#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 19:45:28 2023

@author: dutta26
"""

import os,shutil,sys
import gzip
import subprocess
import numpy as np
import helper
lut1 = np.load('/scratch/bell/dutta26/abell_2390/calib_final.npy')
lut2 = np.load('/scratch/bell/dutta26/abell_2390/calib_final_inf.npy')
from astropy.io import fits
#idNo = sys.argv[1]
#filt = sys.argv[2]
mPhaseArr =[33.35772459, 26.7205951 , 25.77512806, 30.05350671, 31.78522048,
       31.09810409, 37.73820338, 20.7312538 , 35.72899021, 30.64984321,
       29.20009561, 26.51197345, 23.43860268, 21.824414  , 29.15073163,
       31.468955  , 27.89793349, 24.24843731, 40.91277897, 32.83367139,
       29.50347037, 36.71988508, 36.8128183 , 25.7248748 , 30.74386349,
       25.02254098, 29.28575491, 28.36867282, 27.09437222, 30.41092742,
       35.73001405, 38.9784023 , 26.42276918, 33.63151988, 28.18139246,
       24.14327864, 30.10117092, 36.27491867, 24.78616508, 35.29939172,
       30.58399574, 22.40408645, 33.9418769 , 29.08393139, 26.94532914,
       33.76117447, 32.44659933, 37.1253191 , 29.4366403 , 33.19092135,
       32.29346874]

seeingArr=[0.63334611, 0.71109473, 0.95335662, 0.62568536, 0.48102599,
       1.13072337, 1.13142272, 1.1198956 , 1.12947914, 1.14467267,
       0.72036593, 0.71717661, 0.65721645, 0.49186362, 0.83675912,
       1.17046014, 1.19465896, 0.97975643, 1.246128  , 1.10839956,
       0.51901844, 0.69873197, 0.88001856, 0.6122773 , 0.85245784,
       1.04564263, 1.13452122, 1.15994773, 1.14068855, 1.23156282,
       0.5886329 , 0.73767931, 0.78611733, 0.76808936, 0.61768548,
       1.15714109, 1.18778451, 1.03879792, 1.08454227, 1.23645456,
       0.55487526, 1.37729617, 0.96640203, 1.03834031, 0.85559792,
       1.4539305 , 1.43466999, 1.17326914, 1.32429733, 1.19282821]

mjdArr =[70405.39424678, 73777.23429378, 42733.22544178, 45165.29123478,
       94901.20456778, 82147.31337178, 41890.29425878, 26571.22555078,
       36713.39157078, 78048.26300978, 94206.28341578, 93324.42899278,
       49892.40648378, 95632.20294878, 79471.31019478, 48103.26044978,
       53574.37269978, 31249.43030878, 77014.36773378, 36757.34576178,
       23666.37695378, 60555.40513178, 38604.42004478, 54731.21436278,
       76587.23111478, 35448.13149478, 24332.30026078, 31697.23975278,
       82142.24519878, 84620.35122678, 82584.13725578, 55786.34567878,
       26938.42623778, 70406.39424678, 73778.23429378, 42734.22544178,
       45166.29123478, 94900.20456778, 82146.31337178, 41891.29425878,
       26570.22555078, 36714.39157078, 78049.26300978, 94205.28341578,
       93324.42899278, 26938.42623778, 70405.39424678, 73777.23429378, 
       42733.22544178, 78048.26300978]

CM0= [ 0.      , -0.071447,  0.      , -0.094951, -0.115654, -0.02237 ,
        0.      ,  0.      ,  0.      ,  0.      , -0.025684, -0.031689,
       -0.087577, -0.061784,  0.      ,  0.      ,  0.      , -0.034899,
       0.      , -0.071447,  0.      , -0.094951, -0.115654, -0.02237 ,
         0.      ,  0.      ,  0.      ,  0.      , -0.025684, -0.031689,
        -0.087577, -0.061784,  0.      ,  0.      ,  0.      , -0.034899,
        0.      , -0.071447,  0.      , -0.094951, -0.115654, -0.02237 ,
         0.      ,  0.      ,  0.      ,  0.      , -0.025684, -0.031689,
        -0.087577, -0.061784,  0.      ,  0.      ,  0.      , -0.034899]


CM1=[-0.11405166, -0.09178872, -0.05910725, -0.07848782, -0.11476005,
       -0.12028583, -0.08571826, -0.11933576, -0.12852008, -0.04785753,
       -0.09971897, -0.12034122, -0.07024934, -0.08110146, -0.13185267,
       -0.06017125, -0.1415485 , -0.11061495, -0.09363298, -0.05569869,
       -0.12225701, -0.11802213, -0.12108769, -0.12281485, -0.10097572,
       -0.06857276, -0.08678413, -0.06798907, -0.05438334, -0.0996816 ,
       -0.08222609, -0.13472133, -0.08250858, -0.09828945, -0.10185858,
       -0.13441046, -0.04564751, -0.07797755, -0.08063801, -0.0890657 ,
       -0.06782027, -0.12788098, -0.09002002, -0.15509989, -0.09570761,
       -0.10292826, -0.07416492, -0.05994316, -0.08078506, -0.03181791]



CM2=[-0.02547816, -0.07191708, -0.07348027, -0.06220854, -0.04414912,
       -0.04431081, -0.0717718 , -0.10463825, -0.04533428, -0.05071832,
       -0.06132511, -0.04594583, -0.01158064, -0.02466955, -0.0631787 ,
        0.        , -0.0778924 , -0.07580699, -0.06875635, -0.07209737,
       -0.04882157, -0.09754993, -0.02100425, -0.05032188, -0.05269102,
       -0.02577654, -0.07248159, -0.05560495, -0.00855514, -0.02956359,
       -0.09500759, -0.05154435, -0.09307183, -0.05902866, -0.05964318,
       -0.04108874, -0.0244424 , -0.05083271, -0.06400428, -0.04534899,
        0.        , -0.0253987 ,  0.        , -0.03331768, -0.09692519,
       -0.0264804 , -0.10671376, -0.06555297, -0.08106977, -0.07878393]



def clean(folder, catalog):
    
    os.remove(catalog)
    swarpLoc = '/home/dutta26/apps/bin/bin/'
    
    #Remove CELL files
    for files in os.listdir(folder):
        if('CELL' in files):
            os.remove(folder+files)
            
    #Unpack gz files
    for files in os.listdir(folder):
        if('.gz' not in files):
            continue
        
        with gzip.open(folder+files, 'rb') as f_in:
            with open(folder+files[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(folder+files)
        
    #Coadd the files to make one
    swarpfileName = folder+'temp.txt'
    f=open(swarpfileName, 'w+')
    for files in os.listdir(folder):
        if('txt' in files):
            continue
        f.write(folder+files+'\n')
    f.close()    
    swarpCommand = './swarp @'+folder+'temp.txt'+' -c /home/dutta26/codes/wiyn_wl_sim/default.swarp -IMAGEOUT_NAME '+folder+'output.fits -WEIGHTOUT_NAME '+folder+'output.weight.fits -SUBTRACT_BACK N  -WEIGHT_TYPE NONE '
    process = subprocess.Popen(swarpCommand.split(), stdout=subprocess.PIPE, cwd=swarpLoc)
    output, error = process.communicate()
    
    for files in os.listdir(folder):
        if('OTA44' in files):
            continue
        if('output' not in files):
            os.remove(folder+files)
    return 


def make_phosim_imputCatalog(idNo, filt):
    
    dict1={'u':'0', 'g':'1', 'r':'2', 'i':'3', 'z':'4'}

    f=open('/scratch/bell/dutta26/wiyn_sim/sheared_uniform.cat')
    content = f.readlines()
    f.close()
    
    uid = filt+'_'+idNo
    
    #For dithering
    ra = 328.3941 
    dec =17.6697
    delta = 1/60
    if(int(idNo)%4 == 0):
        ra += delta
        dec += delta
    elif(int(idNo)%4 == 1):
        ra -= delta
        dec += delta
    elif(int(idNo)%4 == 2):
        ra -= delta
        dec -= delta
    elif(int(idNo)%4 == 3):
        ra += delta
        dec -= delta
    
    
    f=open('/scratch/bell/dutta26/wiyn_sim/'+str(uid)+'.txt', 'w+')
    f.write('rightascension '+str(ra)+'  \n')    
    f.write('declination '+str(dec)+'  \n')
    f.write('seed '+str(idNo)+ '\n')
    f.write('obshistid '+str(idNo)+ '\n')
    f.write('vistime 60 \n')
    f.write('filter '+dict1[filt]+' \n')
    f.write('moonphase '+str(mPhaseArr[int(idNo)- 1000])+' \n')
    f.write('mjd '+str(mjdArr[int(idNo)- 1000])+' \n')
    f.write('seeing '+str(seeingArr[int(idNo)- 1000])+' \n')
    
    
    for j in range(len(content)):
        f.write(content[j])
    f.close()
    return '/scratch/bell/dutta26/wiyn_sim/'+str(uid)+'.txt'
    
    
def check_n_remove_existing(folder):
    if os.path.isdir(folder):
        for files in os.listdir(folder):
            os.remove(folder+files)
    else:
        os.mkdir(folder)
        
    return
        
        
def make_phosim_star_imputCatalog(idNo, filt):
    dict1={'u':'0', 'g':'1', 'r':'2', 'i':'3', 'z':'4'}
    cent_ra = 328.3941
    cent_dec = 17.6697
    
    raArr =np.arange(cent_ra-(4*0.00555), cent_ra+(4*0.00555), 0.00555 )
    decArr = np.arange(cent_dec -(4*0.00555) ,  cent_dec + (4*0.00555) , 0.005551)    
    
    #For dithering
    ra = 328.3941 
    dec =17.6697
    delta = 1/60
    if(int(idNo)%4 == 0):
        ra += delta
        dec += delta
    elif(int(idNo)%4 == 1):
        ra -= delta
        dec += delta
    elif(int(idNo)%4 == 2):
        ra -= delta
        dec -= delta
    elif(int(idNo)%4 == 3):
        ra += delta
        dec -= delta
    mag = 15
    count = 1
    uid = filt+'_'+idNo
    f=open('/scratch/bell/dutta26/wiyn_sim/star_'+str(uid)+'.txt', 'w+')
    f.write('rightascension '+str(ra)+'  \n')    
    f.write('declination '+str(dec)+'  \n')
    f.write('seed '+str(idNo)+ '\n')
    f.write('obshistid '+str(idNo)+ '\n')
    f.write('vistime 60.0 \n')
    f.write('filter '+dict1[filt]+' \n')
    f.write('moonphase '+str(mPhaseArr[int(idNo)- 1000])+' \n')
    f.write('mjd '+str(mjdArr[int(idNo)- 1000])+' \n')
    f.write('seeing '+str(seeingArr[int(idNo)- 1000])+' \n')
    
    for j in range(len(raArr)):
        mag = mag+1
        for k in range(len(decArr)):
            temp_str = str(count) + ' '+str(raArr[j])[0:9]+ ' '+ str(decArr[k])[0:9] + ' '+str(mag)[0:5]
            f.write('object ' +temp_str + ' ../sky/sed_flat.txt 0.0 0.0 0.0 0.0 0.0 0.0 star none none \n')
            count += 1
    f.close()
    return '/scratch/bell/dutta26/wiyn_sim/star_'+str(uid)+'.txt'


def makeWeight(idNo, filt, loc_img):
    ref_catalog_star = '/scratch/bell/dutta26/wiyn_sim/wcs_reference_star.cat'
    folder_star = loc_img+filt+'/'+str(idNo)+"_star/"
    folder = loc_img+filt+'/'+str(idNo)+"/"
    if not os.path.isdir(folder):
        print ('Could not find folder')
        return -1
    if not os.path.isdir(folder_star):
        print ('Could not find star folder')
        return -1
    
    f=fits.open(loc_img+filt+'/'+str(idNo)+"_star/output.fits")
    img = f[0].data
    f.close()
    bkg = np.median(img[img>0])
    print (bkg)
# =============================================================================
#     cent_ra = 328.3941
#     cent_dec = 17.6697
#     raArr =np.arange(cent_ra-(4*0.00555), cent_ra+(4*0.00555), 0.00555 )
#     decArr = np.arange(cent_dec -(4*0.00555) ,  cent_dec + (4*0.00555) , 0.005551)
# =============================================================================
    
    #New to take into account weird shifts for WIYN
    f=open(ref_catalog_star)
    cat = f.readlines()
    f.close()
    raArr =[]
    decArr = []
    ref_flux=[]
    magArr=[]
    mag = 23
    #Read the list of ra dec and dec into arrays
    for j in range(len(cat)):
        if((cat[j].split()[0]) == '#'):
         continue
         
        raArr.append(float(cat[j].split()[5])) 
        decArr.append(float(cat[j].split()[6]))
        magArr.append(mag)
        
        if(len(raArr)%8 ==0):
            mag -= 1
    
    
    
    fluxArr =[]
    sizeArr=[]
    
    xList, yList = helper.convertToXY(decArr, raArr, loc_img+filt+'/'+str(idNo)+"_star/output.fits") 
    for j in range(len(xList)):
        x=xList[j]
        y=yList[j]
        cut = img[int(y-25):int(y+25), int(x-25):int(x+25)]
        flux_m, mux, muy, e1, e2, bkg_m, psf_m, sigxx, sigyy, sigxy = helper.measure_new(cut, lut1, lut2)
        if(flux_m == None or psf_m == None):
            fluxArr.append(0)
            sizeArr.append(0)
        else:
            fluxArr.append(flux_m)
            sizeArr.append(psf_m)
    
    #print (magArr, sizeArr, fluxArr)
    fluxArr=np.array(fluxArr)
    sizeArr = np.array(sizeArr)
    magArr = np.array(magArr)
    
    #Find avg zp      
    loc = np.where(fluxArr>0)
    #print (fluxArr, loc, magArr)
    zp= np.median( magArr[loc] + 2.5*np.log10(fluxArr[loc]/60))
    #print  (magArr[loc] + 2.5*np.log10(fluxArr[loc]/60))
    med_size = np.median(sizeArr[loc])
    
    wt = 100*np.power(10,(zp-25)/2.5)/((med_size)**2*bkg)
    #print (wt, zp, med_size, bkg)
    #Flux scale assuming ZP= 25
    scaleArr = 10**((25 - magArr[loc] - 2.5*np.log10(fluxArr[loc]/60))/2.5)
    scale = np.median(scaleArr)
    #print (scale, scaleArr, wt)
    #Update weight
    f=fits.open(loc_img+filt+'/'+str(idNo)+"/output.weight.fits", mode ='update')
    img = f[0].data
    loc = np.where(img>0)
    f[0].data[loc] = wt
    f.flush()
    
    #Update fluxscale 
    print (loc_img+filt+'/'+str(idNo)+"/output.fits")
    f=fits.open(loc_img+filt+'/'+str(idNo)+"/output.fits", mode ='update')
    hdr = f[0].header
    hdr['FLXSCALE'] = scale
    hdr['ZP'] = zp
    hdr['SIZE'] = med_size
    hdr['BKG'] = bkg
    f.flush()
    print (wt)
    return scale, zp, med_size

