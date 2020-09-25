#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 15:12:54 2019

@author: altsai
"""

import os
import sys
import shutil
#import re
import numpy as np
from astropy.io import fits
#from astropy import wcs
#from astropy.wcs import WCS
from astropy.io import ascii
import csv
#import pandas as pd
#import pyfits
#import matplotlib.pyplot as plt
#import scipy.ndimage as snd
#import glob
#import subprocess
#from scipy import interpolate
#from scipy import stats
#from scipy.interpolate import griddata
#from time import gmtime, strftime
#import pandas as pd
#from datetime import datetime
#import asciitable


#now=str(datetime.now())

#dir_root='/home/altsai/project/20190801.NCU.Prof.Chen/data/eden/data/'
#os.chdir(dir_root)

#dir_root='/home/altsai/project/20190801.NCU.Prof.Chen/data/eden/data/'
#os.chdir(dir_root)


#folder='LOT20200916'
folder=sys.argv[1]
#folder=input("Enter Folder (ex: LOT2019xxxx): ")
#idx_obj=int(19)
idx_obj=int(sys.argv[2])-1

dir_mod=folder+'_mod'
logfile=dir_mod+'.log'
#sys.stdout=open(logfile,'w')

if os.path.exists(dir_mod):
    shutil.rmtree(dir_mod)
os.makedirs(dir_mod,exist_ok=True)


if os.path.exists(logfile):
    os.remove(logfile)
    
#csv.register_dialect('myDialect',delimiter=",")

#orig_stdout=sys.stdout
#info=csv.write(sys.stdout)
f=open(logfile,'w')
#sys.stdout=f


print('-----------------')
print('input parameters')
print('-----------------')
f.write('-----------------\n')
f.write('input parameters\n')
f.write('-----------------\n')




print(sys.argv)
f.write(str(sys.argv)+'\n')




print("Which Object you are going to process? ")
#filename='/home/altsai/project/20190801.NCU.Prof.Chen/data/eden/data/object_list.txt'
filename='object.list'
file_info=ascii.read(filename)
#print(file_info)
head_column=file_info.colnames
print(head_column)
#sys.exit(0)

column_html=file_info['html']
column_filter=file_info['filter']
column_ra=file_info['RA_hhmmss']
column_dec=file_info['DEC_ddmmss']
column_idx=file_info['Index']
column_obj=file_info['Object']


print('')
show_info=file_info['Index','Object','RA_hhmmss','DEC_ddmmss']
print(show_info)
#idx_obj=int(input("Enter Index of Object: "))-1
obj_name=column_obj[idx_obj]
obj_ra=column_ra[idx_obj]
obj_dec=column_dec[idx_obj]
filter_name=column_filter[idx_obj]
info_choice='Your will modify header for: '+str(obj_name)+' [RA] '+str(obj_ra)+' [Dec] '+str(obj_dec)+' [Filter] '+str(filter_name)
print(info_choice)
print('---------------')
f.write(info_choice+'\n')
f.write('---------------\n')

#sys.exit()

print('folder:',folder)
print('will update fits in:',dir_mod)
print('index number in',filename,': ', idx_obj)
f.write('folder: '+str(folder)+'\n')
f.write('will update fits in: '+str(dir_mod)+'\n')
#f.write('index number in '+str(filename)+': '+str(idx_obj)+'\n')
#f.write('will update header of: '+str(obj_name)+'\n')





'''
nth_line=np.int(idx_obj)+1

with open(filename,'r') as f:
    obj_info=f.readlines()[nth_line:nth_line+1]
#    print(obj_info)
   
n_tab=5
#print(n_tab)

obj=obj_info[0].split('\t',n_tab)[1]
ra=obj_info[0].split('\t',n_tab)[2]
dec=obj_info[0].split('\t',n_tab)[3]
filter_name=obj_info[0].split('\t',n_tab)[4]
print('[OBJ]',obj,'[RA]',ra,'[DEC]',dec,'[FIL]',filter_name)
print()
'''
#sys.exit(0)




print(' ---------------------------')
print(' Exame Fits Header for Bias')
print(' ---------------------------')
f.write(' ---------------------------\n')
f.write(' Exame Fits Header for Bias\n')
f.write(' ---------------------------\n')



cmd_search_file_biasdark="find ./"+folder+"/bias-dark/ | grep 'fits\|fts'|sort "
print(cmd_search_file_biasdark)
f.write(str(cmd_search_file_biasdark)+'\n')
list_file_biasdark=os.popen(cmd_search_file_biasdark,"r").read().splitlines()
#print(list_file_biasdark)
#print('found',len(list_file_biasdark),'Bias')

n_biasdark=len(list_file_biasdark)
bias_info='found '+str(n_biasdark)+' Bias/Dark'
print(bias_info)
f.write(bias_info+'\n')
                      


if n_biasdark > 0:
    k=0
    for i in list_file_biasdark:
        k=k+1
        #    print(k)
        filename_biasdark=[i.split('/',-1)[-1]][0]
#        print(filename_biasdark)

        hdu=fits.open(i)[0]
        imhead=hdu.header
        imdata=hdu.data
        exptime=imhead['EXPTIME']
        #    filter_name=imhead['FILTER']
        imgtype=imhead['IMAGETYP']
        if exptime==0:
            if imgtype=='BIAS':
                info_bias=str(k)+'/'+str(n_biasdark)+' '+str(filename_biasdark)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)
                print(info_bias)
                f.write(info_bias+'\n')
            else:
                info_bias=str(k)+'/'+str(n_biasdark)+' '+str(filename_biasdark)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)+' NOT Bias'
                print(info_bias)
                f.write(info_bias+'\n')
        else:
            if imgtype=='DARK':
                info_dark=str(k)+'/'+str(n_biasdark)+' '+str(filename_biasdark)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)
                print(info_dark)
                f.write(info_dark+'\n')
            else:
                info_dark=str(k)+'/'+str(n_biasdark)+' '+str(filename_biasdark)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)+' NOT Dark'
                print(info_dark)
                f.write(info_dark+'\n')
#        time_calib=str(datetime.now())  
#        imhead.add_history('Header corrected at '+time_calib+' UTC+8 by An-Li Tsai')
        #    fitsname_update=sci_name+'_calib.fits'
        #hdu=fits.PrimaryHDU(calib_sci[i])
        fits.writeto(dir_mod+'/'+filename_biasdark,data=imdata,header=imhead,overwrite=True)
else:
    info_biasdark="no files in ./bias-dark/"
    print(info_biasdark)
    f.write(info_biasdark+'\n')
        
#sys.exit(0)

print(' ---------------------------')
print(' Exame Fits Header for Flat')
print(' ---------------------------')
f.write(' ---------------------------\n')
f.write(' Exame Fits Header for Flat\n')
f.write(' ---------------------------\n')

    
#cmd_search_file_flat="find ./"+folder+"/ | grep flat/ |grep eden| grep -E 'fits|fts'|sort -t'-' -k2.4,2.6 "
cmd_search_file_flat="find ./"+folder+"/ | grep flat/ | grep -E 'fits|fts'|sort -t'-' -k2.4,2.6 "
print(cmd_search_file_flat)
f.write(str(cmd_search_file_flat)+'\n')
list_file_flat=os.popen(cmd_search_file_flat,"r").read().splitlines()
#print(list_file_flat)
n_flat=len(list_file_flat)
info_flat='found '+str(n_flat)+' Flat'
print(info_flat)
f.write(info_flat+'\n')





if len(list_file_flat) > 0:
#    filename_flat=[i.split('/',-1)[-1] for i in list_file_flat]
#    print(filename_flat)
    k=0
    for i in list_file_flat:
        k=k+1
        #    print(k)
        filename_flat=[i.split('/',-1)[-1]][0]
#        print(filename_flat)
#        flat_i1=[filename_flat.split('-',1)[1]][0]
#        print('flat_i1 = ',flat_i1)
#        flat_i2=[flat_i1.split('.',1)[0]][0]
#        print('flat_i2 = ',flat_i2)
#        flat_i3=flat_i2[3:]
#        print('flat_i3 = ',flat_i3)
        hdu=fits.open(i)[0]
        imhead=hdu.header
        imdata=hdu.data
        exptime=imhead['EXPTIME']       
        #    filter_name=imhead['FILTER']
        imgtype=imhead['IMAGETYP']
#        print(k,filename_flat,'[TYPE]',imgtype,'[EXPTIME]',exptime)
        if imgtype=='FLAT':
            if imhead['FILTER']=='Eden':
               imhead['FILTER']=filter_name
            flat_filter=imhead['FILTER']
            info_flat=str(k)+'/'+str(n_flat)+' '+str(filename_flat)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)+' [FIL] '+flat_filter
            print(info_flat)
            f.write(info_flat+'\n')
        else:
            info_flat=str(k)+'/'+str(n_flat)+' '+str(filename_flat)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)
            print(info_flat)
            f.write(info_flat+'\n')
#        time_calib=str(datetime.now())  
#        imhead.add_history('Header corrected at '+time_calib+' UTC+8 by An-Li Tsai')
        #    fitsname_update=sci_name+'_calib.fits'
        #hdu=fits.PrimaryHDU(calib_sci[i])
        fits.writeto(dir_mod+'/'+filename_flat,data=imdata,header=imhead,overwrite=True)

else:
    info_flat="no files in ./flat/"
    print(info_flat)
    f.write(info_flat+'\n')



print(' ---------------------------')
print(' Exame Fits Header for Additional Dark for Flat')
print(' ---------------------------')
f.write(' ---------------------------\n')
f.write(' Exame Fits Header for Additional Dark for Flat\n')
f.write(' ---------------------------\n')


cmd_search_file_darkflat="find ./"+folder+"/ | grep flat/ |grep dark| grep -E 'fits|fts'|sort -t'-' -k2.4,2.6 "
print(cmd_search_file_darkflat)
f.write(str(cmd_search_file_darkflat)+'\n')
list_file_darkflat=os.popen(cmd_search_file_darkflat,"r").read().splitlines()
#print(list_file_darkflat)
#f.write(str(list_file_darkflat)+'\n')
n_darkflat=len(list_file_darkflat)
info_darkflat='found '+str(n_darkflat)+' Additional Dark for Flat'
print(info_darkflat)
f.write(info_darkflat+'\n')


if len(list_file_darkflat) > 0:
    k=0
    for i in list_file_darkflat:
        k=k+1
        #    print(k)
        filename_flat=[i.split('/',-1)[-1]][0]
        hdu=fits.open(i)[0]
        imhead=hdu.header
        imdata=hdu.data
        exptime=imhead['EXPTIME']
        #    filter_name=imhead['FILTER']
        imgtype=imhead['IMAGETYP']
        info_darkflat=str(k)+'/'+str(n_flat)+' '+str(filename_flat)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)+' [FIL] '+flat_filter
        print(info_darkflat)
        f.write(info_darkflat+'\n')
#        print(k,filename_flat,'[TYPE]',imgtype,'[EXPTIME]',exptime)
#        if imgtype=='FLAT':
#            if imhead['FILTER']=='Eden':
#               imhead['FILTER']=filter_name
#            flat_filter=imhead['FILTER']
#            info_flat=str(k)+'/'+str(n_flat)+' '+str(filename_flat)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)+' [FIL] '+flat_filter
#            print(info_flat)
#            f.write(info_flat+'\n')
#        else:
#            info_flat=str(k)+'/'+str(n_flat)+' '+str(filename_flat)+' [TYPE] '+str(imgtype)+' [EXPTIME] '+str(exptime)
#            print(info_flat)
#            f.write(info_flat+'\n')
#        time_calib=str(datetime.now())  
#        imhead.add_history('Header corrected at '+time_calib+' UTC+8 by An-Li Tsai')
        #    fitsname_update=sci_name+'_calib.fits'
        #hdu=fits.PrimaryHDU(calib_sci[i])
#        fits.writeto(dir_mod+'/'+filename_flat,data=imdata,header=imhead,overwrite=True)

else:
    info_darkflat="no additional dark for flat in ./flat/"
    print(info_darkflat)
    f.write(info_darkflat+'\n')



#sys.exit(0)


print(' ---------------------------')
print(' Modify filename for Science Target ')
print(' ---------------------------')
f.write(' ---------------------------\n')
f.write(' Modify filename for Science Target\n')
f.write(' ---------------------------\n')

cmd_search_file_sci_modname="find ./"+folder+"/wchen/ | cut -d / -f1-4| egrep 'fits|fts'  |sort "
print(cmd_search_file_sci_modname)
#f.write(cmd_search_file_sci_modname+'\n')
list_file_sci_modname=os.popen(cmd_search_file_sci_modname,"r").read().splitlines()
print('find files:',list_file_sci_modname)
n_target_modname=len(list_file_sci_modname)
info_target_modname='found '+str(n_target_modname)+' Target\n'
print(info_target_modname)
#f.write(info_target_modname+'\n')

list_filename_sci_modname=[i.split('/',-1)[-1] for i in list_file_sci_modname]
print(list_filename_sci_modname)


#sys.exit(0)


info_update='... update header ...'
print(info_update)
f.write(info_update+'\n')

k=0
for i in list_file_sci_modname:
#    print(i)
    k=k+1
    #filename=i.split('/',-1)[-1]
#    print(k)
    hdu=fits.open(i)[0]
    imhead=hdu.header
    imdata=hdu.data
#    print(imdata.shape)
    exptime=imhead['EXPTIME']
    camera=imhead['CAMERA']
    date=imhead['DATE-OBS']
    time=imhead['TIME-OBS']
    yyyymmdd=date[0:4]+date[5:7]+date[8:10]
#    print(date,"yyyymmdd",yyyymmdd)
    hhmmss=time[0:2]+time[3:5]+time[6:8]
#    print(time,"hhmmss",hhmmss)
    filename=obj_name+'-'+yyyymmdd+'@'+hhmmss+'-Eden.fts'
#    print(filename)
#    idx_time=str(int(exptime))+'S'
#    print(idx_time)
#    print(exptime)
#    naxis=imhead['NAXIS']
#    print(naxis)
#    jd=imhead['JD']
#    obj=imhead['OBJECT']
#    fwhm=imhead['FWHM']
#    zmag=imhead['ZMAG']
#    ra=imhead['RA']
#    dec=imhead['Dec']
#    filter_name=imhead['FILTER']
    imhead.set('FILTER',filter_name)
    imgtype=imhead['IMAGETYP']
#    print('[TYPE] ',imgtype,'[OBJECT] ',obj, '[RA] ',ra, '[DEC] ',dec)    
#    print('[TYPE]',imgtype,'[OBJECT]',obj)    
    imhead.set('RA',obj_ra)
    imhead.set('Dec',obj_dec)
    imhead.set('OBJCTRA',obj_ra)
    imhead.set('OBJCTDec',obj_dec)
#    imhead.set('OBJECT',sci_name)
    imhead.set('OBJECT',obj_name)
#    print(k,'[TYPE]',imgtype,'[OBJ]',obj_name,'[RA]',obj_ra,'[Dec]', obj_dec,'[FIL]',filter_name,'[EXPTIME]',exptime, '[Date]',date,'[Time]', time)   
    info_target_modname=str(k)+'/'+str(n_target_modname)+' '+str(filename)+' [TYPE] '+str(imgtype)+' [OBJ] '+str(obj_name)+' [RA] '+str(obj_ra)+' [Dec] '+str(obj_dec)+' [FIL] '+str(filter_name)+' [EXPTIME] '+str(exptime)+' [Date] '+str(date)+' [Time] '+str(time)
    print(info_target_modname)
    f.write(info_target_modname+'\n')
#    spamwriter=csv.write(f,delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
#    spamwriter.writerow(k,'[TYPE]',imgtype,'[OBJ]',obj_name,'[RA]',obj_ra,'[Dec]', obj_dec,'[FIL]',filter_name,'[EXPTIME]',exptime, '[Date]',date,'[Time]', time)    
#    row=print(k,'[TYPE]',imgtype,'[OBJ]',obj_name,'[RA]',obj_ra,'[Dec]', obj_dec,'[FIL]',filter_name,'[EXPTIME]',exptime, '[Date]',date,'[Time]', time)
#    f.write(row)
#    time_calib=str(datetime.now())  
#    imhead.add_history('RA, Dec, Object added by An-Li Tsai at '+time_calib+' UTC+8')
#    print('filename',filename)
#    filename_root=filename.split('.',1)[0]
#    print(filename_root)
#    fitsname_update=filename_root+'_update.fits'
#    print(fitsname_update)
#    hdu=fits.PrimaryHDU(calib_sci[i])
    fits.writeto(dir_mod+'/'+filename,data=imdata,header=imhead,overwrite=True)

#sys.stdout=orig_stdout
#csvfile.close()

#info_move='... move updated fits files to '+dir_mod+'/ ...'
#print(info_move)
#f.write(info_move+'\n')
info_write='... write log file to '+dir_mod+'.log ...'
print(info_write)
f.write(info_write+'\n')

cmd_search_updated_fits="ls |grep 'fts\|fits'"
list_updated_fits=os.popen(cmd_search_updated_fits,"r").read().splitlines()
#print(list_updated_fits)
#print(len(list_updated_fits))



del list_file_biasdark
del list_file_flat
#del list_file_sci
#del list_sci_1
#del list_sci_2
#del list_filename_sci
#del list_sci_name
del list_filename_sci_modname

info_finish='... finish ...'
print(info_finish)
f.write(info_finish+'\n')

f.close()
