#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 17:56:03 2022

@author: kddazac
"""
################################################################################################################################################
#Libraries
################################################################################################################################################

import numpy as np
from math import log, sin, cos, pow
import math
import astropy
from astropy.io import ascii
from astropy.table import Table, Column
import sys
import re
import csv
import matplotlib
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import collections
from collections import OrderedDict
global debug
debug = 0
import os
import string
import sys
from astropy.io import ascii
from scipy import ndimage
import matplotlib.path as mplPath
#import pyfits
import scipy
from scipy import stats
from scipy.stats import norm as norm
from matplotlib import mlab
#import ephem
from matplotlib.ticker import NullFormatter
import matplotlib.lines as mlines
from scipy.stats import norm as norm
from mpl_toolkits.mplot3d import Axes3D
import random
import datetime
from datetime import date, timedelta
#%matplotlib qt
################################################################################################################################################
#Functions
################################################################################################################################################
def read_file(filename):
    with open(filename, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]
    return data

def extract_data(file_nmnsh):
    exp=[]
    data_dummy=read_file(file_nmnsh)
    data_dummy.remove(data_dummy[0])
    data_dummy.remove(data_dummy[0])  
    for i in range(len(data_dummy)):
        tp=data_dummy[i][1]
        tf=data_dummy[i][2]
        rh=data_dummy[i][3]
        if tp!='' and tp!='' and rh!='':
            my_date = datetime.datetime.strptime(tp, format)
            my_tf=float(tf)
            my_rh=float(rh)
            exp.append([my_date,my_tf,my_rh])
    return exp

def extract_dtr(data_room):
    dtt=[]
    dttf=[]
    dttc=[]
    dtrh=[]
    for i in range(len(data_room)):
        dtt.append(data_room[i][0])
        dttf.append(data_room[i][1])        
        dttc.append((data_room[i][1]-32.)*(5.0/9.0))  
        dtrh.append(data_room[i][2])
    return dtt,dttf,dttc,dtrh

def extract_all_weather(path_files2,filess_weather):
    timess=[]
    tempe_c=[]
    tempe_f=[]
    dew=[]
    humidi=[]
    wind=[]
    pressure=[]
    precipi=[]
    for i in range(len(filess_weather)):
        TT,tt,dd,hu,ws,ps,pr=extract_weather_data(path_files2,filess_weather[i])        
        timess.append(TT)
        tempe_f.append(extract_central_avg(tt))
        tempe_c.append(transform_to_C(extract_central_avg(tt)))
        dew.append(extract_central_avg(dd))
        humidi.append(extract_central_avg(hu))
        wind.append(extract_central_avg(ws))
        pressure.append(extract_central_avg(ps))
        precipi.append(pr)
    return timess,tempe_c,tempe_f,dew,humidi,wind,pressure,precipi

def extract_weather_data(path_files2,file_weather):
    data=[]
    reader = open(path_files2+file_weather, 'rU')
    for row in reader:
        data.append(row)
    reader.close()
    indices = [i for i, x in enumerate(data) if x == 'Max\tAvg\tMin\n']
    last = [i for i, x in enumerate(data) if x == 'Total\n']        
    start_date = date(int(file_weather[0:4]),int(file_weather[-6:][0:2]),1)
    if int(file_weather[-6:][0:2])==12:
        end_date = date(int(file_weather[0:4])+1,1,1)
    else:
        end_date = date(int(file_weather[0:4]),int(file_weather[-6:][0:2])+1,1)
    TT=[]
    for single_date in daterange(start_date, end_date):
        TT.append(single_date)
    new_TT=[]
    for i in range(len(TT)):
        aa=TT[i]
        new_TT.append(datetime.datetime(aa.year,aa.month,aa.day))        
    tt=extract_measurement(data,indices[0]+1,indices[1])
    dd=extract_measurement(data,indices[1]+1,indices[2])
    hu=extract_measurement(data,indices[2]+1,indices[3])
    ws=extract_measurement(data,indices[3]+1,indices[4])
    ps=extract_measurement(data,indices[4]+1,last[0])
    pr=extract_measurement2(data,last[0]+1,len(data))
    return new_TT,tt,dd,hu,ws,ps,pr

def extract_central_avg(tt):
    new_tt=[]
    for i in range(len(tt)):
        new_tt.append(tt[i][1])
    return new_tt      

def transform_to_C(Temperature_F):
    Temperature_C=[]
    for i in range(len(Temperature_F)):
        Temperature_C.append((Temperature_F[i]-32.)*(5./9.))
    return Temperature_C

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
def extract_measurement(data,indi_0,indi_1):
    t_indi=np.arange(indi_0,indi_1,1)
    tt=[]
    for i in t_indi:
        mi,av,ma=split_points(data[i])
        tt.append([mi,av,ma])
    return tt

def split_points(pointt):
    indices = [i for i, x in enumerate(pointt) if x == '\t']
    _1=float(pointt[0:indices[0]])
    _2=float(pointt[indices[0]+1:indices[1]])
    _3=float(pointt[indices[1]+1:-1])   
    return _1,_2,_3

def extract_measurement2(data,indi_0,indi_1):
    t_indi=np.arange(indi_0,indi_1,1)
    tt=[]
    for i in t_indi:
        tt.append(float(data[i][:-1]))
    return tt


################################################################################################################################################
#Importing files
################################################################################################################################################
folder_nmnsh='/Users/lquiroganunez/Documents/Familia/Katherin/northumbria/tesis/NMNSH_data/last/90622/'
plot_fol='/Users/lquiroganunez/Documents/Familia/Katherin/northumbria/tesis/NMNSH_data/plots/'
files_nmnsh=os.listdir(folder_nmnsh) 
#files_nmnsh.remove('.DS_Store')
################################################################################################################################################
#Organizing data per room
################################################################################################################################################
rooms=[]
for i in range(len(files_nmnsh)):
    aa=files_nmnsh[i]
    bb=aa.index('-22.')
    rooms.append(aa[bb+4:-4])


data_rooms={}    
format = '%m/%d/%y %I:%M:%S %p'   
for i in range(len(rooms)):
    data_rooms[rooms[i]]=extract_data(folder_nmnsh+files_nmnsh[i])



dt_TempGallery,dtf_TempGallery,dtc_TempGallery,drh_TempGallery=extract_dtr(data_rooms[rooms[7]])
dt_SouthGallery,dtf_SouthGallery,dtc_SouthGallery,drh_SouthGallery=extract_dtr(data_rooms[rooms[1]])
dt_NuclearMedicine,dtf_NuclearMedicine,dtc_NuclearMedicine,drh_NuclearMedicine=extract_dtr(data_rooms[rooms[4]])
dt_CollectionRoomLeft,dtf_CollectionRoomLeft,dtc_CollectionRoomLeft,drh_CollectionRoomLeft=extract_dtr(data_rooms[rooms[0]])    
dt_ColdWar,dtf_ColdWar,dtc_ColdWar,drh_ColdWar=extract_dtr(data_rooms[rooms[6]])
dt_UraniumGallery,dtf_UraniumGallery,dtc_UraniumGallery,drh_UraniumGallery=extract_dtr(data_rooms[rooms[5]])
dt_TrinityGallery,dtf_TrinityGallery,dtc_TrinityGallery,drh_TrinityGallery=extract_dtr(data_rooms[rooms[3]])
dt_CollectionRoomRight,dtf_CollectionRoomRight,dtc_CollectionRoomRight,drh_CollectionRoomRight=extract_dtr(data_rooms[rooms[2]])

################################################################################################################################################
#Weather
################################################################################################################################################
weather_path='/Users/lquiroganunez/Documents/Familia/Katherin/northumbria/tesis/NMNSH_data/last/Weather/'
filess_weather=os.listdir(weather_path)
filess_weather.sort()
if filess_weather[0]=='.DS_Store':
    filess_weather.remove('.DS_Store')
################################################################################################################################################
#Plotting
################################################################################################################################################
def make_plot_T_RH(weather_path,filess_weather,plot_fol,plot_save='no'):
    '''
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()        
    '''    
    timess_ABQ,tempe_c_ABQ,tempe_f_ABQ,dew_ABQ,humidi_ABQ,wind_ABQ,pressure_ABQ,precipi_ABQ=extract_all_weather(weather_path,filess_weather)        

    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    #ax_c = ax.twinx()
    #ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    #for i in range(len(datetime_2010)):
                 
        

        
    ax.scatter(dt_TempGallery,dtc_TempGallery,color='y',marker='o',label='Temp Gallery')                
    ax.scatter(dt_SouthGallery,dtf_SouthGallery,color='b',marker='^',label='South Gallery')                
    ax.scatter(dt_NuclearMedicine,dtc_NuclearMedicine,color='g',marker='+',label='Temp Gallery')                
    ax.scatter(dt_CollectionRoomLeft,dtc_CollectionRoomLeft,color='k',marker='s',label='South Gallery') 
    ax.scatter(dt_ColdWar,dtc_ColdWar,color='m',marker='P',label='Temp Gallery')                
    ax.scatter(dt_UraniumGallery,dtc_UraniumGallery,color='orange',marker='*',label='South Gallery')                
    ax.scatter(dt_TrinityGallery,dtc_TrinityGallery,color='brown',marker='X',label='Temp Gallery')                
    ax.scatter(dt_CollectionRoomRight,dtc_CollectionRoomRight,color='pink',marker='D',label='South Gallery')    
    for i in range(len(timess_ABQ)):
        if i==0:
            ax.scatter(timess_ABQ[i], tempe_c_ABQ[i],color='grey',marker='p',label='Weather at ABQ')    
        else:
            ax.scatter(timess_ABQ[i], tempe_c_ABQ[i],color='grey',marker='p')        
    
    tit='National Museum of Nuclear Science and History'
    ax.set_title(tit)
    #ax.legend(ncol=7,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    #ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), visible=False)
        
    ax1=fig.add_subplot(212,sharex=ax)
    #for i in range(len(datetime_2010)):
    
    ax1.scatter(dt_TempGallery,drh_TempGallery,color='y',marker='o',label='Temp Gallery')                
    ax1.scatter(dt_SouthGallery,drh_SouthGallery,color='b',marker='^',label='South Gallery')                
    ax1.scatter(dt_NuclearMedicine,drh_NuclearMedicine,color='g',marker='+',label='Nuclear Medicine')                
    ax1.scatter(dt_CollectionRoomLeft,drh_CollectionRoomLeft,color='k',marker='s',label='Collection Room Left') 
    ax1.scatter(dt_ColdWar,drh_ColdWar,color='m',marker='P',label='Cold War')                
    ax1.scatter(dt_UraniumGallery,drh_UraniumGallery,color='orange',marker='*',label='Uranium Gallery')                
    ax1.scatter(dt_TrinityGallery,drh_TrinityGallery,color='brown',marker='X',label='Trinity Gallery')                
    ax1.scatter(dt_CollectionRoomRight,drh_CollectionRoomRight,color='pink',marker='D',label='Collection Room Right')      #ax1.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    for i in range(len(timess_ABQ)):
        if i==0:
            ax1.scatter(timess_ABQ[i], humidi_ABQ[i],color='grey',marker='p',label='Weather at ABQ')    
        else:
            ax1.scatter(timess_ABQ[i], humidi_ABQ[i],color='grey',marker='p')
                     
    #ax1.legend(ncol=1,loc=1,handletextpad=0.1,prop={'size': 20})
    if plot_save=='yes':
        plt.savefig(plot_fol+'NMNSH_6.pdf',bbox_inches='tight')    
    return

make_plot_T_RH(weather_path,filess_weather,plot_fol)

np.median(dtc_TempGallery+dtf_SouthGallery+dtc_NuclearMedicine+dtc_CollectionRoomLeft+dtc_ColdWar+dtc_UraniumGallery+dtc_TrinityGallery+dtc_CollectionRoomRight)

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

np.std(dtc_TempGallery+dtf_SouthGallery+dtc_NuclearMedicine+dtc_CollectionRoomLeft+dtc_ColdWar+dtc_UraniumGallery+dtc_TrinityGallery+dtc_CollectionRoomRight)

