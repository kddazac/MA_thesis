#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 11:21:54 2022

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
import time
import matplotlib.path as mplPath
import datetime
from time import strptime
from datetime import date, timedelta
import math
import pandas
#%matplotlib qt
################################################################################################################################################
#Functions
################################################################################################################################################
def read_file(filename):
    with open(filename, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]
    return data

def extract_date(file_sam):
    aa=file_sam.index('T')
    bb=file_sam.index('.csv')
    ts=file_sam[aa+1:bb]
    if ts[0]==' ':
        ts=ts[1:]  
    dts=datetime.datetime.strptime(ts, "%m-%d-%Y")
    data=read_file(file_sam)
    dts2=data[-1][1]
    dts3=datetime.datetime.strptime(dts2, "%Y-%m-%d")
    dtf=[dts,dts3]
    return dtf
    
def extract_dates(files_sam):
    dates=[]
    for i in range(len(files_sam)):
        #print(i)
        dd=extract_date(files_sam[i])
        dates.append(dd)
    return dates

def extract_data(file_sam):
    data=read_file(file_sam)
    return data

def extract_datas(files_sam):
    datas=[]
    for i in range(len(files_sam)):
        #print(i)
        dd=extract_data(files_sam[i])
        datas.append(dd)
    return datas

def extract_y(data_year,column):
    aa=[]
    for i in range(len(data_year)):
        aa.append(data_year[i][column])
    return aa
def celsius2fahrenheit(temp):
    """
    Returns temperature in Fahrenheit given Celsius temperature.
    """
    return ((9./5.* temp) + 32)

def make_plot_T_RH(data_sam,year,plot_fol,plot_save='no'):
    '''
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()        
    '''    
    weather_path='/Users/kddazac/Documents/Familia/Katherin/northumbria/tesis/sam_Data/weather/'+str(year)+'/'
    filess_weather=os.listdir(weather_path)
    filess_weather.sort()
    if filess_weather[0]=='.DS_Store':
        filess_weather.remove('.DS_Store')
    timess_SEA,tempe_c_SEA,tempe_f_SEA,dew_SEA,humidi_SEA,wind_SEA,pressure_SEA,precipi_SEA=extract_all_weather(weather_path,filess_weather)    
    datetimes=extract_y(data_sam,0)
    #tf_2010=extract_y(data_sam_2010,1)
    tc=extract_y(data_sam,2)
    rh=extract_y(data_sam,3) 

    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    #ax_c = ax.twinx()
    #ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    #for i in range(len(datetime_2010)):
    ax.scatter(datetimes, tc,color='r',marker='o')   
    for i in range(len(timess_SEA)):
        if i==0:
            ax.scatter(timess_SEA[i], tempe_c_SEA[i],color='k',marker='x',label='Weather at SEA')    
        else:
            ax.scatter(timess_SEA[i], tempe_c_SEA[i],color='k',marker='x')             
    tit='Seattle Art Museum '+str(year)
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
    ax1.scatter(datetimes, rh,color='b',marker='o')
    for i in range(len(timess_SEA)):
        if i==0:
            ax1.scatter(timess_SEA[i], humidi_SEA[i],color='k',marker='x',label='Weather at SEA')    
        else:
            ax1.scatter(timess_SEA[i], humidi_SEA[i],color='k',marker='x')
    ax1.legend(loc=8,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    
    if plot_save=='yes':
        plt.savefig(plot_fol+'/SAM_'+str(year)+'.pdf',bbox_inches='tight')    
    return

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
folder_sam='/Users/kddazac/Documents/Familia/Katherin/northumbria/tesis/sam_Data/MA_Thesis_SAM_2'
plot_fol='/Users/kddazac/Documents/Familia/Katherin/northumbria/tesis/sam_Data/plots'
files_sam=os.listdir(folder_sam) 
files_sam.remove('.DS_Store')
################################################################################################################################################
#Clasifiying data
################################################################################################################################################
'''        
dds=extract_dates(files_sam)
dds.sort()
'''
data_sam=extract_datas(files_sam)
result = sum(data_sam, [])
datetime_sam=[]
Tf_sam=[]
Tc_sam=[]
RH_sam=[]
for i in range(len(result)):
    dtts=result[i][1]+'-'+result[i][2]
    try:
        dttsf=datetime.datetime.strptime(dtts, "%Y-%m-%d-%H:%M:%S")    
    except ValueError:
        dttsf=datetime.datetime.strptime(dtts, "%m/%d/%Y-%H:%M:%S")    
    T_f=float(result[i][3])
    T_C=(T_f-32.)/1.8    
    RH=float(result[i][4])
    datetime_sam.append(dttsf)
    Tf_sam.append(T_f)
    Tc_sam.append(T_C)
    RH_sam.append(RH)
################################################################################################################################################
#Splitting files per year
################################################################################################################################################
data_sam_2010=[]
data_sam_2011=[]
data_sam_2012=[]
data_sam_2013=[]
data_sam_2014=[]
data_sam_2015=[]
data_sam_2016=[]
data_sam_2017=[]
data_sam_2018=[]
data_sam_2019=[]
data_sam_2020=[]
data_sam_2021=[]

for i in range(len(datetime_sam)):
    if datetime_sam[i].year==2010:
        data_sam_2010.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])
    elif datetime_sam[i].year==2011:
        data_sam_2011.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])    
    elif datetime_sam[i].year==2012:
        data_sam_2012.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])
    elif datetime_sam[i].year==2013:
        data_sam_2013.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])    
    elif datetime_sam[i].year==2014:
        data_sam_2014.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]]) 
    elif datetime_sam[i].year==2015:
        data_sam_2015.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])    
    elif datetime_sam[i].year==2016:
        data_sam_2016.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])
    elif datetime_sam[i].year==2017:
        data_sam_2017.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])    
    elif datetime_sam[i].year==2018:
        data_sam_2018.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])  
    elif datetime_sam[i].year==2019:
        data_sam_2019.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])
    elif datetime_sam[i].year==2020:
        data_sam_2020.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])    
    elif datetime_sam[i].year==2021:
        data_sam_2021.append([datetime_sam[i],Tf_sam[i],Tc_sam[i],RH_sam[i]])  
    else:
        print('Somthing is  wrong')
        
################################################################################################################################################
#Plots
################################################################################################################################################
make_plot_T_RH(data_sam_2010,2010,plot_fol,'yes')
make_plot_T_RH(data_sam_2011,2011,plot_fol,'yes')
make_plot_T_RH(data_sam_2012,2012,plot_fol,'yes')
make_plot_T_RH(data_sam_2013,2013,plot_fol,'yes')
make_plot_T_RH(data_sam_2014,2014,plot_fol,'yes')
make_plot_T_RH(data_sam_2015,2015,plot_fol,'yes')
make_plot_T_RH(data_sam_2016,2016,plot_fol,'yes')
make_plot_T_RH(data_sam_2017,2017,plot_fol,'yes')
make_plot_T_RH(data_sam_2018,2018,plot_fol,'yes')
make_plot_T_RH(data_sam_2019,2019,plot_fol,'yes')
make_plot_T_RH(data_sam_2020,2020,plot_fol,'yes')
make_plot_T_RH(data_sam_2021,2021,plot_fol,'yes')



def meadian_VAl(data_sam):
    datetimes=extract_y(data_sam,0)     
    tt=[]
    rr=[]
    for i in range(len(data_sam)):
        tc=data_sam[i][2]
        rh=data_sam[i][3]
        tt.append(tc)
        rr.append(rh)
    return tt,rr
    
    
    
t,r=meadian_VAl(data_sam_2010)   