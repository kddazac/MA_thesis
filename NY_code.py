#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  2 11:54:48 2022

@author: kddazac
"""
################################################################################################################################################
#Libraries
################################################################################################################################################
import numpy as np
from math import log, sin, cos, pow
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
import collections
import time
import matplotlib.path as mplPath
import datetime
from time import strptime
from datetime import date, timedelta
import math
import pandas

#%matplotlib qt
################################################################################################################################################
origin_path=os.getcwd()
plot_fol=out_fol=origin_path+'/plots_NY/'
################################################################################################################################################
#Functions
################################################################################################################################################
def celsius2fahrenheit(temp):
    """
    Returns temperature in Fahrenheit given Celsius temperature.
    """
    return ((9./5.* temp) + 32)

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

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
def split_points(pointt):
    indices = [i for i, x in enumerate(pointt) if x == '\t']
    _1=float(pointt[0:indices[0]])
    _2=float(pointt[indices[0]+1:indices[1]])
    _3=float(pointt[indices[1]+1:-1])   
    return _1,_2,_3

def extract_measurement(data,indi_0,indi_1):
    t_indi=np.arange(indi_0,indi_1,1)
    tt=[]
    for i in t_indi:
        mi,av,ma=split_points(data[i])
        tt.append([mi,av,ma])
    return tt

def extract_measurement2(data,indi_0,indi_1):
    t_indi=np.arange(indi_0,indi_1,1)
    tt=[]
    for i in t_indi:
        tt.append(float(data[i][:-1]))
    return tt

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

def plot_CP_weather(timess_DCA,tempe_c_DCA,tempe_f_DCA,dew_DCA,humidi_DCA,wind_DCA,pressure_DCA,precipi_DCA):        
    matplotlib.rcParams.update({'font.size': 18})    
    fig = plt.figure(figsize=(18,18))

    ax7=fig.add_subplot(611)
    #ax_c7 = ax7.twinx()
    #ax7.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess_DCA)):
        ax7.scatter(timess_DCA[i],tempe_c_DCA[i],color='c',marker='s')           
    #ax.set_title('Renwick Gallery')
    ax7.grid()
    ax7.set_ylabel('Temperature \n ($^o$C)')
    ax7.set_title('Weather Station at Central Park NYC 2022')
    #ax.set_ylim([19.8,24.2])
    #ax_c7.set_ylabel('Temperature ($F$)')    
    #plt.xticks(rotation=90)
    plt.setp(ax7.get_xticklabels(), visible=False)
    
    ax8=fig.add_subplot(612,sharex=ax7)    
    #ax_c8 = ax8.twinx()
    #ax8.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()    
    dew_DCA_C=[transform_to_C(dew_DCA[i]) for i in range(len(dew_DCA))]    
    for i in range(len(timess_DCA)):
        ax8.scatter(timess_DCA[i],dew_DCA_C[i],color='m',marker='s')          
    #ax.set_title('Renwick Gallery')
    ax8.grid()
    ax8.set_ylabel('Dew \n Point \n ($^o$C)')
    #ax.set_ylim([19.8,24.2])
    #ax_c8.set_ylabel('Dew Point ($F$)')    
    #plt.xticks(rotation=90)
    plt.setp(ax8.get_xticklabels(), visible=False)
    
    ax9=fig.add_subplot(613,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax9.scatter(timess_DCA[i],humidi_DCA[i],color='g',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax9.grid()
    ax9.set_ylabel('Relative \n Humidity \n (%)')
    plt.setp(ax9.get_xticklabels(), visible=False)
    
    ax10=fig.add_subplot(614,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax10.scatter(timess_DCA[i],wind_DCA[i],color='r',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax10.grid()
    ax10.set_ylabel('Wind \n Speed \n (mph)')
    plt.setp(ax10.get_xticklabels(), visible=False)
    
    ax11=fig.add_subplot(6,1,5,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax11.scatter(timess_DCA[i],pressure_DCA[i],color='k',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax11.grid()
    ax11.set_ylabel('Pressure \n (Hg)')
    plt.setp(ax11.get_xticklabels(), visible=False)

    ax12=fig.add_subplot(6,1,6,sharex=ax7)
    for i in range(len(timess_DCA)):
        ax12.scatter(timess_DCA[i],precipi_DCA[i],color='b',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax12.grid()
    ax12.set_ylabel('Precipitation \n (in)')
    #plt.xticks(rotation=90)
    return

def extract_all_datalogs(path_files,renwick_files,time_up,time_do):
    timess=[]
    tempe_cs=[]
    tempe_fs=[]
    humidis=[]
    tts=[]
    Dataloggers=[]
    luxs=[]
    for i in range(len(renwick_files)):    
        times,tempe_c,tempe_f,humidi,tt,Datalogger,lux=extract_data_NY(path_files,renwick_files[i],time_up,time_do)
        timess.append(times)
        tempe_cs.append(tempe_c)
        tempe_fs.append(tempe_f)
        humidis.append(humidi)
        tts.append(tt)
        Dataloggers.append(Datalogger)
        luxs.append(lux)
    return timess,tempe_cs,tempe_fs,humidis,tts,Dataloggers,luxs

def extract_data_NY(path_files,renwick_file,time_up,time_do):
    print(renwick_file)
    dummy_renwick=pandas.read_csv(path_files+renwick_file,header=1)
    Time=dummy_renwick['Date Time, GMT-05:00']
    times,times2,mm=convert_times(Time)
    RHH=dummy_renwick.keys()[3]
    Humidity=dummy_renwick[RHH][0:mm]
    temperatura=dummy_renwick.keys()[2]
    Temperature_F=dummy_renwick[temperatura][0:mm]
    Temperature_C=transform_to_C(Temperature_F)
    ll=dummy_renwick.keys()[4]
    lums=dummy_renwick[ll][0:mm]
    #hhmm,datess=split_hhmm_date(Time)
    #print(renwick_file)
    timess=[]
    tempe_f=[]
    tempe_c=[]
    humidi=[]
    lux=[]
    for i in range(len(times)):
        if times[i]<time_up and times[i]>time_do :
            timess.append(times[i])
            tempe_f.append(Temperature_F[i])
            tempe_c.append(Temperature_C[i])    #hhmm_=hhmm[4:52]
            humidi.append(Humidity[i])
            lux.append(lums[i])
    tt=set_ticks_times(times2)      
    Datalogger=renwick_file[:-4]
    return timess,tempe_c,tempe_f,humidi,tt,Datalogger,lux

def convert_times(timesss):
    new_times=[]
    mm=0
    for i in range(len(timesss)):
        aa=timesss[i][-2:]
        bb=timesss[i][9:11]
        if aa=='AM':
            change=0
        elif aa=='PM' and bb=='12':
            change=0
        elif aa=='PM' and bb!='12':
            change=12
        new_times.append(datetime.datetime(int(float(timesss[i][6:8])+2000),int(timesss[i][0:2]),int(timesss[i][3:5]),int(timesss[i][9:11])+change,int(timesss[i][12:14]),int(timesss[i][15:17])))       
        mm+=1
    return new_times,timesss,mm

def set_ticks_times(times):
    tick=[]
    fecha_ref=0.0
    zz=np.arange(times.index.start,times.index.stop,times.index.step)
    for i in zz:        
        fecha=times[i][7:]
        hora=times[i][0:5]
        if fecha!=fecha_ref:                        
            tick.append(fecha+'\n'+hora)
            fecha_ref=fecha
        else:
            tick.append(hora)
    return tick

def make_plot_T_RH_2022(colorss):
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(411)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    for i in range(len(timess)):
        ax.scatter(timess[i], tempe_cs[i],color=colorss[i],marker='o',label=Dataloggers[i])           
    for i in range(len(timess_CP)):
        ax.scatter(timess_CP[i], tempe_c_CP[i],color='b',marker='x')
    ax.set_title('XXX Gallery at NY')
    #ax.legend(ncol=7,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature \n ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), visible=False)
        
    ax1=fig.add_subplot(412,sharex=ax)
    for i in range(len(timess)):
            ax1.scatter(timess[i], humidis[i],color=colorss[i],marker='o',label=Dataloggers[i])
    for i in range(len(timess_CP)):
        if i==0:
            ax1.scatter(timess_CP[i], humidi_CP[i],color='b',marker='x',label='Weather at CP')    
        else:
            ax1.scatter(timess_CP[i], humidi_CP[i],color='b',marker='x')
    ax1.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    plt.setp(ax1.get_xticklabels(), visible=False)
    
    ax2=fig.add_subplot(413,sharex=ax)
    for i in range(len(timess)):
            ax2.scatter(timess[i], luxs[i],color=colorss[i],marker='o',label=Dataloggers[i])
    #ax1.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax2.grid()
    ax2.set_ylabel('Uncalibrated \n Intensity (Lux)')
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3=fig.add_subplot(414,sharex=ax)
    for i in range(len(timess)):
            ax3.scatter(timess[i], luxs_cal[i],color=colorss[i],marker='o',label=Dataloggers[i])
    #ax1.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax3.grid()
    ax3.set_ylabel('Calibrated \n Intensity (Lux)')
    return

def calibrated_lux(luxs,Dataloggers):
    cali=[]
    for i in range(len(Dataloggers)):
        aa=equation_calib(luxs[i],Dataloggers[i])            
        cali.append(aa)
    return cali

def equation_calib(data_lux,logger):
    new_data=[]
    if logger=='DSR_4':
        #1.6228x - 1.5432
        m=1.6228
        b=-1.5432
    elif logger=='DSR_6':
        #1.4745x + 14.401
        m=1.4745
        b=14.401
    elif logger=='DSR_7':
        #1.1671x - 16.849
        m=1.1671
        b=-16.849
    elif logger=='DSR_8':
        #1.4983x - 1.503
        m=1.4983
        b=-1.503
    else:
        'SoNYhing is wrong'
    for i in range(len(data_lux)):
        new=data_lux[i]*m+b
        new_data.append(new)  
    return new_data

################################################################################################################################################
#Reading data
################################################################################################################################################
path_files='/Users/kddazac/Documents/Familia/Katherin/NY/data/'
filess_datalogs=os.listdir(path_files)
filess_datalogs.sort()
time_start_2022=datetime.datetime(2022,1,1,00,00,00)
time_end_2022=datetime.datetime(2022,12,31,23,59,59)
timess,tempe_cs,tempe_fs,humidis,tts,Dataloggers,luxs=extract_all_datalogs(path_files,filess_datalogs,time_end_2022,time_start_2022)
luxs_cal=calibrated_lux(luxs,Dataloggers)
colorss=['k','grey','g','orange','cyan', 'r']
################################################################################################################################################
#Reading Weather
################################################################################################################################################
#Weather NYC 20122
path_files='/Users/kddazac/Documents/Familia/Katherin/NY/weather/'
filess_weather=os.listdir(path_files)
filess_weather.sort()
if filess_weather[0]=='.DS_Store':
    filess_weather.remove('.DS_Store')
timess_CP,tempe_c_CP,tempe_f_CP,dew_CP,humidi_CP,wind_CP,pressure_CP,precipi_CP=extract_all_weather(path_files,filess_weather)

#Plot CP weather
plot_CP_weather(timess_CP,tempe_c_CP,tempe_f_CP,dew_CP,humidi_CP,wind_CP,pressure_CP,precipi_CP)
#plt.savefig(plot_fol+'CP_waether_2022.pdf',bbox_inches='tight') 
make_plot_T_RH_2022(colorss)
#plt.savefig(plot_fol+'NY_2022_2.pdf',bbox_inches='tight') 

