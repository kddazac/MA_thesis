#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 21:04:32 2021

@author: kddazac
"""
################################################################################################################################################
#Libraries
################################################################################################################################################
import numpy as np
from math import log, sin, cos, pow
import astropy
from astropy.io import ascii
from astropy.table import Table, Column, vstack
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
import collections
#from aux.functions_v8 import *
#from aux.func_plot_v8 import *
#from aux.vectorastrometry31 import *
#from aux.coordinates31 import *
import time
import matplotlib.path as mplPath
#import pyfits
import scipy
from scipy.stats import norm as norm
from scipy import integrate
from scipy.optimize import fsolve
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column, MaskedColumn
from astropy.io import fits
from scipy.stats import halfnorm
import pandas
import matplotlib.dates as mdates
import datetime
from time import strptime
from datetime import date, timedelta
import math
#%matplotlib qt
################################################################################################################################################
origin_path=os.getcwd()
plot_fol=out_fol=origin_path+'/plots_saam/'
################################################################################################################################################
#Functions
################################################################################################################################################
def celsius2fahrenheit(temp):
    """
    Returns temperature in Fahrenheit given Celsius temperature.
    """
    return ((9./5.* temp) + 32)

def make_plot(colorss):
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
 
    def convert_ax_c_to_celsius2(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c2.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c2.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess)):
        #ax.scatter(timess[i], tempe_cs[i],s=10,color=colorss[i])
        #ax.set_xticklabels(tts[i])
        if Pisos[i]=='1st':
            ax.scatter(timess[i], tempe_cs[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers[i])
        else:
            ax.scatter(timess[i], tempe_cs[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers[i])                
    
    ax.set_title('SAAM at DWRC')
    ax.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)

    ax2=fig.add_subplot(212)
    ax_c2 = ax2.twinx()
    ax2.callbacks.connect("ylim_changed", convert_ax_c_to_celsius2)
    #fig.autofmt_xdate()
    for i in range(len(timess_20)):
        #ax.scatter(timess[i], tempe_cs[i],s=10,color=colorss[i])
        #ax.set_xticklabels(tts[i])
        if Pisos_20[i]=='1st':
            ax2.scatter(timess_20[i], tempe_cs_20[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers_20[i])
        else:
            ax2.scatter(timess_20[i], tempe_cs_20[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers_20[i])                
    
    #ax2.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax2.grid()
    ax2.set_ylabel('Temperature ($^oC$)')
    #ax2.set_ylim([19.8,24.2])    
    ax_c2.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)   
    return

def make_plot_T_RH_2019(colorss):
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess)):
        #ax.scatter(timess[i], tempe_cs[i],s=10,color=colorss[i])
        #ax.set_xticklabels(tts[i])
        if Pisos[i]=='1st':
            ax.scatter(timess[i], tempe_cs[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers[i])
        else:
            ax.scatter(timess[i], tempe_cs[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers[i])                
    ax.set_title('Renwick Gallery')
    #ax.legend(ncol=7,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), visible=False)
        
    ax1=fig.add_subplot(212,sharex=ax)
    for i in range(len(timess)):
        if Pisos[i]=='1st':
            ax1.scatter(timess[i], humidis[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers[i])
        else:
            ax1.scatter(timess[i], humidis[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers[i])
    ax1.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    return

def make_plot_T_RH_2019_DCA():
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    fig.autofmt_xdate()
    #for i in range(len(timess_med)):
    ax.scatter(timess_med, temp_med,color='g',marker='o',label='Renwick Gallery')         

    for i in range(len(timess_DCA)):
        ax.scatter(timess_DCA[i],tempe_c_DCA[i],color='k',marker='s')
    ax.scatter(timess_DCA[0],tempe_c_DCA[0],color='k',marker='s',label='Washington Weather')
    #ax.set_title('Renwick Gallery')
    ax.legend(ncol=4,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), visible=False)
        
    ax1=fig.add_subplot(212,sharex=ax)
    #for i in range(len(timess_meds)):
    ax1.scatter(timess_med_h, humids_med_h,color='g',marker='o')
    ax1.scatter(timess_med_h_saam, humids_med_h_saam,color='r',marker='o')
    ax1.scatter(timess_med_h_vic, humids_med_h_vic,color='b',marker='o')
    for i in range(len(timess_DCA)):
        ax1.scatter(timess_DCA[i],humidi_DCA[i],color='k',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    return

def make_plot_T_RH_2020_DCA():
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    ax.scatter(timess_med_20, temp_med_20,color='g',marker='o',label='Renwick Gallery')

    for i in range(len(timess_DCA20)):
        ax.scatter(timess_DCA20[i],tempe_c_DCA20[i],color='k',marker='s')           
    ax.scatter(timess_DCA20[0],tempe_c_DCA20[0],color='k',marker='s',label='Washington Weather')
    #ax.set_title('Renwick Gallery')
    ax.legend(ncol=4,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)

    ax1=fig.add_subplot(212,sharex=ax)
    ax1.scatter(timess_med_h20, humids_med_h20,color='g',marker='o')
    ax1.scatter(timess_med_h_20_saam, humids_med_h_20_saam,color='r',marker='o')
    ax1.scatter(timess_med_h_20_vic, humids_med_h_20_vic,color='b',marker='o')
    
    for i in range(len(timess_DCA20)):
        ax1.scatter(timess_DCA20[i],humidi_DCA20[i],color='k',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    return

def make_plot_T_RH_2020(colorss):
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    ax_c = ax.twinx()
    ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess_20)):
        #ax.scatter(timess[i], tempe_cs[i],s=10,color=colorss[i])
        #ax.set_xticklabels(tts[i])
        if Pisos_20[i]=='1st':
            ax.scatter(timess_20[i], tempe_cs_20[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers_20[i])
        else:
            ax.scatter(timess_20[i], tempe_cs_20[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers_20[i])                
    
    #ax.set_title('Graphic Art Lab (Victor Building)')
    ax.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    ax_c.set_ylabel('Temperature ($F$)')    
    plt.xticks(rotation=90)

    ax1=fig.add_subplot(212,sharex=ax)
    for i in range(len(timess_20)):
        if Pisos_20[i]=='1st':
            ax1.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers_20[i])
        else:
            ax1.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers_20[i])
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax1.grid()
    ax1.set_ylabel('Relative \n Humidity (%)')
    return

def make_plot_humidity(colorss):
    matplotlib.rcParams.update({'font.size': 23})    
    fig = plt.figure(figsize=(24,8))
    ax=fig.add_subplot(211)
    for i in range(len(timess)):
        if Pisos[i]=='1st':
            ax.scatter(timess[i], humidis[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers[i])
        else:
            ax.scatter(timess[i], humidis[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers[i])
    #ax.set_title('Graphic Art Lab (Victor Building)')
    ax.set_title('SAAM at DWRC')
    #ax.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax.grid()
    ax.set_ylabel('Relative \n Humidity (%)')
    #plt.xticks(rotation=180)

    ax2=fig.add_subplot(212)
    for i in range(len(timess_20)):
        if Pisos_20[i]=='1st':
            ax2.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='o',label='Datalogger #'+Dataloggers_20[i])
        elif Pisos_20[i]=='2nd':
            ax2.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='s', label='Datalogger #'+Dataloggers_20[i])
        elif Pisos_20[i]=='3rd':
            ax2.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers_20[i])
        else:
            ax2.scatter(timess_20[i], humidis_20[i],color=colorss[i],marker='+', label='Datalogger #'+Dataloggers_20[i])                    
    ax2.legend(ncol=6,loc=8,handletextpad=0.1,prop={'size': 20})
    ax2.grid()
    ax2.set_ylabel('Relative \n Humidity (%)')   
    #plt.xticks(rotation=180)   
    return
            
def extract_all_datalogs(renwick_files,time_up,time_do):
    timess=[]
    tempe_cs=[]
    tempe_fs=[]
    humidis=[]
    tts=[]
    Dataloggers=[]
    Pisos=[]
    for i in range(len(renwick_files)):    
        times,tempe_c,tempe_f,humidi,tt,Datalogger,Piso=extract_data_saam(path_files,renwick_files[i],time_up,time_do)
        timess.append(times)
        tempe_cs.append(tempe_c)
        tempe_fs.append(tempe_f)
        humidis.append(humidi)
        tts.append(tt)
        Dataloggers.append(Datalogger)
        Pisos.append(Piso)
    return timess,tempe_cs,tempe_fs,humidis,tts,Dataloggers,Pisos

def extract_data_saam(path_files,renwick_file,time_up,time_do):
    print(renwick_file)
    dummy_renwick=pandas.read_csv(path_files+renwick_file,header=7,skipfooter=11)
    Time=dummy_renwick['Time']
    times,times2,mm=convert_times(Time)
    Humidity=dummy_renwick['Humidity (RH(%))'][0:mm]
    Temperature_F=dummy_renwick['Temperature (F)'][0:mm]
    Temperature_C=transform_to_C(Temperature_F)
    #hhmm,datess=split_hhmm_date(Time)
    #print(renwick_file)
    timess=[]
    tempe_f=[]
    tempe_c=[]
    humidi=[]
    for i in range(len(times)):
        if times[i]<time_up and times[i]>time_do :
            timess.append(times[i])
            tempe_f.append(Temperature_F[i])
            tempe_c.append(Temperature_C[i])    #hhmm_=hhmm[4:52]
            humidi.append(Humidity[i])
    tt=set_ticks_times(times2)      
    kk=renwick_file.index('-')
    try:
        kk2=renwick_file.index('ick')
    except:
        try:
            kk2=renwick_file.index('DWRC')+1
        except:
            kk2=renwick_file.index('tor')
    Datalogger=renwick_file[kk2+3:kk]
    Piso=floor_renwick(Datalogger)    
    return timess,tempe_c,tempe_f,humidi,tt,Datalogger,Piso

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

def extract_central_avg(tt):
    new_tt=[]
    for i in range(len(tt)):
        new_tt.append(tt[i][1])
    return new_tt        

def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)

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

def split_points(pointt):
    indices = [i for i, x in enumerate(pointt) if x == '\t']
    _1=float(pointt[0:indices[0]])
    _2=float(pointt[indices[0]+1:indices[1]])
    _3=float(pointt[indices[1]+1:-1])   
    return _1,_2,_3

def floor_renwick(Datalogger):
    if Datalogger=='7' or Datalogger=='8' or Datalogger=='10' or Datalogger=='13' or Datalogger=='86' or Datalogger=='29':
        piso='2nd'
    elif Datalogger=='9' or Datalogger=='84' or Datalogger=='85' or Datalogger=='79':
        piso='1st'
    elif Datalogger=='87':
        piso='Basement'
    elif Datalogger=='26':
        piso='3rd'
    elif Datalogger=='18' or Datalogger=='24':
        piso='4th'
    else:
        piso='Graphic Arts Lab'
    return piso        

def convert_times(timesss):
    new_times=[]
    mm=0
    for i in range(len(timesss)):
        #print(i)
        try:
            new_times.append(datetime.datetime(int(timesss[i][-4:]),strptime(timesss[i][-8:-5],'%b').tm_mon,int(timesss[i][-11:-9]),int(timesss[i][0:2]),int(timesss[i][3:5])))        
            mm+=1
        except:
            break
    return new_times,timesss,mm

def split_files_buildings(filess):
    renwick_files,victor_files,dwrc_files=[],[],[]
    for i in range(len(filess)):
        aa=filess[i]
        if aa[0]=='R':
            renwick_files.append(aa)
        elif aa[0]=='V':
            victor_files.append(aa)
        elif aa[0]=='D':
            dwrc_files.append(aa)
        else:
            print('Something is wrong please check')
   
    return renwick_files,victor_files,dwrc_files

def transform_to_C(Temperature_F):
    Temperature_C=[]
    for i in range(len(Temperature_F)):
        Temperature_C.append((Temperature_F[i]-32.)*(5./9.))
    return Temperature_C

def split_hhmm_date(Time):
    hhmm=[]
    date=[]
    for i in range(len(Time)):
        aa=Time[i]
        bb=aa.index(' ')
        hhmm.append(aa[0:bb])
        date.append(aa[bb+1:])
    return hhmm,date

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

def zzz(timess,tempe_cs,months,logger):
    maps_month=OrderedDict()
    maps_month_avg=OrderedDict()
    for i in range(31):
        ii=i+1
        maps_month[str(ii)]=[]
    for j in range(len(timess[logger])):
        aa=timess[logger][j]
        if aa.month==months:                  
            maps_month[str(aa.day)].append(tempe_cs[logger][j])        
    for i in range(31):
        ii=i+1
        maps_month_avg[str(ii)]=np.median(maps_month[str(ii)])
    return maps_month_avg

def per_datalogger_avg(timess,tempe_cs,logger):
    maps_avg=OrderedDict()
    for z in range(12):
        maps_avg[str(z+1)]=[]
        maps_avg[str(z+1)].append(zzz(timess,tempe_cs,z+1,logger))            
    return maps_avg

def extract_time_median_T(maps_avgs_dataloggers,month,day,year):
    all_dataloggers=[]
    for gg in range(len(maps_avgs_dataloggers)):
        all_dataloggers.append(maps_avgs_dataloggers['datalogger_'+str(gg)][0][str(month)][0][str(day)])
    try:
        time=datetime.datetime(year,month,day)
        val=np.nanmedian(all_dataloggers)
    except:
        pass
    return time,val

def all_data_(timess,tempe_cs):
    maps_avgs_dataloggers=OrderedDict()
    for pp in range(len(timess)):
        maps_avgs_dataloggers['datalogger_'+str(pp)]=[]
        maps_avgs_dataloggers['datalogger_'+str(pp)].append(per_datalogger_avg(timess,tempe_cs,pp))
    return maps_avgs_dataloggers

def finally_data_(maps_avgs_dataloggers,year):
    finally_data=[]
    for j in range(12):
        for i in range(31):
            try:
                time,val=extract_time_median_T(maps_avgs_dataloggers,j+1,i+1,year)
                finally_data.append([time,val])
            except:
                pass
    return finally_data

def clean_finally_data_(finally_data):
    bad_list_ind=[]
    for i in range(len(finally_data)):
        aa=finally_data[i][-1]
        kk=math.isnan(aa)
        if kk==True:
            bad_list_ind.append(finally_data[i])
    
    for j in bad_list_ind:
        finally_data.remove(j)
    return finally_data

def plot_DCA_weather(timess_DCA,tempe_c_DCA,tempe_f_DCA,dew_DCA,humidi_DCA,wind_DCA,pressure_DCA,precipi_DCA,timess_DCA20,tempe_c_DCA20,tempe_f_DCA20,dew_DCA20,humidi_DCA20,wind_DCA20,pressure_DCA20,precipi_DCA20):
    def convert_ax_c_to_celsius(ax_f):
        """
        Update second axis according with first axis.
        """
        y1, y2 = ax_f.get_ylim()
        ax_c.set_ylim(celsius2fahrenheit(y1), celsius2fahrenheit(y2))
        ax_c.figure.canvas.draw()
        
    matplotlib.rcParams.update({'font.size': 18})    
    fig = plt.figure(figsize=(18,18))

    ax7=fig.add_subplot(621)
    #ax_c7 = ax7.twinx()
    #ax7.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess_DCA)):
        ax7.scatter(timess_DCA[i],tempe_c_DCA[i],color='c',marker='s')           
    #ax.set_title('Renwick Gallery')
    ax7.grid()
    ax7.set_ylabel('Temperature \n ($^oC$)')
    ax7.set_title('2019')
    #ax.set_ylim([19.8,24.2])
    #ax_c7.set_ylabel('Temperature ($F$)')    
    #plt.xticks(rotation=90)
    plt.setp(ax7.get_xticklabels(), visible=False)
    
    ax8=fig.add_subplot(623,sharex=ax7)    
    #ax_c8 = ax8.twinx()
    #ax8.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()    
    dew_DCA_C=[transform_to_C(dew_DCA[i]) for i in range(len(dew_DCA))]    
    for i in range(len(timess_DCA)):
        ax8.scatter(timess_DCA[i],dew_DCA_C[i],color='m',marker='s')          
    #ax.set_title('Renwick Gallery')
    ax8.grid()
    ax8.set_ylabel('Dew \n Point \n ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    #ax_c8.set_ylabel('Dew Point ($F$)')    
    #plt.xticks(rotation=90)
    plt.setp(ax8.get_xticklabels(), visible=False)
    
    ax9=fig.add_subplot(625,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax9.scatter(timess_DCA[i],humidi_DCA[i],color='g',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax9.grid()
    ax9.set_ylabel('Relative \n Humidity \n (%)')
    plt.setp(ax9.get_xticklabels(), visible=False)
    
    ax10=fig.add_subplot(627,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax10.scatter(timess_DCA[i],wind_DCA[i],color='r',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax10.grid()
    ax10.set_ylabel('Wind \n Speed \n (mph)')
    plt.setp(ax10.get_xticklabels(), visible=False)
    
    ax11=fig.add_subplot(6,2,9,sharex=ax7)    
    for i in range(len(timess_DCA)):
        ax11.scatter(timess_DCA[i],pressure_DCA[i],color='k',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax11.grid()
    ax11.set_ylabel('Pressure \n (Hg)')
    plt.setp(ax11.get_xticklabels(), visible=False)

    ax12=fig.add_subplot(6,2,11,sharex=ax7)
    for i in range(len(timess_DCA)):
        ax12.scatter(timess_DCA[i],precipi_DCA[i],color='b',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax12.grid()
    ax12.set_ylabel('Precipitation \n (in)')
    #plt.xticks(rotation=90)

    
    ax=fig.add_subplot(622,sharey=ax7)
    #ax_c = ax.twinx()
    #ax.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()
    for i in range(len(timess_DCA20)):
        ax.scatter(timess_DCA20[i],tempe_c_DCA20[i],color='c',marker='s')           
    #ax.set_title('Renwick Gallery')
    ax.grid()
    #ax.set_ylabel('Temperature ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    #ax_c.set_ylabel('Temperature ($F$)')    
    ax.set_title('2020')    
    plt.xticks(rotation=90)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)

    
    ax2=fig.add_subplot(624,sharey=ax8)    
    #ax_c2 = ax2.twinx()
    #ax2.callbacks.connect("ylim_changed", convert_ax_c_to_celsius)
    #fig.autofmt_xdate()    
    dew_DCA20_C=[transform_to_C(dew_DCA20[i]) for i in range(len(dew_DCA20))]    
    for i in range(len(timess_DCA20)):
        ax2.scatter(timess_DCA20[i],dew_DCA20_C[i],color='m',marker='s')          
    #ax.set_title('Renwick Gallery')
    ax2.grid()
    #ax2.set_ylabel('Dew Point ($^oC$)')
    #ax.set_ylim([19.8,24.2])
    #ax_c2.set_ylabel('Dew Point ($F$)')    
    #plt.xticks(rotation=90)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    
    
    ax3=fig.add_subplot(626,sharex=ax,sharey=ax9)    
    for i in range(len(timess_DCA20)):
        ax3.scatter(timess_DCA20[i],humidi_DCA20[i],color='g',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax3.grid()
    #ax3.set_ylabel('Relative \n Humidity (%)')
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)    
    
    ax4=fig.add_subplot(628,sharex=ax,sharey=ax10)    
    for i in range(len(timess_DCA20)):
        ax4.scatter(timess_DCA20[i],wind_DCA20[i],color='r',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax4.grid()
    #ax4.set_ylabel('Wind \n Speed (mph)')
    plt.setp(ax4.get_xticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    
    ax5=fig.add_subplot(6,2,10,sharex=ax,sharey=ax11)    
    for i in range(len(timess_DCA20)):
        ax5.scatter(timess_DCA20[i],pressure_DCA20[i],color='k',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax5.grid()
    #ax5.set_ylabel('Pressure (Hg)')
    plt.setp(ax5.get_xticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)

    ax6=fig.add_subplot(6,2,12,sharex=ax,sharey=ax12)    
    for i in range(len(timess_DCA20)):
        ax6.scatter(timess_DCA20[i],precipi_DCA20[i],color='b',marker='s')
    #ax1.legend(ncol=6,loc=9,handletextpad=0.1,prop={'size': 20})
    ax6.grid()
    #ax6.set_ylabel('Precipitation (in)')
    plt.setp(ax6.get_yticklabels(), visible=False)
    #plt.xticks(rotation=90)
    return

    
################################################################################################################################################
#Reading data
################################################################################################################################################
path_files='/Users/kddazac/Documents/Familia/Katherin/saam/dataloggers_data/Datalogger/'
filess_datalogs=os.listdir(path_files)
renwick_files,victor_files,dwrc_files=split_files_buildings(filess_datalogs)
time_start_2019=datetime.datetime(2019,1,1,00,00,00)
time_end_2019=datetime.datetime(2019,12,31,23,59,59)
time_start_2020=datetime.datetime(2020,1,1,00,00,00)
time_end_2020=datetime.datetime(2020,12,31,23,59,59)
#Weather DCA 2019
path_files19='/Users/kddazac/Documents/Familia/Katherin/saam/weather_DCA/2019/'
filess_weather=os.listdir(path_files19)
filess_weather.sort()
if filess_weather[0]=='.DS_Store':
    filess_weather.remove('.DS_Store')
timess_DCA,tempe_c_DCA,tempe_f_DCA,dew_DCA,humidi_DCA,wind_DCA,pressure_DCA,precipi_DCA=extract_all_weather(path_files19,filess_weather)
#Weather DCA 2020
path_files20='/Users/kddazac/Documents/Familia/Katherin/saam/weather_DCA/2020/'
filess_weather_20=os.listdir(path_files20)
filess_weather_20.sort()
if filess_weather_20[0]=='.DS_Store':
    filess_weather_20.remove('.DS_Store')
timess_DCA20,tempe_c_DCA20,tempe_f_DCA20,dew_DCA20,humidi_DCA20,wind_DCA20,pressure_DCA20,precipi_DCA20=extract_all_weather(path_files20,filess_weather_20)
#Plot DCA weather
plot_DCA_weather(timess_DCA,tempe_c_DCA,tempe_f_DCA,dew_DCA,humidi_DCA,wind_DCA,pressure_DCA,precipi_DCA,timess_DCA20,tempe_c_DCA20,tempe_f_DCA20,dew_DCA20,humidi_DCA20,wind_DCA20,pressure_DCA20,precipi_DCA20)
#plt.savefig(plot_fol+'DCA_waether_2019_2020.pdf',bbox_inches='tight') 
##############################################################################################################################
#Renwick
##############################################################################################################################
#read_csv_saam(renwick_files)
renwick_files_2019=[renwick_files[i] for i in range(len(renwick_files)) if renwick_files[i][-8:-4]=='2019']
renwick_files_2020=[renwick_files[i] for i in range(len(renwick_files)) if renwick_files[i][-8:-4]=='2020']
renwick_files_2019.sort()
renwick_files_2020.sort()
#2019
timess,tempe_cs,tempe_fs,humidis,tts,Dataloggers,Pisos=extract_all_datalogs(renwick_files_2019,time_end_2019,time_start_2019)
timess_20,tempe_cs_20,tempe_fs_20,humidis_20,tts_20,Dataloggers_20,Pisos_20=extract_all_datalogs(renwick_files_2020,time_end_2020,time_start_2020)

'''
#2019
#Medians Temp
maps_avgs_dataloggers=all_data_(timess,tempe_cs)
finally_data=finally_data_(maps_avgs_dataloggers,2019)
clean_finally_data=clean_finally_data_(finally_data)
timess_med=[clean_finally_data[i][0] for i in range(len(clean_finally_data))]
temp_med=[clean_finally_data[i][1] for i in range(len(clean_finally_data))]
#Medians Humidity
maps_avgs_dataloggers_h=all_data_(timess,humidis)
finally_data_h=finally_data_(maps_avgs_dataloggers_h,2019)
clean_finally_data_h=clean_finally_data_(finally_data_h)
timess_med_h=[clean_finally_data_h[i][0] for i in range(len(clean_finally_data_h))]
humids_med_h=[clean_finally_data_h[i][1] for i in range(len(clean_finally_data_h))]
'''
#2020
'''
maps_avgs_dataloggers20=all_data_(timess_20,tempe_cs_20)
finally_data20=finally_data_(maps_avgs_dataloggers20,2020)
clean_finally_data_20=clean_finally_data_(finally_data20)
timess_med_20=[clean_finally_data_20[i][0] for i in range(len(clean_finally_data_20))]
temp_med_20=[clean_finally_data_20[i][1] for i in range(len(clean_finally_data_20))]
#Medians Temp
maps_avgs_dataloggers_h20=all_data_(timess_20,humidis_20)
finally_data_h20=finally_data_(maps_avgs_dataloggers_h20,2020)
clean_finally_data_h20=clean_finally_data_(finally_data_h20)
timess_med_h20=[clean_finally_data_h20[i][0] for i in range(len(clean_finally_data_h20))]
humids_med_h20=[clean_finally_data_h20[i][1] for i in range(len(clean_finally_data_h20))]                
'''
#Values internally
time_start_2019=datetime.datetime(2019,1,1,00,00,00)
time_end_2019=datetime.datetime(2019,4,30,23,59,59)
timess,tempe_cs,tempe_fs,humidis,tts,Dataloggers,Pisos=extract_all_datalogs(renwick_files_2019,time_end_2019,time_start_2019)
#nanmedian
aa=tempe_cs
print(np.nanmedian(aa[0]))
print(np.nanmedian(aa[1]))
print(np.nanmedian(aa[2]))
print(np.nanmedian(aa[3]))
print(np.nanmedian(aa[4]))
print(np.nanmedian(aa[5]))    
print(np.min(aa[0]))
print(np.min(aa[1]))
print(np.min(aa[2]))
print(np.min(aa[3]))
print(np.min(aa[4]))
print(np.min(aa[5]))    
print(np.max(aa[0]))
print(np.max(aa[1]))
print(np.max(aa[2]))
print(np.max(aa[3]))
print(np.max(aa[4]))
print(np.max(aa[5]))
