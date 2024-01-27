import numpy as np
import matplotlib.pyplot as plt
import ROOT
from root_numpy import root2array, tree2array
import pandas as pd
from datetime import datetime
import matplotlib.dates as md
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import arrow

filename=ROOT.TFile("gap_merged.root")
intree=filename.Get("tree")
array=tree2array(intree)
'''
 uncal_time-toftime-dista/(29.98):run>>gg1
 toftime=0.5*(Tusedtof0+Tusedtof1)”)
 TCut cleanti = "acchi2<30 && npairs>0 && adc0>0 && adc1>0";
 TCut trdacc = "abs(phimiss-phibar)<15 && abs(AntiCrossZ)<35";
 TCut zacc = "abs(AntiCrossZ-cal_zeta)<15";
'''

'''
Yaxis: acczeta-AntiCrossZ

 TCut cleanti = "acchi2<30 && npairs>0 && adc0>0 && adc1>0";
 TCut trdacc = "abs(phimiss-phibar)<15 && abs(AntiCrossZ)<35";
 TCut tacc = "abs(uncal_time-toftime-dista/(29.98))<2";
'''

'''
Charge one Cut:0.8<0.49*(0.5*(Qusedtof0+Qusedtof1)+1.15*QTRD)<1.2

uncal_qside1 uncal_qside0 

 TCut cleanti = "acchi2<30 && npairs>0 && adc0>0 && adc1>0";
 TCut trdacc = "abs(phimiss-phibar)<15 && abs(AntiCrossZ)<35";
 TCut tacc = "abs(uncal_time-toftime-dista/(29.98))<2";
 TCut zacc = "abs(uncal_zeta-AntiCrossZ)<15";
'''
#TofTRD_charge = 0.5*(Qusedtof0+Qusedtof1)-QTRD

dt= array.dtype

reference_time_before_cuts = 0.5*(array['fdata_Tusedtof0'] + array['fdata_Tusedtof1']) + array['fdata_dista']/(29.98)

array_Time = array[np.where((array['fdata_acchi2']<30) & (array['idata_npairs']>0) & (array['idata_adc0']>0) & (array['idata_adc1']>0) & (abs(array['fdata_phimiss']-array['fdata_phibar'])<15) & (abs(array['fdata_AntiCrossZ'])<35) & (abs(array['fdata_AntiCrossZ']-array['fdata_cal_zeta'])<15))[0]]

array_Zeta =  array[np.where(  (array['fdata_acchi2']<30) & (array['idata_npairs']>0) & (array['idata_adc0']>0) & (array['idata_adc1']>0) & (abs(array['fdata_phimiss']-array['fdata_phibar'])<15) & (abs(array['fdata_AntiCrossZ'])<35) & (abs(array['fdata_uncal_time']-reference_time_before_cuts)<2.0 ))[0]]

array_Charge = array[np.where(  (array['fdata_acchi2']<30) & (array['idata_npairs']>0) & (array['idata_adc0']>0) & (array['idata_adc1']>0) & (abs(array['fdata_phimiss']-array['fdata_phibar'])<15) & (abs(array['fdata_AntiCrossZ'])<35) & (abs(array['fdata_uncal_time']-reference_time_before_cuts)<2.0 ) & (abs(array['fdata_uncal_zeta']-array['fdata_AntiCrossZ'])<15))[0]]


array_ChargeOne = array_Charge[np.where( (0.8<((array_Charge['fdata_Qusedtof0']+array_Charge['fdata_Qusedtof1'])*0.5+1.15*array_Charge['fdata_QTRD'])*0.49) & (((array_Charge['fdata_Qusedtof0']+array_Charge['fdata_Qusedtof1'])*0.5+1.15*array_Charge['fdata_QTRD'])*0.49<1.2) )[0]]
array_ChargeTwo = array_Charge[np.where( (1.8<((array_Charge['fdata_Qusedtof0']+array_Charge['fdata_Qusedtof1'])*0.5+1.15*array_Charge['fdata_QTRD'])*0.49) & (((array_Charge['fdata_Qusedtof0']+array_Charge['fdata_Qusedtof1'])*0.5+1.15*array_Charge['fdata_QTRD'])*0.49<2.2) )[0]]


for i in range(1,9):
    locals()['array_T%s' %i] = array_Time[np.where(array_Time['idata_sector']==i)[0]]
    locals()['reference_time%s' %i] = 0.5*(locals()['array_T%s' %i]['fdata_Tusedtof0'] + locals()['array_T%s' %i]['fdata_Tusedtof1']) + locals()['array_T%s' %i]['fdata_dista']/(29.98)
    locals()['array_Qone%s' %i] = array_ChargeOne[np.where(array_ChargeOne['idata_sector']==i)[0]]
    locals()['array_Qtwo%s' %i] = array_ChargeTwo[np.where(array_ChargeTwo['idata_sector']==i)[0]]
    locals()['array_Z%s' %i] = array_Zeta[np.where(array_Zeta['idata_sector']==i)[0]]

'''
array = np.zeros((8,array_Time.shape[0]))
array = np.array(array,dtype=dt)
    array[i] = np.resize(array[i],array_Time[np.where(array_Time['idata_sector']==i)[0]].shape[0])
    array[i][:] = array[i][0:array_Time[np.where(array_Time['idata_sector']==i)[0]].shape[0]]
    array[i][:] = array_Time[np.where(array_Time['idata_sector']==i)[0]]
'''

plt.figure(2)
dates=[datetime.fromtimestamp(i) for i in array_Z1['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Z1['fdata_uncal_zeta']-array_Z1['fdata_AntiCrossZ'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-40,40)
plt.savefig("zeta1.png")


plt.figure(3)   
dates=[datetime.fromtimestamp(i) for i in array_T1['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_T1['fdata_uncal_time']-reference_time1,'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-4,4)
plt.savefig("time1.png")


plt.figure(4)   
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, 0.5*(array_Qtwo8['fdata_Qusedtof0']+array_Qtwo8['fdata_Qusedtof1'])-array_Qtwo8['fdata_QTRD'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-4,6)
plt.savefig("chargetwo8_tof_trd.png")

plt.figure(5)   
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qtwo8['fdata_uncal_qside0']-array_Qtwo8['fdata_QTRD'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-4,6)
plt.savefig("chargetwo8_acc0_trd.png")

plt.figure(6)   
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qtwo8['fdata_uncal_qside0']-0.5*(array_Qtwo8['fdata_Qusedtof0']+array_Qtwo8['fdata_Qusedtof1']),'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-4,6)
plt.savefig("chargetwo8_acc0_tof.png")


plt.figure(7)
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qtwo8['fdata_uncal_qside0']-(0.49*0.5*(array_Qtwo8['fdata_Qusedtof0']+array_Qtwo8['fdata_Qusedtof1'])+1.15*array_Qtwo8['fdata_QTRD']),'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(-4,6)
plt.savefig("chargetwo8_acc0_tof_trd.png")



plt.figure(8)
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qtwo8['fdata_uncal_qside0'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(0,6)
plt.savefig("chargetwo8_acc0.png")


plt.figure(81)
dates=[datetime.fromtimestamp(i) for i in array_Qone8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qone8['fdata_uncal_qside0'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(0,6)
plt.savefig("chargeone8_acc0.png")

plt.figure(9)
dates=[datetime.fromtimestamp(i) for i in array_Qtwo8['idata_run']]
plt.xticks( rotation=25 )
ax=plt.gca()
plt.plot(dates, array_Qtwo8['fdata_uncal_qside1'],'bo',markersize=0.05)
major_f = md.DateFormatter('%Y-%m-%d %H')
major_L = md.DayLocator(bymonthday=[24,25,26,27,28])
ax.xaxis.set_major_locator(major_L)
ax.xaxis.set_major_formatter(major_f)
plt.grid(True,axis='y')
plt.ylim(0,6)
plt.savefig("chargetwo8_acc1.png")























