import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors
from root_numpy import fill_hist
from root_numpy import root2array, tree2array
import PythonPlotDefaultParameters
import binning
import os

highpath = os.getenv('HPCHIGHENERGYDATADIR')

Rigiditybin = binning.bins # 32, 14.1-525
Rigidityedge = binning.Newbinnings_525[26:] # 33 14.1-525
Rigiditycenter = (Rigidityedge[:-1] + Rigidityedge[1:]) / 2

#Rigiditybin = Rigiditybin[0:5]

'''
cc_percentage = []
for i in range(Rigiditybin.shape[0]):
    Correct_MC  = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/ISS_anylsis/data/plot_positive_rigidity" + str(Rigiditybin[i]) + ".npy")
    Confused_MC = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/ISS_anylsis/data/plot_negative_rigidity" + str(Rigiditybin[i]) + ".npy")
    Correct_data = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/ISS_anylsis/data/plot_ISS_positive_rigidity" + str(Rigiditybin[i]) + "_pass7.8.npy")
    Confused_data = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/ISS_anylsis/data/plot_ISS_negative_rigidity" + str(Rigiditybin[i]) + "_pass7.8.npy")

    #cc_percentage.append(Confused_MC.shape[0] / (Correct_MC.shape[0] + Confused_MC.shape[0]))
    cclevel = Confused_MC.shape[0] / (Confused_MC.shape[0] + Correct_MC.shape[0])
    cc_percentage.append( cclevel * (Correct_data.shape[0] + Confused_data.shape[0]) / Confused_data.shape[0])
    #cc_percentage.append( Confused_MC.shape[0] / Correct_MC.shape[0] * Correct_data.shape[0] / Confused_data.shape[0] )

plt.figure(figsize=(30,18))
plt.errorbar(Rigiditycenter, cc_percentage, yerr=0, label="", fmt='o-', markersize=15, color="black")
plt.savefig("CCPercentage" + ".pdf")
'''
###############################################
with open (highpath + '/ISS_anylsis/data/templatefit_pass7.8/negative/results_525version' + '/fit_resultscccut0.20CCN20TRDN16' + '.txt') as fileccp:
    ccp_number = fileccp.readlines()
ccp_number = map(lambda s: float(s.split(',')[1]) ,ccp_number)
ccp_number = list(ccp_number)
ccp_number = np.array(ccp_number)

with open (highpath + '/ISS_anylsis/data/templatefit_pass7.8/negative/results_525version_uncertainty' + '/fit_resultscccut0.20CCN9TRDN11' + '.txt') as fileccp_uncertainty:
    ccp_number_uncertainty = fileccp_uncertainty.readlines()
ccp_number_uncertainty = map(lambda s: float(s.split(',')[1]) ,ccp_number_uncertainty)
ccp_number_uncertainty = list(ccp_number_uncertainty)
ccp_number_uncertainty = np.array(ccp_number_uncertainty)

with open (highpath + '/ISS_anylsis/data/templatefit_pass7.8/negative/results_525version' + '/fit_resultscccut0.20CCN20TRDN16' + '.txt') as fileresult:
    a_number = fileresult.readlines()
a_number = map(lambda s: float(s.split(',')[0]) ,a_number)
a_number = list(a_number)
a_number = np.array(a_number)

with open (highpath + '/ISS_anylsis/data/templatefit_pass7.8/negative/results_525version' + '/fit_resultscccut0.20CCN20TRDN16' + '.txt') as fileresult:
    e_number = fileresult.readlines()
e_number = map(lambda s: float(s.split(',')[2]) ,e_number)
e_number = list(e_number)
e_number = np.array(e_number)

ccp_number = ccp_number[10:]
ccp_number_uncertainty = ccp_number_uncertainty[10:]
a_number = a_number[10:]
e_number = e_number[10:]
Rigiditycenter = Rigiditycenter[10:]

alldata = ccp_number + a_number + e_number
cc_percentage = ccp_number / alldata


plt.figure(figsize=(30,18))
ax=plt.gca()
plt.errorbar(Rigiditycenter, cc_percentage*100, yerr=0, label="", fmt='o-', markersize=15, color="black")
plt.fill_between(Rigiditycenter, 0, cc_percentage*100, color='green', alpha=0.3)
plt.fill_between(Rigiditycenter, cc_percentage*100, 100, color='magenta', alpha=0.3)
ax.axes.set_ylim([0,100])
plt.text(0.6, 0.4, 'Charge Confused Proton', fontsize=50, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
plt.xlabel("Rigidity (GV)",fontsize=40)
plt.ylabel("CCproton percentage (%)",fontsize=40)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
ax.tick_params(direction='in',length=10, width=3, axis='both', which='minor')
ax.tick_params(direction='in',length=20, width=6, axis='both', which='major')
plt.savefig(highpath + "/CCpercentageplot/CCPercentage" + ".pdf")


plt.figure(figsize=(30,18))
ax=plt.gca()
plt.errorbar(Rigiditycenter, alldata, yerr=0, label="", fmt='o-', markersize=15, color="black") 
plt.errorbar(Rigiditycenter, ccp_number+e_number, yerr=0, label="", fmt='o-', markersize=15, color="black")
plt.errorbar(Rigiditycenter, ccp_number, yerr=0, label="", fmt='o-', markersize=15, color="black")
plt.fill_between(Rigiditycenter, ccp_number+e_number, ccp_number, color='blue', alpha=0.3)
plt.fill_between(Rigiditycenter, ccp_number+e_number, alldata, color='red', alpha=0.3)
plt.fill_between(Rigiditycenter, ccp_number, 0, color='green', alpha=0.3)
plt.xlabel("Rigidity (GV)",fontsize=40)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
ax.tick_params(direction='in',length=10, width=3, axis='both', which='minor')
ax.tick_params(direction='in',length=20, width=6, axis='both', which='major')
plt.savefig(highpath + "/CCpercentageplot/composition" + ".pdf")


plt.figure(figsize=(30,18))
ax=plt.gca()
#plt.yscale('log')
plt.errorbar(Rigiditycenter, ccp_number, yerr=np.abs(ccp_number_uncertainty-ccp_number), label="", fmt='o-', markersize=15, color="black")
plt.xlabel("Rigidity (GV)",fontsize=40)
plt.ylabel("CCproton Numbers with uncertainties",fontsize=40)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
ax.tick_params(direction='in',length=10, width=3, axis='both', which='minor')
ax.tick_params(direction='in',length=20, width=6, axis='both', which='major')
plt.savefig(highpath + "/CCpercentageplot/CCProtonNumber_Uncertainties" + ".pdf")




