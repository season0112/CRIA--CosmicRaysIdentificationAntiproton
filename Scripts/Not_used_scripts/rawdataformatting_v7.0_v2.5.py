#!/usr/bin/env python
########################################################################
#                                                                      #
#                           Formatting data                            #
#          from root files to formatted dataset with features          #  
#                                                                      #  
########################################################################

########################## packages ####################################
from __future__ import division
from ROOT import TFile
from root_numpy import root2array, tree2array
import numpy as np
import sys
import os
import heapq
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--binrange', choices=["147-1130GeV", "38.9-147GeV","14.1_38.9GeV","1000-2000GeV","2.97-1130GeV","0.8_1130GeV"], help='bin_range.')
parser.add_argument('--dataset', help='which dataset you choose to do a rawdataformatting')
parser.add_argument('--rootfilenumber', help='rootfilenumber.')
parser.add_argument('--filename', help='automatically choose the iterated job.')
parser.add_argument('--cluster', help='which cluster you run your job.')
arguments = parser.parse_args()

if arguments.cluster == "JUAMS":
    workpath = os.getenv('JUAMSHIGHENERGYDATADIR')
elif arguments.cluster == "HPC":
    workpath = "/hpcwork/jara0052/sichen/analysis_7.1_fixwithIGRFcut"
    workpath = os.getenv('HPCHIGHENERGYDATADIR')

begintime = time.time()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path, exist_ok = True)
######################################################################
##MVARESULT = 37 
##PATTERN = 41 
##TRUERIGIDITY = 43 
MVARESULT = 28
PATTERN = 32
TRUERIGIDITY = 34

######################### root to file ################################
treelistname = sys.argv[8]

rootfile = arguments.rootfilenumber
result = None
with open (sys.argv[8],  "rb") as f:
    lines=f.readlines()
for treecount in range(0,int(rootfile)):
    filename=lines[treecount].strip().decode()
    myfile = TFile(filename)
    intree = myfile.Get('ExampleAnalysisTree')
    array = tree2array(intree)
    if result is None:
        result = array
    else:
        result = np.append(result,array)

############################################  features formatting  ###########################################################
## 1.get features for trainingset. ## 
timeStamp = result['TimeStamp']
run = result['Run']
eventNumber = result['EventNumber']
weight = result['Weight']
isEventMC = result['IsEventMC']
mCParticleID = result['MCParticleID']
mCPrimaryMomentum = result['MCPrimaryMomentum']
#rigidityWithoutThisHit = result['RigidityWithoutThisHit']
ecalEnergyElectron = result['EcalEnergyElectron']
electronCCMVABDT = result['ElectronCCMVABDT']

richbeta=result['RichBeta']
rigidity=result['Rigidity']
pattern=result['Pattern']   
protonCCMVABDT=result['ProtonCCMVABDT']
trdLogLikelihoodRatioProtonHeliumTracker=result['TrdLogLikelihoodRatioProtonHeliumTracker']

trdLogLikelihoodRatioElectronProtonTracker=result['TrdLogLikelihoodRatioElectronProtonTracker']
rigidityAsymmetry=result['RigidityAsymmetry']
rigidityAsymmetryL9=result['RigidityAsymmetryL9']
chi2TrackerYAsymmetry=result['Chi2TrackerYAsymmetry']
innerMaxSpanRigidityMatching=result['InnerMaxSpanRigidityMatching']
innerMaxSpanRigidityMatchingv2 = result['InnerMaxSpanRigidityMatchingv2']
l1L9RigidityMatching=result['L1L9RigidityMatching']
l1L9RigidityMatchingv2 = result['L1L9RigidityMatchingv2']
l24L58RigidityMatching=result['L24L58RigidityMatching']
l24L58RigidityMatchingv2 = result['L24L58RigidityMatchingv2']
log10Chi2TrackerXInner=result['Log10Chi2TrackerXInner']
log10Chi2TrackerYInner=result['Log10Chi2TrackerYInner']
log10Chi2TrackerX=result['Log10Chi2TrackerX']
log10Chi2TrackerY=result['Log10Chi2TrackerY']
trackerL58L24ChargeAsymmetry=result['TrackerL58L24ChargeAsymmetry']
trackerL9Charge=result['TrackerL9Charge']
trackerL78Charge=result['TrackerL78Charge']
upperTofCharge=result['UpperTofCharge']
lowerTofCharge=result['LowerTofCharge']
ecalBDT_EnergyD=result['EcalBDT_EnergyD']
ecalBDT_EnergyD_Smoothed=result['EcalBDT_EnergyD_Smoothed']
##hitZCoordinate=result['HitZLayer']

## 2. formatting hitCoordinate X Y Z && Charge X Y && rigidityWithoutThisHit ##
'''
for i in range(1,10):
        locals()['rigidityWithoutThisHitlayer%s' %i] = []
count=0
for count in range(result.size):
    layernumber=1
    layernumberiterate=1
    for layernumberiterate in range(1,10):
        if (layernumber-1) < hitZCoordinate[count].shape[0]:
            if hitZCoordinate[count][layernumber-1]==layernumberiterate:
                locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]=np.append((locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]),rigidityWithoutThisHit[count][layernumber-1])
                layernumber=layernumber+1
                layernumberiterate=layernumberiterate+1
            else:
                locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]=np.append((locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]),0.0)
                layernumberiterate=layernumberiterate+1
        else:
                locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]=np.append((locals()['rigidityWithoutThisHitlayer%s' %layernumberiterate]),0.0)
                layernumberiterate=layernumberiterate+1
    count=count+1
'''
## 3. getting np.array data ##
data_formatted=trdLogLikelihoodRatioElectronProtonTracker
data_formatted=np.c_[data_formatted,rigidityAsymmetry,\
rigidityAsymmetryL9,chi2TrackerYAsymmetry,innerMaxSpanRigidityMatchingv2,l1L9RigidityMatchingv2,\
l24L58RigidityMatchingv2,innerMaxSpanRigidityMatching,l1L9RigidityMatching,l24L58RigidityMatching,\
log10Chi2TrackerXInner,log10Chi2TrackerYInner,log10Chi2TrackerX,log10Chi2TrackerY,trackerL58L24ChargeAsymmetry,\
trackerL9Charge,trackerL78Charge,upperTofCharge,lowerTofCharge,\
# so far index:0-18, 19 in total.  3 mva not used, 16 in total.#
#rigidityWithoutThisHitlayer1,rigidityWithoutThisHitlayer2,rigidityWithoutThisHitlayer3,rigidityWithoutThisHitlayer4,\
#rigidityWithoutThisHitlayer5,rigidityWithoutThisHitlayer6,rigidityWithoutThisHitlayer7,rigidityWithoutThisHitlayer8,\
#rigidityWithoutThisHitlayer9,\
timeStamp, run, eventNumber, weight, isEventMC, mCParticleID, mCPrimaryMomentum, ecalEnergyElectron, electronCCMVABDT,\
protonCCMVABDT, trdLogLikelihoodRatioProtonHeliumTracker, ecalBDT_EnergyD, ecalBDT_EnergyD_Smoothed, pattern, richbeta, rigidity]
#rigidityWithoutThisHit

## 4. applicable & trdheliumproton ##
appindex=[]
heliumindex=[]
for count in range(data_formatted.shape[0]):
        if data_formatted[count,MVARESULT] == -2:
                appindex.append(count)
allindex = appindex
data_formatted = np.delete (data_formatted, allindex, 0)

## 5.  positive and negative ##
data_formatted_positive=np.array([])
data_formatted_negative=np.array([])
for i in range(data_formatted.shape[0]):
        if data_formatted[i,TRUERIGIDITY]>0:
                if data_formatted_positive.size:
                        data_formatted_positive=np.row_stack((data_formatted_positive,data_formatted[i,:]))
                else:
                        data_formatted_positive=data_formatted[i,:]
        else:
                if data_formatted_negative.size:
                        data_formatted_negative=np.row_stack((data_formatted_negative,data_formatted[i,:]))
                else:
                        data_formatted_negative=data_formatted[i,:]

## 6. split ## 
if arguments.binrange == "38.9-147GeV": 
    binnings = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0])
elif arguments.binrange == "147-1130GeV":
    binnings = np.array([147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
elif arguments.binrange == "14.1_38.9GeV":
    binnings = np.array([14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9])
elif arguments.binrange == "1000-2000GeV":
    binnings = np.array([1000,1300,1600,2000])
elif arguments.binrange == "2.97-1130GeV":
    binnings = np.array([2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
elif arguments.binrange == "0.8_1130GeV":
    binnings = np.array([0.8, 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11.0, 12.0, 13.0, 14.1, 15.3, 16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0, 175.0, 211.0, 250.0, 330.0, 525.0, 643.0, 822.0, 1130.0])
else:
    print("You need to choose a energyrange! Potential choises are: 147-1130GeV, 38.9-147GeV")
    os._exit(0)    

for binleft in range(binnings.shape[0]-1):
    locals()['index_p_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])] = []
    if data_formatted_positive.ndim == 2:
        for i in range(data_formatted_positive.shape[0]):
            if data_formatted_positive[i,TRUERIGIDITY] < binnings[binleft] or data_formatted_positive[i,TRUERIGIDITY]> binnings[binleft+1] :
                locals()['index_p_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])].append(i)
                p_tem = np.delete(data_formatted_positive, locals()['index_p_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])], 0)
    elif data_formatted_positive.shape[0] == 0:
        p_tem = np.array([])
    else:
        if data_formatted_positive[TRUERIGIDITY] < binnings[binleft] or data_formatted_positive[TRUERIGIDITY]> binnings[binleft+1] :
            p_tem = np.array([]) 
        else:
            p_tem = data_formatted_positive
    print('p_tem'+str(binnings[binleft])+'to'+str(binnings[binleft+1])+':'+str(p_tem.shape))
    pattern0p = np.array([])
    if p_tem.ndim == 2:
        for j in range(p_tem.shape[0]):
            if p_tem[j,PATTERN] == 0:
                if pattern0p.shape[0] == 0:
                    pattern0p = p_tem[j,:]
                else:
                    pattern0p = np.row_stack((pattern0p, p_tem[j,:]))
    elif p_tem.shape[0] != 0:
        if p_tem[PATTERN] == 0:
            pattern0p = p_tem

    mkdir(workpath+'/' + str(arguments.dataset) + '/' + str(arguments.binrange) + '_v2' + '/results/rawdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0')
    np.save(workpath+'/' + str(arguments.dataset) + '/' + str(arguments.binrange) + '_v2' + '/results/rawdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+ str(binnings[binleft]) +'_'+ str(binnings[binleft+1]) +'_pattern_0_'+treelistname[-5:]+'.npy',pattern0p)
    print('pattern0p:'+ str(pattern0p.shape))

    locals()['index_n_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])] = []
    if data_formatted_negative.ndim == 2:
        for i in range(data_formatted_negative.shape[0]):
            if data_formatted_negative[i,TRUERIGIDITY] > -binnings[binleft] or data_formatted_negative[i,TRUERIGIDITY] < -binnings[binleft+1] :
                locals()['index_n_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])].append(i)
                n_tem = np.delete(data_formatted_negative, locals()['index_n_' + str(binnings[binleft]) + '_' + str(binnings[binleft+1])], 0)
    elif data_formatted_negative.shape[0] == 0:
        n_tem = np.array([])
    else:
        if data_formatted_negative[TRUERIGIDITY] > -binnings[binleft] or data_formatted_negative[TRUERIGIDITY] < -binnings[binleft+1] :
            n_tem = np.array([])
        else:
            n_tem = data_formatted_negative
    print('n_tem'+str(binnings[binleft])+'to'+str(binnings[binleft+1])+':'+str(n_tem.shape))
    pattern0n = np.array([])
    if n_tem.ndim == 2:
        for j in range(n_tem.shape[0]):
            if n_tem[j,PATTERN] == 0:
                if pattern0n.shape[0] == 0:
                    pattern0n = n_tem[j,:]
                else:
                    pattern0n = np.row_stack((pattern0n, n_tem[j,:]))
    elif n_tem.shape[0] != 0:
        if n_tem[PATTERN] == 0:
            pattern0n = n_tem

    mkdir(workpath+'/' + str(arguments.dataset) + '/' + str(arguments.binrange) + '_v2' + '/results/rawdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0')
    np.save(workpath+'/' +str(arguments.dataset) + '/' + str(arguments.binrange) + '_v2' + '/results/rawdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+ str(binnings[binleft]) +'_'+ str(binnings[binleft+1]) +'_pattern_0_'+treelistname[-5:]+'.npy',pattern0n)
    print('pattern0n:'+ str(pattern0n.shape))

#########################################################################################
endtime = time.time()
print ((endtime - begintime)/60)

