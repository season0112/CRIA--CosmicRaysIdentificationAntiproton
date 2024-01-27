#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend

parser = argparse.ArgumentParser()
parser.add_argument('-D','--dataset', help='which_dataset_you_choose')
parser.add_argument('-E','--energy_range', choices=['147-1000GeV_v2','147-2000GeV_v2'], help='which_energy_you_choose')
parser.add_argument('--ecalBDTCut', dest='ecalBDTCut', action='store_true',help='ecalBDTCut. Default:False')
parser.add_argument('--no-ecalBDTCut', dest='ecalBDTCut', action='store_false',help='ecalBDTCut. Default:False')
parser.set_defaults(ecalBDTCut=False)
parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
parser.set_defaults(trdProtonHeliumCut=False)
parser.add_argument('-C','--cluster', help='which_cluster_you_choose')
arguments = parser.parse_args()

## 
if arguments.energy_range == '147-1000GeV_v2':
    binnings = np.array([147, 175, 211, 259, 330, 525, 1000])
elif arguments.energy_range == '147-2000GeV_v2':
    binnings = np.array([147, 175, 211, 259, 330, 525, 1000, 1300, 1600, 2000])

EcalBDTCut = arguments.ecalBDTCut
TrdProtonHeliumCut = arguments.trdProtonHeliumCut

## events sum up
events = np.array([])
for binleft in range(binnings.shape[0]-1):
    if arguments.cluster == 'HPC': 
        events_tem_p = np.load('/hpcwork/jara0052/sichen/analysis_7.0/' + arguments.dataset + '/' + '147-2000GeV_v2/transferdata/positive/' + str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_' + str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy') 
#        events_tem_n = np.load('/hpcwork/jara0052/sichen/analysis_7.0/' + arguments.dataset + '/' + '147-2000GeV_v2/transferdata/negative/' + str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_' + str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
#        events_tem = np.row_stack((events_tem_p,events_tem_n))
        events_tem = events_tem_p
    elif arguments.cluster == 'JUAMS':
        events_tem_p = np.load('/p/scratch/cvsk10/li8/analysis_v7.0/' + arguments.dataset + "/" + "147-2000GeV_v2" + '/results/rawdata/transferdata/' + 'positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0'+'/'+'positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
#        events_tem_n = np.load('/p/scratch/cvsk10/li8/analysis_v7.0/' + arguments.dataset + "/" + "147-2000GeV_v2" + '/results/rawdata/transferdata/' + 'negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0'+'/'+'negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
#        events_tem = np.row_stack((events_tem_p,events_tem_n))
        events_tem = events_tem_p
    else:
        os._exit(0)
    if events.shape[0] == 0:
        events = events_tem
    else:
        events = np.row_stack((events,events_tem))

print(events.shape)

## TH2D and filling
rigiditybin = np.array(range(147,2000,10),dtype='f')
rigiditybin_rec = np.array(range(147,2000,10),dtype='f')
Unfolding_Matrices = TH2D("Unfolding_Matrices","", rigiditybin.shape[0]-1, rigiditybin, rigiditybin_rec.shape[0]-1, rigiditybin_rec)
#rigiditypoint = np.array( np.zeros(rigiditybin.shape[0]-1) )
#for i in range(rigiditybin.shape[0]-1):
#    rigiditypoint[i]=(rigiditybin[i]+rigiditybin[i+1])/2

## generation rigidity cuts
mCPrimaryMomentum = 34
events = events[np.where(events[:,mCPrimaryMomentum]>147)[0],:]
events = events[np.where(events[:,mCPrimaryMomentum]<2000)[0],:]
## 
EcalBDT_EnergyD = 39
EcalBDTCutvalue_select_proton = 0.0
TrdLogLikelihoodRatioProtonHeliumTracker = 38
TrdProtonHeliumCutValue = 0.3
if (EcalBDTCut):
    events = events[np.where((events[:,EcalBDT_EnergyD]>-2) & (events[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
if (TrdProtonHeliumCut):
    events = events[np.where(events[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

print(events.shape)

for entryNum in range(0,events.shape[0]):
    Unfolding_Matrices.Fill(events[entryNum,34],events[entryNum,43]) ## true x,  reconstructed y.

ROOTFile  = TFile("/p/scratch/cvsk10/li8/analysis_v7.0/RigidityResolution.root","RECREATE")
Unfolding_Matrices.Write("RigidityResolution")
ROOTFile.Close()


