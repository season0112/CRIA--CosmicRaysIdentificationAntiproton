import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist
from root_numpy import array2tree, array2root


parser = argparse.ArgumentParser()
parser.add_argument('-D','--dataset', help='which_dataset_you_choose')
parser.add_argument('-E','--energy_range', choices=['16.6-38.9GeV', '38.9-147GeV', '38.9-147GeV_v2','147-1000GeV','147-1000GeV_v2'], help='which_energy_you_choose')
arguments = parser.parse_args()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)

if arguments.energy_range == '147-1000GeV_v2':
    binnings = np.array([147, 175, 211, 250, 330, 525, 1000])
elif arguments.energy_range == '38.9-147GeV_v2':
    binnings = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0])
else:
    print("Energy from 38.9 to 1000")
    binnings2 = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93.0, 108.0, 125.0, 147.0])
    binnings = np.array([147, 175, 211, 250, 330, 525, 1000])

a0 = np.array([])


for binleft in range(binnings.shape[0]-1):
    a_tem = np.load("/p/scratch/cvsk10/li8/analysis_v7.0" + arguments.dataset + "/" + "147-1000GeV_v2" + '/results/rawdata/transferdata/' + 'positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0'+'/'+'positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
    if a0.shape[0] == 0:
        a0 = a_tem
    else:
        a0 = np.row_stack((a0,a_tem))

if ('binnings2' in locals().keys()):
    for binleft in range(binnings2.shape[0]-1):
        a_tem = np.load("/p/scratch/cvsk10/li8/analysis_v7.0" + arguments.dataset + "/" + "38.9-147GeV_v2" + '/results/rawdata/transferdata/' + 'positive/'+str(binnings2[binleft])+'_'+str(binnings2[binleft+1])+'GeV/pattern0'+'/'+'positive_'+str(binnings2[binleft])+'_'+str(binnings2[binleft+1])+'_pattern_0.npy')
        a0 = np.row_stack((a0,a_tem))


a0.dtype = [('trdLogLikelihoodRatioElectronProtonTracker','float32'), ('rigidityAsymmetry','float32'), ('rigidityAsymmetryL9','float32'), ('chi2TrackerYAsymmetry','float32'), ('innerMaxSpanRigidityMatchingv2','float32'),('l1L9RigidityMatchingv2','float32'),('l24L58RigidityMatchingv2','float32'),('innerMaxSpanRigidityMatching','float32'),('l1L9RigidityMatching','float32'),('l24L58RigidityMatching','float32'),('log10Chi2TrackerXInner','float32'),('log10Chi2TrackerYInner','float32'),('log10Chi2TrackerX','float32'),('log10Chi2TrackerY','float32'),('trackerL58L24ChargeAsymmetry','float32'),('trackerL9Charge','float32'),('trackerL78Charge','float32'),('upperTofCharge','float32'),('lowerTofCharge','float32'),('rigidityWithoutThisHitlayer1','float32'),('rigidityWithoutThisHitlayer2','float32'),('rigidityWithoutThisHitlayer3','float32'),('rigidityWithoutThisHitlayer4','float32'),('rigidityWithoutThisHitlayer5','float32'),('rigidityWithoutThisHitlayer6','float32'),('rigidityWithoutThisHitlayer7','float32'),('rigidityWithoutThisHitlayer8','float32'),('rigidityWithoutThisHitlayer9','float32'),('timeStamp','float32'),('run','float32'),('eventNumber','float32'),('weight','float32'),('isEventMC','float32'),('mCParticleID','float32'),('mCPrimaryMomentum','float32'),('ecalEnergyElectron','float32'),('electronCCMVABDT','float32'),('protonCCMVABDT','float32'),('trdLogLikelihoodRatioProtonHeliumTracker','float32'),('ecalBDT_EnergyD','float32'),('ecalBDT_EnergyD_Smoothed','float32'),('pattern','float32'),('richbeta','float32'),('rigidity','float32')]



array2root(a0, '/p/scratch/cvsk10/li8/analysis_v7.0/test.root', 'tree')



