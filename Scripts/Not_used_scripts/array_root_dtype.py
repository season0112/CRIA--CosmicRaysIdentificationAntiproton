import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist
from root_numpy import array2tree, array2root

a0=np.load('negative_147_175_pattern_0.npy')

dt = dtype=[('trdLogLikelihoodRatioElectronProtonTracker','f4'),('rigidityAsymmetry','f4'),('rigidityAsymmetryL9','f4'),('chi2TrackerYAsymmetry','f4'),('innerMaxSpanRigidityMatchingv2','f4'),('l1L9RigidityMatchingv2','f4'),('l24L58RigidityMatchingv2','f4'),('innerMaxSpanRigidityMatching','f4'),('l1L9RigidityMatching','f4'),('l24L58RigidityMatching','f4'),('log10Chi2TrackerXInner','f4'),('log10Chi2TrackerYInner','f4'),('log10Chi2TrackerX','f4'),('log10Chi2TrackerY','f4'),('trackerL58L24ChargeAsymmetry','f4'),('trackerL9Charge','f4'),('trackerL78Charge','f4'),('upperTofCharge','f4'),('lowerTofCharge','f4'),('rigidityWithoutThisHitlayer1','f4'),('rigidityWithoutThisHitlayer2','f4'),('rigidityWithoutThisHitlayer3','f4'),('rigidityWithoutThisHitlayer4','f4'),('rigidityWithoutThisHitlayer5','f4'),('rigidityWithoutThisHitlayer6','f4'),('rigidityWithoutThisHitlayer7','f4'),('rigidityWithoutThisHitlayer8','f4'),('rigidityWithoutThisHitlayer9','f4'),('timeStamp','f4'),('run','f4'),('eventNumber','f4'),('weight','f4'),('isEventMC','f4'),('mCParticleID','f4'),('mCPrimaryMomentum','f4'),('ecalEnergyElectron','f4'),('electronCCMVABDT','f4'),('protonCCMVABDT','f4'),('trdLogLikelihoodRatioProtonHeliumTracker','f4'),('ecalBDT_EnergyD','f4'),('ecalBDT_EnergyD_Smoothed','f4'),('pattern','f4'),('richbeta','f4'),('rigidity','f4')]

b = np.array(a0,dtype=dt)

tree = array2tree(a0, name='tree')
#array2root(a0, 'selected_tree.root', 'tree')



