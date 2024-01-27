import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TFile, TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import array2tree, array2root, tree2array

parser = argparse.ArgumentParser()
parser.add_argument('--valuemin', type=float, help='get bin value in between.')
parser.add_argument('--valuemax', type=float, help='get bin value in between.')
parser.add_argument('--rigidity', help='which rigidity bin')
parser.add_argument('--cluster', help='which cluster you choose')
arguments = parser.parse_args()

valuemin = arguments.valuemin
valuemax = arguments.valuemax
rigidity = arguments.rigidity

if arguments.cluster=="JUAMS":
    intermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
    os.chdir(intermediateworkpath+'/'+str(rigidity))
elif arguments.cluster=="HPC":
    intermediateworkpath = os.getenv('HPCINTERMEDIATEDIR')
    os.chdir(intermediateworkpath+'/'+str(rigidity))

#### bin 1
antiproton = TFile("B1042_antipr.pl1.1800_7.6_all_Tree_all.root")
antiproton_tree = antiproton.Get('AntiprotonIntermediateEnergyTree')
antiproton_array = tree2array(antiproton_tree)
antiproton_array = antiproton_array[np.where(antiproton_array['Rigidity']>-valuemax)[0]]
antiproton_array = antiproton_array[np.where(antiproton_array['Rigidity']<-valuemin)[0]]
array2root(antiproton_array, 'B1042_antipr.pl1.1800_7.6_all_Tree.root', 'AntiprotonIntermediateEnergyTree')

electron = TFile("B1091_el.pl1.0_25200_7.6_all_Tree_all.root")
electron_tree = electron.Get('AntiprotonIntermediateEnergyTree')
electron_array = tree2array(electron_tree)
electron_array = electron_array[np.where(electron_array['Rigidity']>-valuemax)[0]]
electron_array = electron_array[np.where(electron_array['Rigidity']<-valuemin)[0]]
array2root(electron_array, 'B1091_el.pl1.0_25200_7.6_all_Tree.root', 'AntiprotonIntermediateEnergyTree')

negative = TFile("B1130_pass7_7.7_all_Tree_negative_all.root")
negative_tree = negative.Get('AntiprotonIntermediateEnergyTree')
negative_array = tree2array(negative_tree)
negative_array = negative_array[np.where(negative_array['Rigidity']>-valuemax)[0]]
negative_array = negative_array[np.where(negative_array['Rigidity']<-valuemin)[0]]
array2root(negative_array, 'B1130_pass7_7.7_all_Tree_negative.root', 'AntiprotonIntermediateEnergyTree')

positive = TFile("B1130_pass7_7.7_all_Tree_positive_all.root")
positive_tree = positive.Get('AntiprotonIntermediateEnergyTree')
positive_array = tree2array(positive_tree)
positive_array = positive_array[np.where(positive_array['Rigidity']<valuemax)[0]]
positive_array = positive_array[np.where(positive_array['Rigidity']>valuemin)[0]]
array2root(positive_array, 'B1130_pass7_7.7_all_Tree_positive.root', 'AntiprotonIntermediateEnergyTree')

'''
### bin 2
antiproton_h = TFile("B1042_antipr.pl1.1800_7.6_all_Tree_all.root")
antiproton_tree_h = antiproton_h.Get('AntiprotonIntermediateEnergyTree')
antiproton_array_h = tree2array(antiproton_tree_h)
antiproton_array_h = antiproton_array_h[np.where(antiproton_array_h['Rigidity']<-value)[0]]
array2root(antiproton_array_h, 'B1042_antipr.pl1.1800_7.6_all_Tree_h.root', 'AntiprotonIntermediateEnergyTree')

electron_h = TFile("B1091_el.pl1.0_25200_7.6_all_Tree_all.root")
electron_tree_h = electron_h.Get('AntiprotonIntermediateEnergyTree')
electron_array_h = tree2array(electron_tree_h)
electron_array_h = electron_array_h[np.where(electron_array_h['Rigidity']<-value)[0]]
array2root(electron_array_h, 'B1091_el.pl1.0_25200_7.6_all_Tree_h.root', 'AntiprotonIntermediateEnergyTree')

negative_h = TFile("B1130_pass7_7.7_all_Tree_negative_all.root")
negative_tree_h = negative_h.Get('AntiprotonIntermediateEnergyTree')
negative_array_h = tree2array(negative_tree_h)
negative_array_h = negative_array_h[np.where(negative_array_h['Rigidity']<-value)[0]]
array2root(negative_array_h, 'B1130_pass7_7.7_all_Tree_negative_h.root', 'AntiprotonIntermediateEnergyTree')

positive_h = TFile("B1130_pass7_7.7_all_Tree_positive_all.root")
positive_tree_h = positive_h.Get('AntiprotonIntermediateEnergyTree')
positive_array_h = tree2array(positive_tree_h)
positive_array_h = positive_array_h[np.where(positive_array_h['Rigidity']>value)[0]]
array2root(positive_array_h, 'B1130_pass7_7.7_all_Tree_positive_h.root', 'AntiprotonIntermediateEnergyTree')
'''

