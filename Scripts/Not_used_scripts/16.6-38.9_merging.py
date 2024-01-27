#!/usr/bin/env python
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile
from root_numpy import fill_hist
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument('--issversion', help='ISS data version, pass6 or pass7 or pass7ext')
arguments = parser.parse_args()

if (arguments.issversion):
    ISSversion = arguments.issversion
else:
    print("You need to choose a ISS data cersion!")
    os._exit(0)

for Energybin in ["16.6_18.0","18.0_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9"]:
    ##################### data ####################################################
    if ISSversion == 'pass6':
        template_correct = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_positive_rigidity' + Energybin + '.npy')
        template_confused = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_negative_rigidity'+ Energybin +'.npy')
        template_electron = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/electron_negative'+ Energybin +'.npy')
        data_negative = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity'+ Energybin +'.npy')
        data_positive = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity'+ Energybin +'.npy')
    elif ISSversion == 'pass7':
        template_correct = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_positive_rigidity' + Energybin + '.npy')
        template_confused = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_negative_rigidity'+ Energybin +'.npy')
        template_electron = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/electron_negative'+ Energybin +'.npy')
        data_negative = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity'+ Energybin +'_pass7.npy')
        data_positive = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity'+ Energybin +'_pass7.npy')
    elif ISSversion == 'pass7ext':
        template_correct = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_positive_rigidity' + Energybin + '.npy')
        template_confused = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_negative_rigidity'+ Energybin +'.npy')
        template_electron = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/electron_negative'+ Energybin +'.npy')
        data_negative = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity'+ Energybin +'_pass7ext.npy')
        data_positive = np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity'+ Energybin +'_pass7ext.npy')

    if Energybin == "16.6_18.0":
        template_correct_all = template_correct
        template_confused_all = template_confused
        template_electron_all = template_electron
        data_negative_all = data_negative
        data_positive_all = data_positive
    else:
        template_correct_all = np.row_stack((template_correct_all, template_correct))
        template_confused_all = np.row_stack((template_confused_all, template_confused))
        template_electron_all = np.row_stack((template_electron_all, template_electron))
        data_negative_all = np.row_stack((data_negative_all, data_negative))
        data_positive_all = np.row_stack((data_positive_all, data_positive))

np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_positive_rigidity16.6_38.9.npy', template_correct_all)
np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_negative_rigidity16.6_38.9.npy', template_confused_all)
np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/electron_negative16.6_38.9.npy', template_electron_all)
if ISSversion == 'pass6':
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity16.6_38.9.npy', data_negative_all)
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity16.6_38.9.npy', data_positive_all)
elif ISSversion == 'pass7':
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity16.6_38.9_pass7.npy', data_negative_all)
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity16.6_38.9_pass7.npy', data_positive_all)
elif ISSversion == 'pass7ext':
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity16.6_38.9_pass7ext.npy', data_negative_all)
    np.save('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity16.6_38.9_pass7ext.npy', data_positive_all)



























