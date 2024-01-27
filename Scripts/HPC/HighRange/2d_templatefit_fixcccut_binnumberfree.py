#!/usr/bin/env python
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile
from root_numpy import fill_hist
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
plt.switch_backend('agg')

parser = argparse.ArgumentParser()
parser.add_argument('--issversion', help='ISS data version, pass6 or pass7')
parser.add_argument('--cccut', type=float, help='cccut value,usually is 0.2')
parser.add_argument('--trdlow', type=float, help='trdlikelihood low value,usually is 0.2')
parser.add_argument('--trdhigh', type=float, help='trdlikelihood high value,usually is 0.2')
#parser.add_argument('--CCbinnumber',  type=int, help='CCbinnumbers. Default:20')
#parser.add_argument('--TRDbinnumber',  type=int, help='TRDbinnumbers. Default:12')
arguments = parser.parse_args()

if (arguments.issversion):
    ISSversion = arguments.issversion
else:
    print("You need to choose a ISS data cersion!")
    os._exit(0)

if (arguments.trdlow is not None):
    trdlow_value = arguments.trdlow
else:
    print("you need to give a trd low cut value for TF")
    os._exit(0)

if (arguments.trdhigh):
    trdhigh_value = arguments.trdhigh
else:
    print("you need to give a trd high cut value for TF")
    os._exit(0)

if (arguments.cccut):
    cccutvalue = arguments.cccut
else:
    print("you need to give a charge confusion cut value for TF, otherwise too much charge confused events will destory the TF.")
    os._exit(0)    
'''
if (arguments.CCbinnumber):
    CCbinningnumber = arguments.CCbinnumber
else:
    CCbinningnumber = 20

if (arguments.TRDbinnumber):
    TRDbinningnumber = arguments.TRDbinnumber
else:
    TRDbinningnumber = 12
'''
for Energybin in ["38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93.0","93.0_108.0","108.0_125.0","125.0_147.0","147_175", "175_211", "211_259", "259_330", "330_525",]:
    for TRDbinningnumber in np.arange(10,14):
        for CCbinningnumber in np.arange(9,19,2):
            ##################### data ####################################################
            template_correct=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity' + Energybin + '.npy')
            template_confused=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_negative_rigidity'+ Energybin +'.npy')
            template_electron=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/electron_negative'+ Energybin +'.npy')
            if ISSversion == 'pass6':
                data_negative=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity'+ Energybin +'.npy')
                data_positive=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity'+ Energybin +'.npy')
            elif ISSversion == 'pass7':
                data_negative=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_negative_rigidity'+ Energybin +'_pass7.npy')
                data_positive=np.load('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/plot_ISS_positive_rigidity'+ Energybin +'_pass7.npy')

            ######## write TH2D  ############################### 
            myfile = TFile ('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_' + str(ISSversion) + '/negative/rootfiles/histo_'+ Energybin +'cccut'+str(cccutvalue)+'CCN'+str(CCbinningnumber)+'TRDN'+str(TRDbinningnumber)+'.root','recreate')
            TH_template_correct = TH2D("template_correct","Antiproton", int(CCbinningnumber),float(cccutvalue),1,int(TRDbinningnumber), trdlow_value, trdhigh_value)
            TH_template_confused = TH2D("template_confused","Charge Confused Proton", int(CCbinningnumber),float(cccutvalue),1,int(TRDbinningnumber), trdlow_value, trdhigh_value)
            TH_template_electron = TH2D("template_electron","Electron", int(CCbinningnumber),float(cccutvalue),1,int(TRDbinningnumber), trdlow_value, trdhigh_value)
            TH_data_negative = TH2D("data_negative","Negative Rigidity Events", int(CCbinningnumber),float(cccutvalue),1,int(TRDbinningnumber), trdlow_value, trdhigh_value)
            TH_data_positive = TH2D("data_positive","Positive Rigidity Events", int(CCbinningnumber),float(cccutvalue),1,int(TRDbinningnumber), trdlow_value, trdhigh_value)

            TH_template_correct.Sumw2()
            TH_template_confused.Sumw2()
            TH_template_electron.Sumw2()
            TH_data_negative.Sumw2()
            TH_data_positive.Sumw2()

            fill_hist(TH_template_correct, template_correct)
            fill_hist(TH_template_confused, template_confused)
            fill_hist(TH_template_electron, template_electron)
            fill_hist(TH_data_negative, data_negative)
            fill_hist(TH_data_positive, data_positive)

            scale = 1/TH_template_correct.Integral()
            TH_template_correct.Scale(scale)
            scale = 1/TH_template_confused.Integral()
            TH_template_confused.Scale(scale)
            scale = 1/TH_template_electron.Integral()
            TH_template_electron.Scale(scale)


            TH_template_correct.Write()
            TH_template_confused.Write()
            TH_template_electron.Write()
            TH_data_negative.Write()
            TH_data_positive.Write()

'''
        ######## create graphs ###############################
        c1 = TCanvas()
        gPad.SetGrid()
        gPad.SetFrameFillColor(0)
        TH_template_correct.Draw('COLZ')
        TH_template_correct.SetStats(0)
        TH_template_correct.GetXaxis().SetTitle("#bf{Charge Confusion Estimator}")
        TH_template_correct.GetYaxis().SetTitle("#bf{TrdLikelihood}")
        c1.Update()
        c1.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_" + str(ISSversion) + "/negative/plots/correct"+ Energybin +'cccut'+str(cccutvalue)+'CCN'+str(CCbinningnumber)+'TRDN'+str(TRDbinningnumber)+".pdf")

        c1 = TCanvas()
        gPad.SetGrid()
        gPad.SetFrameFillColor(0)
        TH_template_confused.Draw('COLZ')
        TH_template_confused.SetStats(0)
        TH_template_confused.GetXaxis().SetTitle("#bf{Charge Confusion Estimator}")
        TH_template_confused.GetYaxis().SetTitle("#bf{TrdLikelihood}")
        c1.Update()
        c1.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_" + str(ISSversion) + "/negative/plots/confused"+ Energybin +'cccut'+str(cccutvalue)+'CCN'+str(CCbinningnumber)+'TRDN'+str(TRDbinningnumber)+".pdf")

        c1 = TCanvas()
        gPad.SetGrid()
        gPad.SetFrameFillColor(0)
        TH_template_electron.Draw('COLZ')
        TH_template_electron.SetStats(0)
        TH_template_electron.GetXaxis().SetTitle("#bf{Charge Confusion Estimator}")
        TH_template_electron.GetYaxis().SetTitle("#bf{TrdLikelihood}")
        c1.Update()
        c1.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_" + str(ISSversion) + "/negative/plots/electron"+ Energybin + 'cccut'+str(cccutvalue)+'CCN'+str(CCbinningnumber)+'TRDN'+str(TRDbinningnumber)+ ".pdf")

        c1 = TCanvas()
        gPad.SetGrid()
        gPad.SetFrameFillColor(0)
        TH_data_negative.Draw('COLZ')
        TH_data_negative.SetStats(0)
        TH_data_negative.GetXaxis().SetTitle("#bf{Charge Confusion Estimator}")
        TH_data_negative.GetYaxis().SetTitle("#bf{TrdLikelihood}")
        c1.Update()
        c1.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_" + str(ISSversion) + "/negative/plots/data_negative"+ Energybin +'cccut'+str(cccutvalue)+'CCN'+str(CCbinningnumber)+'TRDN'+str(TRDbinningnumber)+ ".pdf")
'''


