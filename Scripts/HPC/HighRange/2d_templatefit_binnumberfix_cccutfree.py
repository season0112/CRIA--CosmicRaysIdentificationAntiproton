#!/usr/bin/env python
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile
from root_numpy import fill_hist
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import binning
plt.switch_backend('agg')


def main():
    print("IfVGGP0 is " + str(IfVGGP0) )
    if IfVGGP0 == "No":
        vggSuffix = ''
    elif IfVGGP0 == "Yes":
        vggSuffix = '_VGG16NN'

    file_ProtonTemplate_Data             = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ProtonTemplate_Data' + ISSversion   + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')
    file_ChargeConfusedProtomTemplate_MC = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ChargeConfusedProtomTemplate_MC'    + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')
    file_ElectronTemplate_MC             = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ElectronTemplate_MC'                + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')
    file_ISS_negative                    = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ISS_negative' + ISSversion          + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')
    file_ISS_positive                    = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ISS_positive' + ISSversion          + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')
    file_ElectronTemplate_Data           = TFile (workpath + '/templatefit/negative/rootfiles/Histo_' + 'ElectronTemplate_Data' + ISSversion + '_Pattern_' + str(Pattern) + vggSuffix + '.root', 'recreate')

    #### L3,5 FINISHED, L0124 too long. (20 * 12 * 12 = 2880)
    #CCcutvalue_all       = np.arange(0.0, 1, 0.05)
    #CCbinningnumber_all  = np.append( np.arange(9, 31, 2), 20 )
    #TRDbinningnumber_all = np.append( np.arange(9, 31, 2), 20 )

    #### New test (8 * 9 * 9 = 648) (still take several hours)
    #CCcutvalue_all = np.array([0, 0.2, 0.35, 0.4, 0.65, 0.7, 0.8, 0.9])
    #CCbinningnumber_all = np.array([9, 11, 13, 19, 20, 21, 25, 27, 29])
    #TRDbinningnumber_all = np.array([9, 11, 13, 19, 20, 21, 25, 27, 29])

    #for CCbinningnumber in CCbinningnumber_all:
    #    for TRDbinningnumber in TRDbinningnumber_all:

    # official setting (For PhD thesis)
    #CCcutvalue_all = np.array([0, 0.2, 0.35, 0.4, 0.65, 0.7, 0.8, 0.9])
    #for (CCbinningnumber, TRDbinningnumber) in [(30,20), (20,20), (20,16), (20,12), (9,11)]:

    # CC test setting
    CCcutvalue_all = np.array([0, 0.2, 0.4, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
    for (CCbinningnumber, TRDbinningnumber) in [(9,11)]:

            for RigidityBin in ["14.1_15.3", "15.3_16.6", "16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175", "175_211", "211_259", "259_450", "211_250","250_330","330_525"]:

                for cccutvalue in CCcutvalue_all:

                    #### Load data 
                    MergedRigidityBin = binning.bins[0:13] #14.1-38.9 
                    if RigidityBin in MergedRigidityBin:
                        ChargeConfusedProtomTemplate_MC = np.load(workpath + '/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_' + str("14.1_38.9") + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    else:
                        ChargeConfusedProtomTemplate_MC = np.load(workpath + '/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_' + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    ElectronTemplate_Data = np.load(workpath + '/ElectronTemplate_Data__' + suffix + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    ProtonTemplate_Data   = np.load(workpath + '/ISS_positive_' + suffix + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    ElectronTemplate_MC   = np.load(workpath + '/ElectronTemplate_MC_B1091_el.pl1.0_25_2000_7.6_all_Tree_' + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    ISS_positive          = np.load(workpath + '/ISS_positive_' + suffix + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')
                    ISS_negative          = np.load(workpath + '/ISS_negative_' + suffix + str(RigidityBin) + '_Pattern_' + str(Pattern) + vggSuffix + '.npy')

                    #### write TH2D 
                    cccutvalue = "{:.2f}".format(cccutvalue)

                    TH_ProtonTemplate_Data_name             = "ProtonTemplate_Data_"             + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)
                    TH_ChargeConfusedProtomTemplate_MC_name = "ChargeConfusedProtomTemplate_MC_" + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)
                    TH_ElectronTemplate_MC_name             = "ElectronTemplate_MC_"             + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)
                    TH_ISS_negative_name                    = "ISS_negative_"                    + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)
                    TH_ISS_positive_name                    = "ISS_positive_"                    + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)
                    TH_ElectronTemplate_Data_name           = "ElectronTemplate_Data_"           + RigidityBin + '_cccut_' + cccutvalue + '_CCN_' + str(CCbinningnumber) + '_TRDN_' + str(TRDbinningnumber)

                    TH_ProtonTemplate_Data             = TH2D(TH_ProtonTemplate_Data_name            , "Antiproton"              , int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)
                    TH_ChargeConfusedProtomTemplate_MC = TH2D(TH_ChargeConfusedProtomTemplate_MC_name, "Charge Confused Proton"  , int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)
                    TH_ElectronTemplate_MC             = TH2D(TH_ElectronTemplate_MC_name            , "Electron"                , int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)
                    TH_ISS_negative                    = TH2D(TH_ISS_negative_name                   , "Negative Rigidity Events", int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)
                    TH_ISS_positive                    = TH2D(TH_ISS_positive_name                   , "Positive Rigidity Events", int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)
                    TH_ElectronTemplate_Data           = TH2D(TH_ElectronTemplate_Data_name          , "Electron"                , int(CCbinningnumber), float(cccutvalue), 1, int(TRDbinningnumber), trdlow_value, trdhigh_value)

                    TH_ProtonTemplate_Data.Sumw2()
                    TH_ChargeConfusedProtomTemplate_MC.Sumw2()
                    TH_ElectronTemplate_MC.Sumw2()
                    TH_ISS_negative.Sumw2()
                    TH_ISS_positive.Sumw2()
                    TH_ElectronTemplate_Data.Sumw2()

                    #### Fill TH2D
                    fill_hist(TH_ProtonTemplate_Data, ProtonTemplate_Data[:,0:2]) 
                    fill_hist(TH_ChargeConfusedProtomTemplate_MC, ChargeConfusedProtomTemplate_MC[:,0:2])  #No Reweight
                    #fill_hist(TH_ChargeConfusedProtomTemplate_MC, ChargeConfusedProtomTemplate_MC[:,0:2], ChargeConfusedProtomTemplate_MC[:,2])  #Reweight with weight 
                    #fill_hist(TH_ChargeConfusedProtomTemplate_MC, ChargeConfusedProtomTemplate_MC[:,0:2], ChargeConfusedProtomTemplate_MC[:,2]*1.5)  #Reweight with weight * 1.5 (NOTE:same as original reweighting)
                    #fill_hist(TH_ChargeConfusedProtomTemplate_MC, ChargeConfusedProtomTemplate_MC[:,0:2], np.ones(ChargeConfusedProtomTemplate_MC.shape[0]) ) #Reweight with 1.0
                    #fill_hist(TH_ChargeConfusedProtomTemplate_MC, ChargeConfusedProtomTemplate_MC[:,0:2], np.ones(ChargeConfusedProtomTemplate_MC.shape[0])*100 ) #Reweight with 100 (NOTE: same as reweight 1.0)         
                    fill_hist(TH_ElectronTemplate_MC, ElectronTemplate_MC[:,0:2], )
                    fill_hist(TH_ISS_negative, ISS_negative[:,0:2])
                    fill_hist(TH_ISS_positive, ISS_positive[:,0:2])
                    fill_hist(TH_ElectronTemplate_Data, ElectronTemplate_Data[:,0:2])

                    #### Nomalization
                    if TH_ProtonTemplate_Data.Integral() != 0:  
                        scale = 1/TH_ProtonTemplate_Data.Integral()
                        TH_ProtonTemplate_Data.Scale(scale)
                    else:
                        TH_ProtonTemplate_Data.Scale(0)

                    if TH_ChargeConfusedProtomTemplate_MC.Integral() != 0:
                        scale = 1/TH_ChargeConfusedProtomTemplate_MC.Integral()
                        TH_ChargeConfusedProtomTemplate_MC.Scale(scale)
                        #scale = 1/TH_ChargeConfusedProtomTemplate_MC.Integral() # test:nomalization, seems that for binned likelihood fit, nomalization make no differences.
                        #TH_ChargeConfusedProtomTemplate_MC.Scale(scale*5)       # test:nomalization
                    else:
                        TH_ChargeConfusedProtomTemplate_MC.Scale(0)

                    if TH_ElectronTemplate_MC.Integral() != 0:
                        scale = 1/TH_ElectronTemplate_MC.Integral()
                        TH_ElectronTemplate_MC.Scale(scale)
                    else:
                        TH_ElectronTemplate_MC.Scale(0)

                    if TH_ElectronTemplate_Data.Integral() != 0:
                        scale = 1/TH_ElectronTemplate_Data.Integral()
                        TH_ElectronTemplate_Data.Scale(scale)
                    else:
                        TH_ElectronTemplate_Data.Scale(0)


                    #### Write 
                    file_ProtonTemplate_Data.cd()
                    TH_ProtonTemplate_Data.Write()

                    file_ChargeConfusedProtomTemplate_MC.cd()
                    TH_ChargeConfusedProtomTemplate_MC.Write()

                    file_ElectronTemplate_MC.cd()
                    TH_ElectronTemplate_MC.Write()

                    file_ISS_negative.cd()
                    TH_ISS_negative.Write()

                    file_ISS_positive.cd()
                    TH_ISS_positive.Write()

                    file_ElectronTemplate_Data.cd()
                    TH_ElectronTemplate_Data.Write()

    file_ProtonTemplate_Data.Close()
    file_ChargeConfusedProtomTemplate_MC.Close()
    file_ElectronTemplate_MC.Close()
    file_ISS_negative.Close()
    file_ISS_positive.Close()
    file_ElectronTemplate_Data.Close()


if __name__ == "__main__":

    #### Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--trdlow', type=float, help='trdlikelihood low value,usually is 0.2')
    parser.add_argument('--trdhigh', type=float, help='trdlikelihood high value,usually is 1.2')
    parser.add_argument('--pattern', help='which tracker pattern you choose.')
    parser.add_argument('--ifVGGP0', type=str, help='For pattern 0, if turn to Neural Network.')
    arguments = parser.parse_args()

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data cersion!")
        os._exit(0)

    if (arguments.trdlow is not None):
        trdlow_value = arguments.trdlow
    else:
        print("you need to give a trd low cut value for TF.")
        os._exit(0)

    if (arguments.trdhigh):
        trdhigh_value = arguments.trdhigh
    else:
        print("you need to give a trd high cut value for TF.")
        os._exit(0)

    if (arguments.pattern):
        Pattern = arguments.pattern
    else:
        print("you need to give a tracker pattern.")
        os._exit(0)

    if (arguments.ifVGGP0):
        IfVGGP0 = arguments.ifVGGP0
    else:
        print("You need to choose TMVA or Neural Network")
        os._exit(0)

    if ISSversion == "pass7.8":
        suffix = ""
    elif ISSversion == "published2016":
        suffix = "May2015_"
    elif ISSversion == "PhyRep2021":
        suffix = "Nov2017_"

    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'

    main()


