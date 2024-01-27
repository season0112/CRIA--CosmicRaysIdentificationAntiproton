#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors
from root_numpy import fill_hist, root2array, tree2array
import os
import binning
import mvaplots_tool 

def main():

    # For AntiprotonHighEnergy_v2.0
    #MC_Name  = "B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree"
    MC_Name = "B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree"

    # For AntiprotonHighEnergy_v3.0
    #MC_Name  = "B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree"
    #MC_Name  = "B1220_pr.pl1ph.021000_7.8_all_Tree"
    #MC_Name  = "B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree"
    #MC_Name  = "B1042_pr.pl1.1800_7.6_all_Tree"
    #MC_Name  = "B1042_antipr.pl1.1800_7.6_all_Tree"

    Data_Name = "B1130_pass7_7.8_all_Tree"

    #DataPeriod = "May2015_"
    #DataPeriod = "Nov2017_"
    DataPeriod = ""


    #### Rigidity Loop
    #for i in range(0, binning.bins.shape[0], binmerge):  #binning.bins.shape[0]=32
    for i in range(31, 32, binmerge):  # 330_525
    #for i in range(27, 28, binmerge):  # 147_175
    #for i in range(27, 32, binmerge):  


        # Load MC and ISS data
        Correct_MC    = np.array([])
        Confused_MC   = np.array([])
        Correct_data  = np.array([])
        Confused_data = np.array([])

        for filesmerge_index in range(binmerge):

            print( "Now start rigidity bin:" + str(binning.bins[i+filesmerge_index]) )

            Correct_MC_tem    = root2array(highpath + "/" + MC_Name + "_positive_" + binning.bins[i+filesmerge_index] + ".root", selection="Pattern=="+Pattern)
            if Pattern == '0':
                Confused_MC_tem   = root2array(highpath + "/" + MC_Name + "_negative_" + binning.bins[i+filesmerge_index] + ".root", selection="-999<EcalBDT_EnergyD & EcalBDT_EnergyD<-0.7 & Pattern=="+Pattern)
            else:
                Confused_MC_tem   = root2array(highpath + "/" + MC_Name + "_negative_" + binning.bins[i+filesmerge_index] + ".root", selection="Pattern=="+Pattern)
 
            Correct_data_tem  = root2array(highpath + "/" + Data_Name + "_positive_" + DataPeriod + binning.bins[i+filesmerge_index] + ".root", selection="Pattern=="+Pattern)
            if Pattern == '0':
                Confused_data_tem = root2array(highpath + "/" + Data_Name + "_negative_" + DataPeriod + binning.bins[i+filesmerge_index] + ".root", selection="-999<EcalBDT_EnergyD & EcalBDT_EnergyD<-0.7 & Pattern=="+Pattern)
            else:
                Confused_data_tem = root2array(highpath + "/" + Data_Name + "_negative_" + DataPeriod + binning.bins[i+filesmerge_index] + ".root", selection="Pattern=="+Pattern)

            if Correct_MC.shape[0] == 0:
                Correct_MC = Correct_MC_tem
                Confused_MC = Confused_MC_tem
                Correct_data = Correct_data_tem
                Confused_data = Confused_data_tem
            else:
                Correct_MC = np.append(Correct_MC, Correct_MC_tem)
                Confused_MC = np.append(Confused_MC, Confused_MC_tem)
                Correct_data = np.append(Correct_data, Correct_data_tem)
                Confused_data = np.append(Confused_data, Confused_data_tem)


        # Get the rigidity sign 
        sign_Correct_MC    = Correct_MC['Rigidity']   / np.abs(Correct_MC['Rigidity'])
        sign_Confused_MC   = Confused_MC['Rigidity']  / np.abs(Confused_MC['Rigidity'])
        sign_Correct_data  = Correct_data['Rigidity'] / np.abs(Correct_data['Rigidity'])
        sign_Confused_data = Confused_data['Rigidity']/ np.abs(Confused_data['Rigidity'])

        # Get MVA Variables
        Correct_MC    = mvaplots_tool.GetMVAvariables(Correct_MC   , sign_Correct_MC)
        Confused_MC   = mvaplots_tool.GetMVAvariables(Confused_MC  , sign_Confused_MC)
        Correct_data  = mvaplots_tool.GetMVAvariables(Correct_data , sign_Correct_data)
        Confused_data = mvaplots_tool.GetMVAvariables(Confused_data, sign_Confused_data)
        namelabel = ['TRDLikelihood', 'RigidityAsymmetry', 'RigidityAsymmetryL9', 'Chi2TrackerYAsymmetry', 'InnerMaxSpanRigidityMatching', 'L1L9RigidityMatching', 'InnerRigidityMatch', 'InnerTrackerChi2X', 'InnerTrackerChi2Y', 'TrackerChi2X', 'TrackerChi2Y', 'TrackerChargeAsymmetry', 'TrackerL9Charge', 'TrackerL78Charge', 'UpperTofCharge', 'LowerTofCharge']
        namelabel_symbol = [r'\Lambda_{TRD}', r'\it{\delta}_{R}', r'\it{\delta}_{R_{L9}}', r'\it{\delta}_{\chi^{2}y}', r'\Gamma_{Inner}', '\Gamma_{L1L9}', r'\Gamma_{Central}', r'log \chi^{2}_{X_{Inner}}', r'log \chi^{2}_{X_{Inner}}', r'log \chi^{2}_{X}', r'log \chi^{2}_{Y}', r'\it{\delta}_{Q}', r'\it{Q}_{L9}', r'\it{Q}_{L78}', r'\it{Q}_{UppTOF}', r'\it{Q}_{LowTOF}']


        # Define the range of histogram plot
        # index=10: draw ccproton first, because highest value is higher.
        # index=12,13: draw legend on right top corner.
        #low_all  = [0,  -3, -5,  -1,  -5, -3, -10, -4, -4, -4, -4, -3, 0, 0, 0, 0] #147-175
        #high_all = [1.5, 4,  5, 0.2,   5,  3,  10,  2,  2,  1,  1,  3, 7, 5, 2, 2] #147-175

        # index=1,4: draw legend on right top corner.
        # index=3: legend on central
        low_all       = [0.0,  -10,  -6,  -1.2,   -3, -2.5, -10, -4, -4.3, -2.4, -2, -1.7, 0.1, 0.1, 0.75, 0.75] #330-525GV  
        high_all      = [1.2,   80,   7,   0.7,  3.5,  2.6,  10,  2,    2,    1,  1,  2.6,   7, 3.2, 1.55,    2] #330-525GV  
        binnumber_all = [ 80,  100,  80,   100,   80,   90,  80, 80,   80,   80, 80,   80,  80,  80,   80,   80] #330-525GV

        legend_xlow_all  = [0.20, 0.45, 0.55, 0.30, 0.30, 0.30, 0.30, 0.20, 0.18, 0.20, 0.49, 0.54, 0.45, 0.49, 0.30, 0.49]
        legend_xhigh_all = [0.40, 0.85, 0.75, 0.70, 0.70, 0.70, 0.70, 0.40, 0.38, 0.40, 0.89, 0.84, 0.85, 0.89, 0.70, 0.89]
        legend_ylow_all  = [0.68, 0.65, 0.65, 0.25, 0.25, 0.25, 0.20, 0.65, 0.68, 0.68, 0.20, 0.66, 0.65, 0.65, 0.25, 0.69]
        legend_yhigh_all = [0.88, 0.85, 0.85, 0.45, 0.45, 0.45, 0.40, 0.85, 0.88, 0.88, 0.40, 0.86, 0.85, 0.85, 0.45, 0.89]

        #### MVA variables Loop
        for index in range(16):
        #for index in [11]:
            low          = low_all[index]
            high         = high_all[index]
            legend_xlow  = legend_xlow_all[index]
            legend_xhigh = legend_xhigh_all[index]
            legend_ylow  = legend_ylow_all[index]
            legend_yhigh = legend_yhigh_all[index]
            binnumber    = binnumber_all[index]

            TH_Correct_MC    = TH1D("Correct_MC"   , "", binnumber, low, high)
            TH_Confused_MC   = TH1D("Confused_MC"  , "", binnumber, low, high)
            TH_Correct_data  = TH1D("Correct_data" , "", binnumber, low, high)
            TH_Confused_data = TH1D("Confused_data", "", binnumber, low, high)
            '''
            fill_hist(TH_Correct_MC   , Correct_MC[index, :]   , weights=Correct_MC[-1, :]   )
            fill_hist(TH_Confused_MC  , Confused_MC[index, :]  , weights=Confused_MC[-1, :]  )
            fill_hist(TH_Correct_data , Correct_data[index, :] , weights=Correct_data[-1, :] )
            fill_hist(TH_Confused_data, Confused_data[index, :], weights=Confused_data[-1, :])
            '''
            fill_hist(TH_Correct_MC   , Correct_MC[index, :]   )
            fill_hist(TH_Confused_MC  , Confused_MC[index, :]  )
            fill_hist(TH_Correct_data , Correct_data[index, :] )
            fill_hist(TH_Confused_data, Confused_data[index, :])

            TH_Correct_MC.Sumw2();
            TH_Confused_MC.Sumw2();
            TH_Correct_data.Sumw2();
            TH_Confused_data.Sumw2();

            scale = 150.0/TH_Correct_MC.Integral()
            TH_Correct_MC.Scale(scale)
            scale = 150.0/TH_Confused_MC.Integral()
            TH_Confused_MC.Scale(scale)
            scale = 150.0/TH_Correct_data.Integral()
            TH_Correct_data.Scale(scale)
            scale = 150.0/TH_Confused_data.Integral()
            TH_Confused_data.Scale(scale)
            
            '''
            TH_Correct_MC.Scale(TH_Correct_data.GetEntries()/TH_Correct_MC.GetEntries())
            TH_Confused_MC.Scale(TH_Confused_data.GetEntries()/TH_Confused_MC.GetEntries())
            '''

            g_Correct_data  = TGraphErrors(TH_Correct_data)
            g_Confused_data = TGraphErrors(TH_Confused_data)

            for index_Correct_data in range(g_Correct_data.GetN()):
                g_Correct_data.SetPointError(index_Correct_data,0, g_Correct_data.GetErrorY(index_Correct_data) );
                g_Confused_data.SetPointError(index_Correct_data,0, g_Confused_data.GetErrorY(index_Correct_data) );

            # Plot the MVA variables
            c1 = mvaplots_tool.PlotMVAvariables(highpath, i, binmerge, namelabel, namelabel_symbol, index, TH_Correct_MC, TH_Confused_MC, g_Correct_data, g_Confused_data, MC_Name, legend_xlow, legend_xhigh, legend_ylow, legend_yhigh, Pattern)


            ''' 
            # Save 
            ROOTFile  = TFile(highpath + "/mvaplots/proton_" + binning.bins[i] + "_to_" + binning.bins[i+binmerge-1] + str("_") + namelabel[index] + str("_") + str(index) + str("_") + MC_Name + ".root","RECREATE")
            TH_Correct_MC.Write()
            TH_Confused_MC.Write()
            TH_Correct_data.Write()
            TH_Confused_data.Write()
            g_Correct_data.Write()
            g_Confused_data.Write()
            c1.Write()
            ROOTFile.Close()
            '''             

            TH_Correct_MC.Delete()
            TH_Confused_MC.Delete()
            TH_Correct_data.Delete()
            TH_Confused_data.Delete()



if __name__ == '__main__':

    #highpath = os.getenv('HPCHIGHENERGYDATADIR')
    highpath = "/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0"

    binmerge = 1

    Pattern='0'

    main()









