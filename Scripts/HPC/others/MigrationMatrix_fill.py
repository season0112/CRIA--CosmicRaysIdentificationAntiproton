#!/usr/bin/env python
import numpy as np
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import root2array, tree2array
import binning

def main():

    #### Events Sum Up
    if RigidityRange == 'HIGHENERGYRANGE':
        events = np.array([])
        for binleft in range(binnings.shape[0]-1):
            events_tem = root2array( os.getenv('HPCHIGHENERGYDATADIR') + "/" + MCDataSet + "_Tree_" + "{:g}".format(binnings[binleft])  +'_' + "{:g}".format(binnings[binleft+1]) + ".root", "ExampleAnalysisTree", selection='Pattern=='+str(pattern))
            if events.shape[0] == 0:
                events = events_tem
            else:
                events = np.append(events,events_tem)

    elif RigidityRange == 'INTERMEDIATEENERGYRANGE':
            #events = root2array( resultpath + "/" + MCDataSet + ".root", "AntiprotonIntermediateEnergyTree" )
            events = root2array( resultpath + "/" + MCDataSet + ".root", "AntiprotonLowEnergyTree" )
    elif RigidityRange == 'LOWENERGYRANGE':
            events = root2array( resultpath + "/" + MCDataSet + "_Tree.root", "AntiprotonLowEnergyTree" )

    else:
        os._exit(0)


    #### Cut on MCmomentum; Apply Furthrer Cut 
    if energy_range == '14.1_450GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>14.1)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<450)[0]]
    elif energy_range == '14.1_525GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>14.1)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<525)[0]]
    elif energy_range == '2.97_18GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>2.97)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<18.0)[0]]
    elif energy_range == '1.33_5.9GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>1.33)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<5.9)[0]]
    elif energy_range == '0.8_18GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>0.8)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<18)[0]]
    elif energy_range == '1.0_5.9GeV':
        events = events[np.where(events["MCPrimaryMomentum"]>1.0)[0]]
        events = events[np.where(events["MCPrimaryMomentum"]<5.9)[0]]
        events = events[np.where( abs(events["Rigidity"])>1.0 )[0]]
        events = events[np.where( abs(events["Rigidity"])<5.9)[0]]
    
    if pattern == "0":
        if (EcalBDTCut):
            events = events[np.where((events["EcalBDT_EnergyD"]>-2) & (events["EcalBDT_EnergyD"] < EcalBDTCutvalue_select_proton))[0]]

    if (TrdProtonHeliumCut):
        events = events[np.where(events["TrdLogLikelihoodRatioProtonHeliumTracker"] < TrdProtonHeliumCutValue)[0]]

    if (PhysicstriggerCut):
        events = events[np.where((events["TriggerFlags"]&0x3e)>0)[0]]


    #### Fill in TH2D 
    for entryNum in range(0,events.shape[0]):
        Unfolding_Matrices.Fill(np.abs(events["Rigidity"][entryNum]), events["MCPrimaryMomentum"][entryNum]) ## reconstructed x, true y.


    #### Save Unfolding_Matrices root 
    if Cluster == 'HPC':
        if RigidityRange == 'HIGHENERGYRANGE':
            ROOTFile  = TFile( "/hpcwork/jara0052/sichen/Unfolding_Matrices/high/Unfolding_MatricesTH2D_fill_Pattern_" + str(pattern) + str("_") + str(arguments.binningversion) + ".root" , "RECREATE")
            Unfolding_Matrices.Write("Unfolding_Matrices")
            ROOTFile.Close()
        elif RigidityRange == 'INTERMEDIATEENERGYRANGE':
            ROOTFile  = TFile( "/hpcwork/jara0052/sichen/Unfolding_Matrices/intermediate/Unfolding_MatricesTH2D_fill.root", "RECREATE")
            Unfolding_Matrices.Write("Unfolding_Matrices")
            ROOTFile.Close()
        elif RigidityRange == 'LOWENERGYRANGE':
            ROOTFile  = TFile( "/hpcwork/jara0052/sichen/Unfolding_Matrices/low/Unfolding_MatricesTH2D_fill.root", "RECREATE")
            Unfolding_Matrices.Write("Unfolding_Matrices")
            ROOTFile.Close()

    elif Cluster == 'JUAMS':
        ROOTFile  = TFile("/p/scratch/cvsk10/li8/analysis_v7.0/Unfolding_MatricesTH2D_fill.root","RECREATE")
        Unfolding_Matrices.Write("Unfolding_Matrices")
        ROOTFile.Close()


if __name__ == '__main__':

    #### Parse Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern'     , help='which tracker patterns you choose')  #### For high rigidity range, it is a must. For low rigidity, it's ok not to give a value.
    parser.add_argument('--dataset', help='which dataset you choose')
    parser.add_argument('--energy' , help='which Rigidity Range you choose')
    parser.add_argument('--cluster', help='which cluster you are working on')

    parser.add_argument('--binningversion') #### Only for high rigidity range, for other ranges it's ok not to give a value.

    parser.add_argument('--ecalBDTCut', dest='ecalBDTCut', action='store_true',help='ecalBDTCut. Default:False')
    parser.add_argument('--no-ecalBDTCut', dest='ecalBDTCut', action='store_false',help='ecalBDTCut. Default:False')
    parser.set_defaults(ecalBDTCut=False)
    parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
    parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
    parser.set_defaults(trdProtonHeliumCut=False)
    parser.add_argument('--physicstriggerCut', dest='physicstriggerCut', action='store_true',help='physicstriggerCut. Default:False')
    parser.add_argument('--no-physicstriggerCut', dest='physicstriggerCut', action='store_false',help='physicstriggerCut. Default:False')
    parser.set_defaults(physicstriggerCut=False)

    arguments = parser.parse_args()

    if ( arguments.dataset != None and arguments.energy != None and arguments.cluster!= None ):
        pattern       = arguments.pattern
        MCDataSet     = arguments.dataset
        RigidityRange = arguments.energy
        Cluster       = arguments.cluster
    else:
        print("Some necessary arguments are missing, plese check !")
        os._exit(0)

    EcalBDTCut = arguments.ecalBDTCut
    TrdProtonHeliumCut = arguments.trdProtonHeliumCut
    PhysicstriggerCut = arguments.physicstriggerCut

    # WorkPath
    if RigidityRange == 'HIGHENERGYRANGE':
        resultpath = os.getenv('HPCHIGHENERGYDATADIR') + "/ISS_anylsis"
        highpath = os.getenv('HPCHIGHENERGYDATADIR')
    elif RigidityRange == 'INTERMEDIATEENERGYRANGE':
        #resultpath = os.getenv('HPCINTERMEDIATEDIR') + "/total"
        resultpath = os.getenv('HPCLOWENERGYDIR') + "/totalall"
    elif RigidityRange == 'LOWENERGYRANGE':
        resultpath = os.getenv('HPCLOWENERGYDIR') + "/totalall"

    # Cut Value
    TrdProtonHeliumCutValue = 0.3
    EcalBDTCutvalue_select_proton = 0.0

    # Binning
    if RigidityRange == 'HIGHENERGYRANGE':
        if arguments.binningversion == '450version':
            binnings = binning.published2016binnings[26:]
            energy_range = '14.1_450GeV'
            Unfolding_Matrices = TH2D("Unfolding_Matrices","", 31, binnings, 31,  binnings)
        elif arguments.binningversion == '525version':
            binnings = binning.Newbinnings_525_zhili[26:]
            energy_range = '14.1_525GeV'
            Unfolding_Matrices = TH2D("Unfolding_Matrices","", 32, binnings, 32, binnings)
    elif RigidityRange == 'INTERMEDIATEENERGYRANGE':
            binnings = binning.published2016binnings[9:30]
            energy_range = '2.97_18GeV'
            Unfolding_Matrices = TH2D("Unfolding_Matrices", "", 20, binnings, 20,  binnings)
    elif RigidityRange == 'LOWENERGYRANGE':
            binnings = binning.published2016binnings[:17]
            energy_range = '1.0_5.9GeV'
            Unfolding_Matrices = TH2D("Unfolding_Matrices", "", 16, binnings, 16,  binnings)

    main()




