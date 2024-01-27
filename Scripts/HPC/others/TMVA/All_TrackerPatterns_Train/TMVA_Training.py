#!/usr/bin/env python
import ROOT
from ROOT import TMVA, TFile, TCut, gROOT, TChain
import matplotlib.pyplot as plt
import numpy as np
from root_numpy import array2tree, tree2array
import os


def RawToReducedTree(RawTree, pattern):

    ## Cut on Tracker Patterns
    branchlist = ['TrdLogLikelihoodRatioElectronProtonTracker', 'RigidityAsymmetry', 'RigidityAsymmetryL9', 'Chi2TrackerYAsymmetry', 'InnerMaxSpanRigidityMatching', 'L1L9RigidityMatching', 'L24L58RigidityMatching', 'Log10Chi2TrackerXInner', 'Log10Chi2TrackerYInner', 'Log10Chi2TrackerX', 'Log10Chi2TrackerY', 'TrackerL58L24ChargeAsymmetry', 'TrackerL9Charge', 'TrackerL78Charge', 'UpperTofCharge', 'LowerTofCharge']
    ReducedArray = tree2array(RawTree, selection='Pattern=='+str(pattern), branches=branchlist)

    ## Multiply sign of rigidity to three MVA variables: "InnerMaxSpanRigidityMatching", "L1L9RigidityMatching", "L24L58RigidityMatching"
    RawArray_Rigidity = tree2array(RawTree, selection='Pattern=='+str(pattern), branches=['Rigidity'])
    RigiditySign = RawArray_Rigidity['Rigidity']/np.abs(RawArray_Rigidity['Rigidity'])
    ReducedArray['InnerMaxSpanRigidityMatching'] = ReducedArray['InnerMaxSpanRigidityMatching'] * RigiditySign
    ReducedArray['L1L9RigidityMatching']         = ReducedArray['L1L9RigidityMatching']         * RigiditySign
    ReducedArray['L24L58RigidityMatching']       = ReducedArray['L24L58RigidityMatching']       * RigiditySign

    ReducedTree = array2tree(ReducedArray, name='ExampleAnalysisTree')
    return ReducedTree


def main():

    ## Trainning Process
    for pattern in Pattern_all:

        for Binning_Index, Binning_Edge in enumerate(Binning_Edge_all):

            ## Load Signal and Background Events
            chain_Proton   = TChain('ExampleAnalysisTree', '')
            chain_CCproton = TChain('ExampleAnalysisTree', '')
            chain_Electron = TChain('ExampleAnalysisTree', '')
            chain_Proton  .AddFile( datapath + "/" + Correct_MC_name  + "_" + str(Binning_Edge) + ".root")
            chain_CCproton.AddFile( datapath + "/" + Confused_MC_name + "_" + str(Binning_Edge) + ".root")
            chain_Electron.AddFile( datapath + "/" + Electron_MC_name + "_" + str(Binning_Edge) + ".root")
            print("Now the Binning_Edge is " + str(Binning_Edge))
            print("Proton Number:"   + str(chain_Proton.GetEntries()))
            print("CCProton Number:" + str(chain_CCproton.GetEntries()))
            print("Electron Number:" + str(chain_Electron.GetEntries()))
            print('\n')
            gROOT.cd() 
            SignalTree_Proton       = chain_Proton.GetTree()
            SignalTree_Electron     = chain_Electron.GetTree()
            BackgroundTree_CCProton = chain_CCproton.GetTree()

            ## Get MVA Variables and Cut on Tracker Patterns 
            SignalReducedTree_Proton       = RawToReducedTree(SignalTree_Proton      , pattern)
            SignalReducedTree_Electron     = RawToReducedTree(SignalTree_Electron    , pattern)
            BackgroundReducedTree_CCProton = RawToReducedTree(BackgroundTree_CCProton, pattern)

            ## Add MVA variables to TMVA dataloader
            if Binning_Edge == "14.1_80.5":
                dataloader = ROOT.TMVA.DataLoader( 'dataset_pymva_' + str("14.1_80.5") + "_Pattern_" + str(pattern)  )
            else:
                dataloader = ROOT.TMVA.DataLoader( 'dataset_pymva_' + str(Binning_Edge) + "_Pattern_" + str(pattern)  )
            for branch in SignalReducedTree_Proton.GetListOfBranches():
                dataloader.AddVariable(branch.GetName())

            ## Add Training data to TMVA dataloader
            dataloader.AddSignalTree(SignalReducedTree_Proton          , 0.8 * 1/SignalReducedTree_Proton.GetEntries()      ) #(TTree, weight)
            dataloader.AddSignalTree(SignalReducedTree_Electron        , 0.2 * 1/SignalReducedTree_Electron.GetEntries()    )
            dataloader.AddBackgroundTree(BackgroundReducedTree_CCProton, 1.0 * 1/BackgroundReducedTree_CCProton.GetEntries())

            ## Prepare
            trainTestSplit = 0.7
            dataloader.PrepareTrainingAndTestTree(ROOT.TCut(''),
                    'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
                    'TrainTestSplit_Background={}:'.format(trainTestSplit)+\
                    'SplitMode=Random')

            ## Setup TMVA Factory
            ROOT.TMVA.Tools.Instance()
            if Binning_Edge == "14.1_80.5":
                outputFile = ROOT.TFile.Open('TMVAOutputPyMVA_' + str("14.1_80.5") + "_Pattern_" + str(pattern) + '.root', 'RECREATE')
            else:
                outputFile = ROOT.TFile.Open('TMVAOutputPyMVA_' + str(Binning_Edge) + "_Pattern_" + str(pattern) + '.root', 'RECREATE')
            factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
                    '!V:!Silent:Color:DrawProgressBar:Transformations=I,G:'+\
                    'AnalysisType=Classification')

            ## Book method
            factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=50:MaxDepth=3:IgnoreNegWeightsInTraining")

            ## Training
            factory.TrainAllMethods()
            factory.TestAllMethods()
            factory.EvaluateAllMethods()

            ## Print ROC plot (Enable Javascript for ROOT so that we can draw the canvas)
            #%jsroot on
            canvas = factory.GetROCCurve(dataloader)
            canvas.Draw()
            if Binning_Edge == "14.1_80.5":
                canvas.Print("ROC_Curve_" + str("14.1_80.5") + "_Pattern_" + str(pattern) + ".png")
            else:
                canvas.Print("ROC_Curve_" + str(Binning_Edge) + "_Pattern_" + str(pattern) + ".png")


if __name__ == "__main__":

    # Rigidity Binnings
    Binning_Edge_all = ["14.1_80.5", '80.5_93', '93_108', '108_125', '125_147', '147_175', '175_211',
       '211_259', "259_450"]  # Binning Used in 2016 AMS Pbar PRL Paper.
    Binning_Edge_all = np.append(Binning_Edge_all, ["211_250", "250_330", "330_525"]) # Updated binning used in 2021 AMS Physics Report.

    # Tracker Patterns
    Pattern_all = [-1, 0, 1, 2, 3, 4, 5]

    # MC dataset For Training
    Correct_MC_name  = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive'
    Confused_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative'
    Electron_MC_name = 'B1091_el.pl1.0_25_2000_7.6_all_Tree'

    datapath = os.getenv('HPCHIGHENERGYDATADIR')

    main()



