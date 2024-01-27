#!/usr/bin/env python
import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString, TBranch, gROOT, TChain
import matplotlib.pyplot as plt
import numpy as np
import binning
from root_numpy import array2tree, array2root, tree2array, root2array
import os
import TMVA_Application_model

def RawToReducedTree(RawTree, pattern):

    RawArray_Rigidity = tree2array(RawTree, selection='Pattern=='+str(pattern), branches=['Rigidity'])
    Sign = RawArray_Rigidity['Rigidity']/np.abs(RawArray_Rigidity['Rigidity'])

    branchlist = ['TrdLogLikelihoodRatioElectronProtonTracker', 'RigidityAsymmetry', 'RigidityAsymmetryL9', 'Chi2TrackerYAsymmetry', 'InnerMaxSpanRigidityMatching', 'L1L9RigidityMatching', 'L24L58RigidityMatching', 'Log10Chi2TrackerXInner', 'Log10Chi2TrackerYInner', 'Log10Chi2TrackerX', 'Log10Chi2TrackerY', 'TrackerL58L24ChargeAsymmetry', 'TrackerL9Charge', 'TrackerL78Charge', 'UpperTofCharge', 'LowerTofCharge']
    ReducedArray = tree2array(RawTree, selection='Pattern=='+str(pattern), branches=branchlist)

    ReducedArray['InnerMaxSpanRigidityMatching'] = ReducedArray['InnerMaxSpanRigidityMatching'] * Sign
    ReducedArray['L1L9RigidityMatching'] = ReducedArray['L1L9RigidityMatching'] * Sign
    ReducedArray['L24L58RigidityMatching'] = ReducedArray['L24L58RigidityMatching'] * Sign

    ReducedTree = array2tree(ReducedArray, name='ExampleAnalysisTree')
    return ReducedTree


def main():
    #### Parameters
    Binning_Edge_all = ["14.1_80.5"]
    Binning_Edge_all = np.append(Binning_Edge_all, binning.bins[23:])
    Binning_Edge_all = np.append(Binning_Edge_all, ["211_250", "250_330", "259_450"])
    #Binning_Edge_all = ["211_250", "250_330", "259_450"]

    #Pattern_all = [-1, 0, 1, 2, 3, 4, 5]
    Pattern_all = [3]

    datapath = os.getenv('HPCHIGHENERGYDATADIR')

    #Correct_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive'                    # Large
    Correct_MC_name = 'B1042_antipr.pl1.1800_7.6_all_Tree'                                       # Default
    #Correct_MC_name = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_positive'          # Test
    #Confused_MC_name = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_negative'          # Test
    Confused_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative'                     # Default
    Electron_MC_name = 'B1091_el.pl1.0_25_2000_7.6_all_Tree'                                      # Default

    #### Trainning Process
    for pattern in Pattern_all:
        for Binning_Index, Binning_Edge in enumerate(Binning_Edge_all):

            chain_Proton   = TChain('ExampleAnalysisTree','')
            chain_CCproton = TChain('ExampleAnalysisTree','')
            chain_Electron = TChain('ExampleAnalysisTree','')

            chain_Proton.AddFile( datapath + "/" + Correct_MC_name + "_" + str(Binning_Edge) + ".root")
            chain_CCproton.AddFile( datapath + "/" + Confused_MC_name + "_" + str(Binning_Edge) + ".root")
            chain_Electron.AddFile( datapath + "/" + Electron_MC_name + "_" + str(Binning_Edge) + ".root")

            print("Now the Binning_Edge is " + str(Binning_Edge))
            print(chain_Proton.GetEntries())
            print(chain_CCproton.GetEntries())
            print(chain_Electron.GetEntries())
            print('\n')

            #Proton   = TFile.Open("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_" + str(Binning_Edge) + ".root") 
            #CCproton = TFile.Open("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + str(Binning_Edge) + ".root") 
            #Electron = TFile.Open("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/B1091_el.pl1.0_25_2000_7.6_all_Tree_" + str(Binning_Edge) + ".root") 
            gROOT.cd() 

            # Get signal and background trees from file
            #background1RawTree = CCproton.Get('ExampleAnalysisTree') 
            #signal1RawTree     = Proton.Get('ExampleAnalysisTree') 
            #signal2RawTree     = Electron.Get('ExampleAnalysisTree') 
            signal1RawTree     = chain_Proton.GetTree()
            background1RawTree = chain_CCproton.GetTree()
            signal2RawTree     = chain_Electron.GetTree()

            # Get MVA Variables and Cut on Tracker Patterens 
            signal1ReducedTree     = RawToReducedTree(signal1RawTree    , pattern)
            signal2ReducedTree     = RawToReducedTree(signal2RawTree    , pattern)
            background1ReducedTree = RawToReducedTree(background1RawTree, pattern)

            print("Check to make sure the merged trees are looded:")
            print(signal1ReducedTree.GetEntries())
            print(signal2ReducedTree.GetEntries())
            print(background1ReducedTree.GetEntries())

            # Take same numbers for all samples
            #MinNumber = min(signal1ReducedTree.GetEntries(), signal2ReducedTree.GetEntries(), background1ReducedTree.GetEntries())
            #signal1ReducedTree.SetEntries(MinNumber)
            #signal2ReducedTree.SetEntries(MinNumber)
            #background1ReducedTree.SetEntries(MinNumber)

            # Add variables to dataloader
            if Binning_Edge == "14.1_80.5":
                dataloader = ROOT.TMVA.DataLoader( 'dataset_pymva_' + str("14.1_80.5") + "_Pattern_" + str(pattern)  )
            else:
                dataloader = ROOT.TMVA.DataLoader( 'dataset_pymva_' + str(Binning_Edge) + "_Pattern_" + str(pattern)  )

            for branch in signal1ReducedTree.GetListOfBranches():
                dataloader.AddVariable(branch.GetName())

            # Add trees to dataloader
            dataloader.AddSignalTree(signal1ReducedTree, 0.8 * 1/signal1ReducedTree.GetEntries()) #(TTree, weight)
            dataloader.AddSignalTree(signal2ReducedTree, 0.2 * 1/signal2ReducedTree.GetEntries())
            dataloader.AddBackgroundTree(background1ReducedTree, 1.0 * 1/background1ReducedTree.GetEntries())

            trainTestSplit = 0.8
            dataloader.PrepareTrainingAndTestTree(ROOT.TCut(''),
                    'TrainTestSplit_Signal={}:'.format(trainTestSplit)+\
                    'TrainTestSplit_Background={}:'.format(trainTestSplit)+\
                    'SplitMode=Random')

            # Setup TMVA
            ROOT.TMVA.Tools.Instance()
            if Binning_Edge == "14.1_80.5":
                outputFile = ROOT.TFile.Open('TMVAOutputPyMVA_' + str("14.1_80.5") + "_Pattern_" + str(pattern) + '.root', 'RECREATE')
            else:
                outputFile = ROOT.TFile.Open('TMVAOutputPyMVA_' + str(Binning_Edge) + "_Pattern_" + str(pattern) + '.root', 'RECREATE')
            factory = ROOT.TMVA.Factory('TMVAClassification', outputFile,
                    '!V:!Silent:Color:DrawProgressBar:Transformations=I,G:'+\
                    'AnalysisType=Classification')

            # Book method
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" )
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG","!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" )
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTB","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" )
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTD","!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" )
            ## Used
            factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=50:MaxDepth=3:IgnoreNegWeightsInTraining")
            #factory.BookMethod(dataloader, TMVA.Types.kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" )
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDT","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" )
            #factory.BookMethod(dataloader, TMVA.Types.kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" )
            #factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "Likelihood","H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
            #factory.BookMethod( dataloader, TMVA.Types.kLikelihood, "LikelihoodKDE","!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" )
            #factory.BookMethod( dataloader, TMVA.Types.kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" )
            ## Fabian
            #factory.BookMethod(dataloader,TMVA.Types.kBDT, "BDT", "!H:!V:NTrees=300:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:IgnoreNegWeightsInTraining")
            #factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTG", "!H:!V:NTrees=200:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=20:MaxDepth=2:IgnoreNegWeightsInTraining")

            factory.TrainAllMethods()
            factory.TestAllMethods()
            factory.EvaluateAllMethods()

            # Print ROC plot (Enable Javascript for ROOT so that we can draw the canvas)
            #%jsroot on
            canvas = factory.GetROCCurve(dataloader)
            canvas.Draw()
            if Binning_Edge == "14.1_80.5":
                canvas.Print("ROC_Curve_" + str("14.1_80.5") + "_Pattern_" + str(pattern) + ".png")
            else:
                canvas.Print("ROC_Curve_" + str(Binning_Edge) + "_Pattern_" + str(pattern) + ".png")


if __name__ == "__main__":
    main()



