import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString
import matplotlib.pyplot as plt
import numpy as np
import array 
from root_numpy import array2tree, array2root, tree2array
import binning
import CNN_models


def TMVAEvaluate(InputArray, InputArray_Weight, InputArray_MCPrimaryMomentum, Pattern, RigidityBin, InputSpecies, datapath):

    v1 = array.array('f',[0])
    v2 = array.array('f',[0])
    v3 = array.array('f',[0])
    v4 = array.array('f',[0])
    v5 = array.array('f',[0])
    v6 = array.array('f',[0])
    v7 = array.array('f',[0])
    v8 = array.array('f',[0])
    v9 = array.array('f',[0])
    v10 = array.array('f',[0])
    v11 = array.array('f',[0])
    v12 = array.array('f',[0])
    v13 = array.array('f',[0])
    v14 = array.array('f',[0])
    v15 = array.array('f',[0])
    v16 = array.array('f',[0])

    reader = ROOT.TMVA.Reader()

    reader.AddVariable("TrdLogLikelihoodRatioElectronProtonTracker", v1)
    reader.AddVariable("RigidityAsymmetry"                         , v2)
    reader.AddVariable("RigidityAsymmetryL9"                       , v3)
    reader.AddVariable("Chi2TrackerYAsymmetry"                     , v4)
    reader.AddVariable("InnerMaxSpanRigidityMatching"              , v5)
    reader.AddVariable("L1L9RigidityMatching"                      , v6)
    reader.AddVariable("L24L58RigidityMatching"                    , v7)
    reader.AddVariable("Log10Chi2TrackerXInner"                    , v8)
    reader.AddVariable("Log10Chi2TrackerYInner"                    , v9)
    reader.AddVariable("Log10Chi2TrackerX"                         , v10)
    reader.AddVariable("Log10Chi2TrackerY"                         , v11)
    reader.AddVariable("TrackerL58L24ChargeAsymmetry"              , v12)
    reader.AddVariable("TrackerL9Charge"                           , v13)
    reader.AddVariable("TrackerL78Charge"                          , v14)
    reader.AddVariable("UpperTofCharge"                            , v15)
    reader.AddVariable("LowerTofCharge"                            , v16)

    LowMergedList = binning.bins[0:23]
    if RigidityBin in LowMergedList:
        reader.BookMVA("BDTG","/home/bo791269/CCTMVA_Weights/dataset_pymva_" + str("14.1_80.5") + "_Pattern_" + str(Pattern) + "/weights/TMVAClassification_BDTG.weights.xml")
    else:
        reader.BookMVA("BDTG","/home/bo791269/CCTMVA_Weights/dataset_pymva_" + str(RigidityBin) + "_Pattern_" + str(Pattern) + "/weights/TMVAClassification_BDTG.weights.xml")

    CCEstimator = []

    #### Evaluate MVA
    for ievt in range(InputArray.shape[0]):
        v1[0]  = InputArray['TrdLogLikelihoodRatioElectronProtonTracker'][ievt][0]
        v2[0]  = InputArray['RigidityAsymmetry'][ievt][0]
        v3[0]  = InputArray['RigidityAsymmetryL9'][ievt][0]
        v4[0]  = InputArray['Chi2TrackerYAsymmetry'][ievt][0]
        v5[0]  = InputArray['InnerMaxSpanRigidityMatching'][ievt][0]
        v6[0]  = InputArray['L1L9RigidityMatching'][ievt][0]
        v7[0]  = InputArray['L24L58RigidityMatching'][ievt][0]
        v8[0]  = InputArray['Log10Chi2TrackerXInner'][ievt][0]
        v9[0]  = InputArray['Log10Chi2TrackerYInner'][ievt][0]
        v10[0] = InputArray['Log10Chi2TrackerX'][ievt][0]
        v11[0] = InputArray['Log10Chi2TrackerY'][ievt][0]
        v12[0] = InputArray['TrackerL58L24ChargeAsymmetry'][ievt][0]
        v13[0] = InputArray['TrackerL9Charge'][ievt][0]
        v14[0] = InputArray['TrackerL78Charge'][ievt][0]
        v15[0] = InputArray['UpperTofCharge'][ievt][0]
        v16[0] = InputArray['LowerTofCharge'][ievt][0]

        mva_value = reader.EvaluateMVA("BDTG")
        mva_value = (mva_value + 1)/2 #change from -1 to 1 to 0 to 1.

        CCEstimator.append(mva_value)
    
  
    CCEstimator = np.array(CCEstimator)
    
    np.save(datapath + '/ISS_anylsis/data/' + InputSpecies + str(RigidityBin) + "_Pattern_" + str(Pattern) + ".npy", np.column_stack(( CCEstimator, InputArray['TrdLogLikelihoodRatioElectronProtonTracker'][:,0], InputArray_Weight, InputArray_MCPrimaryMomentum )) )


def VGG16NeuralNetworkEvaluate(InputArray, InputArray_Weight, InputArray_MCPrimaryMomentum, RigidityBin, InputSpecies, datapath):

        #### Transpose and add a dimention
        InputArray = np.transpose(InputArray)
        InputArray = np.expand_dims(InputArray, axis=2)

        ####  model prediction
        import CNN_models
        model = CNN_models.VGG16(16)

        if RigidityBin == "80.5_93.0" or RigidityBin == "93.0_108.0" or RigidityBin == "108.0_125.0" or RigidityBin == "125.0_147.0" or RigidityBin == "147_175" or RigidityBin == "175_211" or RigidityBin == "211_259" or RigidityBin == "259_330" or RigidityBin == "330_525":
            model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+RigidityBin+'_'+ 'mva_variables' +'.h5')
        elif RigidityBin == "259_450":
            model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+'259_330'+'_'+ 'mva_variables' +'.h5')
        elif RigidityBin == "211_250":
            model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+'211_259'+'_'+ 'mva_variables' +'.h5')
        elif RigidityBin == "250_330":
            model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+'259_330'+'_'+ 'mva_variables' +'.h5')
        else:
            model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+'16.6_80.5'+'_'+ 'mva_variables' +'.h5')

        y_pred = model.predict(InputArray)

        #### Fill histogram with weights
        #Confused_MC_hist = TH2D("Confused_MC_hist","Confused_MC_hist", 31, binnings, 31,  binnings)
        #Confused_MC_hist = TH1D("Confused_MC_hist","Confused_MC_hist", 31, 0, 1 )
        #fill_hist( Confused_MC_hist, y_pred[Correct_MC.shape[0]:(Correct_MC.shape[0]+Confused_MC.shape[0]),1], weights=Confused_MC['Weight'])
        #### Convert histogram back to numpy array
        #Confused_MC_cc = hist2array( Confused_MC_hist)

        #print("D1:")
        #print(y_pred.shape)
        #print(y_pred)
        #print("D2:")
        #print(InputArray.shape)
        #print(InputArray[:,0,0])

        ''' 
        #if y_pred[:,1].shape[0]>1:
        print("y_pred[:,1]:")
        print(y_pred)
        print(y_pred[:,1])
        print(y_pred[:,1].shape)

        print("InputArray[:,0,0]:")
        print(InputArray[:,0,0])
        print(InputArray[:,0,0].shape)

        print("InputArray_Weight:")
        print(InputArray_Weight.shape)
        print(InputArray_Weight)

        print("InputArray_MCPrimaryMomentum")     
        print(InputArray_MCPrimaryMomentum.shape) 
        print(InputArray_MCPrimaryMomentum)
        '''


        np.save(datapath + '/ISS_anylsis/data/' + InputSpecies + str(RigidityBin) + "_Pattern_0" + "_VGG16NN" + ".npy", np.transpose( [y_pred[:,1], InputArray[:,0,0], InputArray_Weight, InputArray_MCPrimaryMomentum ] ))

        '''
        #### save Template and ISS data (TRDLikelihood and CCestimator)
        np.save(datapath + '/ISS_anylsis/data/plot_ISS_positive' + suffix + 'rigidity' + RigidityBin + '_pass7.8.npy', np.transpose( [y_pred_ISS[0:ISSpositive.shape[0],1], ISSpositive['TrdLogLikelihoodRatioElectronProtonTracker'], ISSpositive['Weight'] ] ))
        np.save(datapath + '/ISS_anylsis/data/plot_ISS_negative' + suffix + 'rigidity'+ RigidityBin + '_pass7.8.npy', np.transpose( [y_pred_ISS[ISSpositive.shape[0]:(ISSpositive.shape[0]+ISSnegative.shape[0]),1], ISSnegative['TrdLogLikelihoodRatioElectronProtonTracker'], ISSnegative['Weight'] ] ))
        np.save(datapath + '/ISS_anylsis/data/ElectronTemplateData'+ suffix + RigidityBin + '.npy', np.transpose( [y_pred_ElectronTemplateData[:,1], Electron_Template_Data['TrdLogLikelihoodRatioElectronProtonTracker'], Electron_Template_Data['Weight'] ] ))

        np.save(datapath + '/ISS_anylsis/data/' + 'MCChargeCorrect_' + Correct_MC_name + '_' + RigidityBin + '.npy', np.transpose( [y_pred[0:Correct_MC.shape[0],1], Correct_MC['TrdLogLikelihoodRatioElectronProtonTracker'], Correct_MC['Weight'] ] ))
        np.save(datapath + '/ISS_anylsis/data/' + 'MCChargeConfused_' + Confused_MC_name + '_' + RigidityBin + '.npy', np.transpose( [ y_pred[Correct_MC.shape[0]:(Correct_MC.shape[0]+Confused_MC.shape[0]),1], Confused_MC['TrdLogLikelihoodRatioElectronProtonTracker'], Confused_MC['Weight'] ] ))
        np.save(datapath + '/ISS_anylsis/data/MCelectron_negative' + RigidityBin + '.npy', np.transpose( [y_pred[(Correct_MC.shape[0]+Confused_MC.shape[0]):,1], Electron_MC['TrdLogLikelihoodRatioElectronProtonTracker'], Electron_MC['Weight']] ))
        '''

