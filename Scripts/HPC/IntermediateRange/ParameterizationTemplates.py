#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F
import uproot
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import binning 
import argparse
import math
from iminuit import Minuit
from iminuit.util import describe, make_func_code
from iminuit.cost import LeastSquares
import uncertainties
import ParameterizationTemplates_function
import ParameterizationTemplates_tool
import ParameterizationTemplates_InitialValues
import ParameterizationTemplates_Plot 
import ParameterizationTemplates_Fit
import os 

def main():

    ParameterError_TRD_all = []
    ParameterError_TOF_all = []

    #### Loop in Rigidity
    for index, binname in enumerate(UsedBinningEdge): ## Don't change to part loop, otherwise the index would be wrong. or start from beginning.
    #for index, binname in enumerate(UsedBinningEdge[0:1]):

        print("\033[1;35mNow Rigidity is \033[0m" + str(binname))
        print("index is:" + str(index))

        for splitindex in range(0, SplitTotal, mergestep):  ## If time dependent, loop in time index, for time averaged, it's 1.
            print("splitindex" + str(splitindex))
            #### Load data
            if TimeMode == "TimeAveraged": 
                if RigidityRange == 'intermediate':
                    f_ROOT_ParametrilisedTemplates = TFile( IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/RigidityRootFiles/averaged_ratio_" + binname + "_" + ISSversion + "_new.root", "RECREATE")
                    f_uproot_OriginalTemplate = uproot.open( IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/RigidityRootFiles/averaged_ratio_" + binname + "_" + ISSversion + ".root") 
                    dic = {'data_pass7_positive':'Antiproton', 'template_electron_Data':'Electron'}
                elif RigidityRange == 'low':
                    f_ROOT_ParametrilisedTemplates = TFile( LowPath  + "/totalall/Time_Averaged_ratio_Low/binmerge1/tof/averaged_ratio_" + binname + "_" + ISSversion + "_new.root", "RECREATE")
                    f_uproot_OriginalTemplate = uproot.open( LowPath + "/totalall/Time_Averaged_ratio_Low/binmerge1/tof/averaged_ratio_" + binname + "_" + ISSversion + ".root")
                    dic = {'TofTRD_data_pass7_positive_template':'Antiproton', 'TofTRD_template_ElectronData':'Electron', 'TofTRD_template_PionData':'Pion'}
            elif TimeMode == "TimeDependent":
                if RigidityRange == 'intermediate':
                    dic = {'data_pass7_positive':'Antiproton', 'template_electron_Data':'Electron'}
                    print("in progress")
                elif RigidityRange == 'low':
                    if TimeSplitMode == "3BartalRotation": 
                        f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_3BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + "_new.root", "RECREATE")
                        f_uproot_OriginalTemplate = uproot.open(LowPath + "/totalall/rootfiles_3BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root")
                    elif TimeSplitMode == "6BartalRotation":
                        f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex)+ "_new.root", "RECREATE")
                        f_uproot_OriginalTemplate = uproot.open(LowPath + "/totalall/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root")
                        print('load data')
                    elif TimeSplitMode == "6Months":
                        f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_6months_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + "_new.root", "RECREATE")
                        f_uproot_OriginalTemplate = uproot.open(LowPath + "/totalall/rootfiles_6months_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root")
                    dic = {'TofTRD_data_pass7_positive_template':'Antiproton', 'TofTRD_template_ElectronData':'Electron', 'TofTRD_template_PionData':'Pion'}


            #### Loop for all templates
            for key, value in dic.items():
                print("\033[1;36mLoad Template is :\033[0m"   + str(value))
                print("\033[1;36mLoad Template Histogram is :\033[0m" + str(key))


                #### Get Original Values
                if RigidityRange == 'intermediate':
                    x_Edge_trd, x_Center_trd, y_Value_trd, y_Error_trd, totalnumbers, x_More_trd = ParameterizationTemplates_tool.GetValues(RigidityRange, f_uproot_OriginalTemplate, key)
                    # Naive Fix for bin content is 0
                    ParameterizationTemplates_InitialValues.FixForNaf( {"variable":[y_Value_trd], "error":[y_Error_trd]} )
                elif RigidityRange == 'low':
                    x_Edge_trd, x_Center_trd, y_Value_trd, y_Error_trd, totalnumbers, x_More_trd, x_Edge_tof, x_Center_tof, y_Value_tof, y_Error_tof, x_More_tof = ParameterizationTemplates_tool.GetValues(RigidityRange, f_uproot_OriginalTemplate, key)
                    # Naive Fix for bin content is 0
                    ParameterizationTemplates_InitialValues.FixForNaf( {"variable":[y_Value_trd, y_Value_tof], "error":[y_Error_trd, y_Error_tof]} )


                #### Fit the histogram
                ## (1). TRD Part
                ## Used fit function choice
                TrdFitFunction = ParameterizationTemplates_tool.Trd_UsedFitFunction(value, RigidityRange)
                Dimension = "TRD"
                ## Initial Values and Fit
                if TrdFitFunction == "OneGaussOneNovo":
                    UsedFunction_TRD = ParameterizationTemplates_function.FitFunction
                    a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0 = ParameterizationTemplates_InitialValues.TRD_InitialValues_OneGaussOneNovo(value, index, RigidityRange)
                    if RigidityRange == 'low':
                        if value == 'Pion':
                            Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.MinuitFit_OneGaussOneNovo(Dimension, RigidityRange, value, TrdFitFunction, x_Center_trd, y_Value_trd, y_Error_trd, a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index)
                    else:
                        Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.CurveFit_OneGaussOneNovo(Dimension, RigidityRange, value, x_Center_trd, y_Value_trd, y_Error_trd, a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index)
                elif TrdFitFunction == "OneGauss":
                    UsedFunction_TRD = ParameterizationTemplates_function.Gauss_function
                    a_0, mean_gauss_0, sigma_gauss_0 = ParameterizationTemplates_InitialValues.TRD_InitialValues_OneGauss(value, index, RigidityRange)
                    Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.CurveFit_OneGauss       (Dimension, RigidityRange, value, x_Center_trd, y_Value_trd, y_Error_trd, a_0, mean_gauss_0, sigma_gauss_0, index)
                    #Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.MinuitFit_OneGauss      (Dimension, RigidityRange, value, TrdFitFunction, x_Center_trd, y_Value_trd, y_Error_trd, a_0, mean_gauss_0, sigma_gauss_0, index)
                elif TrdFitFunction == "ExponentialFun":
                    UsedFunction_TRD = ParameterizationTemplates_function.ExponentialFun
                    k1_0, k2_0, b_0, c_0 = ParameterizationTemplates_InitialValues.TRD_InitialValues_ExponentialFun(value, index, RigidityRange)
                    #Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.CurveFit_ExponentialFun (Dimension, RigidityRange, value, TrdFitFunction, x_Center_trd, y_Value_trd, y_Error_trd, k1_0, k2_0, b_0, c_0, index)
                    Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.MinuitFit_ExponentialFun (Dimension, RigidityRange, value, TrdFitFunction, x_Center_trd, y_Value_trd, y_Error_trd, k1_0, k2_0, b_0, c_0, index)
                elif TrdFitFunction == "OneNovo":
                    UsedFunction_TRD = ParameterizationTemplates_function.Novosibirsk_function
                    b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0 = ParameterizationTemplates_InitialValues.TRD_InitialValues_OneNovo(value, index, RigidityRange)
                    Para_TRD, Covan_TRD, ParameterError_TRD = ParameterizationTemplates_Fit.CurveFit_OneNovo(Dimension, RigidityRange, value, x_Center_trd, y_Value_trd, y_Error_trd, b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index)
                ## (2). TOF Part
                if RigidityRange == 'low':
                    ## Used fit function choice
                    TofFitFunction = ParameterizationTemplates_tool.Tof_UsedFitFunction(value, RigidityRange) 
                    Dimension = "TOF"
                    ## Initial Values and Fit
                    if TofFitFunction == "OneGauss":
                        UsedFunction_TOF = ParameterizationTemplates_function.Gauss_function
                        a_0, mean_gauss_0, sigma_gauss_0 = ParameterizationTemplates_InitialValues.TOF_InitialValues_OneGauss(value, index, RigidityRange)
                        Para_TOF, Covan_TOF, ParameterError_TOF = ParameterizationTemplates_Fit.CurveFit_OneGauss       (Dimension, RigidityRange, value, x_Center_tof, y_Value_tof, y_Error_tof, a_0, mean_gauss_0, sigma_gauss_0, index)
                    elif TofFitFunction == "OneNovo":
                        UsedFunction_TOF = ParameterizationTemplates_function.Novosibirsk_function
                        b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0 = ParameterizationTemplates_InitialValues.TOF_InitialValues_OneNovo(value, index, RigidityRange)
                        Para_TOF, Covan_TOF, ParameterError_TOF = ParameterizationTemplates_Fit.CurveFit_OneNovo(Dimension, RigidityRange, value, x_Center_tof, y_Value_tof, y_Error_tof, b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index)


                ## Uncertainty
                # (1). TRD Part
                Para_TRD_uncertainty = Para_TRD - ParameterError_TRD
                print("Para_TRD:" + str(Para_TRD))
                print("Para_TRD (new) is " + str(Para_TRD_uncertainty))
                # (2). TOF Part
                if RigidityRange == 'low':
                    Para_TOF_uncertainty = Para_TOF - ParameterError_TOF
                    print("Para_TOF:" + str(Para_TOF))
                    print("Para_TOF (new) is " + str(Para_TOF_uncertainty))

                
                #### Gaussion Distribution
                if ErrorParametrizationMode == "Yes":
                    NewParameterError_TRDSet = np.load('ParameterError_TRD' + str(value) + '_new.npy')
                    ParameterError_TRD = NewParameterError_TRDSet[index,:]  
                    if RigidityRange == 'low':
                        NewParameterError_TOFSet = np.load('ParameterError_TOF' + str(value) + '_new.npy')
                        ParameterError_TOF = NewParameterError_TOFSet[index,:]
                # (1). Generate parameters
                parameter_set_trd      = np.zeros((len(Para_TRD), GenerateNumber))
                for i in range(len(Para_TRD)):
                    parameter_set_trd[i,:] = np.random.normal(Para_TRD[i], ParameterError_TRD[i], GenerateNumber)
                if RigidityRange == 'low':
                    parameter_set_tof = np.zeros((len(Para_TOF), GenerateNumber))
                    for i in range(len(Para_TOF)):
                        parameter_set_tof[i,:] = np.random.normal(Para_TOF[i], ParameterError_TOF[i], GenerateNumber)
                # (2). Plot Generated Parameters with Gaussion assumption
                if TimeMode == "TimeAveraged":
                    ParameterizationTemplates_Plot.Plot_GeneratedParameters(parameter_set_trd, "parameterset_TRD", value, binname, RigidityRange)
                    if RigidityRange == 'low':
                        ParameterizationTemplates_Plot.Plot_GeneratedParameters(parameter_set_tof, "parameterset_TOF", value, binname, RigidityRange)

                for j in range(GenerateNumber):
                    # (3). Get Fit Curve for randomlized template
                    y_trdfit_value_random, y_trdfit_error_random, y_trdfit_value_uncer_random, y_trdfit_error_uncer_random, y_trdfit_value_more_random, y_trdfit_error_more_random, y_trdfit_value_uncer_more_random, y_trdfit_error_uncer_more_random = ParameterizationTemplates_tool.GetFitCurve_TRD(RigidityRange, UsedFunction_TRD, x_Center_trd, Para_TRD, parameter_set_trd[:,j], x_More_trd, y_Value_trd, y_Error_trd, index)
                    if RigidityRange == 'low':
                        y_toffit_value_random, y_toffit_error_random, y_toffit_value_uncer_random, y_toffit_error_uncer_random, y_toffit_value_more_random, y_toffit_error_more_random, y_toffit_value_uncer_more_random, y_toffit_error_uncer_more_random = ParameterizationTemplates_tool.GetFitCurve_TOF(RigidityRange, UsedFunction_TOF, x_Center_tof, Para_TOF, parameter_set_tof[:,j], x_More_tof, y_Value_tof, y_Error_tof) 
                    # (4). Convert to histogram and Save HistogramTemplate into ROOT file
                    if RigidityRange == 'intermediate':
                        ParameterizationTemplates_tool.MakeAndSaveHistogramTemplate_Intermediate(RigidityRange, index, key, x_Center_trd, x_Edge_trd, y_Value_trd, y_Error_trd, y_trdfit_value_uncer_random, y_trdfit_error_uncer_random, totalnumbers, "randomlized", j )
                    elif RigidityRange == 'low':
                        ParameterizationTemplates_tool.MakeAndSaveHistogramTemplate_Low         (RigidityRange, index, key, x_Center_trd, x_Edge_trd, y_Value_trd, y_Error_trd, y_trdfit_value_uncer_random, y_trdfit_error_uncer_random, totalnumbers, x_Center_tof, x_Edge_tof, y_Value_tof, y_Error_tof, y_toffit_value_uncer_random, y_toffit_error_uncer_random, "randomlized", j )
                #(5). Plot Genereted Templates
                if TimeMode == "TimeAveraged":
                    ParameterizationTemplates_Plot.Plot_GeneratedMultiTemplates(x_Center_trd, x_More_trd, y_Value_trd, parameter_set_trd, value, binname, RigidityRange, UsedFunction_TRD, GenerateNumber, "TRD")
                    if RigidityRange == 'low':
                        ParameterizationTemplates_Plot.Plot_GeneratedMultiTemplates(x_Center_tof, x_More_tof, y_Value_tof, parameter_set_tof, value, binname, RigidityRange, UsedFunction_TOF, GenerateNumber, "TOF")
                

                #### Get Fit Curve
                ## (1). TRD dimension:
                y_trdfit_value, y_trdfit_error, y_trdfit_value_uncer, y_trdfit_error_uncer, y_trdfit_value_more, y_trdfit_error_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more = ParameterizationTemplates_tool.GetFitCurve_TRD(RigidityRange, UsedFunction_TRD, x_Center_trd, Para_TRD, Para_TRD_uncertainty, x_More_trd, y_Value_trd, y_Error_trd, index)
                ## (2). TOF dimension:
                if RigidityRange == 'low':
                    y_toffit_value, y_toffit_error, y_toffit_value_uncer, y_toffit_error_uncer, y_toffit_value_more, y_toffit_error_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more = ParameterizationTemplates_tool.GetFitCurve_TOF(RigidityRange, UsedFunction_TOF, x_Center_tof, Para_TOF, Para_TOF_uncertainty, x_More_tof, y_Value_tof, y_Error_tof)


                #### Calculate Chi2
                ## (1). TRD dimension:
                reduced_chi_squared_TRD = ParameterizationTemplates_tool.CalculateChi2(x_Center_trd, Para_TRD, y_trdfit_value, y_Value_trd, y_Error_trd)
                if value == 'Antiproton':
                    chi2_AntiprotonFit_TRD.append(reduced_chi_squared_TRD)
                    Para_TRD_Anti_all.append(Para_TRD)
                elif value == 'Electron':
                    chi2_ElectronFit_TRD.append(reduced_chi_squared_TRD)
                    Para_TRD_Elec_all.append(Para_TRD)
                elif value == 'Pion':
                    chi2_PionFit_TRD.append(reduced_chi_squared_TRD) 
                    Para_TRD_Pion_all.append(Para_TRD)
                print("chi2/dof(TRD): " + str(reduced_chi_squared_TRD) )
                ## (2). TOF dimension:
                if RigidityRange == 'low':
                    reduced_chi_squared_TOF = ParameterizationTemplates_tool.CalculateChi2(x_Center_tof, Para_TOF, y_toffit_value, y_Value_tof, y_Error_tof)
                    if value == 'Antiproton':
                        chi2_AntiprotonFit_TOF.append(reduced_chi_squared_TOF)
                        Para_TOF_Anti_all.append(Para_TOF)
                    elif value == 'Electron':
                        chi2_ElectronFit_TOF.append(reduced_chi_squared_TOF)
                        Para_TOF_Elec_all.append(Para_TOF)
                    elif value == 'Pion':
                        chi2_PionFit_TOF.append(reduced_chi_squared_TOF)
                        Para_TOF_Pion_all.append(Para_TOF)
                    print("chi2/dof(TOF): " + str(reduced_chi_squared_TOF) )



                #### Convert to histogram and Save HistogramTemplate into ROOT file
                if RigidityRange == 'intermediate':
                    ParameterizationTemplates_tool.MakeAndSaveHistogramTemplate_Intermediate(RigidityRange, index, key, x_Center_trd, x_Edge_trd, y_Value_trd, y_Error_trd, y_trdfit_value_uncer, y_trdfit_error_uncer, totalnumbers, "uncertainty", -999)
                elif RigidityRange == 'low':
                    ParameterizationTemplates_tool.MakeAndSaveHistogramTemplate_Low         (RigidityRange, index, key, x_Center_trd, x_Edge_trd, y_Value_trd, y_Error_trd, y_trdfit_value_uncer, y_trdfit_error_uncer, totalnumbers, x_Center_tof, x_Edge_tof, y_Value_tof, y_Error_tof, y_toffit_value_uncer, y_toffit_error_uncer, "uncertainty", -999)

                
                #### Plot TemplateFit
                if TimeMode == "TimeAveraged":
                    ## (1). TRD Part
                    if TrdFitFunction == "OneGaussOneNovo":
                        ParameterizationTemplates_Plot.Plot_Trd_OneGaussOneNovo(x_Center_trd, x_More_trd, y_Value_trd, y_trdfit_value, y_trdfit_value_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more, value, binname, reduced_chi_squared_TRD, RigidityRange, Para_TRD, Para_TRD_uncertainty)
                    elif TrdFitFunction == "OneGauss":
                        ParameterizationTemplates_Plot.Plot_Trd_OneGauss       (x_Center_trd, x_More_trd, y_Value_trd, y_trdfit_value, y_trdfit_value_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more, value, binname, reduced_chi_squared_TRD, RigidityRange)
                    elif TrdFitFunction == "ExponentialFun":
                        ParameterizationTemplates_Plot.Plot_Trd_ExponentialFun (x_Center_trd, x_More_trd, y_Value_trd, y_trdfit_value, y_trdfit_value_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more, value, binname, reduced_chi_squared_TRD, RigidityRange)
                    elif TrdFitFunction == "OneNovo":
                        ParameterizationTemplates_Plot.Plot_Trd_OneNovo        (x_Center_trd, x_More_trd, y_Value_trd, y_trdfit_value, y_trdfit_value_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more, value, binname, reduced_chi_squared_TRD, RigidityRange)
                    ## (2). TOF Part
                    if RigidityRange == 'low':
                        if TofFitFunction == "OneGauss":
                            ParameterizationTemplates_Plot.Plot_Tof_OneGauss(x_Center_tof, x_More_tof, y_Value_tof, y_toffit_value, y_toffit_value_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more, value, binname, reduced_chi_squared_TOF, RigidityRange)
                        elif TofFitFunction == "OneNovo":
                            ParameterizationTemplates_Plot.Plot_Tof_OneNovo (x_Center_tof, x_More_tof, y_Value_tof, y_toffit_value, y_toffit_value_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more, value, binname, reduced_chi_squared_TOF, RigidityRange)     


            f_ROOT_ParametrilisedTemplates.Close()

            '''
            #### make ParameterError_TRD_all
            ParameterError_TRD_all.append(ParameterError_TRD)
            if RigidityRange == 'low':
                ParameterError_TOF_all.append(ParameterError_TOF)
                print("ParameterError_TOF:"+ str(ParameterError_TOF))
            '''

            #### Save ISS Data
            ParameterizationTemplates_tool.SaveISSData(RigidityRange, IntermediatePath, LowPath, binname, ISSversion, TimeMode, TimeSplitMode, splitindex) 

    if TimeMode == "TimeAveraged": 
        #### Plot Parameter vs Rigidity
        ParameterizationTemplates_Plot.Plot_ParavsRigidty(BinningCenterAll, Para_TRD_Anti_all, "Antiproton", "TRD", RigidityRange)
        ParameterizationTemplates_Plot.Plot_ParavsRigidty(BinningCenterAll, Para_TRD_Elec_all, "Electron"  , "TRD", RigidityRange)
        if RigidityRange == 'low':
            ParameterizationTemplates_Plot.Plot_ParavsRigidty(BinningCenterAll, Para_TOF_Anti_all, "Antiproton", "TOF", RigidityRange)
            ParameterizationTemplates_Plot.Plot_ParavsRigidty(BinningCenterAll, Para_TOF_Elec_all, "Electron"  , "TOF", RigidityRange)
            ParameterizationTemplates_Plot.Plot_ParavsRigidty(BinningCenterAll, Para_TOF_Pion_all, "Pion"      , "TOF", RigidityRange)



        #### Plot Chi2/dof
        ParameterizationTemplates_Plot.Plot_Chi2(BinningCenterAll, chi2_AntiprotonFit_TRD, "chi2_AntiprotonFit_TRD", RigidityRange)
        ParameterizationTemplates_Plot.Plot_Chi2(BinningCenterAll, chi2_ElectronFit_TRD,   "chi2_ElectronFit_TRD"  , RigidityRange)
        if RigidityRange == 'low':
            ParameterizationTemplates_Plot.Plot_Chi2(BinningCenterAll, chi2_AntiprotonFit_TOF, "chi2_AntiprotonFit_TOF", RigidityRange)
            ParameterizationTemplates_Plot.Plot_Chi2(BinningCenterAll, chi2_ElectronFit_TOF,   "chi2_ElectronFit_TOF"  , RigidityRange)
            ParameterizationTemplates_Plot.Plot_Chi2(BinningCenterAll, chi2_PionFit_TOF,       "chi2_PionFit_TOF"      , RigidityRange)
        
        '''
        #### Save ParameterError only for TRD side (Will rewrite old file)
        np.save("ParameterError_TRD" + value + '.npy' , ParameterError_TRD_all)
        if RigidityRange == 'low':
            np.save("ParameterError_TOF" + value + '.npy' , ParameterError_TOF_all)
        '''

if __name__ == '__main__':

    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--Range'         , help='Rigidity Range, low, intermediate or high')
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or 2016paper')
    parser.add_argument('--GenerateNumber', type=int, help='how many tempalte are generated to calculate sys err due to template shape')
    parser.add_argument('--ErrorParametrizationMode', help='if you want to use Parametrilised Error')
    parser.add_argument('--TimeMode',       help='You want to use time averaged or which time dependent data')
    parser.add_argument('--TimeSplitMode',  help='which time slit mode you choose')
    arguments = parser.parse_args()

    if (arguments.issversion and arguments.issversion and arguments.GenerateNumber and arguments.ErrorParametrizationMode and arguments.TimeMode):
        ISSversion               = arguments.issversion
        RigidityRange            = arguments.Range
        GenerateNumber           = arguments.GenerateNumber
        ErrorParametrizationMode = arguments.ErrorParametrizationMode
        TimeMode                 = arguments.TimeMode
    else:
        print("You need to provide all parameters!")
        os._exit(0)
    TimeSplitMode = arguments.TimeSplitMode

    IntermediatePath = os.getenv('HPCINTERMEDIATEDIR')
    LowPath          = os.getenv('HPCLOWENERGYDIR')

    if RigidityRange == 'intermediate':
        BinningAll       = binning.binslow[10:]
        BinningAll[19]   = '16.6_18'
        BinningCenterAll = binning.published2016binnings_center[9:29]
    elif RigidityRange == 'low':
        BinningAll       = binning.binslow[1:17]
        BinningAll[0]    = '1_1.16'
        BinningCenterAll = binning.published2016binnings_center[0:16]
    BinsEdgeName_TimeDependent_Low = binning.BinsEdgeName_TimeDependent_Low

    if TimeMode == "TimeAveraged":
        UsedBinningEdge = BinningAll
        SplitTotal = 1
        mergestep = 1
    elif TimeMode == "TimeDependent":
        if RigidityRange == 'intermediate':
            print("in progress")     
            os._exit(0)
        elif RigidityRange == 'low':        
            UsedBinningEdge = BinsEdgeName_TimeDependent_Low

        if TimeSplitMode == "3BartalRotation":
            SplitTotal = 123
            mergestep = 3
        elif TimeSplitMode == "6BartalRotation":
            SplitTotal = 121
            mergestep = 6
        elif TimeSplitMode == "6Months":
            print("has problem!")
            os._exit(0) 
            SplitTotal = 20
            mergestep = 1 ##???!!!

    chi2_AntiprotonFit_TRD = []
    chi2_ElectronFit_TRD   = []
    chi2_PionFit_TRD       = []
    chi2_AntiprotonFit_TOF = []
    chi2_ElectronFit_TOF   = []
    chi2_PionFit_TOF       = []

    Para_TRD_Anti_all = []
    Para_TRD_Elec_all = []
    Para_TRD_Pion_all = []
    Para_TOF_Anti_all = []
    Para_TOF_Elec_all = []
    Para_TOF_Pion_all = []

    main()

