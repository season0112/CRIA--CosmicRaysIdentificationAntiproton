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
import ParameterizationTemplates_InitialValues
import inspect, re



def GetValues(RigidityRange, f_uproot_OriginalTemplate, key):
    if RigidityRange == 'intermediate':
        x_edge = f_uproot_OriginalTemplate[key].axis().edges()
        x_center = []
        for i in range(x_edge.shape[0]-1):
            x_center.append( (x_edge[i] + x_edge[i+1])/2 )
        x_center = np.array(x_center)
        y_value = f_uproot_OriginalTemplate[key].values()
        y_error = f_uproot_OriginalTemplate[key].errors()
        totalnumbers = sum(y_value)
        x_more = np.arange(x_edge[0], x_edge[-1], 0.001)

        return x_edge, x_center, y_value, y_error, totalnumbers, x_more

    elif RigidityRange == 'low':
        # TRD
        x_edge = f_uproot_OriginalTemplate[key].axis(0).edges() # index=0:TRD, index=1:TOF
        x_center = []
        for i in range(x_edge.shape[0]-1):
            x_center.append( (x_edge[i] + x_edge[i+1])/2 )
        x_center = np.array(x_center)
        y_value = []
        y_error = []
        for index in range(f_uproot_OriginalTemplate[key].values().shape[0]):
            y_value.append( sum(f_uproot_OriginalTemplate[key].values()[index, :]) )
            y_error.append( np.sqrt(sum(f_uproot_OriginalTemplate[key].values()[index, :])) )
        totalnumbers = sum(y_value)
        x_more = np.arange(x_edge[0], x_edge[-1], 0.001)
        # TOF
        x_edge_tof = f_uproot_OriginalTemplate[key].axis(1).edges() # index=0:TRD, index=1:TOF
        x_center_tof = []
        for i in range(x_edge_tof.shape[0]-1):
            x_center_tof.append( (x_edge_tof[i] + x_edge_tof[i+1])/2 )
        x_center_tof = np.array(x_center_tof)
        y_value_tof = []
        y_error_tof = []
        for index in range(f_uproot_OriginalTemplate[key].values().shape[1]):
            y_value_tof.append( sum(f_uproot_OriginalTemplate[key].values()[:, index]) )
            y_error_tof.append( np.sqrt(sum(f_uproot_OriginalTemplate[key].values()[:, index])) )
        x_more_tof = np.arange(x_edge_tof[0], x_edge_tof[-1], 0.001)
        return x_edge, x_center, y_value, y_error, totalnumbers, x_more, x_edge_tof, x_center_tof, y_value_tof, y_error_tof, x_more_tof 


def Trd_UsedFitFunction(TemplateName, RigidityRange):

    if RigidityRange == 'intermediate':
        if TemplateName == "Antiproton":
            FitFunction    = "OneGaussOneNovo"
        elif TemplateName == "Electron":
            FitFunction    = "OneGaussOneNovo"
        elif TemplateName == "Pion":
            FitFunction    = "OneGaussOneNovo"

    elif RigidityRange == 'low':
        if TemplateName == "Antiproton":
            FitFunction    = "OneGauss"
        elif TemplateName == "Electron":
            FitFunction    = "ExponentialFun"
        elif TemplateName == "Pion":
            #FitFunction    = "OneNovo"
            FitFunction    = "OneGaussOneNovo"

    return FitFunction


def Tof_UsedFitFunction(TemplateName, RigidityRange):

    if RigidityRange == 'low':
        if TemplateName == "Antiproton":
            FitFunction    = "OneGauss"
        elif TemplateName == "Electron":
            FitFunction    = "OneNovo"
        elif TemplateName == "Pion":
            FitFunction    = "OneNovo"

    return FitFunction


def GetFitCurve_TRD(RigidityRange, UsedFunction_TRD, x_Center_trd, Para_TRD, Para_TRD_uncertainty, x_More_trd, y_Value_trd, y_Error_trd, index):

    if RigidityRange == 'low' or RigidityRange == 'intermediate':
        y_trdfit_value            = UsedFunction_TRD(x_Center_trd, *Para_TRD)
        y_trdfit_error            = np.sqrt(y_trdfit_value)
        y_trdfit_value_uncer      = UsedFunction_TRD(x_Center_trd, *Para_TRD_uncertainty)
        y_trdfit_error_uncer      = np.sqrt(y_trdfit_value_uncer)
        y_trdfit_value_more       = UsedFunction_TRD(x_More_trd, *Para_TRD)
        y_trdfit_error_more       = np.sqrt(y_trdfit_value_more)
        y_trdfit_value_uncer_more = UsedFunction_TRD(x_More_trd, *Para_TRD_uncertainty)
        y_trdfit_error_uncer_more = np.sqrt(y_trdfit_value_uncer_more)
        # Naive Fix for bin content is 0
        ParameterizationTemplates_InitialValues.FixForNaf( {"variable":[y_trdfit_value, y_trdfit_value_uncer, y_trdfit_value_more, y_trdfit_value_uncer_more], "error":[y_trdfit_error, y_trdfit_error_uncer, y_trdfit_error_more, y_trdfit_error_uncer_more] })
        # Fix for 1.33-1.51 GV
        if index == 2:
            #if value == "Electron" or value == "Pion":
            y_trdfit_value_uncer[0] = y_Value_trd[0]
            y_trdfit_error_uncer[0] = y_Error_trd[0]
            y_trdfit_value_uncer[1] = y_Value_trd[1]
            y_trdfit_error_uncer[1] = y_Error_trd[1]
            y_trdfit_value_uncer[2] = y_Value_trd[2]
            y_trdfit_error_uncer[2] = y_Error_trd[2]
    return y_trdfit_value, y_trdfit_error, y_trdfit_value_uncer, y_trdfit_error_uncer, y_trdfit_value_more, y_trdfit_error_more, y_trdfit_value_uncer_more, y_trdfit_error_uncer_more 


def GetFitCurve_TOF(RigidityRange, UsedFunction_TOF, x_Center_tof, Para_TOF, Para_TOF_uncertainty, x_More_tof, y_Value_tof, y_Error_tof):

    if RigidityRange == 'low':
        y_toffit_value            = UsedFunction_TOF(x_Center_tof, *Para_TOF)
        y_toffit_error            = np.sqrt(y_toffit_value)
        y_toffit_value_uncer      = UsedFunction_TOF(x_Center_tof, *Para_TOF_uncertainty)
        y_toffit_error_uncer      = np.sqrt(y_toffit_value_uncer)
        y_toffit_value_more       = UsedFunction_TOF(x_More_tof, *Para_TOF)
        y_toffit_error_more       = np.sqrt(y_toffit_value_more)
        y_toffit_value_uncer_more = UsedFunction_TOF(x_More_tof, *Para_TOF_uncertainty)
        y_toffit_error_uncer_more = np.sqrt(y_toffit_value_uncer_more)
        # Naive Fix for bin content is 0
        ParameterizationTemplates_InitialValues.FixForNaf( {"variable":[y_toffit_value, y_toffit_value_uncer, y_toffit_value_more, y_toffit_value_uncer_more], "error":[y_toffit_error, y_toffit_error_uncer, y_toffit_error_more, y_toffit_error_uncer_more]} )

    return y_toffit_value, y_toffit_error, y_toffit_value_uncer, y_toffit_error_uncer, y_toffit_value_more, y_toffit_error_more, y_toffit_value_uncer_more, y_toffit_error_uncer_more


def CalculateChi2(x_center, popt, y_fit_value, y_value, y_error):
    #pulls = ( y_value - y_fit_value ) / y_error
    pulls = ( y_value - y_fit_value ) / np.sqrt(y_fit_value)
    #print("pulls:" + str(pulls))

    chi_squared = np.sum(pulls**2)
    reduced_chi_squared = chi_squared/( x_center.shape[0] - len(popt) )

    return reduced_chi_squared


def MakeAndSaveHistogramTemplate_Intermediate(RigidityRange, Rigidityindex, key, x_center, x_edge, y_Value_trd, y_Error_trd, y_trdfit_value_uncer, y_trdfit_error_uncer, totalnumbers, templatemode, RandomTemplateindex):

    h_template = TH1F("h_template", "", x_center.shape[0], x_edge)
    for i in range(h_template.GetNbinsX()):
        h_template.SetBinContent(i+1, y_trdfit_value_uncer[i])
        h_template.SetBinError  (i+1, y_trdfit_error_uncer[i])
    ## Nomalization of template
    #h_template.Scale(totalnumbers)

    ## Save
    if templatemode == "uncertainty":
        h_template.Write( key + "_fit" )
    elif templatemode == "randomlized":
        h_template.Write( key + "_RandomIndex_" + str(RandomTemplateindex) )
    h_template.Delete()


def MakeAndSaveHistogramTemplate_Low(RigidityRange, Rigidityindex, key, x_Center_trd, x_Edge_trd, y_Value_trd, y_Error_trd, y_trdfit_value_uncer, y_trdfit_error_uncer, totalnumbers, x_Center_tof, x_Edge_tof, y_Value_tof, y_Error_tof, y_toffit_value_uncer, y_toffit_error_uncer, templatemode, RandomTemplateindex):

    h_trd = TH1F("h_trd", "", x_Center_trd.shape[0], x_Edge_trd)
    h_tof = TH1F("h_tof", "", x_Center_tof.shape[0], x_Edge_tof)
    for i in range(h_trd.GetNbinsX()):
        h_trd.SetBinContent(i+1, y_trdfit_value_uncer[i])
        h_trd.SetBinError  (i+1, y_trdfit_error_uncer[i])
    for i in range(h_tof.GetNbinsX()):
        h_tof.SetBinContent(i+1, y_Value_tof[i])
        h_tof.SetBinError  (i+1, y_Error_tof[i])

    h_template2D = TH2F("h_template2D", "", x_Center_trd.shape[0], x_Edge_trd, x_Center_tof.shape[0], x_Edge_tof)
    for i in range(h_trd.GetNbinsX()):
        for j in range(h_tof.GetNbinsX()):
            h_template2D.SetBinContent(i+1, j+1, y_trdfit_value_uncer[i] * y_toffit_value_uncer[j])
            h_template2D.SetBinError  (i+1, j+1, np.sqrt(y_trdfit_value_uncer[i] * y_toffit_value_uncer[j]) )

    ## Nomalization of template
    totalnumbers2 = h_template2D.Integral()
    h_template2D.Scale(1/totalnumbers2*totalnumbers)

    ## Save
    if templatemode == "uncertainty":
        h_template2D.Write( key + "_fit" )
    elif templatemode == "randomlized":
        h_template2D.Write( key + "_RandomIndex_" + str(RandomTemplateindex) )

    ## Plot
    #c2 = TCanvas("c2", "c2", 1000, 500);
    #h_template2D.Draw("COLZ")
    #c2.SaveAs(str(key) + "_RigidityIndex_"+ str(Rigidityindex) + "_fit.pdf");

    h_template2D.Delete()


def SaveISSData(RigidityRange, IntermediatePath, LowPath, binname, ISSversion, TimeMode, TimeSplitMode, splitindex):
    if TimeMode == "TimeAveraged":
        if RigidityRange == 'intermediate':
            ## Open Original ROOT File
            f_ROOT_OriginalTemplate = TFile( IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/RigidityRootFiles/averaged_ratio_" + binname + "_" + ISSversion + ".root", "READ")
            ISSNegativeData = f_ROOT_OriginalTemplate.Get("data_pass7_negative")
            ISSPositiveData = f_ROOT_OriginalTemplate.Get("data_pass7_positive")
            ## Save ISS Data Template
            f_ROOT_ParametrilisedTemplates = TFile( IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/RigidityRootFiles/averaged_ratio_" + binname + "_" + ISSversion + "_new.root", "UPDATE")
            ISSNegativeData.Write("data_pass7_negative")
            ISSPositiveData.Write("data_pass7_positive")
            f_ROOT_ParametrilisedTemplates.Close()

            f_ROOT_OriginalTemplate.Close()

        elif RigidityRange == 'low':
            ## Open Original ROOT File
            f_ROOT_OriginalTemplate = TFile( LowPath + "/totalall/Time_Averaged_ratio_Low/binmerge1/tof/averaged_ratio_" + binname + "_" + ISSversion + ".root", "READ")
            ISSNegativeData      = f_ROOT_OriginalTemplate.Get("TofTRD_data_pass7_negative")
            ISSPositiveData      = f_ROOT_OriginalTemplate.Get("TofTRD_data_pass7_positive_data")
            ISSPositiveData_test = f_ROOT_OriginalTemplate.Get("TofTRD_data_pass7_positive_data_test")
            ## Save ISS Data Template
            f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/Time_Averaged_ratio_Low/binmerge1/tof/averaged_ratio_" + binname + "_" + ISSversion + "_new.root", "UPDATE")
            ISSNegativeData.Write("TofTRD_data_pass7_negative")
            ISSPositiveData.Write("TofTRD_data_pass7_positive_data")
            ISSPositiveData_test.Write("TofTRD_data_pass7_positive_data_test")

            f_ROOT_ParametrilisedTemplates.Close()

            f_ROOT_OriginalTemplate.Close()

    elif TimeMode == "TimeDependent":
        ## Open Original ROOT File
        if RigidityRange == 'intermediate':
            print("in progess")
        elif RigidityRange == 'low':
            if TimeSplitMode == "3BartalRotation": 
                f_ROOT_OriginalTemplate = TFile( LowPath + "/totalall/rootfiles_3BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root", "READ")
            elif TimeSplitMode == "6BartalRotation":
                f_ROOT_OriginalTemplate = TFile( LowPath + "/totalall/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root", "READ")
            elif TimeSplitMode == "6Months":
                f_ROOT_OriginalTemplate = TFile( LowPath + "/totalall/rootfiles_6months_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + ".root", "READ")
            ISSNegativeData      = f_ROOT_OriginalTemplate.Get("TofTRD_data_pass7_negative")
            ISSPositiveData      = f_ROOT_OriginalTemplate.Get("TofTRD_data_pass7_positive_data")
        ## Save ISS Data Template
        if RigidityRange == 'intermediate':
            print("in progess")
        elif RigidityRange == 'low':
            if TimeSplitMode == "3BartalRotation":
                f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_3BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + "_new.root", "UPDATE")
            elif TimeSplitMode == "6BartalRotation":
                f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + "_new.root", "UPDATE")
            elif TimeSplitMode == "6Months":
                f_ROOT_ParametrilisedTemplates = TFile( LowPath + "/totalall/rootfiles_6months_template/TimeDependentTemplatesAndData_"+ binname +"_" + str(splitindex) + "_new.root", "UPDATE")
            ISSNegativeData.Write("TofTRD_data_pass7_negative")
            ISSPositiveData.Write("TofTRD_data_pass7_positive_data")
        f_ROOT_ParametrilisedTemplates.Close()
        f_ROOT_OriginalTemplate.Close()











