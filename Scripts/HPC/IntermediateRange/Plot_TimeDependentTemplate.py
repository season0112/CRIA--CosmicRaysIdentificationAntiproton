import numpy as np
import matplotlib.pyplot as plt
import argparse
import binning
import uproot
import TimeDependentTimeRange
import uncertainties
import ParameterizationTemplates_function
import ParameterizationTemplates_tool
import ParameterizationTemplates_InitialValues
import ParameterizationTemplates_Plot
import ParameterizationTemplates_Fit
import os
from datetime import datetime
import PythonPlotDefaultParameters
from uncertainties import ufloat
from uncertainties import unumpy
import Plot_TimeDependentTemplate_tool
import scipy.stats as stats
from scipy.signal import savgol_filter


def main():

    for value in value_all:
    
        print("Now the template is: " + str(value))
        
        for (R_index, RigidityBin) in enumerate(BinningEdgeName):

            #if R_index != 5: #low
            if R_index != 4: #intermediate
                continue

            print('\n'); print("R_index:" + str(R_index) + ", RigidityBin:" + str(RigidityBin))

            miu_all   = np.array(())
            sigma_all = np.array(())
            chi2_all  = []
            
            for time_index in range(0, 132+1, 6):
            #for time_index in range(0, 3, 6):


                if rigidityrange == 'intermediate':
                    file = uproot.open("/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v6.0/total/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_fullRange_" + str(RigidityBin) + "_" + str(time_index) + ".root")
                    if value == "Antiproton":
                        UsedFunction = 'Novosibirsk'
                        FitHalfRange = 7
                    elif value == "Electron":
                        UsedFunction = 'Novosibirsk'
                        FitHalfRange = 7
                elif rigidityrange == 'low':
                    file = uproot.open("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/rootfiles_6BartalRotation_template/TimeDependentTemplatesAndData_" + str(RigidityBin) + "_" + str(time_index) + ".root")
                    if value == "Antiproton":
                        UsedFunction = 'Gaussion'
                        FitHalfRange = 7
                    elif value == "Electron":
                        UsedFunction = 'Gaussion'
                        FitHalfRange = 7
                    elif value == "Pion":
                        #UsedFunction = 'Novosibirsk'
                        UsedFunction = 'Gaussion'
                        FitHalfRange = 9
                
                if UsedFunction == 'Novosibirsk':
                    InitialValue = [0, 0, 0.2, 0.1]
                elif UsedFunction == 'Gaussion':
                    InitialValue = [0, 0, 0.2]

                # Get Fit Range, and Initial Value
                FileName, Dimension = Plot_TimeDependentTemplate_tool.GetFileName(rigidityrange, value)
                # Get Value
                if rigidityrange == 'intermediate':
                    BinCenter    = (file[FileName].to_numpy()[1][0:-1] + file[FileName].to_numpy()[1][1:])/2
                    BinContent   = file[FileName].to_numpy()[0]
                elif rigidityrange == 'low':
                    BinCenter    = (file[FileName].to_numpy()[2][0:-1] + file[FileName].to_numpy()[2][1:])/2
                    BinContent   = file[FileName].to_numpy()[0].sum(axis=0)
                # Smooth the content
                #if rigidityrange == 'low':
                #    BinContent = savgol_filter(BinContent, 3, 1) 

                # [0]:peak, [1]:miu, [2]:sigma
                print('\n')
                InitialValue[0] = np.max(BinContent)
                InitialValue[1] = BinCenter[ np.where(np.array(BinContent)==np.max(BinContent))[0][0] ]
                maxvalue_index = np.where(np.array(BinContent)==np.max(BinContent))[0]
                #lowedge  = int(InitialValue[1] - np.std(np.array(BinContent))/2 - 3)
                #highedge = int(InitialValue[1] + np.std(np.array(BinContent))/2 + 3)
                lowedge  = int(maxvalue_index[0] - FitHalfRange)
                highedge = int(maxvalue_index[0] + FitHalfRange)
                if lowedge < 0:
                    fit_start = 0
                else: 
                    fit_start = int(lowedge)
                if highedge > BinCenter.shape[0]:
                    fit_end = int(BinCenter.shape[0])
                else:
                    fit_end = int(highedge)
                #print("time_index:" + str(time_index))
                #print("BinContent:" + str(BinContent))
                #print("max value:" + str(np.max(BinContent)))
                #print("max value index:" + str(maxvalue_index) )
                #print("np.std(np.array(BinContent)):" + str(np.std(np.array(BinContent))))
                #print("InitialValue[0]: " + str(InitialValue[0]))
                #print("InitialValue[1]: " + str(InitialValue[1]))
                #print("fit_start:" + str(fit_start)) 
                #print("fit_end:" + str(fit_end))                  

                
                InitialValue, fit_start, fit_end = Plot_TimeDependentTemplate_tool.TunningForFitSettings(rigidityrange, value, UsedFunction, R_index, time_index, InitialValue, fit_start, fit_end)
                #FileName, Dimension, InitialValue, fit_start, fit_end = Plot_TimeDependentTemplate_tool.InitialValue(rigidityrange, value, UsedFunction, R_index)


                # Fit
                if UsedFunction == 'Gaussion':
                    Para, Covan, ParameterError = ParameterizationTemplates_Fit.CurveFit_OneGauss       (Dimension, rigidityrange, value, BinCenter[fit_start:fit_end], BinContent[fit_start:fit_end], np.sqrt(BinContent[fit_start:fit_end]), InitialValue[0], InitialValue[1], InitialValue[2], time_index)
                elif UsedFunction == 'Novosibirsk':
                    Para, Covan, ParameterError = ParameterizationTemplates_Fit.CurveFit_OneNovo (Dimension, rigidityrange, value, BinCenter[fit_start:fit_end], BinContent[fit_start:fit_end], np.sqrt(BinContent[fit_start:fit_end]), InitialValue[0], InitialValue[1], InitialValue[2], InitialValue[3], time_index)
                #print("Para" + str(Para)); print("Covan" + str(Covan)); print("ParameterError" + str(ParameterError))
                miu_all   = np.append(miu_all  , ufloat(Para[1], ParameterError[1]) ) # Gaussion and Novosibirsk
                sigma_all = np.append(sigma_all, ufloat(Para[2], ParameterError[2]) ) # Gaussion and Novosibirsk

                # Get Fit value
                if UsedFunction == 'Gaussion':
                    y_fit = ParameterizationTemplates_function.Gauss_function(BinCenter, Para[0], Para[1], Para[2]) # Gaussion
                elif UsedFunction == 'Novosibirsk':
                    y_fit = ParameterizationTemplates_function.Novosibirsk_function(BinCenter, Para[0], Para[1], Para[2], Para[3]) # Novosibirsk
                if rigidityrange == 'intermediate':
                    x_more     = np.arange(-2.0, 0, 0.005) #TRD
                elif rigidityrange == 'low':
                    x_more     = BinCenter   #TOF
                if UsedFunction == 'Gaussion':
                    y_fit_more = ParameterizationTemplates_function.Gauss_function(x_more, Para[0], Para[1], Para[2]) 
                elif UsedFunction == 'Novosibirsk':
                    y_fit_more = ParameterizationTemplates_function.Novosibirsk_function(x_more, Para[0], Para[1], Para[2], Para[3]) 
                
                # Calculate Chi2
                reduced_chi_squared_TRD = ParameterizationTemplates_tool.CalculateChi2(BinCenter[fit_start:fit_end], Para, y_fit[fit_start:fit_end], BinContent[fit_start:fit_end], np.sqrt(BinContent[fit_start:fit_end]))  
                #print(stats.chisquare(f_obs=BinContent, f_exp=y_fit)[0]/(BinCenter.shape[0]-len(Para)))
                #print("reduced_chi_squared_TRD: " + str(reduced_chi_squared_TRD))
                chi2_all.append(reduced_chi_squared_TRD)

                # Plot Template
                #Plot_TimeDependentTemplate_tool.Plot_Template(BinCenter, y_fit, reduced_chi_squared_TRD, BinContent, value, RigidityBin, time_index, Dimension, fit_start, fit_end)
                
                if rigidityrange == 'intermediate':
                    if value_all == ["Antiproton", "Electron"]:
                        if time_index == TemplatesFitPlotTimeIndex:
                            if RigidityBin == TemplatesFitPlotRigidityBin:
                                if value == "Antiproton":
                                    y_value_antiproton      = BinContent
                                    y_fit_antiproton        = y_fit_more
                                    reduced_chi2_antiproton = reduced_chi_squared_TRD
                                elif value == "Electron":
                                    y_value_electron      = BinContent
                                    y_fit_electron        = y_fit_more
                                    reduced_chi2_electron = reduced_chi_squared_TRD

                
            print("Chi2: " + str(chi2_all))
            Plot_TimeDependentTemplate_tool.Plot_Miu_of_Gaussion  (Time6B, unumpy.nominal_values(miu_all)  , unumpy.std_devs(miu_all)  , value, RigidityBin, Dimension, yaxisfactor, rigidityrange)
            Plot_TimeDependentTemplate_tool.Plot_Sigma_of_Gaussion(Time6B, unumpy.nominal_values(sigma_all), unumpy.std_devs(sigma_all), value, RigidityBin, Dimension, yaxisfactor, rigidityrange)
            #Plot_TimeDependentTemplate_tool.Plot_Chi2(Time6B, chi2_all, "chi2_" + value + "Fit_" + Dimension, rigidityrange, RigidityBin)

            
    
    if value_all == ["Antiproton", "Electron"] and rigidityrange == 'intermediate':
        Plot_TimeDependentTemplate_tool.Plot_TemplatetFitPlot(BinCenter, x_more, y_value_antiproton, y_value_electron, y_fit_antiproton, y_fit_electron, reduced_chi2_antiproton, reduced_chi2_electron, TemplatesFitPlotRigidityBin, TemplatesFitPlotTimeIndex, Dimension)


if __name__ == "__main__":


    rigidityrange = 'intermediate'
    #rigidityrange = 'low'

    if rigidityrange == 'intermediate':
        BinningEdgeName = binning.BinsEdgeName_TimeDependent_Intermediate
        BinningEdgeName[6] = "9.26_11.0"
        BinningEdgeName[7] = "11.0_13.0"
        BinningEdgeName[8] = "13.0_15.3"
        yaxisfactor = 1.01
    elif rigidityrange == 'low':
        BinningEdgeName = binning.BinsEdgeName_TimeDependent_Low
        yaxisfactor = 1.1

    Time6B = [datetime.fromtimestamp(i) for i in binning.Bartals6Unixtime[0:23]]

    #value_all = ["Antiproton", "Electron", "Pion"]
    value_all = ["Antiproton", "Electron"]
    #value_all = ["Antiproton"]
    #value_all = ["Electron"]
    #value_all = ["Pion"]
   
    TemplatesFitPlotTimeIndex   = 78  # 0 or 78
    TemplatesFitPlotRigidityBin = "6.47_7.76" 

    main()



