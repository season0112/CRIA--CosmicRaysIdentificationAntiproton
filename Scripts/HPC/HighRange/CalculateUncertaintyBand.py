#!/usr/bin/env python
from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
import argparse
import os
import binning
from root_numpy import root2array, tree2array
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import ROOT
import uncertainties
from uncertainties import unumpy
import scipy as sp
import scipy.stats as stats
from pylab import *
from kapteyn import kmpfit
import CalculateUncertaintyBand_function
from iminuit import Minuit
from iminuit.cost import LeastSquares
from jacobi import propagate


def main():

    #RatioChoice = 'CCProtonOverProtonRatio' 
    RatioChoice = 'CCLevelISSToMc'

    if RatioChoice == 'CCProtonOverProtonRatio':
        Name = 'CCProtonOverProton'
    elif RatioChoice == 'CCLevelISSToMc':
        Name = 'CCLevel'

    #### Calculate CI as uncertainties 

    for SignalEfficiency in SignalEfficiency_all:
        print('\n')
        print("SignalEfficiency: " + str(SignalEfficiency))
        for CCcut_TF in CCcut_TF_all:
            print("CCcut_TF: " + str(CCcut_TF))

            #### Load Result and the Ratio, used to Calculate Uncertainty Band.
            ValueA_AllPointsForNoFurtherRemoved, ValueAError_AllPointsForNoFurtherRemoved, ValueB_AllPointsForNoFurtherRemoved, ValueBError_AllPointsForNoFurtherRemoved, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, RigidityBinningCenter_AllPointsForNoFurtherRemoved = CalculateUncertaintyBand_function.LoadResult(workpath, pattern, NNsuffix, ISSversion, Binningversion, SignalEfficiency, CCcut_TF, RatioChoice)
            print("First fit index: " + str(PointsForFit_Begin))
            #### Remove first two points due to statistics for no matther CC cut value
            ValueA_AllPointsForNoFurtherRemoved                = ValueA_AllPointsForNoFurtherRemoved[2:]
            ValueAError_AllPointsForNoFurtherRemoved           = ValueAError_AllPointsForNoFurtherRemoved[2:]
            ValueB_AllPointsForNoFurtherRemoved                = ValueB_AllPointsForNoFurtherRemoved[2:]
            ValueBError_AllPointsForNoFurtherRemoved           = ValueBError_AllPointsForNoFurtherRemoved[2:]
            RatioUserd_AllPointsForNoFurtherRemoved            = RatioUserd_AllPointsForNoFurtherRemoved[2:]
            RatioErrorUserd_AllPointsForNoFurtherRemoved       = RatioErrorUserd_AllPointsForNoFurtherRemoved[2:]
            RigidityBinningCenter_AllPointsForNoFurtherRemoved = RigidityBinningCenter_AllPointsForNoFurtherRemoved[2:]


            #### Scaler of data error 
            #Scaler = 2.57 # (For CC cut = 0)
            Scaler = 1.7231928866332435 # (For CC cut = 0.2)
            #Scaler = 1
            RatioErrorUserd_AllPointsForNoFurtherRemoved = RatioErrorUserd_AllPointsForNoFurtherRemoved * Scaler        

            #### Further Remove Points For fit in different CC Cut values. 
            if CCcut_TF == "0.00":
                FurtherPointRemoved = 0
            elif CCcut_TF == "0.20":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.40":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.60":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.65":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.70":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.75":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.80":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.85":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.90":
                FurtherPointRemoved = 2
            elif CCcut_TF == "0.95":
                FurtherPointRemoved = 2
            ValueA                = ValueA_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            ValueAError           = ValueAError_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            ValueB                = ValueB_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            ValueBError           = ValueBError_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            RatioUserd            = RatioUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            RatioErrorUserd       = RatioErrorUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]
            RigidityBinningCenter = RigidityBinningCenter_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]


            #### Option 1: Linear fit
            print('\n')
            print("Option1: ")
            y_model, y_more, ConfidenceInterval_show_67, ConfidenceInterval_Option1, PredictionInterval_show_67, PredictionInterval_show_95 = CalculateUncertaintyBand_function.LinearFit(NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, RigidityBinningCenter)
            print(ConfidenceInterval_Option1)

            #### Option 2: kmpfit (FIXME: When estimating the parameter error, the data error is not used in this package?) 
            print('\n')
            print("Option2 KMPfit:")
            ## 2.1 Poly0 Fit   
            kmpCL_67_P0, kmpCL_67_save_P0, yhat_67_P0, upper_67_P0, lower_67_P0 = CalculateUncertaintyBand_function.KMPFit(CalculateUncertaintyBand_function.PolyFit0Model, [1]            , NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, CalculateUncertaintyBand_function.dfdpPolyFit0, RigidityBinningCenter)
            yhat_67_P0 = np.ones(kmpCL_67_P0.size) * yhat_67_P0
            print('\n')
            ## 2.2 Poly1 Fit
            #kmpCL_67_P1, kmpCL_67_save_P1, yhat_67_P1, upper_67_P1, lower_67_P1 = CalculateUncertaintyBand_function.KMPFit(CalculateUncertaintyBand_function.LinearFitModel, [1, 1]        , NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, CalculateUncertaintyBand_function.dfdpLinearFit, RigidityBinningCenter)
            #print('\n')
            ## 2.3 Poly2 Fit
            #kmpCL_67_P3, kmpCL_67_save_P3, yhat_67_P3, upper_67_P3, lower_67_P3 = CalculateUncertaintyBand_function.KMPFit(CalculateUncertaintyBand_function.PolyFit3Model, [1, 1, 1, 1], NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, CalculateUncertaintyBand_function.dfdpPolyFit3, RigidityBinningCenter)
            #print('\n')
            ## 2.4 Poly3 Fit
            #kmpCL_67_P4, kmpCL_67_save_P4, yhat_67_P4, upper_67_P4, lower_67_P4 = CalculateUncertaintyBand_function.KMPFit(CalculateUncertaintyBand_function.PolyFit4Model, [1, 1, 1, 1, 1], NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, CalculateUncertaintyBand_function.dfdpPolyFit4, RigidityBinningCenter)
            #print('\n')

            #### Option 3: curve_fit (Only For Const Fit)
            print("Option3 curve_fit:")
            print("Data: " + str(RatioUserd))
            print("Data Error: " + str(RatioErrorUserd))
            popt, pcov = sp.optimize.curve_fit( CalculateUncertaintyBand_function.PolyFit0Model_ForCurveFit, RigidityBinningCenter, RatioUserd, 1.0, RatioErrorUserd, True )
            y_model_option3 = popt
            print("Parameters: " + str(y_model_option3))
            ParameterError = np.sqrt(np.diag(pcov))
            print("Parameters Error: " + str(ParameterError))
            residual = RatioUserd - CalculateUncertaintyBand_function.PolyFit0Model_ForCurveFit( RigidityBinningCenter, *popt )
            print("residual:" + str(residual))
            #print(RatioErrorUserd)
            chi2 = (np.sum(residual**2/RatioErrorUserd**2))
            print("chi2:" + str((chi2)))
            dof = len(RigidityBinningCenter)-1
            print("dof:" + str(dof))
            chi2dof = chi2/dof
            print("chi2/dof:" + str(chi2dof))
            y_more_option3 = CalculateUncertaintyBand_function.PolyFit0Model_ForCurveFit(x_more, *popt) 
            ConfidenceInterval_Option3 = [ParameterError[0]*2] * 32 
            print('\n')

            #### Option 4: Manually BandWidth For 68% Points included. (4,5,6th points (126.25, 161.0, 193.0 GV) are out of band for CC cut=0.0.)
            print("Option 4: BandWidth For 68% Points included. (Currently only for CC cut = 0.00)")
            print("Residual:" + str(abs(RatioUserd-y_model_option3)))
            BandWidthFor68PercentPoints_Option4 = [0.10416928]*32
            print("BandWidthFor68PercentPoints_Option4:" + str(BandWidthFor68PercentPoints_Option4))
            BandWidthFor68PercentPoints_more = [0.10416928]*len(x_more)
            print('\n')

            #### Option 4.2: BandWidth For Points including a bump
            print("Option 4.2: BandWidth For Points included with a bump")
            ## Fit for bump
            if CCcut_TF == "0.00":
                BumpStart_in_ReducedCCLevelPoint = 2
                BumpEnd_in_ReducedCCLevelPoint   = 5
                BumpRangeStart_in_x_more         = 13
                BumpRangeEnd_in_x_more           = 29
                BumpRange_PeakValue_in_x_more    = 21
                popt_bump, pcov_bump = sp.optimize.curve_fit( CalculateUncertaintyBand_function.PolyFit2Model_ForCurveFit, RigidityBinningCenter[BumpStart_in_ReducedCCLevelPoint:BumpEnd_in_ReducedCCLevelPoint], RatioUserd[BumpStart_in_ReducedCCLevelPoint:BumpEnd_in_ReducedCCLevelPoint], [1, 1, 2], RatioErrorUserd[BumpStart_in_ReducedCCLevelPoint:BumpEnd_in_ReducedCCLevelPoint], True )
                BumpRangeStart_in_NNpoint = 23
                BumpRangeEnd_in_NNpoint   = 27
            elif CCcut_TF == "0.20":
                BumpStart_in_ReducedCCLevelPoint = 0
                BumpEnd_in_ReducedCCLevelPoint   = 6
                BumpRangeStart_in_x_more         = 13
                BumpRangeEnd_in_x_more           = 52
                BumpRange_PeakValue_in_x_more    = 32
                popt_bump, pcov_bump = sp.optimize.curve_fit( CalculateUncertaintyBand_function.PolyFit2Model_ForCurveFit, RigidityBinningCenter[0:3].tolist()+RigidityBinningCenter[4:6].tolist(), RatioUserd[0:3].tolist()+RatioUserd[4:6].tolist(), [1, 1, 2], RatioErrorUserd[0:3].tolist()+RatioErrorUserd[4:6].tolist(), True )
                print("popt_bump:" + str(popt_bump))
                #popt_bump = [ 2.66864225e-05 -9.35731526e-03  1.58098460e+00]
                #print("popt_bump:" + str(popt_bump))
                BumpRangeStart_in_NNpoint = 23
                BumpRangeEnd_in_NNpoint   = 29
            
            y_more_Bump = CalculateUncertaintyBand_function.PolyFit2Model_ForCurveFit( x_more, *popt_bump )
            print("3rd point:" + str(RigidityBinningCenter[2]))
            print("4th point:" + str(RigidityBinningCenter[3]))
            print("5th point:" + str(RigidityBinningCenter[4]))
            print("6th point:" + str(RigidityBinningCenter[5]))
            print("x_more:" + str(x_more[BumpRangeStart_in_x_more]))
            print("x_more:" + str(x_more[BumpRangeEnd_in_x_more]))
            print('NNpoint:' + str(NNpoint))
            print('NNpoint size:' + str(len(NNpoint)))
            print('NNpoint_BumpStart:' + str(NNpoint[BumpRangeStart_in_NNpoint]))
            print('NNpoint_BumpEnd:' + str(NNpoint[BumpRangeEnd_in_NNpoint]))
            print('NNpoint in Bump: ' + str(NNpoint[BumpRangeStart_in_NNpoint:BumpRangeEnd_in_NNpoint+1]))
 
            if CCcut_TF == "0.00": 
                ## BandWith For 100 Percent Points
                BandWidthFor100PercentPoints_Option4_2_UnBumpRange      = [0.1113668]*32
                BandWidthFor100PercentPoints_Option4_2_UnBumpRange_more = [0.1113668]*len(x_more)
                print(BandWidthFor100PercentPoints_Option4_2_UnBumpRange[0:BumpRangeStart_in_NNpoint])
                print(BandWidthFor100PercentPoints_Option4_2_UnBumpRange[BumpRangeEnd_in_NNpoint+1:])
            elif CCcut_TF == "0.20":
                ## BandWith For 68 sigma band
                BandWidthForScaled1SigmaBand_Option4_2_UnBumpRange      = [0.12162104119447018]*32
                BandWidthForScaled1SigmaBand_Option4_2_UnBumpRange_more = [0.12162104119447018]*len(x_more)

            y_Bump = CalculateUncertaintyBand_function.PolyFit2Model_ForCurveFit( np.array(NNpoint[BumpRangeStart_in_NNpoint:BumpRangeEnd_in_NNpoint+1]), *popt_bump )
            print("y_Bump:" + str(y_Bump))
            print("y_more_Bump:" + str(y_more_Bump[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more]))
            y_BumpWidth      = y_model_option3 - y_Bump
            y_BumpWidth_more = y_model_option3 - y_more_Bump[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more]         
            print("y_BumpWidth:"          + str(y_BumpWidth))
            print("y_BumpWidth_more:"     + str(y_BumpWidth_more))
            print("y_BumpWidth_more[18]:" + str(y_BumpWidth_more[18]))

            if CCcut_TF == "0.00":
                BandWidth_Option4_2_Final = BandWidthFor100PercentPoints_Option4_2_UnBumpRange[0:BumpRangeStart_in_NNpoint] + y_BumpWidth.tolist() + BandWidthFor100PercentPoints_Option4_2_UnBumpRange[BumpRangeEnd_in_NNpoint+1:]
            elif CCcut_TF == "0.20":
                BandWidth_Option4_2_Final = BandWidthForScaled1SigmaBand_Option4_2_UnBumpRange[0:BumpRangeStart_in_NNpoint] + y_BumpWidth.tolist() + BandWidthForScaled1SigmaBand_Option4_2_UnBumpRange[BumpRangeEnd_in_NNpoint+1:] 
            print('\n')

            #### Option 4.3: BandWidth For Points including a bump
            # Original version
            #BandWidth_Option4_3_FirstConst_more = [0.27444715002994047] * BumpRange_PeakValue_in_x_more  
            #BandWidth_Option4_3_Final           = [0.27444715002994047] * 28 + [0.27157744444156495, 0.19866301960469146, 0.12162104119447018, 0.12162104119447018]  
            ## Test v1
            ShiftValue = 0.05
            BandWidth_Option4_3_FirstConst_more = [0.27444715002994047 - ShiftValue] * BumpRange_PeakValue_in_x_more
            BandWidth_Option4_3_Final           = [0.27444715002994047 - ShiftValue] * 28 + [0.27157744444156495-ShiftValue, 0.19866301960469146-ShiftValue, 0.12162104119447018, 0.12162104119447018] # test v1



            #### Option 5: iMinuit
            popt, pcov = sp.optimize.curve_fit( CalculateUncertaintyBand_function.PolyFit0Model_ForCurveFit, RigidityBinningCenter, RatioUserd, 1.0, RatioErrorUserd, True )
            print("Option4 iMinuit:")
            def PolyFit0Model_ForForMinuit(x, par):
                return par
            least_squares = LeastSquares(RigidityBinningCenter, RatioUserd, RatioErrorUserd, PolyFit0Model_ForForMinuit)
            m = Minuit(least_squares, par=1)
            m.migrad()  # run optimiser
            m.hesse()   # run covariance estimator
            print(m.values)  
            print(m.errors) 
            print("chi2/dof:" + str(m.fval/(len(RigidityBinningCenter) - m.nfit)))
            y, ycov = propagate(lambda p: PolyFit0Model_ForForMinuit(RigidityBinningCenter, p), m.values, m.covariance)
            yerr_prop = np.diag(ycov) ** 0.5
            print("1 Sigma Band:" + str(yerr_prop))

            #### Plotting
            print('\n')  #_AllPointsForNoFurtherRemoved
            print("Plot Section:")
            CalculateUncertaintyBand_function.Plot_ValueAandB(ValueA_AllPointsForNoFurtherRemoved, ValueAError_AllPointsForNoFurtherRemoved, ValueB_AllPointsForNoFurtherRemoved, ValueBError_AllPointsForNoFurtherRemoved, RigidityBinningCenter_AllPointsForNoFurtherRemoved, pattern, NNsuffix, Binningversion, suffix, workpath, Name, CCcut_TF)
            CalculateUncertaintyBand_function.Plot_Ratio     (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, CCcut_TF)
            ## Option 1
            CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option1_LinearFit(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, y_model, x_more, y_more, ConfidenceInterval_show_67, PredictionInterval_show_67, PredictionInterval_show_95, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, CCcut_TF)
            ## Option 2
            #CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option2_PolyFit  (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, yhat_67_P0, upper_67_P0, lower_67_P0, "P0", RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, CCcut_TF)
            #CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option2_PolyFit  (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, yhat_67_P1, upper_67_P1, lower_67_P1, "P1", RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, CCcut_TF)
            #CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option2_PolyFit  (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, yhat_67_P3, upper_67_P3, lower_67_P3, "P3", RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, CCcut_TF)
            #CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option2_PolyFit  (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, yhat_67_P4, upper_67_P4, lower_67_P4, "P4", RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, CCcut_TF)
            ## Option 3
            CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option3_curvefit (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, y_model_option3, x_more, y_more_option3, ParameterError, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, chi2dof, Scaler, True, '1sigma', CCcut_TF)
            CalculateUncertaintyBand_function.Plot_ConfidenceInterval_Option3_curvefit (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, y_model_option3, x_more, y_more_option3, ParameterError, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, chi2dof, Scaler, False, '1sigma', CCcut_TF)
            ## Option 4 (Only for CC cut = 0)
            CalculateUncertaintyBand_function.Plot_Option4_BandWidthFor68PercentPoints          (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidthFor68PercentPoints_more, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, True, y_model_option3, y_more_option3, CCcut_TF)
            CalculateUncertaintyBand_function.Plot_Option4_BandWidthFor68PercentPoints          (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidthFor68PercentPoints_more, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, False, y_model_option3, y_more_option3, CCcut_TF)
            ## Option 4.2
            if CCcut_TF == "0.00":
                BandWidth_UnBumpRange = BandWidthFor100PercentPoints_Option4_2_UnBumpRange_more
            elif CCcut_TF == "0.20": 
                BandWidth_UnBumpRange = BandWidthForScaled1SigmaBand_Option4_2_UnBumpRange
            CalculateUncertaintyBand_function.Plot_Option4_2_BandWidth_WithBump(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidth_UnBumpRange, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, False, y_model_option3, y_more_option3, y_more_Bump, BumpRangeStart_in_x_more, BumpRangeEnd_in_x_more, CCcut_TF)
            CalculateUncertaintyBand_function.Plot_Option4_2_BandWidth_WithBump(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidth_UnBumpRange, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, True, y_model_option3, y_more_option3, y_more_Bump, BumpRangeStart_in_x_more, BumpRangeEnd_in_x_more, CCcut_TF)
            ## Option 4.3
            CalculateUncertaintyBand_function.Plot_Option4_3_SmoothedBandWidth(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidth_UnBumpRange, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, True, y_model_option3, y_more_option3, y_more_Bump, BumpRangeStart_in_x_more, BumpRangeEnd_in_x_more, CCcut_TF, BandWidth_Option4_3_Final, BandWidth_Option4_3_FirstConst_more, BumpRange_PeakValue_in_x_more, FurtherPointRemoved, ShiftValue)


            #### Save Uncertainty Band
            print("\n")
            #ConfidenceInterval = ConfidenceInterval_Option3
            #ConfidenceInterval = BandWidthFor68PercentPoints_Option4 
            #ConfidenceInterval = BandWidth_Option4_2_Final
            ConfidenceInterval = BandWidth_Option4_3_Final
            print("Saving values:")
            print(ConfidenceInterval)
            f_UncertaintyBand = TFile.Open(workpath + "/templatefit/negative/FixedCC/UncertaintyBand_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "RECREATE")
            f_UncertaintyBand.WriteObject( CalculateUncertaintyBand_function.ListToVector(ConfidenceInterval), "ConfidenceInterval")  # Linear Fit 
            #f_UncertaintyBand.WriteObject( CalculateUncertaintyBand_function.ListToVector(kmpCL_67_save_P4)  , "KMPConfidenceInterval_67_P4")  # KMPfit
            f_UncertaintyBand.Close()

if __name__ == '__main__':

    #### Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    parser.add_argument('--pattern'       , help='Which tracker patterns you choose')
    parser.add_argument('--ifVGGNN'       , help='if you want to turn to VGGNN, please turn to Yes.')
    arguments = parser.parse_args()

    if (arguments.ifVGGNN):
        IfVGGNN = arguments.ifVGGNN
    else:
        print("You need to choose which CC estimator you want to use!")
        os._exit(0)
    if IfVGGNN == "No":
        NNsuffix = ""
    elif IfVGGNN == "Yes":
        NNsuffix = "_VGG16NN"

    if (arguments.binningversion):
        Binningversion = arguments.binningversion
    else:
        print("You need to choose a binning version!")
        os._exit(0)

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data version!")
        os._exit(0)

    if ISSversion == "pass7.8":
        suffix = ""
    elif ISSversion == "published2016":
        suffix = "_May2015"
    elif ISSversion == "PhyRep2021":
        suffix = "_Nov2017"

    pattern = arguments.pattern

    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'

    #### To Do: Currently only one element, otherwise the output file (like root and pdf) need to change name.
    #SignalEfficiency_all = ["0.7", "0.8", "0.9", "0.95", "1.0"]
    SignalEfficiency_all = ["1.0"]
    #CCcut_TF_all = ["0.00", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45"];
    #CCcut_TF_all = ["0.00"]
    CCcut_TF_all = ["0.20"]
    #CCcut_TF_all = ["0.20", "0.65", "0.70", "0.80"]
    #CCcut_TF_all = ["0.00", "0.20", "0.40", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95"]

    ## All Range: PointsForFit_Begin = 0, PointsForFit_End = 32  (NNpoint.shape[0]=32, 32 points in total for 525 version )
    # official:
    #PointsForFit_Begin = 20
    #PointsForFit_End   = 32
    # test: (for constant fit y=1.)
    PointsForFit_Begin = 0
    PointsForFit_End   = 32

    if Binningversion == "525version":
        NNbinnings = binning.Newbinnings_525_zhili[26:59] # 14.1-525:published2016binnings[26:59]
        resultdir = "results_525version"
        protonnumberindex = "525version_"
        binname = binning.bins_525_zhili
    elif Binningversion == "450version":
        NNbinnings = binning.published2016binnings[26:58] #  14.1-450:published2016binnings[26:58]
        resultdir = "results_450version"
        protonnumberindex = "450version_"
        binname = binning.bins_450

    NNpoint=np.array( np.zeros(NNbinnings.shape[0]-1) )
    for i in range(NNpoint.shape[0]):
        NNpoint[i]=(NNbinnings[i]+NNbinnings[i+1])/2

    # More points for illustration
    #x_more = np.linspace(np.min(NNpoint), np.max(NNpoint), 100)
    x_more = np.linspace(14, 525, 100)
    main()



