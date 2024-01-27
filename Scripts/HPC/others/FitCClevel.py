#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
import ROOT
from root_numpy import root2array, tree2array
import PythonPlotDefaultParameters
from scipy.optimize import curve_fit
import binning
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.optimize import minimize
import uncertainties
from uncertainties import unumpy
import FitCClevel_functions
import PythonPlotDefaultParameters


def main():
    #### Load CCLevelISS and RigidityResolution
    f_CCLevelISS    = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
    CCLevelISS      = np.asarray(f_CCLevelISS.Get("CCLevelISS_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    CCLevelISSError = np.asarray(f_CCLevelISS.Get("CCLevelISS_Error_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    f_CCLevelISS.Close()

    f_RigidityResolution = TFile.Open(workpath + "/templatefit/negative/FixedCC/RigidityResolution_Pattern_" + str(pattern) + "_" + str(Binningversion) + ".root", "READ")
    Resolution           = np.asarray(f_RigidityResolution.Get("Sigma"))
    ResolutionError      = np.asarray(f_RigidityResolution.Get("SigmaError"))
    f_RigidityResolution.Close()


    #### Fix the first a few points due to limited statistics.
    if str(pattern) == "0":
        print("CCLevelISS:")
        print(CCLevelISS)
        CCLevelISS[0]      = CCLevelISS[2] * 0.6
        CCLevelISSError[0] = CCLevelISSError[2] * 1.8


    #### Linear Fit for ISSCCLevel and Rigidity Model 
    par_Res        = np.polyfit(np.log10(BinningCenter[0:]), np.log10(Resolution), 1)
    par_CCLevel    = np.polyfit(np.log10(BinningCenter[0:]), np.log10(CCLevelISS), 1)
    p_ResModel     = np.poly1d(par_Res)
    p_CCLevelModel = np.poly1d(par_CCLevel)


    #### Residual Fit
    def ResidualModel(alpha, BinningCenter, par_Res, par_CCLevel, CCLevelISS):
        Residual = 0
        for i in range(BinningCenter[1:].shape[0]):
            addterm = ( ( (alpha * par_Res[0] * np.log10(BinningCenter[1:][i]) + par_Res[1]) - par_Res[1] + par_CCLevel[1]) -  np.log10(CCLevelISS[i]) )
            Residual = Residual + np.abs(addterm)
        return Residual

    Result = minimize(ResidualModel, 12.0, args=(BinningCenter, par_Res, par_CCLevel, CCLevelISS))
    print('\n')
    print("#######################################")
    print("Fit Result:")
    print(Result)


    #### Rescale the Resolution model accordiing to the fit
    ResolutionFitLine_Rescale = 10** ( ( Result.x[0] * par_Res[0] * np.log10(BinningCenter[0:]) + par_Res[1] ) - par_Res[1] + par_CCLevel[1] )
    ResolutionFitLine         = 10** (                 par_Res[0] * np.log10(BinningCenter[0:]) + par_Res[1] )

    Resolution_Rescale        = ResolutionFitLine_Rescale / ResolutionFitLine * Resolution
    ResolutionError_Rescale   = ResolutionFitLine_Rescale / ResolutionFitLine * ResolutionError 

    print("\n")
    print("Rigidity Model Error check: (Percentage)")
    print(ResolutionError_Rescale/Resolution_Rescale)


    #### Plot
    FitCClevel_functions.Plot_RigidityResolution              (workpath, BinningCenter, Resolution, ResolutionError, p_ResModel, pattern, ISSversion, Binningversion, NNsuffix)
    FitCClevel_functions.Plot_CCLevelISS                      (workpath, BinningCenter, CCLevelISS, CCLevelISSError, p_CCLevelModel, pattern, ISSversion, Binningversion, NNsuffix)
    FitCClevel_functions.Plot_CCLevelISSOverRigidityResolution(workpath, pattern, ISSversion, Binningversion, BinningCenter, CCLevelISS, CCLevelISSError, Resolution_Rescale, ResolutionError_Rescale, ResolutionFitLine_Rescale, p_CCLevelModel, p_ResModel, NNsuffix)
    FitCClevel_functions.Plot_CCLevelISS_AllPatterns          (workpath, BinningCenter, ISSversion, Binningversion, NNsuffix)


    #### Save Ratio
    CCLevelISS_uncer = unumpy.uarray(CCLevelISS, CCLevelISSError)
    Resolution_uncer = unumpy.uarray(Resolution_Rescale, ResolutionError_Rescale) 
    CCLevelISS_Resolution2_Ratio_uncer = CCLevelISS_uncer/Resolution_uncer

    CCLevelISS_Resolution2_Ratio      = []
    CCLevelISS_Resolution2_RatioError = []
    for i in range(CCLevelISS_Resolution2_Ratio_uncer.shape[0]):
        CCLevelISS_Resolution2_Ratio.append(CCLevelISS_Resolution2_Ratio_uncer[i].n)
        CCLevelISS_Resolution2_RatioError.append(CCLevelISS_Resolution2_Ratio_uncer[i].s)

    #CCLevelISS_Resolution2_RatioError = np.sqrt(CCLevelISS_Resolution2_Ratio) ### FIXME

    f_CCLevelISS_Over_RigidityResolution = TFile.Open(workpath + "/templatefit/negative/FixedCC/ISSCCLevelOverRigidityResolution_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "RECREATE")
    f_CCLevelISS_Over_RigidityResolution.WriteObject( FitCClevel_functions.ListToVector( CCLevelISS_Resolution2_Ratio),      "CCLevelISSOverResolution")
    f_CCLevelISS_Over_RigidityResolution.WriteObject( FitCClevel_functions.ListToVector( CCLevelISS_Resolution2_RatioError), "CCLevelISSOverResolution_Error")
    f_CCLevelISS_Over_RigidityResolution.Close()


if __name__ == '__main__':
    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern'       , help='Which tracker patterns you choose')
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    parser.add_argument('--ifVGGNN', help='if you want to turn to VGGNN, please turn to Yes.')
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

    pattern = arguments.pattern
    ISSversion = arguments.issversion
    Binningversion = arguments.binningversion

    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'
    BinningCenter = binning.Newbinnings_525_center_zhili[26:]

    main()


'''
#Obsolete:
#### Resolution Fit
def ResolutionModel(R, alpha):
    Resolution_Rescale = []
    for i in range(len(R)):
        Resolution_Rescale.append( (alpha * par_Res[0] * np.log10(BinningCenter[i]) + par_Res[1]) - par_Res[1] + par_CCLevel[1] )
    return Resolution_Rescale

least_squares = LeastSquares(BinningCenter[0:], np.log10(CCLevelISS), 0, ResolutionModel, loss="soft_l1")
minuit = Minuit(least_squares, alpha=12)
minuit.errordef = 1.0
minuit.migrad()
minuit.hesse()
print("fit result: ")
print(minuit.values)
'''

