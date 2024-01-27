import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F, gStyle
import binning
import argparse
import os
import matplotlib.pyplot as plt
import PythonPlotDefaultParameters
import CalculateSysErrorFromShape_Plot

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def FillinHist(rootfile, BinsNumber, Ratio_min, Ratio_max, generatenumber, RigidityRange, SysErrorShape_original, RIndex, BinningAll, Timemode):
        hist_ratio = TH1D("", "", BinsNumber, Ratio_min, Ratio_max)
        ## Fill ratio histogram
        for numberindex in range(generatenumber):
            if RigidityRange == 'low':
                if Timemode == "TimeAveraged":
                    g_Ratio = rootfile.Get( "ratio_tof_with_effective" + str(numberindex) )
                elif Timemode == "TimeDependent":
                    g_Ratio = rootfile.Get( "g_ratio_tof" + str(numberindex) )
            elif RigidityRange == 'intermediate':
                g_Ratio = rootfile.Get( "g_ratio_with_effective_acceptance" + str(numberindex) )

            if RigidityRange == 'low':
                hist_ratio.Fill(g_Ratio.GetY()[RIndex])
            elif RigidityRange == 'intermediate':
                hist_ratio.Fill(g_Ratio.GetY()[RIndex]*100000)
        ## Plot
        #CalculateSysErrorFromShape_Plot.Plot_ratio_histogram(hist_ratio, BinningAll, RIndex, RigidityRange)
        ## Save RMS
        if RigidityRange == 'low':
            SysErrorShape_original.append(hist_ratio.GetRMS())
        elif RigidityRange == 'intermediate':
            SysErrorShape_original.append(hist_ratio.GetRMS()/100000)
        ## Clear
        hist_ratio.Reset()



