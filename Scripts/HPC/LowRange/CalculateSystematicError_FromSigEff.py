import numpy as np
import uproot
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F, gStyle
import binning
import TimeDependentTimeRange
import argparse
import os
import matplotlib.pyplot as plt
import CalculateSystematicError_FromSigEff_Plot
import ROOT
import math
import PythonPlotDefaultParameters

def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector


def main():

    #### Loop in TimeIndex
    for TimeIndex in range(0, SplitTotal, mergestep):  ## If time dependent, loop in time index, for time averaged, it's 1.

        print('\n')
        print("\033[1;31;40mTimeIndex: \033[0m" + str(TimeIndex))

        #### For 3Bartels, 123-126 is empty, remove this point.
        if TimeIndex == 123 and TimeMode == "TimeDependent" and TimeSplitMode == "3BartalRotation": 
            if RigidityRange == 'low':
                f_SysErr = TFile(LowPath          + "/totalall/SysError_Eff/SysErr_SigEff_" + TimeSplitMode + "_" + str(TimeIndex) + ".root", "RECREATE")
                f_SysErr.Close()
            elif RigidityRange == 'intermediate':
                f_SysErr = TFile(IntermediatePath + "/total/SysError_Eff/SysErr_SigEff_" + TimeSplitMode + "_" + str(TimeIndex) + ".root", "RECREATE")
                f_SysErr.Close()
            continue

        SystemError_SigEff = []
        OfficialTRDEff_All = []
        OfficialTOFEff_All = []


        #### Load RootFile for Ratio
        if TimeMode == "TimeAveraged":
            if RigidityRange == 'low':
                RatioRootFile = uproot.open(LowPath          + "/totalall/Time_Averaged_ratio_Low/binmerge1/plots/Ratio_" + ISSversion + ".root")
            elif RigidityRange == 'intermediate':
                RatioRootFile = uproot.open(IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/intermediate_0124_free_" + ISSversion + "binmerge1.root")
        elif TimeMode == "TimeDependent":
            if RigidityRange == 'low':
                RatioRootFile = uproot.open(LowPath          + "/totalall/SysError_Eff/Ratio_" + TimeSplitMode + "_TimeIndex" + str(TimeIndex) + ".root")
            elif RigidityRange == 'intermediate':
                RatioRootFile = uproot.open(IntermediatePath + "/total/SysError_Eff/Ratio_"    + TimeSplitMode + "_TimeIndex" + str(TimeIndex) + ".root")

        #### Loop Rigidity
        for rigidityindex in [6]:
        #for rigidityindex in range( UsedBinningCenter.shape[0] ):
            print("\n")
            print("\033[1;35mNow Rigidityindex is \033[0m" + str(rigidityindex)) 
            print("Rigidity is " +str(UsedBinningEdge[rigidityindex]))

            #### Fill list of ratio and efficiency
            l_TotalEff   = []
            l_TrdEff     = []
            l_TofEff     = []
            ratioall     = []


            if RigidityRange == 'low':
                if rigidityindex < 5:
                    TrdEfficiencyAll = np.arange(0.78, 0.96, 0.02)
                    TofEfficiencyAll = np.arange(0.83, 0.98, 0.02)
                else:
                    TrdEfficiencyAll = np.arange(0.60, 0.96, 0.02)
                    TofEfficiencyAll = np.arange(0.61, 0.99, 0.02) 
            elif RigidityRange == 'intermediate':
                TrdEfficiencyAll = np.arange(0.60, 1.00, 0.01)


            #### Loop in Trd Efficiency
            for trdeff in TrdEfficiencyAll:
                #print("trdeff:" + str("{:.2f}".format(trdeff)) )
                if RigidityRange == 'low':
                    #### Loop in Tof Efficiency
                    for tofeff in TofEfficiencyAll:
                        #print("tofeff:" + str("{:.2f}".format(tofeff)) )
                        if TimeMode == "TimeAveraged":
                            ratio = RatioRootFile["ratio_tof_with_effective_TRDeff_" + str("{:.2f}".format(trdeff)) + "_TOFeff_" + str("{:.2f}".format(tofeff)) ].member('fY')[rigidityindex]
                        elif TimeMode == "TimeDependent":
                            ratio = RatioRootFile["g_ratio_tof_TRDeff_" + str("{:.2f}".format(trdeff)) + "_TOFeff_" + str("{:.2f}".format(tofeff)) ].member('fY')[rigidityindex]

                        #l_TotalEff.append(trdeff * tofeff)
                        l_TotalEff.append( (trdeff+0.02) * (tofeff+0.03) ) ## FIXME: In order to reach 100%.

                        l_TrdEff.append(trdeff)
                        l_TofEff.append(tofeff)
                        ratioall.append(ratio)
                elif RigidityRange == 'intermediate':
                    if TimeMode == "TimeAveraged":
                        ratio = RatioRootFile[ "g_ratio_with_effective_acceptance_TRDeff_" + str("{:.2f}".format(trdeff)) ].member('fY')[rigidityindex]
                    elif TimeMode == "TimeDependent":
                        ratio = RatioRootFile[ "g_ratio_TRDeff_" + str("{:.2f}".format(trdeff)) ].member('fY')[rigidityindex]
                    l_TotalEff.append(trdeff)
                    l_TrdEff.append(trdeff)
                    ratioall.append(ratio)


            #### Remove Nan result ratio And Assign 0 if all results are Nan.
            print("ratioall lens: " + str(len(ratioall)))
            ratioall = [ x for x in ratioall if math.isnan(x) == False ]
            print("ratioall lens: " + str(len(ratioall)))
            if len(ratioall) == 0:
                ratioall = [0]
            #print('min:' + str(np.min(np.array(ratioall))))
            #print('max:' + str(np.max(np.array(ratioall))))


            #### Convert list into Histogram
            hist_ratio = TH1D("", "", BinsNumber, min(ratioall), max(ratioall))
            for i in ratioall:
                hist_ratio.Fill(i)

            #print("mean is " + str(hist_ratio.GetMean()))
            #print("all are :" + str(ratioall))

            mean = hist_ratio.GetMean()
            OfficianlValue      = min(ratioall, key=lambda x:abs(x-mean))
            OfficianlValueIndex = ratioall.index(OfficianlValue)
            #print("OfficianlValue is " + str(OfficianlValue))

            #### Find Settings for best fit
            OfficialTRDEff     = l_TrdEff[OfficianlValueIndex]
            OfficialTRDEff_All.append(OfficialTRDEff)
            print("OfficialTRDEff:" + str(OfficialTRDEff))
            if RigidityRange == 'low':
                OfficialTOFEff     = l_TofEff[OfficianlValueIndex]
                OfficialTOFEff_All.append(OfficialTOFEff)
                print("OfficialTOFEff:" + str(OfficialTOFEff))
                print("OfficialEff: " + str( float(OfficialTRDEff) * float(OfficialTOFEff) ))

            ## check
            #print("check:")
            #if RigidityRange == 'low':
            #    print( RatioRootFile["ratio_tof_with_effective_TRDeff_" + str("{:.2f}".format(OfficialTRDEff)) + "_TOFeff_" + str("{:.2f}".format(OfficialTOFEff)) ].member('fY')[rigidityindex] )

            #### Calculate RMS
            SystemError_SigEff.append(hist_ratio.GetRMS())

            #### Plot ratio vs signal efficiency (Only have no Nan result in original ratioall, otherwise dimensions don't match)
            CalculateSystematicError_FromSigEff_Plot.Plot_RatiovsSignalEfficiency(l_TotalEff, ratioall, UsedBinningEdge[rigidityindex], RigidityRange, TimeMode, TimeSplitMode, TimeIndex)

            #### Plot RMS as Systematic Error 
            CalculateSystematicError_FromSigEff_Plot.Plot_RMS( hist_ratio, UsedBinningEdge[rigidityindex], RigidityRange, TimeMode, TimeSplitMode, TimeIndex)

            #### Clear
            hist_ratio.Reset()


        ## Make TGraph for Systematic Error
        #print(SystemError_SigEff)
        SystemError_SigEff = np.array(SystemError_SigEff)
        g_SysError_SigEff  = TGraph(UsedBinningCenter.shape[0], UsedBinningCenter, SystemError_SigEff)


        #### Plot RMS vs Rigidity
        #CalculateSystematicError_FromSigEff_Plot.Plot_RMSvsRigidity(UsedBinningCenter, SystemError_SigEff, RigidityRange, TimeMode, TimeSplitMode, TimeIndex)

        '''
        #### Save Systematic Error, and Official TRD & TOF Efficiency
        if TimeMode == "TimeAveraged":
            if RigidityRange == 'low':
                f_SysErr = TFile(LowPath          + "/totalall/Time_Averaged_ratio_Low/binmerge1/plots/SysErr_SigEff_" + ISSversion + ".root", "RECREATE")
                f_SysErr.WriteObject( ListToVector(OfficialTRDEff_All), "OfficialTRDEff")
                f_SysErr.WriteObject( ListToVector(OfficialTOFEff_All), "OfficialTOFEff")
            elif RigidityRange == 'intermediate':
                f_SysErr = TFile(IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/SysErr_SigEff_" + ISSversion + ".root", "RECREATE")
                f_SysErr.WriteObject( ListToVector(OfficialTRDEff_All), "OfficialTRDEff")
        elif TimeMode == "TimeDependent":
            if RigidityRange == 'low':
                f_SysErr = TFile(LowPath          + "/totalall/SysError_Eff/SysErr_SigEff_" + TimeSplitMode + "_" + str(TimeIndex) + ".root", "RECREATE")
                f_SysErr.WriteObject( ListToVector(OfficialTRDEff_All), "OfficialTRDEff")
                f_SysErr.WriteObject( ListToVector(OfficialTOFEff_All), "OfficialTOFEff")
            elif RigidityRange == 'intermediate':        
                f_SysErr = TFile(IntermediatePath + "/total/SysError_Eff/SysErr_SigEff_" + TimeSplitMode + "_" + str(TimeIndex) + ".root", "RECREATE")
                f_SysErr.WriteObject( ListToVector(OfficialTRDEff_All), "OfficialTRDEff")
        g_SysError_SigEff.Write("g_SysError_SigEff")
 
        f_SysErr.Close()
        '''


if __name__ == '__main__':

    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--Range'         , help='Rigidity Range, low, intermediate or high')
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or 2016paper')
    parser.add_argument('--TimeMode'      , help='You want to use time averaged or which time dependent data')
    parser.add_argument('--TimeSplitMode' , help='which time slit mode you choose')
    arguments = parser.parse_args()

    if (arguments.Range and arguments.TimeMode):
        RigidityRange  = arguments.Range
        TimeMode       = arguments.TimeMode
    else:
        print("You need to provide all parameters!")
        os._exit(0)

    TimeSplitMode = arguments.TimeSplitMode
    ISSversion     = arguments.issversion
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

    BinsNumber = 50

    if TimeMode == "TimeAveraged":
        UsedBinningEdge   = BinningAll
        UsedBinningCenter = BinningCenterAll
        SplitTotal = 1
        mergestep = 1

    elif TimeMode == "TimeDependent":
        if RigidityRange == 'intermediate':
            UsedBinningEdge   = binning.BinsEdgeName_TimeDependent_Intermediate
            UsedBinningCenter = binning.BinsCenter_TimeDependent_Intermediate
        elif RigidityRange == 'low':
            UsedBinningEdge   = binning.BinsEdgeName_TimeDependent_Low
            UsedBinningCenter = binning.BinsCenter_TimeDependent_Low

        if TimeSplitMode == "3BartalRotation":
            SplitTotal = TimeDependentTimeRange.SplitTotal_3B
            mergestep  = TimeDependentTimeRange.mergestep_3B 
        elif TimeSplitMode == "6BartalRotation":
            SplitTotal = TimeDependentTimeRange.SplitTotal_6B
            mergestep  = TimeDependentTimeRange.mergestep_6B
        elif TimeSplitMode == "6months":
            SplitTotal = TimeDependentTimeRange.SplitTotal_6M
            mergestep  = TimeDependentTimeRange.mergestep_6M


    ## Set1: (0.70-0.96, step:0.02)
    ## use this one:
    #TrdEfficiencyAll = np.arange(0.7, 0.95, 0.02)
    #TofEfficiencyAll = np.arange(0.90, 0.95, 0.02)
    ## others:
    #TofEfficiencyAll = np.arange(0.7, 0.95, 0.02)
    #TofEfficiencyAll = np.arange(0.8, 0.95, 0.02)
    #TofEfficiencyAll = np.arange(0.88, 0.95, 0.02)
    #TofEfficiencyAll = np.arange(0.90, 0.95, 0.02)
    #TofEfficiencyAll = np.arange(0.94, 0.95, 0.02)

    ## Set2:
    # For intermediate:
    #TrdEfficiencyAll = np.arange(0.6, 1.00, 0.01)
    # For low: 
    #TrdEfficiencyAll = np.arange(0.60, 1.00, 0.02)
    #TofEfficiencyAll = np.arange(0.61, 0.99, 0.02)
    # To remove bad result
    #TrdEfficiencyAll = np.arange(0.60, 0.96, 0.02)
    #TofEfficiencyAll = np.arange(0.61, 0.99, 0.02)


    main()


