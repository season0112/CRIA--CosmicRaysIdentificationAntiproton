import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import binning
from ROOT import TFile, TH1D, TCanvas, gPad, gStyle, TPad, TLegend
import PythonPlotDefaultParameters
import TimeDependentMeasuringTime_tool

def main():

    MeasuringTime      = []
    Measuringtime_allR = []

    for i in range(23):
        ROOTFile  = TFile( "/hpcwork/jara0052/sichen/Measuringtime/TimeDependent/MeasuringTime_pass7.8_11_2021_GEOMETRIC25_1.2_" + str(Time6B_Edge[i]) + "_" + str(Time6B_Edge[i+1]) + ".root", "OPEN")
        h1 = ROOTFile.Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff") # 62 bins

        MeasuringTime.append(h1.GetBinContent(50))

        h1.Rebin(2) # 31 bins
        #print(h1.GetBinContent(10))  # h1->GetBinLowEdge(10) = 6.47 GV; h1->GetBinLowEdge(25) = 69.7 GV
        tem = []
        for j in range(h1.GetNbinsX()):
            tem.append(h1.GetBinContent(j+1))
        Measuringtime_allR.append(tem)


    ## Plot
    TimeDependentMeasuringTime_tool.Plot_TimeDependentMeasuringTime(Time6B, MeasuringTime)

    ## Save
    np.save('/rwthfs/rz/cluster/hpcwork/jara0052/sichen/Measuringtime/TimeDependent/TimeDependentMeasuringTime_AllRigidity_6B.npy', Measuringtime_allR)
    


if __name__ == '__main__':

    Time6B_Edge = binning.Bartals6Unixtime_edge
    Time6B = [datetime.fromtimestamp(i) for i in binning.Bartals6Unixtime[0:23]]

    main()



