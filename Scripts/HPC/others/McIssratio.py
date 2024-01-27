import numpy as np
import matplotlib.pyplot as plt
from root_numpy import array2tree, array2root, fill_hist, root2array
from ROOT import TFile, TH1D, TCanvas, gPad, gStyle, TPad, TLegend
import argparse
import McIssratio_tool
import binning
import uproot

def main():

    #### MC data
    PbarNumber_MC = np.array([])
    ## Low and Intermediate
    for Energybin in binning.binslow[1:-2]:
        print("Rigidity:" + str(Energybin))
        Number_MC    = root2array('/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/' + mcdataset + '_Tree_' + str(Energybin) + '.root', 'AntiprotonLowEnergyTree').shape[0]
        PbarNumber_MC = np.append(PbarNumber_MC, Number_MC)
    print(PbarNumber_MC.shape)
    ## High
    for Energybin in binning.bins[1:]:
        print("Rigidity:" + str(Energybin))
        Number_MC     = root2array('/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/' + mcdataset + '_Tree_' + str(Energybin) + '.root', 'ExampleAnalysisTree').shape[0]
        PbarNumber_MC = np.append(PbarNumber_MC, Number_MC)
    print(PbarNumber_MC.shape)


    #### ISS data
    PbarNumber_ISS = np.array([])
    file = uproot.open("/home/bo791269/Software/AntiprotonAnalysis/Macros/CheckResult.root")    
    PbarNumber_Low          = file["g_Antiproton_number_unfolded_Low"]         .values()[1]
    PbarNumber_Intermediate = file["g_Antiproton_number_unfolded_Intermediate"].values()[1]
    PbarNumber_High         = file["g_Antiproton_number_unfolded_High"]        .values()[1]
    PbarNumber_ISS = np.append(PbarNumber_ISS, PbarNumber_Low)
    PbarNumber_ISS = np.append(PbarNumber_ISS, PbarNumber_Intermediate)
    PbarNumber_ISS = np.append(PbarNumber_ISS, PbarNumber_High)
    print("PbarNumber_ISS:" + str(PbarNumber_ISS))
    #PbarNumber_ISS = np.array([ 3428, 3031, 2659, 2362, 2167, 1864, 1615, 1438, 1182, 1052, 982, 793, 645, 622, 543, 465, 403, 374, 329, 284, 230, 401, 304, 214, 189, 180, 104, 108, 85, 94]) #pass7
    #CCPNumber_ISS  = np.array([ 26.1449, 9.84125, 25.1866, 3.91486, 34.7692, 42.5148, 56.1257, 62.6984, 66.1572, 66.8996, 87.4114, 85.8961, 86.1828, 88.9212, 134.723, 138.889, 143.871, 145.723, 138.386, 165.07, 157.301, 484.817, 614.611, 918.961, 1367.93, 2113.02, 2917.31, 4122.84, 5913.62, 12257.4]) #pass7ext no cccut.


    #### Calculate Ratio
    if mcdataset == "B1042_pr.pl1.flux.l1a9.2016000_7.6_all":
        MCccproton = PbarNumber_MC
        ratio      = MCccproton/CCPNumber_ISS
    else:
        ratio      = PbarNumber_MC/PbarNumber_ISS

    
    #### plot
    NNbinnings     = binning.Newbinnings_525_zhili
    NNbinningWidth = binning.NewbinningWidth_525_zhili

    TH_ratio = TH1D("ratio", "", 58, NNbinnings)
    TH_iss   = TH1D("iss"  , "", 58, NNbinnings)
    TH_mc    = TH1D("mc"   , "", 58, NNbinnings)

    for i in range(1,59):
        if mcdataset == "B1042_pr.pl1.flux.l1a9.2016000_7.6_all":
            TH_iss.SetBinContent(i, CCPNumber_ISS[i-1]/NNbinningWidth[i-1])
            TH_mc .SetBinContent(i, MCccproton[i-1]   /NNbinningWidth[i-1])
        else:
            TH_iss.SetBinContent(i, PbarNumber_ISS[i-1]/NNbinningWidth[i-1])
            TH_mc .SetBinContent(i, PbarNumber_MC[i-1] /NNbinningWidth[i-1])
        TH_ratio.SetBinContent(i,ratio[i-1])

    McIssratio_tool.plot(TH_ratio, TH_iss, TH_mc, mcdataset)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--mcdata', help='which mc dataset')
    arguments = parser.parse_args()

    mcdataset = arguments.mcdata

    main()













