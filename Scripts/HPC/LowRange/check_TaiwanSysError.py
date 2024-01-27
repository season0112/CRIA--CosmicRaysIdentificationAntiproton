#### check_TaiwanSysError.py
## check different terms of systematic errors in Taiwan Low energy antiproton to proton ratio study.

from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import numpy as np

Taiwan_result_6 = TFile("$MY_ANALYSIS/ReferenceFiles/apflux_tme6M.root")

hTErrTmpl_6 = Taiwan_result_6.Get("hTErrTmpl") #p5
hTErrNcnt_6 = Taiwan_result_6.Get("hTErrNcnt") #p6
hTErrRfnc_6 = Taiwan_result_6.Get("hTErrRfnc") #p7
hTErrAccp_6 = Taiwan_result_6.Get("hTErrAccp") #p8
hTErrCC_6   = Taiwan_result_6.Get("hTErrCC")   #p9
hTErrSyst_6 = Taiwan_result_6.Get("hTErrSyst")



#GetNbinsX()=19, GetNbinsY()=18

for i in range(19):
    for j in range(18):
        p5 = hTErrTmpl_6.GetBinContent(i,j)
        p6 = hTErrNcnt_6.GetBinContent(i,j)
        p7 = hTErrRfnc_6.GetBinContent(i,j)
        p8 = hTErrAccp_6.GetBinContent(i,j)
        p9 = hTErrCC_6.GetBinContent(i,j)
        pSys = hTErrSyst_6.GetBinContent(i,j)

        pSys_cal = np.sqrt(p5**2 + p6**2 + p7**2 + p8**2 + p9**2)
        print("i,j=" + str(i) + str(j)) 
        print(pSys-pSys_cal)

