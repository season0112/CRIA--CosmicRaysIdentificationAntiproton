import draw
from draw import *

def drawPosOverEleRatio(lepton_result, lepton_result_Maura, lowworkpath, intermediateworkpath, rigidityrange, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, ElectronPositronIndex, binmerge, PosiOverElec, PosiOverElec_error, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura):
    import asireader
    reader = asireader.ASIReader_E()
    LeptonLowEnergyBin = binning.LeptonLowEnergyBin

    MauraTimeRange_all          = binning.Bartals1Unixtime[1:103] #102 bins in total in Maura result, including 2 empty bins.
    MauraTimeRange_WithoutEmpty = np.delete(np.array(MauraTimeRange_all),[46,47]) #100 bins have results
    dates_1Bartel_Maura         = [datetime.fromtimestamp(i) for i in MauraTimeRange_WithoutEmpty]

    Positron_Maura   = np.array(Positron_Maura)
    Electron_Maura   = np.array(Electron_Maura)
    PosOverEle_Maura = np.array(Positron_Maura)/np.array(Electron_Maura)

    ### load path
    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"

    # Load PRL result
    pos_ele = reader.read('/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/TimeDependentAnalysis/TimeDepedentEleAndPosPRL/pos_ele_AMS_PRL2018_ekin_036.xml')

    # time averaged (Niko)
    lepton_result = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/LeptonFluxes_Zimmermann_04_09_2018_2D.root")
    grElectron_RWTH_FinalFlux = lepton_result.Get("SxFluxTools_Electron_RWTH").Get("grElectron_RWTH_FinalFlux") # TGraphAsymmErrors; grElectron_RWTH_FinalFlux.GetN()=74
    grPositron_RWTH_FinalFlux = lepton_result.Get("SxFluxTools_Positron_RWTH").Get("grPositron_RWTH_FinalFlux") # TGraphAsymmErrors;
    grPosiOverElec_RWTH = lepton_result.Get("SxRatioTools_PosiOverElec_RWTH").Get("grPosiOverElec_RWTH_PosiOverElec")
    grElectron_RWTH_RelErrorTot = lepton_result.Get("SxFluxTools_Electron_RWTH").Get("grElectron_RWTH_RelErrorTot") 
    grPositron_RWTH_RelErrorTot = lepton_result.Get("SxFluxTools_Positron_RWTH").Get("grPositron_RWTH_RelErrorTot") 
    grPosiOverElec_RWTH_RelErrorTot = lepton_result.Get("SxRatioTools_PosiOverElec_RWTH").Get("grPosiOverElec_RWTH_RelErrorTot")

    electron_to_positron_error = np.sqrt( ((Electron_error_relative*Electron)/Positron)**2 + (Electron/Positron**2*(Positron_error_relative*Positron))**2)
    positron_to_electron_error = np.sqrt( ((Positron_error_relative*Positron)/Electron)**2 + (Positron/Electron**2*(Electron_error_relative*Electron))**2)

    electron_to_positron_error_Maura = np.sqrt( ((ElectronError_Maura)/Positron_Maura)**2 + (Electron_Maura/Positron_Maura**2*(PositronError_Maura))**2)
    positron_to_electron_error_Maura = np.sqrt( ((PositronError_Maura)/Electron_Maura)**2 + (Positron_Maura/Electron_Maura**2*(ElectronError_Maura))**2)

    #### GetErrorY method
    Pos = []
    Pos_Error = []
    Pos_Error_Rel = []
    Ele = []
    Ele_Error = []
    Ele_Error_Rel = []
    for i in range(grPositron_RWTH_FinalFlux.GetN()):
        Pos.append(grPositron_RWTH_FinalFlux.GetY()[i])
        Pos_Error.append(grPositron_RWTH_FinalFlux.GetErrorY(i))
        Pos_Error_Rel.append(grPositron_RWTH_RelErrorTot.GetY()[i])
    for i in range(grElectron_RWTH_FinalFlux.GetN()):
        Ele.append(grElectron_RWTH_FinalFlux.GetY()[i])
        Ele_Error.append(grElectron_RWTH_FinalFlux.GetErrorY(i))
        Ele_Error_Rel.append(grElectron_RWTH_RelErrorTot.GetY()[i])
    PosOverEle_Error = np.sqrt( np.power(Pos_Error,2)/np.power(Ele,2) + np.power(Ele_Error,2) * np.power(Pos,2) / np.power(Ele,4) )

    #### GetRelativeError method
    PosOverEle_Error_FromRelative = []
    for i in range(grPosiOverElec_RWTH_RelErrorTot.GetN()):
        PosOverEle_Error_FromRelative.append( grPosiOverElec_RWTH_RelErrorTot.GetY()[i] * grPosiOverElec_RWTH.GetY()[i] )

    # Time averaged e+/e
    c1 = TCanvas()
    grPosiOverElec_RWTH.Draw("")
    xaxis_ratio = grPosiOverElec_RWTH.GetXaxis();
    yaxis_ratio = grPosiOverElec_RWTH.GetYaxis();
    xaxis_ratio.SetLimits(0.5,50);
    yaxis_ratio.SetRangeUser(0.04,0.2);
    c1.SetLogx()
    c1.SetLogy()
    c1.SaveAs( path + "Time_Averaged_PosiOverElec.png")

    # Time dependent e+/e
    fig_PosOverEle, ax_PosOverEle, = plt.subplots(figsize=(30,18))
    plt.errorbar(np.delete(draw.dates_1Bartel,[46,47]), np.delete(Positron/Electron,[46,47]), yerr=np.delete(positron_to_electron_error,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Niko flux analysis)")
    plt.errorbar(np.delete(draw.dates_1Bartel,[46,47]), np.delete(PosiOverElec,[46,47]), yerr=np.delete(PosiOverElec_error,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="b",ecolor="b", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Niko ratio analysis)")
    plt.errorbar(np.delete(draw.dates_1Bartel,[46,47])[0:79], pos_ele.fluxratio, yerr=pos_ele.err_total, marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(PRL ratio analysis)")
    plt.errorbar(dates_1Bartel_Maura, PosOverEle_Maura, yerr=positron_to_electron_error_Maura, marker='s', linestyle="None", markersize=12, markerfacecolor="magenta",ecolor="magenta", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Maura flux analysis)")

    plt.legend(loc='best',fontsize=40)
    ax_PosOverEle.axes.set_ylim([0.035, 0.082])
    plt.ylabel(r'$\rm{e^+}$/$\rm{e^-}$',fontsize=70)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    plt.savefig( path + "Positron_Over_Electron_" + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + ".png")

    # Time dependent e+/e relative error
    fig = plt.figure(figsize=(38,18))
    ax_RelErr = fig.add_subplot(111)
    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]),  np.delete(positron_to_electron_error,[46,47])/np.delete(Positron/Electron,[46,47])*100, yerr=0, marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Niko flux analysis)")
    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]), np.delete(PosiOverElec_error,[46,47])/np.delete(PosiOverElec,[46,47])*100, yerr=0, marker='s', linestyle="None", markersize=12, markerfacecolor="b",ecolor="b", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Niko ratio analysis)")
    plt.errorbar(np.delete(draw.dates_1Bartel,[46,47])[0:79], pos_ele.err_total/pos_ele.fluxratio*100, yerr=0, marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(PRL ratio analysis)")
    plt.errorbar(dates_1Bartel_Maura, positron_to_electron_error_Maura/PosOverEle_Maura*100, yerr=0, marker='s', linestyle="None", markersize=12, markerfacecolor="magenta",ecolor="magenta", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Maura flux analysis)")

    ax_RelErr.axes.set_ylim([0, 7])
    plt.legend(loc='best',fontsize=40)
    ax_RelErr.tick_params(axis='y', direction='in', length=12, width=10, colors='black')
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.ylabel(r'$\rm{e^+}$/$\rm{e^-}$' + ' relative error (%)',fontsize=55)
    plt.savefig( path + "Positron_Over_Electron_RelativeError" + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + ".png")





