import draw 
from draw import *


def Plot_old_bin_ElectronPositron(Electron, Positron, electron_to_positron_error, LeptonLowEnergyBin, ElectronPositronIndex, dates_1Bartel_Maura, EleOverPos_Maura, EleOverPosRatioError_Maura, yaxisfactor, path): 

    fig, ax, = plt.subplots(figsize=(30,18))

    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]), np.delete(Electron/Positron,[46,47]), yerr=np.delete(electron_to_positron_error,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="b",ecolor="b", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Original) (Niko)")
    plt.errorbar( dates_1Bartel_Maura                  , EleOverPos_Maura                    , yerr=EleOverPosRatioError_Maura                   , marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Original) (Maura)")
    #plt.errorbar(dates_1Bartel_Maura                  , ratio_Perugia                       , yerr=ratioerror_Perugia                           , marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Maura)")
    #plt.errorbar(np.delete(draw.dates_1Bartel,[46,47]), ratio_Aachen                        , yerr=ratioerror_Aachen                            , marker='s', linestyle="None", markersize=12, markerfacecolor="y",ecolor="y", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Niko)")

    plt.ylabel(r'$\rm{e^-}$/$\rm{e^+}$',fontsize=80)
    plt.legend(loc='best',fontsize=50)
    ymin = plt.ylim()[0]
    ymax = plt.ylim()[1]
    ax.axes.set_ylim([ ymin-((ymax-ymin)/2)*yaxisfactor,  ymax+((ymax-ymin)/2)*yaxisfactor ])
    #ax.set_ylabel('Electron to Positron ratio',fontsize=55, color="r")
    ax.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)

    plt.savefig(path + "Electron_Over_Positron_ElectronPositronIndex_" + str(ElectronPositronIndex) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + ".pdf")

    '''
    ## Fit with fermi function
    params2, params_covariance2 = optimize.curve_fit(Fermi_func,  np.delete(draw.dates_1Bartel_unixsecond,[46,47]),  np.delete(Electron/Positron,[46,47]) , maxfev=10000)
    plt.plot( draw.dates_1Bartel[0:46]+draw.dates_1Bartel[48:] ,  Fermi_func( np.delete(draw.dates_1Bartel_unixsecond,[46,47]), params2[0], params2[1],params2[2], params2[3]) ,linewidth=5, color='m')
    plt.savefig(path + "Electron_Over_Positron_ElectronPositronIndex_" + str(ElectronPositronIndex) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + "_withFit.pdf")
    '''
    plt.close()


def Plot_old_bin_ElectronPositron_RelatvieError(electron_to_positron_error, Electron, Positron, dates_1Bartel_Maura, EleOverPosRatioError_Maura, EleOverPos_Maura, LeptonLowEnergyBin, ElectronPositronIndex, path):

    fig, ax = plt.subplots(figsize=(30,18))

    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]), np.delete(electron_to_positron_error,[46,47])/np.delete(Electron/Positron,[46,47])*100, yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Original) (Niko)" )
    plt.errorbar( dates_1Bartel_Maura                  , EleOverPosRatioError_Maura/(EleOverPos_Maura)*100                                     , yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="g"  , ecolor="g"  , label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Original) (Maura)")

    plt.ylabel(r'$\rm{e^-}$/$\rm{e^+}$' + ' relative error (%)',fontsize=80)
    plt.legend(loc='best', fontsize=50)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)

    plt.savefig(path + "Electron_Over_Positron_RelativeError_ElectronPositronIndex_" + str(ElectronPositronIndex) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + ".pdf")

    plt.close()


def Plot_new_bin_ElectronPositron(ratio_Aachen, ratioerror_Aachen, index, binmerge, dates_1Bartel_Maura, ratio_Perugia, ratioerror_Perugia, path, yaxisfactor):

    fig, ax = plt.subplots(figsize=(30,18))

    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]), ratio_Aachen , yerr=ratioerror_Aachen , marker='s', linestyle="None", markersize=12, markerfacecolor="b",ecolor="b", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Niko)")
    plt.errorbar( dates_1Bartel_Maura                  , ratio_Perugia, yerr=ratioerror_Perugia, marker='s', linestyle="None", markersize=12, markerfacecolor="g",ecolor="g", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Maura)")

    plt.legend(loc='best',fontsize=50)
    ymin = plt.ylim()[0]
    ymax = plt.ylim()[1]

    ax.axes.set_ylim([ ymin-((ymax-ymin)/2)*yaxisfactor,  ymax+((ymax-ymin)/2)*yaxisfactor ])
    #ax.set_ylabel('Electron to Positron ratio',fontsize=55, color="r")
    ax.tick_params(direction='in',length=10, width=3)
    #plt.xticks( rotation=40 )
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)

    plt.savefig( path + "Electron_Over_Positron_" + str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + ".pdf")


def Plot_new_bin_ElectronPositron_RelatvieError(ratioerror_Aachen, ratio_Aachen, index, binmerge, dates_1Bartel_Maura, ratioerror_Perugia, ratio_Perugia, path):
    fig = plt.figure(figsize=(38,18))
    ax  = fig.add_subplot(111)

    plt.errorbar( np.delete(draw.dates_1Bartel,[46,47]), ratioerror_Aachen/ratio_Aachen*100  , yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="b"  , ecolor="b"  , label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Niko)")
    plt.errorbar( dates_1Bartel_Maura                  , ratioerror_Perugia/ratio_Perugia*100, yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Maura)" )

    ax.tick_params(axis='y', direction='in', length=12, width=10, colors='black')
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.ylabel(r'$\rm{e^-}$/$\rm{e^+}$' + ' relative error',fontsize=55)
    plt.legend(loc='best',fontsize=50)
    plt.text(0.7, 0.9, format(binning.published2016binnings[index], '.1f') + "-" + format(binning.published2016binnings[index+int(binmerge)], '.1f') + ' GV', fontsize=80, transform=ax.transAxes)

    plt.savefig( path + "EleOverPos_RelativieError_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    plt.close()


def Plot_old_new_bin_ElectronPositron_RelatvieErrorComparision(dates_1Bartel_Maura, EleOverPosRatioError_Maura, EleOverPos_Maura, ElectronPositronIndex, ratioerror_Perugia, ratio_Perugia, index, binmerge, path, LeptonLowEnergyBin):

    fig = plt.figure(figsize=(38,18))
    ax  = fig.add_subplot(111)

    plt.errorbar( dates_1Bartel_Maura, EleOverPosRatioError_Maura/(EleOverPos_Maura)*100, yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="g"  , ecolor="g"  , label=str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) + " GV" +"(Original) (Maura)")
    plt.errorbar( dates_1Bartel_Maura, ratioerror_Perugia/ratio_Perugia*100             , yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", label=str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)]) + " GV" +"(Fitted) (Maura)" )
    #plt.errorbar(dates_1Bartel_Maura, ratioerror_Perugia/ratio_Perugia*100, yerr=0, marker='s', linestyle="None", markersize=25, markerfacecolor="red",ecolor="red", label="4.54_5.0" + " GV" +"(Fitted) (Maura)" )

    ax.tick_params(axis='y', direction='in', length=12, width=10, colors='black')
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.ylabel(r'$\rm{e^-}$/$\rm{e^+}$' + ' relative error',fontsize=55)
    plt.legend(loc='best',fontsize=50)

    plt.savefig( path + "EleOverPos_RelativieError_old_new_compare_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    plt.close()


#### main function()

def drawEleOverPosRatio(lepton_result, lepton_result_Maura, lowworkpath, intermediateworkpath, rigidityrange, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, ElectronPositronIndex, binmerge, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura):

    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"

    #### Load
    ## Load time axis and old bin leption result error 
    # Niko
    LeptonLowEnergyBin           = binning.LeptonLowEnergyBin
    dates_1Bartel_unixsecond_all = binning.Bartals1Unixtime[0:123]
    dates_1Bartel_all            = [datetime.fromtimestamp(i) for i in dates_1Bartel_unixsecond_all]
    electron_to_positron_error   = np.sqrt( ((Electron_error_relative*Electron)/Positron)**2 + (Electron/Positron**2*(Positron_error_relative*Positron))**2)
    # Maura
    MauraTimeRange_all           = binning.Bartals1Unixtime[1:103] #102 bins in total in Maura result, including 2 empty bins.
    MauraTimeRange_WithoutEmpty  = np.delete(np.array(MauraTimeRange_all),[46,47]) #100 bins have results
    dates_1Bartel_Maura          = [datetime.fromtimestamp(i) for i in MauraTimeRange_WithoutEmpty]
    ElectronError_Maura          = np.array(ElectronError_Maura)
    PositronError_Maura          = np.array(PositronError_Maura)
    Positron_Maura               = np.array(Positron_Maura)
    Electron_Maura               = np.array(Electron_Maura)
    EleOverPos_Maura             = np.array(Electron_Maura)/np.array(Positron_Maura)
    EleOverPosRatioError_Maura   = np.sqrt( ((ElectronError_Maura)/Positron_Maura)**2 + (Electron_Maura/Positron_Maura**2*(PositronError_Maura))**2)
    ## Rebin lepton result
    ratio_Aachen, ratioerror_Aachen   = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Aachen")
    ratio_Perugia, ratioerror_Perugia = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Perugia")

    #### Plot
    Plot_old_bin_ElectronPositron              (Electron, Positron, electron_to_positron_error, LeptonLowEnergyBin, ElectronPositronIndex, dates_1Bartel_Maura, EleOverPos_Maura, EleOverPosRatioError_Maura, yaxisfactor, path)  # (the bin choosen depends on ElectronPositronIndex in mannly input)  
    Plot_old_bin_ElectronPositron_RelatvieError(electron_to_positron_error, Electron, Positron, dates_1Bartel_Maura, EleOverPosRatioError_Maura, EleOverPos_Maura, LeptonLowEnergyBin, ElectronPositronIndex, path)
    Plot_new_bin_ElectronPositron              (ratio_Aachen, ratioerror_Aachen, index, binmerge, dates_1Bartel_Maura, ratio_Perugia, ratioerror_Perugia, path, yaxisfactor)
    Plot_new_bin_ElectronPositron_RelatvieError(ratioerror_Aachen, ratio_Aachen, index, binmerge, dates_1Bartel_Maura, ratioerror_Perugia, ratio_Perugia, path)   
    Plot_old_new_bin_ElectronPositron_RelatvieErrorComparision(dates_1Bartel_Maura, EleOverPosRatioError_Maura, EleOverPos_Maura, ElectronPositronIndex, ratioerror_Perugia, ratio_Perugia, index, binmerge, path, LeptonLowEnergyBin)


