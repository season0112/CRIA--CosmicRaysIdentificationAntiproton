import draw
from draw import *
from matplotlib.ticker import ScalarFormatter
import PythonPlotDefaultParameters

def Plot_ElectronPositronRatio(ax1, dates_1Bartel, EleOverPosRatio_Rebin_Perugia, EleOverPosRatioError_Rebin_Perugia):
    ## (1). original result in electron to positron
    #electron_to_positron_error = np.sqrt( ((Electron_error_relative*Electron)/Positron)**2 + (Electron/Positron**2*(Positron_error_relative*Positron))**2)
    #plt.errorbar(np.delete(dates_1Bartel,[46,47]), np.delete(Electron/Positron,[46,47]), yerr=np.delete(electron_to_positron_error,[46,47]), marker='s', linestyle="None", markersize=15, markerfacecolor="r",ecolor="r", label= format(LeptonLowEnergyBin[ElectronPositronIndex-2+3], '.1f') + "_" + format(LeptonLowEnergyBin[ElectronPositronIndex-1+3], '.1f') + ' GeV') # old electron bin
    ## (2). rebined electron/positron ratio
    plt.errorbar(dates_1Bartel, EleOverPosRatio_Rebin_Perugia, yerr=EleOverPosRatioError_Rebin_Perugia, marker='s', linestyle="None", markersize=0, markerfacecolor="white",ecolor="white", label='_nolegend_' )

    ax1.yaxis.set_ticks_position('right')
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    ax1.tick_params(axis='y', which='both', colors='white', direction='in', length=0, width=0)
    plt.yticks(fontsize=0)
    plt.xticks(fontsize=80)

def Plot_AntiprotonProtonRatio(ax1, Antiprotonratio_Scaler, mode):
    ax2 = ax1.twinx()
    if mode == "6months":
        ax2.errorbar(draw.datesM6                       , draw.M6result * Antiprotonratio_Scaler                   , yerr=draw.Error_Total6M * Antiprotonratio_Scaler, marker='o', linestyle="None", markersize=40, linewidth=7, markerfacecolor="b", ecolor="b", label= r'$\overline{\rm {p}}$/p')
    elif mode == "3Bartels":
        ax2.errorbar(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.B3result,[15,38]) * Antiprotonratio_Scaler, yerr=np.delete(draw.Error_Total3B,[15,38]) * Antiprotonratio_Scaler, marker='o',linestyle="None",markersize=40, linewidth=7, markerfacecolor="b", ecolor="b", label= r'$\overline{\rm {p}}$/p' )
    elif mode == "6Bartels":
        ax2.errorbar(draw.dates6B                       , draw.B6result * Antiprotonratio_Scaler                   , yerr=draw.Error_Total6B * Antiprotonratio_Scaler, marker='o', linestyle="None", markersize=40, linewidth=7, markerfacecolor="b", ecolor="b", label= r'$\overline{\rm {p}}$/p')  #label= r'$\overline{\rm {p}}$/p $\times$ 10$^{-5}$', label= r'$\overline{\rm {p}}$/p $\times$' + str('{:.7f}'.format(1/Antiprotonratio_Scaler))

    #ax2.axes.set_ylim([ ax2.get_ylim()[0]-((ax2.get_ylim()[1]-ax2.get_ylim()[0])/2)*yaxisfactor,  ax2.get_ylim()[1]+((ax2.get_ylim()[1]-ax2.get_ylim()[0])/2)*yaxisfactor ])
    #ax2.set_ylabel('',fontsize=70, color="b")
    #ax2.tick_params(axis='y', direction='in', length=30, width=10, colors='black', pad=40.0)

    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position("left")
    ax1.yaxis.set_ticks_position('right') #Here needs to rewrite it again, to display properly.

    plt.xticks(fontsize=80)
    plt.yticks(fontsize=80)
    #major_f = md.DateFormatter('%Y-%m') #major_f = md.DateFormatter('%Y-%m-%d %H:%M:%S')
    #major_L = md.MonthLocator(bymonth=[1,7])
    #ax2.xaxis.set_major_locator(major_L) # for some reason, it shows ax2 xaxis.
    #ax2.xaxis.set_major_formatter(major_f)

    ax2.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
    ax2.get_yaxis().get_major_formatter().set_scientific(True)
    ax2.yaxis.offsetText.set_fontsize(80)

    #ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

    return ax2


def RePlot_ElectronPositronInSingleYAxis(mode, ax1, ax2, Antiprotonratio_Scaler, dates_1Bartel, EleOverPosRatio_Rebin_Perugia, EleOverPosRatioError_Rebin_Perugia):

    #### Replot electron plot in single y axis
    '''
    #(1). mean of e/e+ is same as pbar/p
    print("Electron_Scaler:")
    if mode == "6months":
        Electron_Scaler = np.mean(draw.M6result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia)
    elif mode == "3Bartels":
        Electron_Scaler = np.mean(draw.B3result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia)
    elif mode == "6Bartels":
        Electron_Scaler = np.mean(draw.B6result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia)
    print(Electron_Scaler)
    print('\n')
    '''
    #(2). mean of e/e+ is same as mean of pbar/p for: first 4 points for pbar/p, first 23 points for e/e+. (2011 to 2013)
    if mode == "6months":
        Electron_Scaler     = np.mean(draw.M6result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia[0:22])
        #Electron_Scaler_ori = np.mean(draw.M6result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPos_Maura[0:22])
    elif mode == "3Bartels":
        Electron_Scaler     = np.mean(draw.B3result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia[0:22])
        #Electron_Scaler_ori = np.mean(draw.B3result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPos_Maura[0:22])
    elif mode == "6Bartels":
        Electron_Scaler     = np.mean(draw.B6result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia[0:22])
        #Electron_Scaler_ori = np.mean(draw.B6result[0:3] * Antiprotonratio_Scaler) / np.mean(EleOverPos_Maura[0:22])

    ax3 = ax1.twinx()
    ax3.errorbar(dates_1Bartel, EleOverPosRatio_Rebin_Perugia * Electron_Scaler, yerr = EleOverPosRatioError_Rebin_Perugia * Electron_Scaler, marker='o', linestyle="None", markersize=13, markerfacecolor="r", ecolor="r", markeredgecolor='r', label=r'$\rm{e^-}$/$\rm{e^+}$' + r'$\times$' + "{:.2}".format(Electron_Scaler*1e6) + r'$\times$' + r'10$^{-6}$') # label=r'$\rm{e^-}$/$\rm{e^+}$' +  r'$\times$' + "{:.2}".format(Electron_Scaler*1e6) + r'$\times$' + r'10$^{-6}$';  label=r'$\rm{e^-}$/$\rm{e^+}$' + r'$\times$' + "{:.2}".format(Electron_Scaler)
    #ax3.errorbar(dates_1Bartel, EleOverPos_Maura * Electron_Scaler_ori        , yerr = EleOverPosRatioError_Maura * Electron_Scaler_ori    , marker='o', linestyle="None", markersize=13, markerfacecolor="r",ecolor="r", label=r'$\rm{e^-}$/$\rm{e^+}$' + r'$\times$' + "{:.2}".format(Electron_Scaler_ori))

    if ax3.get_ylim()[0] < ax2.get_ylim()[0]:
        ax2.axes.set_ylim([ax3.get_ylim()[0], ax2.get_ylim()[1]])
    if ax3.get_ylim()[1] > ax2.get_ylim()[1]:
        ax2.axes.set_ylim([ax2.get_ylim()[0], ax3.get_ylim()[1]])

    ax3.tick_params(axis='y', which='both', direction='in', length=0, width=0, colors='white')
    #plt.yticks(fontsize=80)

    ax3.axes.set_ylim([ax2.get_ylim()[0], ax2.get_ylim()[1]])
    print("electron low (replot):" + str( ax3.get_ylim()[0] ))
    print("electron high (replot):" + str( ax3.get_ylim()[1] ))

    return ax3


#### main function()
def drawCompareWithElectron(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, trackerpattern, richcut, plotrange, binmerge, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, mode, ElectronPositronIndex, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura, group):

    #### Load Maura original result
    ElectronError_Maura        = np.array(ElectronError_Maura) # no use
    PositronError_Maura        = np.array(PositronError_Maura) # no use
    Positron_Maura             = np.array(Positron_Maura)      # used
    Electron_Maura             = np.array(Electron_Maura)      # used
    EleOverPos_Maura           = np.array(Electron_Maura)/np.array(Positron_Maura) # no use
    EleOverPosRatioError_Maura = np.sqrt( ((ElectronError_Maura)/Positron_Maura)**2 + (Electron_Maura/Positron_Maura**2*(PositronError_Maura))**2) # no use

    #### Load Maura time bins
    LeptonLowEnergyBin          = binning.LeptonLowEnergyBin
    MauraTimeRange_all          = binning.Bartals1Unixtime[1:103]                 # 102 bins in total in Maura result, including 2 empty bins.
    MauraTimeRange_WithoutEmpty = np.delete(np.array(MauraTimeRange_all),[46,47]) # 100 bins have results
    dates_1Bartel_Maura         = [datetime.fromtimestamp(i) for i in MauraTimeRange_WithoutEmpty]

    #### Load rebinned result
    if group == "Perugia":
        EleOverPosRatio_Rebin_Perugia, EleOverPosRatioError_Rebin_Perugia = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Perugia")
    elif group == "Aachen":
        EleOverPosRatio_Rebin_Aachen , EleOverPosRatioError_Rebin_Aachen  = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Aachen")

    #### Remove Maura TTCS OFF period
    if group == "Perugia":
        del dates_1Bartel_Maura[-3]
        EleOverPosRatio_Rebin_Perugia      = np.delete(EleOverPosRatio_Rebin_Perugia      , [-3])  ## -3 and -5 are partly ttcsoff period, -3 removed, -5 roughly ok so keeped.
        EleOverPosRatioError_Rebin_Perugia = np.delete( EleOverPosRatioError_Rebin_Perugia, [-3])
        ElectronError_Maura                = np.delete(ElectronError_Maura                , [-3])
        PositronError_Maura                = np.delete(PositronError_Maura                , [-3])
        Positron_Maura                     = np.delete(Positron_Maura                     , [-3])
        Electron_Maura                     = np.delete(Electron_Maura                     , [-3])
        EleOverPos_Maura                   = np.delete(EleOverPos_Maura                   , [-3])
        EleOverPosRatioError_Maura         = np.delete(EleOverPosRatioError_Maura         , [-3])


    #### Scaler
    yaxisfactor = 1.0
    #Antiprotonratio_Scaler = 100000
    Antiprotonratio_Scaler = 1

    #### Choice of electron over positron result
    if group == "Perugia":
        dates_1Bartel              = dates_1Bartel_Maura 
        EleOverPosRatio_Rebin      = EleOverPosRatio_Rebin_Perugia
        EleOverPosRatioError_Rebin = EleOverPosRatioError_Rebin_Perugia
    elif group == "Aachen":
        dates_1Bartel              = np.delete(draw.dates_1Bartel,[46,47])
        EleOverPosRatio_Rebin      = EleOverPosRatio_Rebin_Aachen   
        EleOverPosRatioError_Rebin = EleOverPosRatioError_Rebin_Aachen


    #### Plot
    fig = plt.figure(figsize=(40, 20))
    ax1 = fig.add_subplot(111)

    Plot_ElectronPositronRatio              (ax1, dates_1Bartel, EleOverPosRatio_Rebin, EleOverPosRatioError_Rebin) # x tick drar here
    ax2=Plot_AntiprotonProtonRatio          (ax1, Antiprotonratio_Scaler, mode)  # y tick draw here
    ax3=RePlot_ElectronPositronInSingleYAxis(mode, ax1, ax2, Antiprotonratio_Scaler, dates_1Bartel, EleOverPosRatio_Rebin, EleOverPosRatioError_Rebin)  

    ## plot legend
    leg = fig.legend(loc="upper right", bbox_to_anchor=(0.87,1), bbox_transform=ax3.transAxes, fontsize=60, frameon=False)
    leg.texts[0].set_color('blue')
    leg.texts[1].set_color('red')

    ## plot text for rigidity information
    RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    plt.text(0.05, 0.15, RigidityBinText, transform = ax1.transAxes, fontsize=60)

    ## Save plots
    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"
    if mode == "6months":
        plt.savefig( path + "Antiproton_Electron_" + group + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.pdf")
    elif mode == "3Bartels":
        plt.savefig( path + "Antiproton_Electron_" + group + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3Bartels.pdf")
    elif mode == "6Bartels":
        plt.savefig( path + "Antiproton_Electron_" + group + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6Bartels.pdf")

    plt.close()







