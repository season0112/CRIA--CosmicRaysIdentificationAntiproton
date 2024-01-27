import draw
from draw import *

def drawCompareWithElectron(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, trackerpattern, richcut, plotrange, binmerge, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, mode, ElectronPositronIndex):

    ##Load rebinned result
    LeptonLowEnergyBin = binning.LeptonLowEnergyBin
    EleOverPosRatio_Rebin_Perugia, EleOverPosRatioError_Rebin_Perugia = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Perugia")
    MauraTimeRange_all = binning.Bartals1Unixtime[1:103] #102 bins in total in Maura result, including 2 empty bins.
    MauraTimeRange_WithoutEmpty = np.delete(np.array(MauraTimeRange_all),[46,47]) #100 bins have results
    dates_1Bartel_Maura = [datetime.fromtimestamp(i) for i in MauraTimeRange_WithoutEmpty]

    ## TTCS OFF period
    del dates_1Bartel_Maura[-3]
    EleOverPosRatio_Rebin_Perugia = np.delete(EleOverPosRatio_Rebin_Perugia, [-3])  ## -3 and -5 are partly ttcsoff period, -3 removed, -5 roughly ok so keeped.
    EleOverPosRatioError_Rebin_Perugia = np.delete( EleOverPosRatioError_Rebin_Perugia, [-3])

    ## Scaler
    yaxisfactor=1.0
    Antiprotonratio_Scaler = 100000

    #### Draw electron/positron ratio
    fig = plt.figure(figsize=(38,18))
    ax1 = fig.add_subplot(111)
    ## (1). original result in electron to positron 
    #electron_to_positron_error = np.sqrt( ((Electron_error_relative*Electron)/Positron)**2 + (Electron/Positron**2*(Positron_error_relative*Positron))**2)
    #plt.errorbar(np.delete(dates_1Bartel,[46,47]), np.delete(Electron/Positron,[46,47]), yerr=np.delete(electron_to_positron_error,[46,47]), marker='s', linestyle="None", markersize=15, markerfacecolor="r",ecolor="r", label= format(LeptonLowEnergyBin[ElectronPositronIndex-2+3], '.1f') + "_" + format(LeptonLowEnergyBin[ElectronPositronIndex-1+3], '.1f') + ' GeV') # old electron bin
    ## (2). rebined electron/positron ratio
    plt.errorbar(dates_1Bartel_Maura, EleOverPosRatio_Rebin_Perugia, yerr=EleOverPosRatioError_Rebin_Perugia, marker='s', linestyle="None", markersize=0, markerfacecolor="white",ecolor="white", label='_nolegend_' )  #label=r'$\rm{e^-}$/$\rm{e^+}$' +  r' $\times$ ' + "{:.5}".format(Electron_Scaler), label='_nolegend_'

    #ax1.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    ax1.yaxis.set_ticks_position('right')
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.tick_right()
    #plt.tick_params(axis='y',top=False,bottom=True,left=False,right=False)
    ax1.tick_params(axis='y', labeltop=False,labelbottom=True, labelleft=True, labelright=True, colors='white', direction='in', length=0, width=0)
    ax1.tick_params(axis='x', direction='in', length=30, width=10)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)


    #### Draw antiproton/proton ratio
    ax2 = ax1.twinx()
    if mode == "6months":
        ax2.errorbar(draw.datesM6, draw.M6result * Antiprotonratio_Scaler, yerr=draw.ErrorM6 * Antiprotonratio_Scaler, marker='o', linestyle="None", markersize=25, markerfacecolor="b", ecolor="b", linewidth=5, label= r'$\overline{\rm {p}}$/p $\times$ 10$^{-5}$')  
    elif mode == "3Bartels":
        ax2.errorbar(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.B3result,[15,38]) * Antiprotonratio_Scaler, yerr=np.delete(Error3b,[15,38]) * Antiprotonratio_Scaler, marker='o',linestyle="None",markersize=25, markerfacecolor="b", ecolor="b", label= r'$\overline{\rm {p}}$/p $\times$ 10$^{-5}$' )
    elif mode == "6Bartels":
        ax2.errorbar(draw.dates6B, draw.B6result * Antiprotonratio_Scaler, yerr=draw.Error6b * Antiprotonratio_Scaler, marker='o', linestyle="None", markersize=25, markerfacecolor="b", ecolor="b", linewidth=5, label= r'$\overline{\rm {p}}$/p $\times$ 10$^{-5}$')  #label= r'$\overline{\rm {p}}$/p $\times$ 10$^{-5}$', label= r'$\overline{\rm {p}}$/p $\times$' + str('{:.7f}'.format(1/Antiprotonratio_Scaler))

    '''
    print("ax2.get_ylim()[0]:")
    print(plt.ylim()[0])
    print(plt.ylim()[1])
    print(ax2.get_ylim()[0])
    print(ax2.get_ylim()[1])
    if ax2.get_ylim()[0]<0:
        ax2.axes.set_ylim(0,ax2.get_ylim()[1])
    '''

    ax2.set_ylim([0, plt.ylim()[1]])

    #ax2.axes.set_ylim([ ax2.get_ylim()[0]-((ax2.get_ylim()[1]-ax2.get_ylim()[0])/2)*yaxisfactor,  ax2.get_ylim()[1]+((ax2.get_ylim()[1]-ax2.get_ylim()[0])/2)*yaxisfactor ])
    #ax2.set_ylabel('',fontsize=70, color="b")
    ax2.tick_params(axis='y', direction='in', length=12, width=10, colors='blue')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position("left")
    ax1.yaxis.set_ticks_position('right') #Here needs to rewrite it again, to display properly.
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    #major_f = md.DateFormatter('%Y-%m') #major_f = md.DateFormatter('%Y-%m-%d %H:%M:%S')
    #major_L = md.MonthLocator(bymonth=[1,7])
    #ax2.xaxis.set_major_locator(major_L) # for some reason, it shows ax2 xaxis.
    #ax2.xaxis.set_major_formatter(major_f)


     
    #### Replot electron plot in single y axis
    ymin1 = ax1.get_ylim()[0]
    ymax1 = ax1.get_ylim()[1]
    ymin2 = ax2.get_ylim()[0]
    ymax2 = ax2.get_ylim()[1]

    print('\n')
    print("electron low:" + str(ymin1))
    print("electron high:" + str(ymax1))
    print("antiproton low:" + str(ymin2))
    print("antiproton high:" + str(ymax2))

    Rebinfactor = (ymax2-ymin2)/(ymax1-ymin1)
    print("Rebinfactor:")
    print(Rebinfactor)

    #### Replot electron plot in single y axis
    print("Electron_Scaler:")
    if mode == "6months":
        Electron_Scaler = np.mean(draw.M6result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia * Rebinfactor)
    elif mode == "3Bartels":
        Electron_Scaler = np.mean(draw.B3result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia * Rebinfactor)
    elif mode == "6Bartels":
        Electron_Scaler = np.mean(draw.B6result * Antiprotonratio_Scaler) / np.mean(EleOverPosRatio_Rebin_Perugia * Rebinfactor)
    print(Electron_Scaler)
    print('\n')


    print("Rebinfactor*Electron_Scaler")
    print(Rebinfactor*Electron_Scaler)

    ShiftValue = -ymin1*Rebinfactor+ymin2
    print("shift value:" )
    print(ShiftValue)
    print('\n')
    
    ax3 = ax1.twinx()
    ax3.errorbar(dates_1Bartel_Maura, EleOverPosRatio_Rebin_Perugia * Rebinfactor, yerr=EleOverPosRatioError_Rebin_Perugia * Rebinfactor, marker='s', linestyle="None", markersize=15, markerfacecolor="r",ecolor="r", label='_nolegend_')

    ax3.errorbar(dates_1Bartel_Maura, EleOverPosRatio_Rebin_Perugia * Rebinfactor - ax3.get_ylim()[0] + ymin2, yerr=EleOverPosRatioError_Rebin_Perugia * Rebinfactor, marker='s', linestyle="None", markersize=15, markerfacecolor="r",ecolor="r", label=r'$\rm{e^-}$/$\rm{e^+}$' +  r' $\times$ ' + "{:.2}".format(Rebinfactor) + "-" + str("{:.2f}".format(ax3.get_ylim()[0])))

    #ax3.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    ax3.tick_params(axis='y', direction='in', length=12, width=10, colors='red')
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)

    ax3.axes.set_ylim([ymin2,ymax2])
    ymin3 = ax3.get_ylim()[0]
    ymax3 = ax3.get_ylim()[1]
    print("electron low (replot):" + str(ymin3))
    print("electron high (replot):" + str(ymax3))
    

    '''
    #### Replot electron plot in single y axis and has comparable time variation effect
    '''


    #### plot legend
    leg = fig.legend(loc="upper right", bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes, fontsize=70)
    #leg.texts[0].set_color('red')
    #leg.texts[1].set_color('blue')
    leg.texts[0].set_color('blue')
    leg.texts[1].set_color('red')

    #### plot text for rigidity information
    plt.text(0.05, 0.9, format(binning.published2016binnings[index], '.1f') + "-" + format(binning.published2016binnings[index+int(binmerge)], '.1f') + ' GV', fontsize=80, transform=ax1.transAxes)

    #### Save plots
    if rigidityrange == "low":
        path = lowworkpath +"/totalall/results/plot/"
    elif rigidityrange == "intermediate":
        path = intermediateworkpath +"/total/results/plot/"
    
    if mode == "6months":
        plt.savefig( path + "Antiproton_Electron_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.pdf")
    elif mode == "3Bartels":
        plt.savefig( path + "Antiproton_Electron_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3Bartels.pdf")
    elif mode == "6Bartels":
        plt.savefig( path + "Antiproton_Electron_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6Bartels.pdf")
 
    plt.close()
    



