import draw
from draw import *

def drawRelativeResult(rigidityrange, time_averaged_ratio, mode, yaxisfactor, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, index, lepton_result, lepton_result_Maura):

    #### Load Electron to Positron result
    LeptonLowEnergyBin                                                = binning.LeptonLowEnergyBin
    EleOverPosRatio_Rebin_Perugia, EleOverPosRatioError_Rebin_Perugia = model.binnin_revised(rigidityrange, lowworkpath, intermediateworkpath, lepton_result, lepton_result_Maura, index, binmerge, "Perugia")
    MauraTimeRange_all           = binning.Bartals1Unixtime[1:103] #102 bins in total in Maura result, including 2 empty bins.
    MauraTimeRange_WithoutEmpty  = np.delete(np.array(MauraTimeRange_all),[46,47]) #100 bins have results
    dates_1Bartel_Maura          = [datetime.fromtimestamp(i) for i in MauraTimeRange_WithoutEmpty]
    EleOverPosRatio_TimeAveraged = np.mean(EleOverPosRatio_Rebin_Perugia)


    #### Plot    
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    # Plot Antiproton to proton relative ratio
    if mode == "6months":
        plt.errorbar(draw.datesM6, draw.M6_relative, yerr=draw.M6relative_error, marker='o',linestyle="None",markersize=12, label="Antiproton_to_Proton("+ str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)])+"GV)")
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, draw.B6_relative, yerr=draw.B6relative_error, marker='o',linestyle="None",markersize=12, label="Antiproton_to_Proton("+ str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)])+"GV)")
    elif mode == "3Bartels":
        plt.errorbar(draw.dates, draw.B3_relative, yerr=draw.B3relative_error, marker='o',linestyle="None",markersize=12, label="Antiproton to Proton Ratio")

    # Plot Electron to Positron relative ratio
    plt.errorbar(dates_1Bartel_Maura, EleOverPosRatio_Rebin_Perugia/EleOverPosRatio_TimeAveraged, yerr=0, marker='o',linestyle="None",markersize=12, label="Electron_to_Positron("+ str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)])+"GV)")

    # format
    plt.legend(loc='best',fontsize=55)
    plt.ylabel("Relative ratio",fontsize=55)
    major_f = md.DateFormatter('%Y-%m')
    major_L = md.MonthLocator(bymonth=[1,7])
    ax.xaxis.set_major_locator(major_L)
    ax.xaxis.set_major_formatter(major_f)
    plt.xticks( rotation=40 )
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    #plt.vlines(datetime(2011, 5, 19, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 19.05.2011 ISS data taking
    #plt.vlines(datetime(2014, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2014 Proton Minimum.
    #plt.vlines(datetime(2017, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2017 Proton Maximum.
 
    # Save Plot
    if rigidityrange == "low":
        if mode == "6months":
            plt.savefig(lowworkpath +"/totalall/results/plot/" +"6MonthsRotation_relative_Low_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" +"6BartalRotation_relative_Low_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" +"3BartalRotation_relative_Low_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    elif rigidityrange == "intermediate":
        if mode == "6months":
            plt.savefig(intermediateworkpath +"/total/results/plot/" +"6MonthsRotation_relative_Intermediate_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" +"6BartalRotation_relative_Intermediate_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" +"3BartalRotation_relative_Intermediate_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()


