import draw
from draw import *


def drawProtontoAntiprotonRelativeResult(mode, intermediateworkpath, trackerpattern, richcut, binmerge, index):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    if mode == "6months":
        plt.errorbar(draw.datesM6, 1./M6_relative, yerr=M6relative_error, marker='o',linestyle="None",markersize=12, label="Proton_to_Antiproton("+ str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)])+"GV)") # error to fix.
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, 1./B6_relative, yerr=B6relative_error, marker='o',linestyle="None",markersize=12, label="Proton_to_Antiproton("+ str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+int(binmerge)])+"GV)") # error to fix.
    elif mode == "3Bartels":
        plt.errorbar(draw.dates, 1./draw.B3_relative, yerr=draw.B3_relative_proton_to_antiproton, marker='o',linestyle="None",markersize=12, markerfacecolor="b", ecolor="b", label="Proton to Antiproton Ratio")
    plt.errorbar(np.delete(draw.dates_1Bartel,[46,47]), np.delete(PosiOverElec/E_ave,[46,47]), yerr= np.delete(PosiOverElec_relative_error,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="r", ecolor="r", label="Positron to Electron Ratio") # option: label="Positron_to_Electron("+ str(LeptonLowEnergyBin[ElectronPositronIndex-2+3]) + "_" + str(LeptonLowEnergyBin[ElectronPositronIndex-1+3]) +"GeV)"
    plt.legend(loc='best',fontsize=55)
    plt.ylabel("Relative ratio",fontsize=55)
    major_f = md.DateFormatter('%Y-%m')
    major_L = md.MonthLocator(bymonth=[1,7])
    ax.xaxis.set_major_locator(major_L)
    ax.xaxis.set_major_formatter(major_f)
    ax.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    #ax.axes.set_ylim([lowlimit_relative,highlimit_relative])
    plt.xticks( rotation=40 )
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    #plt.vlines(datetime(2011, 5, 19, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 19.05.2011 ISS data taking
    #plt.vlines(datetime(2014, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2014 Proton Minimum.
    #plt.vlines(datetime(2017, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2017 Proton Maximum.
    plt.savefig(intermediateworkpath +"/total/results/plot/" +"3BartalRotation_relative_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".png")
    plt.close()

