import draw
from draw import *


def drawAllbinsRelative(B3_relative):
    plt.figure(figsize=(30,18))
    dates3B = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime[0:B3_relative_all[0,:].shape[0]]]
    plt.xticks( rotation=25 )
    ax=plt.gca()
    if mode == "6months":
        for i in range(M6_relative_all.shape[0]):
            plt.errorbar(datesM6, M6_relative_all[i,:], yerr=M6relative_error_all[i,:], marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV" )
    elif mode == "3Bartels":
        for i in range(B3_relative_all.shape[0]):
            plt.errorbar(dates3B, B3_relative_all[i,:], yerr=B3relative_error_all[i,:], marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV" )
    elif mode == "6Bartels":
        for i in range(B6_relative_all.shape[0]):
            plt.errorbar(dates6B, B6_relative_all[i,:], yerr=B6relative_error_all[i,:], marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV" )
    major_f = md.DateFormatter('%Y-%m')
    major_L = md.MonthLocator(bymonth=[1,7])
    ax.xaxis.set_major_locator(major_L)
    ax.xaxis.set_major_formatter(major_f)
    plt.grid(True,axis='y')
    ax.axes.set_ylim([lowlimit_relative,highlimit_relative])
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=40)
    plt.legend(loc='best',fontsize=35)
    plt.ylabel("Relative Antiproton to Proton Ratio",fontsize=40)
    plt.vlines(datetime(2011, 5, 19, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 19.05.2011 ISS data taking
    plt.vlines(datetime(2014, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2014 Proton Minimum.
    plt.vlines(datetime(2017, 2, 15, 0, 0), lowlimit_relative,highlimit_relative, linestyles = "solid") ## 15.02.2017 Proton Maximum.
    if mode == "6months":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6MonthsRotation_relative_ratio_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    elif mode == "3Bartels":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "3BartalRotation_relative_ratio_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    elif mode == "6Bartels":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6BartalRotation_relative_ratio_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    plt.close()
