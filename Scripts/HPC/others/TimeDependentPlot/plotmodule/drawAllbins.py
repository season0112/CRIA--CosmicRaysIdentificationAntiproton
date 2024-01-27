import draw
from draw import *


def drawAllbins(mode):
    #dates3B = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime[0:draw.B3_relative_all[0,:].shape[0]]]
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    if mode == "6months":
        for i in range(M6result_all.shape[0]):
            plt.errorbar(datesM6, M6result_all[i,:]*10000, yerr=Error6m_all[i,:]*10000, marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV")
    elif mode == "3Bartels":
        for i in range(draw.B3result_all.shape[0]):
            plt.errorbar(dates3B, draw.B3result_all[i,:]*10000, yerr=Error3b_all[i,:]*10000, marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV")
    elif mode == "6Bartels":
        for i in range(B6result_all.shape[0]):
            plt.errorbar(dates6B, B6result_all[i,:]*10000, yerr=Error6b_all[i,:]*10000, marker='o', linestyle="None", markersize=12, label="{:g}".format(binning.published2016binnings[plotrange[i]]) + "_" + "{:g}".format(binning.published2016binnings[plotrange[i]+int(binmerge)]) + "GV")
    major_f = md.DateFormatter('%Y-%m')
    major_L = md.MonthLocator(bymonth=[1,7])
    ax.xaxis.set_major_locator(major_L)
    ax.xaxis.set_major_formatter(major_f)
    plt.grid(True,axis='y')
    plt.tick_params(direction='in',length=10, width=3)
    ax.axes.set_ylim([lowlimit,highlimit])
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=40)
    plt.legend(loc='best',fontsize=35)
    plt.ylabel("Antiproton to Proton ratio (10$^{-4}$)",fontsize=40)
    plt.vlines(datetime(2011, 5, 19, 0, 0), lowlimit, highlimit, linestyles = "solid") ## 19.05.2011 ISS data taking
    plt.vlines(datetime(2014, 2, 15, 0, 0), lowlimit, highlimit, linestyles = "solid") ## 15.02.2014 Proton Minimum.
    plt.vlines(datetime(2017, 2, 15, 0, 0), lowlimit, highlimit, linestyles = "solid") ## 15.02.2017 Proton Maximum.
    if mode == "6months":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6MonthsRotation_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    elif mode == "3Bartels":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "3BartalRotation_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    elif mode == "6Bartels":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6BartalRotation_all" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + ".png")
    plt.close()
