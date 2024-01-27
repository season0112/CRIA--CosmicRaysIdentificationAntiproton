import draw
from draw import *


def drawStatisticalErrorComparision(mode, rigidityrange, Taiwan_error_all, Ihep_error_all, Taiwan_error_all_6, Taiwan_all_6, Taiwan_all, Ihep_all, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge):

    ## Plot1: plot absolute Statistical error
    plt.figure(figsize=(30,18))
    #dates=[datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime[0:B3result.shape[0]]] # first bartel rotation included.
    dates_unixsecond = binning.Bartals3Unixtime_WithoutFirstBartel[0:draw.B3result.shape[0]] # first bartel rotation removed.
    dates=[datetime.fromtimestamp(i) for i in dates_unixsecond]
    ax=plt.gca()

    if mode == "3Bartels":
        plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])*10000, marker='o',linestyle="dashed",markersize=14,label="StaErr")
        #plt.plot(draw.dates_TaiwanRange[0:15]+draw.dates_TaiwanRange[16:], np.delete(Taiwan_error_all[index-5-int((index-9)/2)-1], 15)*10000, marker='s',linestyle="dashed",markersize=14, label="Taiwan", color="red")
        #plt.plot(draw.dates_IhepRange[0:15]+draw.dates_IhepRange[16:], np.delete(Ihep_error_all[index-6-int((index-9)/2)-1], 15)*10000, marker='*',linestyle="dashed",markersize=14, label="IHEP", color="magenta")
    elif mode == "6Bartels":
        plt.plot(draw.dates6B, draw.Error6b*10000, marker='o',linestyle="None",markersize=15, label="Aachen")
        plt.plot(draw.dates6B_TaiwanRange, np.array( Taiwan_error_all_6[index-5-int((index-9)/2)-1-1] *10000 ), marker='s', linestyle="None", markersize=12, label="Taiwan")

    #major_f = md.DateFormatter('%Y-%m')
    #major_L = md.MonthLocator(bymonth=[1,7])
    #ax.xaxis.set_major_locator(major_L)
    #ax.xaxis.set_major_formatter(major_f)
    plt.legend(loc='best',fontsize=35)
    plt.grid(True,axis='y')
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    #plt.ylabel("Statistical Error",fontsize=30)
    plt.ylabel("Error",fontsize=30)

    if rigidityrange == "low":
        if mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "Statistical_Error_Compare_3Bartels_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "Statistical_Error_Compare_6Bartels_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    elif rigidityrange == "intermediate":
        if mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Statistical_Error_Compare_3Bartels_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Statistical_Error_Compare_6Bartels_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()


    ## Plot2: plot relative Statistical error

    plt.figure(figsize=(30,18))
    ax=plt.gca()

    if mode == "3Bartels":
        plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])/np.delete(draw.B3result,[15,38]), marker='o',linestyle="dashed",markersize=14,label="Aachen")
        plt.plot(draw.dates_TaiwanRange[0:15]+draw.dates_TaiwanRange[16:], np.delete(Taiwan_error_all[index-5-int((index-9)/2)-1], 15) / np.delete(np.array(Taiwan_all[index-5-int((index-9)/2)-1]), 15), marker='s',linestyle="dashed",markersize=14, label="Taiwan", color="red")
        plt.plot(draw.dates_IhepRange[0:15]+draw.dates_IhepRange[16:], np.delete(Ihep_error_all[index-6-int((index-9)/2)-1], 15) / np.delete( np.array(Ihep_all[index-6-int((index-9)/2)-1]), 15), marker='*',linestyle="dashed",markersize=14, label="IHEP", color="magenta")
    elif mode == "6Bartels":
        plt.plot(draw.dates6B, draw.Error6b/draw.B6result, marker='o',linestyle="None",markersize=12, label="Aachen")
        plt.plot(draw.dates6B_TaiwanRange, np.array(Taiwan_error_all_6[index-5-int((index-9)/2)-1-1]) / np.array(Taiwan_all_6[index-5-int((index-9)/2)-1-1]), marker='s', linestyle="None", markersize=12, label="Taiwan")

    #major_f = md.DateFormatter('%Y-%m')
    #major_L = md.MonthLocator(bymonth=[1,7])
    #ax.xaxis.set_major_locator(major_L)
    #ax.xaxis.set_major_formatter(major_f)
    plt.legend(loc='best',fontsize=35)
    plt.grid(True,axis='y')
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.ylabel("Relative Statistical Error",fontsize=30)

    if rigidityrange == "low":
        if mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "RelativeStatistical_Error_Compare_3Bartels_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "RelativeStatistical_Error_Compare_6Bartels_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    elif rigidityrange == "intermediate":
        if mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "RelativeStatistical_Error_Compare_3Bartels_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "RelativeStatistical_Error_Compare_6Bartels_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()




