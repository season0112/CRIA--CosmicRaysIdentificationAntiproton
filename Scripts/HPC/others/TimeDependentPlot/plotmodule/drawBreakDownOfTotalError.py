import draw
from draw import *
import PythonPlotDefaultParameters


def drawStatisticalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, pattern):

    plt.figure()

    dates_unixsecond = binning.Bartals3Unixtime_WithoutFirstBartel[0:draw.B3result.shape[0]] # first bartel rotation removed.
    dates=[datetime.fromtimestamp(i) for i in dates_unixsecond]
    ax=plt.gca()

    MarkersizeValue = 30

    if mode == "3Bartels":
        if pattern == 'absolute':
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])       , marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError")
        if pattern == 'relative':
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])/np.delete(draw.B3result,[15,38])*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError")

    elif mode == "6Bartels":
        if pattern == 'absolute':
            plt.plot(draw.dates6B, draw.Error6b         , color='red'   , marker='o', linestyle="None", markersize=MarkersizeValue, label="Statistical uncertainty")
        if pattern == 'relative':
            plt.plot(draw.dates6B, draw.Error6b/draw.B6result*100         , color='red'   , marker='o', linestyle="None", markersize=MarkersizeValue, label="Relative Statistical uncertainty")

    yaxisfactor = 1.4
    ax.set_ylim(0, plt.ylim()[1]*yaxisfactor)

    if pattern == 'absolute':
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
        ax.yaxis.get_offset_text().set_fontsize(60)

    if pattern == 'absolute':
        plt.ylabel(r'Absolute statistical uncertainty'    , horizontalalignment='right', y=1.0)
    elif pattern == 'relative':
        plt.ylabel(r'Relative statistical uncertainty (%)', horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        if pattern == 'absolute': 
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "StatisticalAbsoluteError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "StatisticalRelativeError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    elif rigidityrange == "intermediate":
        if pattern == 'absolute':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "StatisticalAbsoluteError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "StatisticalRelativeError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()




def drawBreakDownOfTotalSysError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, pattern):

    plt.figure()

    dates_unixsecond = binning.Bartals3Unixtime_WithoutFirstBartel[0:draw.B3result.shape[0]] # first bartel rotation removed.
    dates=[datetime.fromtimestamp(i) for i in dates_unixsecond]

    ax=plt.gca()

    MarkersizeValue = 30

    if mode == "6Bartels":
        if pattern == 'absolute':
            plt.plot(draw.dates6B, draw.Error_Sys6B     , color='magenta', marker='o', linestyle="None", markersize=MarkersizeValue, label="Systematic uncertainty (FitRange)")
            plt.plot(draw.dates6B, draw.Error_SysAcc6B  , color='green'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="Systematic uncertainty (Acceptance)")
            plt.plot(draw.dates6B, draw.Error_TotalSys6B, color='blue'   , marker='o', linestyle="None", markersize=MarkersizeValue, label="Total systematic uncertainty")
        if pattern == 'relative':
            plt.plot(draw.dates6B, draw.Error_Sys6B/draw.B6result     *100, color='magenta', marker='o', linestyle="None", markersize=MarkersizeValue, label="Systematic uncertainty (FitRange)")
            plt.plot(draw.dates6B, draw.Error_SysAcc6B/draw.B6result  *100, color='green'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="Systematic uncertainty (Acceptance)")
            plt.plot(draw.dates6B, draw.Error_TotalSys6B/draw.B6result*100, color='blue'   , marker='o', linestyle="None", markersize=MarkersizeValue, label="Total systematic uncertainty")

    yaxisfactor = 1.4
    ax.set_ylim(0, plt.ylim()[1]*yaxisfactor)

    #RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    #plt.text(0.10, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    plt.legend()

    if pattern == 'absolute':
        ax.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
        ax.yaxis.offsetText.set_fontsize(60)
        plt.ylabel(r'Absolute uncertainty'    , horizontalalignment='right', y=1.0)
    elif pattern == 'relative':
        plt.ylabel(r'Relative uncertainty (%)', horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        if pattern == 'absolute':
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "BreakDownOfTotalSysError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "RelativeBreakDownOfTotalSysError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    elif rigidityrange == "intermediate":
        if pattern == 'absolute':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "BreakDownOfTotalSysError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "RelativeBreakDownOfTotalSysError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()



def drawBreakDownOfTotalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, pattern):

    plt.figure()

    dates_unixsecond = binning.Bartals3Unixtime_WithoutFirstBartel[0:draw.B3result.shape[0]] # first bartel rotation removed.
    dates=[datetime.fromtimestamp(i) for i in dates_unixsecond]
    ax=plt.gca()

    MarkersizeValue = 30

    if mode == "3Bartels":
        if pattern == 'absolute':
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])*10000       , marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_Sys3B,[15,38])*10000   , marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(FitRange)")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_SysAcc3B,[15,38])*10000, marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(Acceptance)")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_Total3B,[15,38])*10000 , marker='o', linestyle="None", markersize=MarkersizeValue, label="TotalError")
        if pattern == 'relative':
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error3b,[15,38])    /np.delete(draw.B3result,[15,38])*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_Sys3B,[15,38])/np.delete(draw.B3result,[15,38])*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(FitRange)")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_SysAcc3B,[15,38])/np.delete(draw.B3result,[15,38])*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(Acceptance)")
            plt.plot(dates[0:15]+dates[16:38]+dates[39:], np.delete(draw.Error_Total3B,[15,38])/np.delete(draw.B3result,[15,38])*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="TotalError")

    elif mode == "6Bartels":
        if pattern == 'absolute':
            plt.plot(draw.dates6B, draw.Error6b*10000         , color='red'   , marker='o', linestyle="None", markersize=MarkersizeValue, label="Statistical uncertainty") 
            #plt.plot(draw.dates6B, draw.Error_Sys6B*10000    , color='orange' , marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(FitRange)")
            #plt.plot(draw.dates6B, draw.Error_SysAcc6B*10000 , color='green'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(Acceptance)")
            plt.plot(draw.dates6B, draw.Error_TotalSys6B*10000, color='blue', marker='o', linestyle="None", markersize=MarkersizeValue, label="Systematic uncertainty")
            plt.plot(draw.dates6B, draw.Error_Total6B*10000   , color='black'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="Total uncertainty")
        if pattern == 'relative':
            plt.plot(draw.dates6B, draw.Error6b/draw.B6result*100       , color='red'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="Statistical uncertainty")
            plt.plot(draw.dates6B, draw.Error_Sys6B/draw.B6result*100   , color='blue' , marker='o', linestyle="None", markersize=MarkersizeValue, label="Time-dependent systematic uncertainty") # FitRange
            plt.plot(draw.dates6B, draw.Error_SysAcc6B/draw.B6result*100, color='green', marker='o', linestyle="None", markersize=MarkersizeValue, label="Time-independent systematic uncertainty") # Acceptance
            #plt.plot(draw.dates6B, draw.Error_TotalSys6B/draw.B6result*100, color='orange'  , marker='o', linestyle="None", markersize=MarkersizeValue, label="Total systematic uncertainty")
            plt.plot(draw.dates6B, draw.Error_Total6B/draw.B6result*100   , color='black' , marker='o', linestyle="None", markersize=MarkersizeValue, label="Total uncertainty")

    elif mode == "6months":
        if pattern == 'absolute':
            plt.plot(draw.datesM6, draw.ErrorM6*10000       , marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError") 
            plt.plot(draw.datesM6, draw.Error_Sys6M*10000   , marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(FitRange)")
            plt.plot(draw.datesM6, draw.Error_SysAcc6M*10000, marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(Acceptance)")
            plt.plot(draw.datesM6, draw.Error_Total6M*10000 , marker='o', linestyle="None", markersize=MarkersizeValue, label="TotalError")
        if pattern == 'relative':
            plt.plot(draw.datesM6, draw.ErrorM6/draw.M6result*100       , marker='o', linestyle="None", markersize=MarkersizeValue, label="StaError")
            plt.plot(draw.datesM6, draw.Error_Sys6M/draw.M6result*100   , marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(FitRange)")
            plt.plot(draw.datesM6, draw.Error_SysAcc6M/draw.M6result*100, marker='o', linestyle="None", markersize=MarkersizeValue, label="SysError(Acceptance)")
            plt.plot(draw.datesM6, draw.Error_Total6M/draw.M6result*100 , marker='o', linestyle="None", markersize=MarkersizeValue, label="TotalError")

    yaxisfactor = 1.4
    ax.set_ylim(0, plt.ylim()[1]*yaxisfactor)

    #RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    #plt.text(0.10, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    plt.legend()

    if pattern == 'absolute':
        plt.ylabel(r'Uncertainty ($10^{-4}$)', fontsize=60, horizontalalignment='right', y=1.0)
    elif pattern == 'relative':
        plt.ylabel(r'Relative uncertainty (%)'   , fontsize=60, horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        if pattern == 'absolute':
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "BreakDownOfTotalError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "RelativeBreakDownOfTotalError_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    elif rigidityrange == "intermediate":
        if pattern == 'absolute':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "BreakDownOfTotalError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif pattern == 'relative':
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "RelativeBreakDownOfTotalError_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()


def drawErrorPersentage(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start):

    plt.figure(figsize=(30,18))
    ax=plt.gca()

    dates_unixsecond = binning.Bartals3Unixtime_WithoutFirstBartel[0:draw.B3result.shape[0]] # first bartel rotation removed.
    dates=[datetime.fromtimestamp(i) for i in dates_unixsecond]

    MarkersizeValue = 35

    if mode == "3Bartels":
        plt.plot( dates[0:15]+dates[16:38]+dates[39:], (np.delete(draw.Error3b,[15,38]))**2/(np.delete(draw.Error_Total3B,[15,38]))**2, marker='o', color='red', linestyle="None", markersize=MarkersizeValue)
    elif mode == "6Bartels":
        plt.plot( draw.dates6B, draw.Error6b**2/draw.Error_Total6B**2, marker='o', color='red', linestyle="None", markersize=MarkersizeValue)
    elif mode == "6months":
        plt.plot(draw.datesM6, draw.ErrorM6**2/draw.Error_Total6M**2, marker='o', color='red', linestyle="None", markersize=MarkersizeValue)

    ax.set_ylim(0,1)
    ax.tick_params(axis='x', labeltop=False, labelbottom=True , labelleft=False, labelright=False, direction='in', length=30, width=10)
    ax.tick_params(axis='y', labeltop=False, labelbottom=False, labelleft=True , labelright=False, direction='in', length=30, width=10)

    RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    #plt.text(0.10, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    plt.grid(True,axis='y')
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=70)
    plt.ylabel("Statistical uncertainty percentage", fontsize=60)

    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "StaErrorPercentage_" + mode + "_binmerge_" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "StaErrorPercentage_" + mode + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()



