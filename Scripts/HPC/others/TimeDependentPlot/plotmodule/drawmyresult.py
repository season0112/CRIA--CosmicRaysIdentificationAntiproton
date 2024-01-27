import draw
from draw import *
import TimeDependentTimeRange
import PythonPlotDefaultParameters 

def drawmyresult_StaErrOnly(lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index, RemovedList_3B):
    ## Plot
    fig_mine, ax = plt.subplots()

    if mode == "6months":  
        plt.errorbar(draw.datesM6, draw.M6result*100000, yerr=draw.ErrorM6*100000, marker='o', linestyle="None", markersize=30, markerfacecolor="blue", ecolor="blue", linewidth=5 )
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, draw.B6result*100000, yerr=draw.Error6b*100000, marker='o', linestyle="None", markersize=30, markerfacecolor="blue", ecolor="blue", linewidth=5 )
    elif mode == "3Bartels":
        plt.errorbar(np.delete(draw.dates, RemovedList_3B),  np.delete(draw.B3result, RemovedList_3B)*100000, yerr=np.delete(draw.Error3b, RemovedList_3B)*100000, marker='o',linestyle="None", markersize=30, markerfacecolor="blue",ecolor="blue")

    #RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    #plt.text(0.60, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    ax.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])

    plt.ylabel(r'$\rm{\overline{p}/p}$ ($\times$ $10^{5}$)', horizontalalignment='right', y=1.0)
    plt.text(0.13, 0.85, "Statistical uncertainty only", transform = ax.transAxes, fontsize=60)

    ## Save Plot
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + mode + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +"_StaErrOnly.pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + mode + "_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_StaErrOnly.pdf")
    plt.close()


def drawmyresult(lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index, RemovedList_3B):
    ## Plot
    fig_mine, ax = plt.subplots(figsize=(30, 18))

    if mode == "6months":
        plt.errorbar(draw.datesM6, draw.M6result*100000, yerr=draw.Error_Total6M*100000, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5 )
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, draw.B6result*100000, yerr=draw.Error_Total6B*100000, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5 )
    elif mode == "3Bartels":
        plt.errorbar(np.delete(draw.dates, RemovedList_3B),  np.delete(draw.B3result, RemovedList_3B)*100000, yerr=np.delete(draw.Error_Total3B, RemovedList_3B)*100000, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue")

    RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    plt.text(0.60, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    ax.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    ax.tick_params(axis='x', labeltop=False, labelbottom=True , labelleft=False, labelright=False, direction='in', length=30, width=10)
    ax.tick_params(axis='y', labeltop=False, labelbottom=False, labelleft=True , labelright=False, direction='in', length=30, width=10)

    plt.xticks(fontsize=50)
    plt.yticks(fontsize=70)
    plt.ylabel("Antiproton ratio (10$^{-5}$)",fontsize=70)

    ## Save Plot
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + mode + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + mode + "_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    plt.close()


def drawmyresult_alltime(rigidityrange, lowworkpath, intermediateworkpath, index, binmerge, RemovedList_3B):

    fig_alltime, ax_alltime = plt.subplots(figsize=(30,18))

    #plt.errorbar(draw.datesM6, draw.M6result*10000, yerr=draw.Error_Total6M*10000, marker='o', linestyle="None", markersize=15, markerfacecolor="red" , ecolor="red" , linewidth=5, label="6 Months" )
    plt.errorbar(draw.dates6B, draw.B6result*10000, yerr=draw.Error_Total6B*10000, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5, label="6 Bartel's Rotation")
    plt.errorbar( np.delete(draw.dates, RemovedList_3B), np.delete(draw.B3result, RemovedList_3B)*10000, yerr=np.delete(draw.Error_Total3B, RemovedList_3B)*10000, marker='o', linestyle="None", markersize=25, markerfacecolor="green", ecolor="green", linewidth=5, label="3 Bartel's Rotation")


    ax_alltime.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.legend(loc='best',fontsize=50)

    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "AllTimeBin_Low_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "AllTimeBin_Low_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
    plt.close()


def drawTimeDependentRatio_vs_Rigditiy(rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, mode, plotrange):

    plotrange_Low          = range(1, 17, int(binmerge))
    plotrange_Intermediate = range(9, 29, int(binmerge))
    timeindexall_6M = int(TimeDependentTimeRange.SplitTotal_6M/TimeDependentTimeRange.mergestep_6M)
    timeindexall_3B = int(TimeDependentTimeRange.SplitTotal_3B/TimeDependentTimeRange.mergestep_3B)
    timeindexall_6B = int(TimeDependentTimeRange.SplitTotal_6B/TimeDependentTimeRange.mergestep_6B)
  
    #RevmovePointInLow          = 1
    #RevmovePointInIntermediate = 3
    RevmovePointInLow          = 4
    RevmovePointInIntermediate = 0 

    #### Load time averaged ratio (low and intermediate range)
    # low
    time_averaged_ratio_low                = np.array([])
    time_averaged_ratio_with_effective_low = np.array([])
    time_averaged_error_low                = np.array([])
    f_averaged_low = TFile(lowworkpath + "/totalall/Time_Averaged_ratio_Low/binmerge" + str(binmerge) + "/plots/Ratio_pass7.8.root")
    g_ratio_low                = f_averaged_low.Get("ratio_tof_TRDeff_0.94_TOFeff_0.95")
    g_ratio_with_effective_low = f_averaged_low.Get("ratio_tof_with_effective_TRDeff_0.94_TOFeff_0.95")
    g_error_low                = f_averaged_low.Get("g_error_TRDeff_0.94_TOFeff_0.95")
    for i in range(g_ratio_low.GetN()):
        time_averaged_ratio_low                = np.append(time_averaged_ratio_low               , g_ratio_low.GetY()[i])
        time_averaged_ratio_with_effective_low = np.append(time_averaged_ratio_with_effective_low, g_ratio_with_effective_low.GetY()[i])
        time_averaged_error_low                = np.append(time_averaged_error_low               , g_error_low.GetY()[i])
    time_averaged_ratio_low                = time_averaged_ratio_low[0:-RevmovePointInLow]
    time_averaged_ratio_with_effective_low = time_averaged_ratio_with_effective_low[0:-RevmovePointInLow]
    time_averaged_error_low                = time_averaged_error_low[0:-RevmovePointInLow]

    # intermediate
    time_averaged_ratio_intermediate                = np.array([])
    time_averaged_ratio_with_effective_intermediate = np.array([])
    time_averaged_error_intermediate                = np.array([])
    f_averaged_intermediate = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/binmerge" + str(binmerge) + "/intermediate_"+ str(trackerpattern) + "_" + str(richcut) + "_" + "pass7.8" + "binmerge" + str(binmerge) + ".root")
    g_ratio_intermediate                = f_averaged_intermediate.Get("g_ratio_TRDeff_0.99")
    g_ratio_with_effective_intermediate = f_averaged_intermediate.Get("g_ratio_with_effective_acceptance_TRDeff_0.99")
    g_error_intermediate                = f_averaged_intermediate.Get("g_error_TRDeff_0.99")
    for i in range(g_ratio_intermediate.GetN()):
        time_averaged_ratio_intermediate                = np.append(time_averaged_ratio_intermediate               , g_ratio_intermediate.GetY()[i])
        time_averaged_ratio_with_effective_intermediate = np.append(time_averaged_ratio_with_effective_intermediate, g_ratio_with_effective_intermediate.GetY()[i])
        time_averaged_error_intermediate                = np.append(time_averaged_error_intermediate               , g_error_intermediate.GetY()[i])
    time_averaged_ratio_intermediate                = time_averaged_ratio_intermediate[RevmovePointInIntermediate:]
    time_averaged_ratio_with_effective_intermediate = time_averaged_ratio_with_effective_intermediate[RevmovePointInIntermediate:]
    time_averaged_error_intermediate                = time_averaged_error_intermediate[RevmovePointInIntermediate:]

    # total error
    f_totalerror = TFile("/home/bo791269/Software/AntiprotonAnalysis/Macros/CheckResult.root")
    g_TotalError_Low          = f_totalerror.Get("g_TotalError_Low")
    g_TotalError_Intermediate = f_totalerror.Get("g_TotalError_Intermediate")
    #low
    g_TotalError_Low.RemovePoint(0) # remove 1.0-1.16GV (not included)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 5.37-5.9GV   (last 1th dependent point) (only have to delete once, since 5.9-6.47GV has been removed.)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 4.88-5.37 GV (last 2th dependent point)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 4.43-4.88 GV (last 2th dependent point)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 4.02-4.43 GV (last 3th dependent point)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 3.64-4.02 GV (last 3th dependent point)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 3.29-3.64 GV (last 4th dependent point)
    g_TotalError_Low.RemovePoint(g_TotalError_Low.GetN()-1) # remove 2.97-3.29 GV (last 4th dependent point)
    time_averaged_totalerror_low = np.array([])  
    for i in range(0, g_TotalError_Low.GetN(), 2): # assume binmerge=2, and RevmovePointInLow = 4, RevmovePointInIntermediate = 0
        time_averaged_totalerror_low = np.append( time_averaged_totalerror_low, (g_TotalError_Low.GetY()[i] + g_TotalError_Low.GetY()[i+1] )/2 )
    #intermediate
    time_averaged_totalerror_intermediate = np.array([])
    for i in range(0, g_TotalError_Intermediate.GetN(), 2): # assume binmerge=2, and RevmovePointInLow = 4, RevmovePointInIntermediate = 0
        time_averaged_totalerror_intermediate = np.append( time_averaged_totalerror_intermediate, (g_TotalError_Intermediate.GetY()[i] + g_TotalError_Intermediate.GetY()[i+1] )/2 )


 
    #### Plot

    #### 6 months
    fig_Ratio_R, ax_Ratio_R = plt.subplots(figsize=(30,18))
    for timeindex in range(timeindexall_6M):
        ## Load Time Dependent result in Low
        Rigiditybincenter_low = []
        M6ratio_vs_R_low = []
        M6ratioerror_vs_R_low = []
        for Rindex in plotrange_Low:
            Rigiditybincenter_low.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_6months.txt") as f6m:
                lines6m = f6m.readlines()
            M6ratio_vs_R_low.append( float(lines6m[timeindex].strip())*100000 )
            with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_6months.txt") as e6m:
                error6m = e6m.readlines()
            M6ratioerror_vs_R_low.append( float(error6m[timeindex].strip())*100000 )

        ## Plot Time Dependent in Low range
        if timeindex == 0:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], M6ratio_vs_R_low[0:-RevmovePointInLow], yerr=M6ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green"   , ecolor="green", label='6 months')
        else:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], M6ratio_vs_R_low[0:-RevmovePointInLow], yerr=M6ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green"   , ecolor="green")
        ## Load Time Dependent result in Intermediate 
        Rigiditybincenter_intermediate = []
        M6ratio_vs_R_intermediate = []
        M6ratioerror_vs_R_intermediate = []
        for Rindex in plotrange_Intermediate:
            Rigiditybincenter_intermediate.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex++int(binmerge)]) + "_6months.txt") as f6m:
                lines6m=f6m.readlines()
            M6ratio_vs_R_intermediate.append( float(lines6m[timeindex].strip())*100000 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex++int(binmerge)]) + "_6months.txt") as e6m:
                error6m = e6m.readlines()
            M6ratioerror_vs_R_intermediate.append( float(error6m[timeindex].strip())*100000 )

        ## Plot Time Dependent in Intermediate range
        plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], M6ratio_vs_R_intermediate[RevmovePointInIntermediate:], yerr=M6ratioerror_vs_R_intermediate[RevmovePointInIntermediate:], marker='o', linestyle="None", markersize=12, markerfacecolor="green"   , ecolor="green")

    ## Plot time averaged ratio
    plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow]                , time_averaged_ratio_low                , yerr=time_averaged_totalerror_low, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8, label='Time Averaged')
    plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], time_averaged_ratio_intermediate*100000, yerr=time_averaged_totalerror_intermediate*100000, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8)

    ax_Ratio_R.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=70)
    plt.legend(loc='best', fontsize=35)
    plt.xlabel("Rigidity (GV)",fontsize=50)
    plt.ylabel("Antiproton ratio (10$^{-5}$)", fontsize=60)
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_6M_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_6M_" + "binmerge" + binmerge +".pdf")
    plt.close()


    #### 3 Bartels
    fig_Ratio_R, ax_Ratio_R = plt.subplots(figsize=(30,18))
    for timeindex in range(timeindexall_3B):    
        ## Load Time Dependent result in Low
        Rigiditybincenter_low = []
        B3ratio_vs_R_low = []
        B3ratioerror_vs_R_low = []
        for Rindex in plotrange_Low:
            Rigiditybincenter_low.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_3BartalRotation.txt") as f3b:
                lines3b = f3b.readlines()
            B3ratio_vs_R_low.append( float(lines3b[timeindex].strip())*100000 )
            with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_3BartalRotation.txt") as e3b:
                error3b = e3b.readlines()
            B3ratioerror_vs_R_low.append( float(error3b[timeindex].strip())*100000 )

        ## Plot Time Dependent in Low range
        if timeindex == 0:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], B3ratio_vs_R_low[0:-RevmovePointInLow], yerr=B3ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green", ecolor="green", label='3 Bartels')
        else:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], B3ratio_vs_R_low[0:-RevmovePointInLow], yerr=B3ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green", ecolor="green")
        ## Load Time Dependent result in Intermediate
        Rigiditybincenter_intermediate = []
        B3ratio_vs_R_intermediate = []
        B3ratioerror_vs_R_intermediate = []
        for Rindex in plotrange_Intermediate:
            Rigiditybincenter_intermediate.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex++int(binmerge)]) + "_3BartalRotation.txt") as f3b:
                lines3b = f3b.readlines()
            B3ratio_vs_R_intermediate.append( float(lines3b[timeindex].strip())*100000 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_3BartalRotation.txt") as e3b:
                error3b = e3b.readlines()
            B3ratioerror_vs_R_intermediate.append( float(error3b[timeindex].strip())*100000 )

        ## Plot Time Dependent in Intermediate range
        plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], B3ratio_vs_R_intermediate[RevmovePointInIntermediate:], yerr=B3ratioerror_vs_R_intermediate[RevmovePointInIntermediate:], marker='o', linestyle="None", markersize=12, markerfacecolor="green", ecolor="green")
    ## Plot time averaged ratio
    plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], time_averaged_ratio_low, yerr=time_averaged_totalerror_low, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8, label='Time Averaged')
    plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], time_averaged_ratio_intermediate*100000, yerr=time_averaged_totalerror_intermediate*100000, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8)

    ax_Ratio_R.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=70)
    plt.legend(loc='best', fontsize=35)
    plt.xlabel("Rigidity (GV)",fontsize=50)
    plt.ylabel("Antiproton ratio (10$^{-5}$)", fontsize=60)
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_3B_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_3B_" + "binmerge" + binmerge +".pdf")
    plt.close()

    ## 6 Bartels
    fig_Ratio_R, ax_Ratio_R = plt.subplots(figsize=(30,18))
    for timeindex in range(timeindexall_6B):
        ## Load Time Dependent result in Low
        Rigiditybincenter_low = []
        B6ratio_vs_R_low = []
        B6ratioerror_vs_R_low = []
        for Rindex in plotrange_Low:
            Rigiditybincenter_low.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_6BartalRotation.txt") as f6b:
                lines6b = f6b.readlines()
            B6ratio_vs_R_low.append( float(lines6b[timeindex].strip())*100000 )
            with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_6BartalRotation.txt") as e6b:
                error6b = e6b.readlines()
            B6ratioerror_vs_R_low.append( float(error6b[timeindex].strip())*100000 )

        ## Plot Time Dependent in Low range
        if timeindex == 0:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], B6ratio_vs_R_low[0:-RevmovePointInLow], yerr=B6ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green"  , ecolor="green", label='6 Bartels Rotations results')
        else:
            plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow], B6ratio_vs_R_low[0:-RevmovePointInLow], yerr=B6ratioerror_vs_R_low[0:-RevmovePointInLow], marker='o', linestyle="None", markersize=12, markerfacecolor="green"  , ecolor="green")
        ## Load Time Dependent result in Intermediate
        Rigiditybincenter_intermediate = []
        B6ratio_vs_R_intermediate = []
        B6ratioerror_vs_R_intermediate = []
        for Rindex in plotrange_Intermediate:
            Rigiditybincenter_intermediate.append( (binning.published2016binnings_center[Rindex] + binning.published2016binnings_center[Rindex+1])/2 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex++int(binmerge)]) + "_6BartalRotation.txt") as f6b:
                lines6b = f6b.readlines()
            B6ratio_vs_R_intermediate.append( float(lines6b[timeindex].strip())*100000 )
            with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[Rindex]) + "_" + "{:g}".format(binning.published2016binnings[Rindex+int(binmerge)]) + "_6BartalRotation.txt") as e6b:
                error6b = e6b.readlines()
            B6ratioerror_vs_R_intermediate.append( float(error6b[timeindex].strip())*100000 )
        ## Plot Time Dependent in Intermediate range
        plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], B6ratio_vs_R_intermediate[RevmovePointInIntermediate:], yerr=B6ratioerror_vs_R_intermediate[RevmovePointInIntermediate:], marker='o', linestyle="None", markersize=12, markerfacecolor="green"  , ecolor="green")
    #### Plot time averaged ratio
    plt.errorbar( Rigiditybincenter_low[0:-RevmovePointInLow]                , time_averaged_ratio_low                , yerr=time_averaged_totalerror_low, marker='o', linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8, label='Time Averaged result')
    plt.errorbar( Rigiditybincenter_intermediate[RevmovePointInIntermediate:], time_averaged_ratio_intermediate*100000, yerr=time_averaged_totalerror_intermediate*100000, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8)

    ax_Ratio_R.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=70)
    plt.legend(loc='best', fontsize=50, frameon=False)
    plt.xlabel("Rigidity (GV)", fontsize=50)
    #plt.ylabel("Antiproton to proton ratio (10$^{-5}$)", fontsize=50)
    plt.ylabel(r'$\rm{\overline{p}/p}$ ($\times$ $10^{-5}$)', fontsize=60)

    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_6B_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "TimeDependentAndAveragedRatio_vs_R_6B_" + "binmerge" + binmerge +".pdf")
    plt.close()    


def draw6Bartel6MonthsCompare(lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index):

    ## Plot
    fig_mine, ax = plt.subplots(figsize=(30,18))

    plt.errorbar(draw.datesM6,  draw.M6result*10000, yerr=draw.Error_Total6M*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue", label='6 Months' )
    plt.errorbar(draw.dates6B,  draw.B6result*10000, yerr=draw.Error_Total6B*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="red" ,ecolor="red",  label='6 Bartels')

    ax.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    ax.tick_params(direction='in', length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.legend(loc='best', fontsize=35)
    plt.ylabel("Antiproton ratio (10$^{-4}$)", fontsize=60)
    plt.title( "{:g}".format(binning.published2016binnings[index])  + "-" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ' GV', fontsize=60)

    #plt.ylabel( r'$\frac{\rm IHEP-\rm Taiwan}{\rm Taiwan}$'+" (%)",fontsize=40)

    ## Save Plot
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "6Bartel6MonthsCompare_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6Bartel6MonthsCompare_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()




