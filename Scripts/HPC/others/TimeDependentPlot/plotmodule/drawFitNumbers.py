import draw
import matplotlib
from draw import *
import uncertainties
from uncertainties import unumpy
from uncertainties import ufloat
import PythonPlotDefaultParameters

def drawFitNumbers(mode, rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, index, time_averaged_proton, time_averaged_antiproton, MeasuringTime_6B, yaxisfactor):

    ## 1. Plot Antiproton number

    ## Print time averaged total numbers
    #print("Time averaged total Pbar numbers:" + str(time_averaged_antiproton))

    # Plot
    fig, ax = plt.subplots()

    if mode == "6months":
        plt.errorbar(draw.datesM6, draw.M6_PbarNumber, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5 )
        #print("6months Total Pbarnumber:" + str(sum(draw.M6_PbarNumber)))
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, draw.B6_PbarNumber / (MeasuringTime_6B[:, int((index+1)/2)]/60/60/24), yerr=0, marker='o', linestyle="None", markersize=30, markerfacecolor="blue", ecolor="blue", linewidth=5)
        #print("6Bartels Total Pbarnumber:" + str(sum(draw.B6_PbarNumber)))
    elif mode == "3Bartels":
        plt.errorbar(draw.dates[0:15]+draw.dates[16:],  np.delete(draw.B3_PbarNumber, 15), yerr=0, marker='o',linestyle="None",markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5)
        #print("3Bartels Total Pbarnumber:" + str(sum(draw.B3_PbarNumber)))

    if mode == "6Bartels":
        plt.ylabel(r'$\rm{N}_{\rm{\overline{p}}}$ / $\Delta \rm{T}$ / (1/day)', horizontalalignment='right', y=1.0) #"Antiproton Numbers / Measuring time(days)"
    else:
        plt.ylabel("Antiproton Numbers")

    # Save Plot
    if rigidityrange == "low":
        if mode == "6months":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "AntiprotonNumber_6months_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "AntiprotonNumber_6BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "AntiprotonNumber_3BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")

    elif rigidityrange == "intermediate":
        if mode == "6months":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "AntiprotonNumber_6months_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "AntiprotonNumber_6BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "AntiprotonNumber_3BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()


    ## 2. Plot Proton number

    ## Print time averaged total numbers
    #print("Time averaged total Proton numbers:" + str(time_averaged_proton))

    # Plot
    fig_proton, ax_proton = plt.subplots(figsize=(30,18))
    if mode == "6months":
        plt.errorbar(draw.datesM6,  draw.M6_ProtonNumber, yerr=0, marker='o',linestyle="None",markersize=25, markerfacecolor="blue",ecolor="blue", linewidth=5)
        #print("6months Total ProtonNumber: " + str(sum(draw.M6_ProtonNumber)))
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B,  draw.B6_ProtonNumber / (MeasuringTime_6B[:, int((index+1)/2)]/60/60/24), yerr=0, marker='o',linestyle="None",markersize=25, markerfacecolor="blue",ecolor="blue", linewidth=5)
        #print("6Bartels Total ProtonNumber: " + str(sum(draw.B6_ProtonNumber)))
    elif mode == "3Bartels":
        plt.errorbar(draw.dates[0:15]+draw.dates[16:],  np.delete(draw.B3_ProtonNumber, 15), yerr=0, marker='o',linestyle="None",markersize=25, markerfacecolor="blue", ecolor="blue", linewidth=5)
        #print("3Bartels Total ProtonNumber: " + str(sum(draw.B3_ProtonNumber)))

    ax_proton.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)

    if mode == "6Bartels":
        plt.ylabel("Proton Numbers / Measuring time(days) ",fontsize=50)
    else:
        plt.ylabel("Proton Numbers",fontsize=50)
    ax_proton.ticklabel_format(style='sci', scilimits=(0,0), axis='y', useMathText=True)
    ax_proton.yaxis.offsetText.set_fontsize(80)

    # Save Plot
    if rigidityrange == "low":
        if mode == "6months":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "ProtonNumber_6months_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "ProtonNumber_6BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "ProtonNumber_3BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")

    elif rigidityrange == "intermediate":
        if mode == "6months":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "ProtonNumber_6months_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "ProtonNumber_6BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "ProtonNumber_3BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")

    plt.close()


def drawBinmergedRatio(time_averaged_antiproton, time_averaged_proton, lowworkpath, intermediateworkpath, binmerge, PbarNumber_Merged_6M, ProtonNumber_Merged_6M, PbarNumber_Merged_6B, ProtonNumber_Merged_6B, PbarNumber_Merged_3B, ProtonNumber_Merged_3B, mode, rigidityrange, trackerpattern, richcut):

    RevmovePointInLow          = 4
    RevmovePointInIntermediate = 0


    #### Load Time Averaged
    time_averaged_ratio_low                = np.array([])
    time_averaged_ratio_with_effective_low = np.array([])
    time_averaged_error_low                = np.array([])
    time_averaged_proton_low               = np.array([])
    time_averaged_antiproton_low           = np.array([])
    time_averaged_ratio_intermediate                = np.array([])
    time_averaged_ratio_with_effective_intermediate = np.array([])
    time_averaged_error_intermediate                = np.array([])
    time_averaged_proton_intermediate               = np.array([])
    time_averaged_antiproton_intermediate           = np.array([])

    f_averaged_low = TFile(lowworkpath + "/totalall/Time_Averaged_ratio_Low/binmerge" + str(binmerge) + "/plots/Ratio_pass7.8.root")
    # TOFEff:0.61-0.99 TRDEff:0.60-1.00
    TRDeff = "0.98"  # 0.94
    TOFeff = "0.97"  # 0.95
    g_ratio_low                = f_averaged_low.Get("ratio_tof_TRDeff_"                + TRDeff + "_TOFeff_" + TOFeff)
    g_ratio_with_effective_low = f_averaged_low.Get("ratio_tof_with_effective_TRDeff_" + TRDeff + "_TOFeff_" + TOFeff)
    g_error_low                = f_averaged_low.Get("g_error_TRDeff_"                  + TRDeff + "_TOFeff_" + TOFeff)
    g_proton_low               = f_averaged_low.Get("g_proton_number_TRDeff_"          + TRDeff + "_TOFeff_" + TOFeff)
    g_antiproton_low           = f_averaged_low.Get("g_antiproton_number_TRDeff_"      + TRDeff + "_TOFeff_" + TOFeff)
    f_averaged_intermediate = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/binmerge" + str(binmerge) + "/intermediate_"+ str(trackerpattern) + "_" + str(richcut) + "_" + "pass7.8" + "binmerge" + str(binmerge) + ".root")
    TRDeff_intermediate = "0.99"
    g_ratio_intermediate                = f_averaged_intermediate.Get("g_ratio_TRDeff_"                           + TRDeff_intermediate)
    g_ratio_with_effective_intermediate = f_averaged_intermediate.Get("g_ratio_with_effective_acceptance_TRDeff_" + TRDeff_intermediate)
    g_error_intermediate                = f_averaged_intermediate.Get("g_error_TRDeff_"                           + TRDeff_intermediate)
    g_proton_intermediate               = f_averaged_intermediate.Get("g_proton_TRDeff_"                          + TRDeff_intermediate)
    g_antiproton_intermediate           = f_averaged_intermediate.Get("g_antiproton_TRDeff_"                      + TRDeff_intermediate)

    for i in range(g_ratio_with_effective_low.GetN()):
        time_averaged_ratio_low                = np.append(time_averaged_ratio_low               , g_ratio_low.GetY()[i])
        time_averaged_ratio_with_effective_low = np.append(time_averaged_ratio_with_effective_low, g_ratio_with_effective_low.GetY()[i])
        time_averaged_error_low                = np.append(time_averaged_error_low               , g_error_low.GetY()[i])
        time_averaged_proton_low               = np.append(time_averaged_proton_low              , g_proton_low.GetY()[i])
        time_averaged_antiproton_low           = np.append(time_averaged_antiproton_low          , g_antiproton_low.GetY()[i])
    for i in range(g_ratio_with_effective_intermediate.GetN()):
        time_averaged_ratio_intermediate                = np.append(time_averaged_ratio_intermediate               , g_ratio_intermediate.GetY()[i])
        time_averaged_ratio_with_effective_intermediate = np.append(time_averaged_ratio_with_effective_intermediate, g_ratio_with_effective_intermediate.GetY()[i])
        time_averaged_error_intermediate                = np.append(time_averaged_error_intermediate               , g_error_intermediate.GetY()[i])
        time_averaged_proton_intermediate               = np.append(time_averaged_proton_intermediate              , g_proton_intermediate.GetY()[i])
        time_averaged_antiproton_intermediate           = np.append(time_averaged_antiproton_intermediate          , g_antiproton_intermediate.GetY()[i])

    #### Load total error
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


    #### Load Merged Numbers
    PbarNumber_Merged_6M_low           = [] 
    ProtonNumber_Merged_6M_low         = [] 
    PbarNumber_Merged_6B_low           = []         
    ProtonNumber_Merged_6B_low         = []
    PbarNumber_Merged_3B_low           = []         
    ProtonNumber_Merged_3B_low         = []
    PbarNumber_Merged_6M_intermediate  = []
    ProtonNumber_Merged_6M_intermediate= []
    PbarNumber_Merged_6B_intermediate  = []
    ProtonNumber_Merged_6B_intermediate= []
    PbarNumber_Merged_3B_intermediate  = []
    ProtonNumber_Merged_3B_intermediate= []

    plotrange_Low          = range(1, 17, int(binmerge))
    plotrange_Intermediate = range(9, 29, int(binmerge))

    ## Load 6 months
    for index in plotrange_Low:
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as pbarNum_6m:
            PbarNum_6m = pbarNum_6m.readlines()
        M6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6m)))
        PbarNumber_Merged_6M_low.append(sum(M6_PbarNumber))
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as protonNum_6m:
            ProtonNum_6m = protonNum_6m.readlines()
        M6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6m)))
        ProtonNumber_Merged_6M_low.append(sum(M6_ProtonNumber))
    for index in plotrange_Intermediate:
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as pbarNum_6m:
            PbarNum_6m = pbarNum_6m.readlines()
        M6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6m)))
        PbarNumber_Merged_6M_intermediate.append(sum(M6_PbarNumber))
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as protonNum_6m:
            ProtonNum_6m = protonNum_6m.readlines()
        M6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6m)))
        ProtonNumber_Merged_6M_intermediate.append(sum(M6_ProtonNumber))

    ## Load 6 Bartels
    for index in plotrange_Low:
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as pbarNum_6b:
            PbarNum_6b = pbarNum_6b.readlines()
        B6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6b)))
        PbarNumber_Merged_6B_low.append(sum(B6_PbarNumber))
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as protonNum_6b:
            ProtonNum_6b = protonNum_6b.readlines()
        B6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6b)))
        ProtonNumber_Merged_6B_low.append(sum(B6_ProtonNumber))
    for index in plotrange_Intermediate:
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as pbarNum_6b:
            PbarNum_6b = pbarNum_6b.readlines()
        B6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6b)))
        PbarNumber_Merged_6B_intermediate.append(sum(B6_PbarNumber))
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as protonNum_6b:
            ProtonNum_6b = protonNum_6b.readlines()
        B6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6b)))
        ProtonNumber_Merged_6B_intermediate.append(sum(B6_ProtonNumber))

    ## Load 3 Bartels
    for index in plotrange_Low:
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as pbarNum_3b:
            PbarNum_3b = pbarNum_3b.readlines()
        B3_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_3b)))
        PbarNumber_Merged_3B_low.append(sum(B3_PbarNumber))
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as protonNum_3b:
            ProtonNum_3b = protonNum_3b.readlines()
        B3_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_3b)))
        ProtonNumber_Merged_3B_low.append(sum(B3_ProtonNumber))
    for index in plotrange_Intermediate:
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as pbarNum_3b:
            PbarNum_3b = pbarNum_3b.readlines()
        B3_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_3b)))
        PbarNumber_Merged_3B_intermediate.append(sum(B3_PbarNumber))
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as protonNum_3b:
            ProtonNum_3b = protonNum_3b.readlines()
        B3_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_3b)))
        ProtonNumber_Merged_3B_intermediate.append(sum(B3_ProtonNumber))

    PbarNumber_Merged_6M_low           = np.array( PbarNumber_Merged_6M_low )
    ProtonNumber_Merged_6M_low         = np.array( ProtonNumber_Merged_6M_low )
    PbarNumber_Merged_6B_low           = np.array( PbarNumber_Merged_6B_low )
    ProtonNumber_Merged_6B_low         = np.array( ProtonNumber_Merged_6B_low )
    PbarNumber_Merged_3B_low           = np.array( PbarNumber_Merged_3B_low )
    ProtonNumber_Merged_3B_low         = np.array( ProtonNumber_Merged_3B_low )
    PbarNumber_Merged_6M_intermediate  = np.array( PbarNumber_Merged_6M_intermediate )
    ProtonNumber_Merged_6M_intermediate= np.array( ProtonNumber_Merged_6M_intermediate )
    PbarNumber_Merged_6B_intermediate  = np.array( PbarNumber_Merged_6B_intermediate )
    ProtonNumber_Merged_6B_intermediate= np.array( ProtonNumber_Merged_6B_intermediate )
    PbarNumber_Merged_3B_intermediate  = np.array( PbarNumber_Merged_3B_intermediate )
    ProtonNumber_Merged_3B_intermediate= np.array( ProtonNumber_Merged_3B_intermediate )


    #### Plot
    ## Define variables
    time_averaged_antiproton_low          = time_averaged_antiproton_low         [0:-RevmovePointInLow] * 100000
    time_averaged_proton_low              = time_averaged_proton_low             [0:-RevmovePointInLow]
    time_averaged_antiproton_intermediate = time_averaged_antiproton_intermediate[RevmovePointInIntermediate:] * 100000
    time_averaged_proton_intermediate     = time_averaged_proton_intermediate    [RevmovePointInIntermediate:] 
    time_averaged_totalerror_low          = time_averaged_totalerror_low          * 100000
    time_averaged_totalerror_intermediate = time_averaged_totalerror_intermediate * 100000

    PbarNumber_Merged_6M_low            = PbarNumber_Merged_6M_low[0:-RevmovePointInLow] 
    ProtonNumber_Merged_6M_low          = ProtonNumber_Merged_6M_low[0:-RevmovePointInLow] 
    PbarNumber_Merged_6M_intermediate   = PbarNumber_Merged_6M_intermediate[RevmovePointInIntermediate:]
    ProtonNumber_Merged_6M_intermediate = ProtonNumber_Merged_6M_intermediate[RevmovePointInIntermediate:] 

    PbarNumber_Merged_6B_low            = PbarNumber_Merged_6B_low[0:-RevmovePointInLow] 
    ProtonNumber_Merged_6B_low          = ProtonNumber_Merged_6B_low[0:-RevmovePointInLow]
    PbarNumber_Merged_6B_intermediate   = PbarNumber_Merged_6B_intermediate[RevmovePointInIntermediate:] 
    ProtonNumber_Merged_6B_intermediate = ProtonNumber_Merged_6B_intermediate[RevmovePointInIntermediate:]

    PbarNumber_Merged_3B_low            = PbarNumber_Merged_3B_low[0:-RevmovePointInLow] 
    ProtonNumber_Merged_3B_low          = ProtonNumber_Merged_3B_low[0:-RevmovePointInLow] 
    PbarNumber_Merged_3B_intermediate   = PbarNumber_Merged_3B_intermediate[RevmovePointInIntermediate:] 
    ProtonNumber_Merged_3B_intermediate = ProtonNumber_Merged_3B_intermediate[RevmovePointInIntermediate:] 


    PbarNumber_Merged_6B_low_uncertainty            = unumpy.uarray(PbarNumber_Merged_6B_low           , np.sqrt(PbarNumber_Merged_6B_low))
    ProtonNumber_Merged_6B_low_uncertainty          = unumpy.uarray(ProtonNumber_Merged_6B_low         , np.sqrt(ProtonNumber_Merged_6B_low))
    PbarNumber_Merged_6B_intermediate_uncertainty   = unumpy.uarray(PbarNumber_Merged_6B_intermediate  , np.sqrt(PbarNumber_Merged_6B_intermediate))
    ProtonNumber_Merged_6B_intermediate_uncertainty = unumpy.uarray(ProtonNumber_Merged_6B_intermediate, np.sqrt(ProtonNumber_Merged_6B_intermediate))
    PbarRatio_Merged_6B_low_uncertainty             = PbarNumber_Merged_6B_low_uncertainty         /ProtonNumber_Merged_6B_low_uncertainty * 100000  
    PbarRatio_Merged_6B_intermediate_uncertainty    = PbarNumber_Merged_6B_intermediate_uncertainty/ProtonNumber_Merged_6B_intermediate_uncertainty * 100000

    PbarNumber_Merged_3B_low_uncertainty            = unumpy.uarray(PbarNumber_Merged_3B_low           , np.sqrt(PbarNumber_Merged_3B_low))
    ProtonNumber_Merged_3B_low_uncertainty          = unumpy.uarray(ProtonNumber_Merged_3B_low         , np.sqrt(ProtonNumber_Merged_3B_low))
    PbarNumber_Merged_3B_intermediate_uncertainty   = unumpy.uarray(PbarNumber_Merged_3B_intermediate  , np.sqrt(PbarNumber_Merged_3B_intermediate))
    ProtonNumber_Merged_3B_intermediate_uncertainty = unumpy.uarray(ProtonNumber_Merged_3B_intermediate, np.sqrt(ProtonNumber_Merged_3B_intermediate))
    PbarRatio_Merged_3B_low_uncertainty             = PbarNumber_Merged_3B_low_uncertainty         /ProtonNumber_Merged_3B_low_uncertainty * 100000
    PbarRatio_Merged_3B_intermediate_uncertainty    = PbarNumber_Merged_3B_intermediate_uncertainty/ProtonNumber_Merged_3B_intermediate_uncertainty * 100000

    PbarNumber_Merged_6M_low_uncertainty            = unumpy.uarray(PbarNumber_Merged_6M_low           , np.sqrt(PbarNumber_Merged_6M_low))
    ProtonNumber_Merged_6M_low_uncertainty          = unumpy.uarray(ProtonNumber_Merged_6M_low         , np.sqrt(ProtonNumber_Merged_6M_low))
    PbarNumber_Merged_6M_intermediate_uncertainty   = unumpy.uarray(PbarNumber_Merged_6M_intermediate  , np.sqrt(PbarNumber_Merged_6M_intermediate))
    ProtonNumber_Merged_6M_intermediate_uncertainty = unumpy.uarray(ProtonNumber_Merged_6M_intermediate, np.sqrt(ProtonNumber_Merged_6M_intermediate))
    PbarRatio_Merged_6M_low_uncertainty             = PbarNumber_Merged_6M_low_uncertainty         /ProtonNumber_Merged_6M_low_uncertainty * 100000
    PbarRatio_Merged_6M_intermediate_uncertainty    = PbarNumber_Merged_6M_intermediate_uncertainty/ProtonNumber_Merged_6M_intermediate_uncertainty * 100000

    ## Plot 6 months
    fig, ax = plt.subplots(figsize=(30,18))

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , time_averaged_antiproton_low/time_averaged_proton_low                  , yerr=time_averaged_totalerror_low, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", markeredgecolor='red', elinewidth=15, capsize=15, capthick=8, label='Time Averaged')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], time_averaged_antiproton_intermediate/time_averaged_proton_intermediate, yerr=time_averaged_totalerror_intermediate, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", markeredgecolor='red', elinewidth=15, capsize=15, capthick=8,)

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , unumpy.nominal_values(PbarRatio_Merged_6M_low_uncertainty)         , yerr=unumpy.std_devs(PbarRatio_Merged_6M_low_uncertainty)         , marker='o',linestyle="None", markersize=20, markerfacecolor="green", ecolor="green", label='Averaged 6 Months results')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], unumpy.nominal_values(PbarRatio_Merged_6M_intermediate_uncertainty), yerr=unumpy.std_devs(PbarRatio_Merged_6M_intermediate_uncertainty), marker='o',linestyle="None", markersize=20, markerfacecolor="green", ecolor="green")

    plt.legend(loc='best', fontsize=35)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=70)
    plt.xlabel("Rigidity (GV)"               , fontsize=60, horizontalalignment='right', x=1.0)
    plt.ylabel("Antiproton ratio (10$^{-5}$)", fontsize=60, horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_6months_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_6months_" + "binmerge" + binmerge +".pdf")

    plt.close()

    ## Plot 6 Bartels 
    fig, ax = plt.subplots(figsize=(30,18))

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , time_averaged_antiproton_low/time_averaged_proton_low                  , yerr=time_averaged_totalerror_low         , marker='o',linestyle="None", markersize=25, markerfacecolor="red", markeredgecolor="red", ecolor="red", elinewidth=10, capsize=15, capthick=8, label='Time Averaged result')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], time_averaged_antiproton_intermediate/time_averaged_proton_intermediate, yerr=time_averaged_totalerror_intermediate, marker='o',linestyle="None", markersize=25, markerfacecolor="red", markeredgecolor="red", ecolor="red", elinewidth=10, capsize=15, capthick=8,)

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , unumpy.nominal_values(PbarRatio_Merged_6B_low_uncertainty)         , yerr=unumpy.std_devs(PbarRatio_Merged_6B_low_uncertainty)         , marker='o',linestyle="None", markersize=20, markerfacecolor="green", markeredgecolor="green", ecolor="green", label='Averaged 6 Bartels Rotations results')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], unumpy.nominal_values(PbarRatio_Merged_6B_intermediate_uncertainty), yerr=unumpy.std_devs(PbarRatio_Merged_6B_intermediate_uncertainty), marker='o',linestyle="None", markersize=20, markerfacecolor="green", markeredgecolor="green", ecolor="green")  


    plt.legend(loc='best', fontsize=50, frameon=False)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    plt.xlabel("Rigidity (GV)"                              , fontsize=50, horizontalalignment='right', x=1.0)
    #plt.ylabel("Antiproton to proton ratio (10$^{-5}$)"    , fontsize=50, horizontalalignment='right', y=1.0)
    plt.ylabel(r'$\rm{\overline{p}/p}$ ($\times$ $10^{-5}$)', fontsize=60, horizontalalignment='right', y=1.0)

    ax.tick_params(axis='both', which='both', direction='in', length=10, width=3)

    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +".pdf")
    ax.set_yscale('log')
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.get_yaxis().get_major_formatter().labelOnlyBase = False
    ax.set_yticks([1.0, 2.0, 5.0, 10.0])
    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_LogY.pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_LogY.pdf")

    plt.close()

    ## Plot 6 Bartels (With Residuall)
    fig = plt.figure(figsize=(30,18))
    grid = plt.GridSpec(4, 1, wspace=0, hspace=0)
    # Upper plot
    ax1 = fig.add_subplot(grid[0:3,0])

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , time_averaged_antiproton_low/time_averaged_proton_low                  , yerr=time_averaged_totalerror_low         , marker='o',linestyle="None", markersize=35, markerfacecolor="red", markeredgecolor="red", ecolor="red", elinewidth=10, capsize=15, capthick=8, label='Time Averaged result')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], time_averaged_antiproton_intermediate/time_averaged_proton_intermediate, yerr=time_averaged_totalerror_intermediate, marker='o',linestyle="None", markersize=35, markerfacecolor="red", markeredgecolor="red", ecolor="red", elinewidth=10, capsize=15, capthick=8,)

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , unumpy.nominal_values(PbarRatio_Merged_6B_low_uncertainty)         , yerr=unumpy.std_devs(PbarRatio_Merged_6B_low_uncertainty)         , marker='o',linestyle="None", markersize=35, markerfacecolor="green", markeredgecolor="green", ecolor="green", label='Averaged 6 Bartels Rotations results')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], unumpy.nominal_values(PbarRatio_Merged_6B_intermediate_uncertainty), yerr=unumpy.std_devs(PbarRatio_Merged_6B_intermediate_uncertainty), marker='o',linestyle="None", markersize=35, markerfacecolor="green", markeredgecolor="green", ecolor="green")

    plt.legend(loc='best', fontsize=50, frameon=False)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    plt.ylabel(r'$\rm{\overline{p}/p}$ ($\times$ $10^{-5}$)', fontsize=60, horizontalalignment='right', y=1.0)
    ax1.tick_params(axis='x', which='both', direction='in', length=0, width=0)
    ax1.tick_params(axis='y', which='both', direction='in', length=30, width=10)

    # Lower plot
    ax2 = fig.add_subplot(grid[3,0])  #error of flux ratio is much less than flux ratio, so residual should be divided by flux ratio.

    residual_low          = ( (time_averaged_antiproton_low/time_averaged_proton_low) - unumpy.nominal_values(PbarRatio_Merged_6B_low_uncertainty) ) / unumpy.nominal_values(PbarRatio_Merged_6B_low_uncertainty)
    #residual_low          = ( (time_averaged_antiproton_low/time_averaged_proton_low) - unumpy.nominal_values(PbarRatio_Merged_6B_low_uncertainty) ) / time_averaged_totalerror_low
    residual_intermediate = ( (time_averaged_antiproton_intermediate/time_averaged_proton_intermediate) - unumpy.nominal_values(PbarRatio_Merged_6B_intermediate_uncertainty) ) / unumpy.nominal_values(PbarRatio_Merged_6B_intermediate_uncertainty)
    #residual_intermediate = ( (time_averaged_antiproton_intermediate/time_averaged_proton_intermediate) - unumpy.nominal_values(PbarRatio_Merged_6B_intermediate_uncertainty) ) / time_averaged_totalerror_intermediate

    ax2.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , residual_low         , yerr=0, marker='s',linestyle="None", markersize=25, markerfacecolor="k", ecolor="k" )
    ax2.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], residual_intermediate, yerr=0, marker='s',linestyle="None", markersize=25, markerfacecolor="k", ecolor="k" ) 
    #plt.axhline(y=0.5, color='r', linestyle='-')
    ax2.axes.set_ylim(-0.3, 0.1)
    ax2.tick_params(axis='both', which='both', direction='in', length=30, width=10)

    plt.legend(loc='best', fontsize=50, frameon=False)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=45)
    plt.xlabel("Rigidity (GV)"  , fontsize=60, horizontalalignment='right', x=1.0)
    #plt.ylabel(r'$\frac{{\rm{Time}} \; {\rm{Averaged}} - {\rm{Averaged}} \; 6 \; {\rm{Bartels}}}{{\rm{Averaged}} \; 6 \; {\rm{Bartels}}}$', fontsize=25)
    plt.ylabel(r'$\frac{{\rm{Residual}}}{{\rm{Averaged}} \; 6 \; {\rm{Bartels}}}$', fontsize=30, horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_residual.pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_residual.pdf")

    ax1.set_yscale('log')
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.get_yaxis().get_major_formatter().labelOnlyBase = False
    ax1.set_yticks([2.0, 5.0, 10.0])
    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_LogY_residual.pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_6Bartels_" + "binmerge" + binmerge +"_LogY_residual.pdf")
    plt.close()


    ## Plot 3 Bartels
    fig, ax = plt.subplots(figsize=(30,18))

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , time_averaged_antiproton_low/time_averaged_proton_low                  , yerr=time_averaged_totalerror_low, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8, label='Time Averaged')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], time_averaged_antiproton_intermediate/time_averaged_proton_intermediate, yerr=time_averaged_totalerror_intermediate, marker='o',linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", elinewidth=15, capsize=15, capthick=8,)

    plt.errorbar( binning.BinsCenter_TimeDependent_Low[0:-RevmovePointInLow]                , unumpy.nominal_values(PbarRatio_Merged_3B_low_uncertainty)         , yerr=unumpy.std_devs(PbarRatio_Merged_3B_low_uncertainty)         , marker='o',linestyle="None", markersize=20, markerfacecolor="green", ecolor="green", label='Averaged 3 Bartels Rotations results')
    plt.errorbar( binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], unumpy.nominal_values(PbarRatio_Merged_3B_intermediate_uncertainty), yerr=unumpy.std_devs(PbarRatio_Merged_3B_intermediate_uncertainty), marker='o',linestyle="None", markersize=20, markerfacecolor="green", ecolor="green")

    plt.legend(loc='best', fontsize=35)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=70)
    plt.xlabel("Rigidity (GV)"               , fontsize=60, horizontalalignment='right', x=1.0)
    plt.ylabel("Antiproton ratio (10$^{-5}$)", fontsize=60, horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        plt.savefig(lowworkpath + "/totalall/results/plot/" + "MergedFitResult_vs_AveragedResult_3Bartels_" + "binmerge" + binmerge +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath + "/total/results/plot/" + "MergedFitResult_vs_AveragedResult_3Bartels_" + "binmerge" + binmerge +".pdf")

    plt.close()


