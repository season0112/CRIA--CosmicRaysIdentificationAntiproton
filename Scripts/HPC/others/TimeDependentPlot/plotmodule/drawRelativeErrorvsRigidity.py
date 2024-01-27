import draw
import numpy as np
import binning
import matplotlib.pyplot as plt
import PythonPlotDefaultParameters


def drawRelativeErrorvsRigidity(lowworkpath, intermediateworkpath, mode, binmerge):

    # Load Numpy file
    if mode == "6months":
        x_data = draw.datesM6
        nameindex = "_6months"
    elif mode == "3Bartels":
        x_data = draw.dates[0:15]+draw.dates[16:38]+draw.dates[39:]
        nameindex = "_3Bartels"
    elif mode == "6Bartels":
        x_data = draw.dates6B
        nameindex = "_6Bartels"

    PbarResult = np.load(lowworkpath + "/totalall/results/TimeDependentResult" + nameindex + ".npz")

    RevmovePointInLow          = 4
    RevmovePointInIntermediate = 0

    MarkersizeValue = 30

    # Remove first point of 1.16-1.51 GV, which is index=0.
    StartIndex = 1

    #### Plot
    ## RelaiveStatistical
    fig, axes = plt.subplots(1,1, sharex=True)

    plt.plot(binning.BinsCenter_TimeDependent_Low[StartIndex:-RevmovePointInLow]                , np.average(PbarResult['RelativeStatisticalErrorLowAll'],axis=1)[StartIndex:-RevmovePointInLow]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="r", markeredgecolor='r', linewidth=5)
    plt.plot(binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], np.average(PbarResult['RelativeStatisticalErrorIntermediateAll'],axis=1)[RevmovePointInIntermediate:]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="r", markeredgecolor='r', linewidth=5)

    plt.xlabel("|R| / (GV)")
    plt.ylabel("Relative Statistical Error (%)")

    plt.savefig( lowworkpath + "/totalall/results/plot/" + "RelativeStatisticalErrorvsRigidity_" + mode + ".pdf")


    ## RelativeSystematic
    fig, axes = plt.subplots(1,1,figsize=(30, 18), sharex=True)

    plt.plot(binning.BinsCenter_TimeDependent_Low[StartIndex:-RevmovePointInLow]                , np.average(PbarResult['RelativeSystematicErrorLowAll'],axis=1)[StartIndex:-RevmovePointInLow]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="b", linewidth=5, markeredgecolor='b')
    plt.plot(binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], np.average(PbarResult['RelativeSystematicErrorIntermediateAll'],axis=1)[RevmovePointInIntermediate:]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="b", linewidth=5, markeredgecolor='b',)

    plt.xlabel("|R| / (GV)")
    plt.ylabel("Relative Systematic Error (%)")

    plt.savefig( lowworkpath + "/totalall/results/plot/" + "RelativeSystematicErrorvsRigidity_" + mode + ".pdf")


    ## Relaive Together
    fig, axes = plt.subplots(1,1, sharex=True)

    plt.plot(binning.BinsCenter_TimeDependent_Low[StartIndex:-RevmovePointInLow]                , np.average(PbarResult['RelativeStatisticalErrorLowAll'],axis=1)[StartIndex:-RevmovePointInLow]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="r", linewidth=5, markeredgecolor='r', label="Statistical uncertainty")
    plt.plot(binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], np.average(PbarResult['RelativeStatisticalErrorIntermediateAll'],axis=1)[RevmovePointInIntermediate:]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="r", markeredgecolor='r', linewidth=5)

    plt.plot(binning.BinsCenter_TimeDependent_Low[StartIndex:-RevmovePointInLow]                , np.average(PbarResult['RelativeSystematicErrorTimeDependentLowAll'],axis=1)[StartIndex:-RevmovePointInLow]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="blue", linewidth=5, markeredgecolor='blue', label="Time-dependent Systematic uncertainty")
    plt.plot(binning.BinsCenter_TimeDependent_Low[StartIndex:-RevmovePointInLow]                , np.average(PbarResult['RelativeSystematicErrorTimeInDependentLowAll'],axis=1)[StartIndex:-RevmovePointInLow]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="green", linewidth=5, markeredgecolor='green', label="Time-independent Systematic uncertainty")

    plt.plot(binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], np.average(PbarResult['RelativeSystematicErrorTimeDependentIntermediateAll'],axis=1)[RevmovePointInIntermediate:]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="blue", linewidth=5, markeredgecolor='blue',)
    plt.plot(binning.BinsCenter_TimeDependent_Intermediate[RevmovePointInIntermediate:], np.average(PbarResult['RelativeSystematicErrorTimeInDependentIntermediateAll'],axis=1)[RevmovePointInIntermediate:]*100, marker='o', linestyle="None", markersize=MarkersizeValue, markerfacecolor="green", linewidth=5, markeredgecolor='green',)

    axes.set_ylim(0, 20)
    plt.xlabel("|R| / (GV)"              , horizontalalignment='right', x=1.0)
    plt.ylabel("Relative uncertainty (%)", horizontalalignment='right', y=1.0)

    plt.legend()

    plt.savefig( lowworkpath + "/totalall/results/plot/" + "RelativeErrorvsRigidity_" + mode + ".pdf")






