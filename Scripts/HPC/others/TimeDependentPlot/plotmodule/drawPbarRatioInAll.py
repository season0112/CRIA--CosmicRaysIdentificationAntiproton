import draw
from draw import *
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import PythonPlotDefaultParameters

def LoadTimeDependentResult(mode, lowworkpath, intermediateworkpath, binmerge, trackerpattern, richcut):

    if mode == "6months":
        x_data    = draw.datesM6
        nameindex = "_6months"
    elif mode == "3Bartels":
        x_data    = dates[0:15]+dates[16:38]+dates[39:]
        nameindex = "_3BartalRotation"
    elif mode == "6Bartels":
        x_data    = draw.dates6B
        nameindex = "_6BartalRotation"

    plotrange_low          = range(1, 17, 2)
    plotrange_intermediate = range(9, 29, 2)

    PbarRatioLowAll          = []
    PbarErrorLowAll          = []
    PbarRatioIntermediateAll = []
    PbarErrorIntermediateAll = []

    for index in plotrange_low:
        with open ( lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + nameindex + ".txt") as RatioLow:
            lineRatioLow = RatioLow.readlines()
        PbarRatio_Low = np.array(list(map(lambda s: float(s.strip()), lineRatioLow)))
        PbarRatioLowAll.append(PbarRatio_Low) 
        with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + nameindex + ".txt") as ErrorLow:
            lineErrorLow = ErrorLow.readlines()
        PbarRatioError_Low = np.array(list(map(lambda s: float(s.strip()), lineErrorLow)))  
        PbarErrorLowAll.append(PbarRatioError_Low) 

    for index in plotrange_intermediate:
        with open ( intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + nameindex + ".txt" ) as RatioIntermediate:
            lineRatioIntermediate = RatioIntermediate.readlines()
        PbarRatio_Intermediate = np.array(list(map(lambda s: float(s.strip()), lineRatioIntermediate)))
        PbarRatioIntermediateAll.append(PbarRatio_Intermediate)
        with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + nameindex + ".txt") as ErrorIntermediate:
            lineErrorIntermediate = ErrorIntermediate.readlines()
        PbarRatioError_Intermediate = np.array(list(map(lambda s: float(s.strip()), lineErrorIntermediate)))
        PbarErrorIntermediateAll.append(PbarRatioError_Intermediate)

    if mode == "3Bartels":
        for i in len(PbarRatioLowAll):
            PbarRatioLowAll          = np.delete(PbarRatioLowAll[i], [15,38])
            PbarErrorLowAll          = np.delete(PbarErrorLowAll[i], [15,38])
            PbarRatioIntermediateAll = np.delete(PbarRatioIntermediateAll[i], [15,38])
            PbarErrorIntermediateAll = np.delete(PbarErrorIntermediateAll[i], [15,38])
    print( "len PbarRatioLowAll: " + str(len(PbarRatioLowAll)))  #8
    print( "len PbarRatioIntermediateAll: " + str(len(PbarRatioIntermediateAll))) #10

    return x_data, PbarRatioLowAll, PbarErrorLowAll, PbarRatioIntermediateAll, PbarErrorIntermediateAll


def drawPbarRatioInAllRigidity(lowworkpath, intermediateworkpath, mode, binmerge):

    if mode == "6months":
        x_data    = draw.datesM6
        nameindex = "_6months"
    elif mode == "3Bartels":
        x_data    = draw.dates[0:15]+draw.dates[16:38]+draw.dates[39:]
        nameindex = "_3Bartels"
    elif mode == "6Bartels":
        x_data    = draw.dates6B
        nameindex = "_6Bartels"

    PbarResult = np.load(lowworkpath + "/totalall/results/TimeDependentResult" + nameindex + ".npz")

    #showindex_low          = range(1, 17, 2) # 1, 3,  5, 7, 9, 11, 13, 15
    #showindex_intermediate = range(9, 29, 2) # 9, 11, 13, 15, 17, 19, 21, 23, 25, 27
    #showindex_low          = range(3, 7, 2)   #3 , 5
    #showindex_intermediate = range(13, 17, 2) #13, 15

    showindex_low          = [5 ,  7, 9]
    showindex_intermediate = [13, 21]

    #showindex_low          = [5, 9]
    #showindex_intermediate = [13, 21]


    showindex_all          = list(showindex_low) + list(showindex_intermediate)

    NormalizedFactor = 100000

    yaxisfactor = 0.16


    #### Plot Absoute Pbar Ratio
    fig, axes = plt.subplots(len(showindex_all), 1, figsize=(20,40), sharex=True)

    for i, row in enumerate(axes[0:len(showindex_low)]):

        row.errorbar(x_data, PbarResult['PbarRatioLowAll'][int((showindex_low[i]-1)/2),:] * NormalizedFactor, PbarResult['PbarTotErrorTimeDependentLowAll'][int((showindex_low[i]-1)/2),:] * NormalizedFactor, marker='o', linestyle="None", markersize=20, markerfacecolor="b", ecolor="b", linewidth=5)

        RigidityBinText = "{:g}".format(binning.published2016binnings[showindex_low[i]]) + " - " + "{:g}".format(binning.published2016binnings[showindex_low[i]+int(binmerge)]) + " GV"
        plt.text(0.65, 0.75, RigidityBinText, transform = row.transAxes, fontsize=35)

        plt.setp(row.get_yticklabels(), fontsize=50)   # this control the y label size in low range.
        row.axes.set_ylim( [row.get_ylim()[0]*(1-yaxisfactor), row.get_ylim()[1]*(1+yaxisfactor)])
        row.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

        if i == 0:
            row.set_yticks([2, 3, 4])
        if i == 1:
            row.set_yticks([3, 4, 5, 6])
        if i == 2:
            row.set_yticks([4, 5, 6, 7])

    for i, row in enumerate(axes [len(showindex_low):len(showindex_all)]):

        row.errorbar(x_data, PbarResult['PbarRatioIntermediateAll'][int((showindex_intermediate[i]-9)/2),:] * NormalizedFactor,  PbarResult['PbarTotErrorTimeDependentIntermediateAll'][int((showindex_intermediate[i]-9)/2),:] * NormalizedFactor, marker='o', linestyle="None", markersize=20, markerfacecolor="b", ecolor="b", linewidth=5)

        RigidityBinText = "{:g}".format(binning.published2016binnings[showindex_intermediate[i]]) + " - " + "{:g}".format(binning.published2016binnings[showindex_intermediate[i]+int(binmerge)]) + " GV"
        plt.text(0.65, 0.75, RigidityBinText, transform = row.transAxes, fontsize=35)

        row.axes.set_ylim( [row.get_ylim()[0]*(1-yaxisfactor), row.get_ylim()[1]*(1+yaxisfactor)])
        #row.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
        plt.setp(row.get_yticklabels(), fontsize=50)   # this control the y label size in intermediate range.

        if i == 0:
            row.set_yticks([7, 8, 9, 10, 11, 12])
        if i == 1:
            row.set_yticks([13, 15, 17, 19, 21])

    plt.xticks(fontsize=40)
    plt.subplots_adjust(wspace=0, hspace=0)

    #plt.text(-0.13, 4.0, r'Antiproton to proton flux ratio ($\times$ $10^{-5}$)', transform = row.transAxes, rotation=90, fontsize=50) ## For 4 bins display
    plt.text(-0.16, 3.0, r'$\rm{\overline{p}/p}$ ($\times$ $10^{5}$)', transform = row.transAxes, rotation=90, fontsize=60)

    plt.savefig( lowworkpath + "/totalall/results/plot/" + "PbarRatioInAllRigidityBins" + mode + ".pdf")


    #### Plot Relatvie Pbar Ratio (Not active)
    fig, ax = plt.subplots(figsize=(50,18))

    for i in range(len(showindex_low)):
        RigidityBinText = "{:g}".format(binning.published2016binnings[showindex_low[i]]) + " - " + "{:g}".format(binning.published2016binnings[showindex_low[i]+int(binmerge)]) + " GV"
        mean = np.mean(PbarResult['PbarRatioLowAll'][int((showindex_low[i]-1)/2),:])
        ax.errorbar(x_data, PbarResult['PbarRatioLowAll'][int((showindex_low[i]-1)/2),:]/mean, PbarResult['PbarTotErrorLowAll'][int((showindex_low[i]-1)/2),:]/mean, marker='o', linestyle="None", markersize=40, linewidth=3, label=RigidityBinText)
    for i in range(len(showindex_intermediate)):
        RigidityBinText = "{:g}".format(binning.published2016binnings[showindex_intermediate[i]]) + " - " + "{:g}".format(binning.published2016binnings[showindex_intermediate[i]+int(binmerge)]) + " GV"
        mean = np.mean(PbarResult['PbarRatioIntermediateAll'][int((showindex_intermediate[i]-9)/2),:])
        ax.errorbar(x_data, PbarResult['PbarRatioIntermediateAll'][int((showindex_intermediate[i]-9)/2),:]/mean, PbarResult['PbarTotErrorIntermediateAll'][int((showindex_intermediate[i]-9)/2),:]/mean, marker='o', linestyle="None", markersize=40, linewidth=3, label=RigidityBinText)

    ax.tick_params(axis='x', labeltop=False, labelbottom=True , labelleft=False, labelright=False, direction='in', length=10, width=5)
    ax.tick_params(axis='y', labeltop=False, labelbottom=False, labelleft=True , labelright=False, direction='in', length=10, width=5)
    #ax.axes.set_ylim( [row.get_ylim()[0]*(1-yaxisfactor), row.get_ylim()[1]*(1+yaxisfactor)])
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    plt.ylabel("antiproton to proton relative flux ratio",fontsize=50)
    plt.legend(loc='best',fontsize=50)
    #plt.text(-0.1, 4.0, r'antiproton to proton relative flux ratio', transform = ax.transAxes, rotation=90, fontsize=27)

    plt.savefig( lowworkpath + "/totalall/results/plot/" + "PbarRelativeRatioInAllRigidityBins" + mode + ".pdf")



