import draw
from draw import *

def drawChi2(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, rigidityrange, binmerge, index, mode, RemovedList_3B):

    fig_mine_chi2, ax = plt.subplots()

    if mode == "6months":
        plt.errorbar(draw.datesM6, draw.M6Chi2, yerr=0, marker='o', linestyle="None", markersize=30, markerfacecolor="black", ecolor="black", linewidth=5 )
    elif mode == "6Bartels":
        plt.errorbar(draw.dates6B, draw.B6Chi2, yerr=0, marker='o', linestyle="None", markersize=30, markerfacecolor="black", ecolor="black", linewidth=5)
    elif mode == "3Bartels":
        plt.errorbar(np.delete(draw.dates, RemovedList_3B),  np.delete(draw.B3Chi2, RemovedList_3B), yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="black", ecolor="black", linewidth=5 )

    ax.set_ylim([0, 2])

    #RigidityBinText = "{:g}".format(binning.published2016binnings[index]) + " - " + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + " GV"
    #plt.text(0.60, 0.85, RigidityBinText, transform = ax.transAxes, fontsize=70)

    #plt.ylabel("chi2/ndf", horizontalalignment='right', y=1.0)
    plt.ylabel(r'$\rm{\chi^2}$/ndf', horizontalalignment='right', y=1.0)

    if rigidityrange == "low":
        if mode == "6months":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "Chi2dof_6months_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "6Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "Chi2dof_6BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
        elif mode == "3Bartels":
            plt.savefig(lowworkpath +"/totalall/results/plot/" + "Chi2dof_3BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +".pdf")
    elif rigidityrange == "intermediate":
        if mode == "6months":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Chi2dof_6MonthsRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "6Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Chi2dof_6BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
        elif mode == "3Bartels":
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Chi2dof_3BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".pdf")
    plt.close()



def drawChi2_6Bartel6MonthsCompare(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, rigidityrange, binmerge, index, mode):

    fig_mine_chi2, ax = plt.subplots(figsize=(30,18))

    plt.scatter(draw.datesM6, draw.M6Chi2, s=180, label='6 Months' )
    plt.scatter(draw.dates6B, draw.B6Chi2, s=180, label='6 Bartels')

    ax.tick_params(direction='in',length=10, width=3)
    plt.ylabel("chi2/dof",fontsize=60)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    ax.axes.set_ylim([0, 2])
    plt.legend(loc='best',fontsize=35)
    plt.title( "{:g}".format(binning.published2016binnings[index])  + "-" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ' GV', fontsize=60)

    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "6Bartel6MonthsCompare_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +"_Chi2.pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "6Bartel6MonthsCompare_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" +"{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_Chi2.pdf")

    plt.close()




