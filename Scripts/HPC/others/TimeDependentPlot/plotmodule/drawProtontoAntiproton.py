import draw
from draw import *

def drawProtontoAntiproton(rigidityrange):
    if rigidityrange == "intermediate":
        if mode == "3Bartels":
            #### Compare with Positron_to_Electron ratio.
            fig = plt.figure(figsize=(30,18))
            ax1 = fig.add_subplot(111)
            #ax1.errorbar(dates[0:15]+dates[16:], 1./np.delete(B3result,15)/10000, yerr=np.delete(B3_P_Anti_Error,15)/10000, marker='o',linestyle="None",markersize=12, markerfacecolor="b", ecolor="b")
            ax1.set_ylabel('Proton to antiproton Ratio (10$^{4}$)',fontsize=55, color="b")
            plt.xticks( rotation=40 )
            plt.xticks(fontsize=40)
            plt.yticks(fontsize=40)
            ax1.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])

            ax2 = ax1.twinx()
            positron_to_electron_error = np.sqrt( ((Positron_error_relative*Positron)/Electron)**2 + (Positron/Electron**2*(Electron_error_relative*Electron))**2)
            plt.errorbar(np.delete(dates_1Bartel,[46,47]), np.delete(Positron/Electron,[46,47]), yerr=np.delete( positron_to_electron_error ,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r") # Fix me: this result is different from PosiOverElec value.
            #plt.errorbar(np.delete(dates_1Bartel,[46,47]), np.delete(PosiOverElec,[46,47]), yerr=np.delete(PosiOverElec_error,[46,47]), marker='s', linestyle="None", markersize=12, markerfacecolor="r",ecolor="r")
            ax2.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
            ax2.set_ylabel('Positron to Electron Ratio',fontsize=55, color="r")
            major_f = md.DateFormatter('%Y-%m') #major_f = md.DateFormatter('%Y-%m-%d %H:%M:%S')
            major_L = md.MonthLocator(bymonth=[1,7])
            ax1.xaxis.set_major_locator(major_L)
            ax1.xaxis.set_major_formatter(major_f)
            plt.xticks( rotation=40 )
            plt.xticks(fontsize=40)
            plt.yticks(fontsize=40)
            plt.savefig(intermediateworkpath +"/total/results/plot/" + "Proton_Antiproton_ratio_3BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + ".png")
            plt.close()


