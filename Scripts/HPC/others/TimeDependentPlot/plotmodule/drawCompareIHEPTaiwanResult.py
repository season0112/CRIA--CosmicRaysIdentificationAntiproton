import draw
from draw import *

def drawCompareIHEPTaiwanResult(rigidityrange, Taiwan_all, Taiwan_error_all, Ihep_all, Ihep_error_all, index, yaxisfactor, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge):

    # Plot
    fig = plt.figure(figsize=(30,18))
    grid = plt.GridSpec(4, 1, wspace=0, hspace=0)

    # Plot ratio
    ax1 = fig.add_subplot(grid[0:3,0])
    ax1.errorbar(np.delete(draw.dates_TaiwanRange, [15,36]), np.delete(np.array(Taiwan_all[index-5-int((index-9)/2)-1]), [15,36])*10000, yerr=np.delete(np.array(Taiwan_error_all[index-5-int((index-9)/2)-1]), [15,36])*10000, marker='s', linestyle="None",markersize=12, markerfacecolor="red",ecolor="red", label="Taiwan")
    ax1.errorbar(np.delete(draw.dates_IhepRange, [15,36]), np.delete(np.array(Ihep_all[index-6-int((index-9)/2)-1]),[15,36])*10000, yerr=np.delete(np.array(Ihep_error_all[index-6-int((index-9)/2)-1]),[15,36])*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="magenta",ecolor="magenta",label="IHEP")

    plt.ylabel("Antiproton ratio (10$^{-4}$)",fontsize=40)
    ax1.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    ax1.tick_params(direction='in',length=10, width=3)
    plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=40)
    plt.legend(loc='best',fontsize=35)

    # Plot Residual
    ax2 = fig.add_subplot(grid[3,0])

    taiwan_ttcson = np.delete(np.array(Taiwan_all[index-5-int((index-9)/2)-1]),[15,36]) *10000
    ihep_ttcson = np.delete( np.array(Ihep_all[index-6-int((index-9)/2)-1]),[15,36]) *10000
    taiwanerror_ttcson = np.delete(np.array(Taiwan_error_all[index-5-int((index-9)/2)-1]), [15,36])*10000
    iheperror_ttcson = np.delete(np.array(Ihep_error_all[index-6-int((index-9)/2)-1]), [15,36])*10000

    residual_Ihep_Taiwan = (np.delete(np.array(Ihep_all[index-6-int((index-9)/2)-1]),[15,36]) - np.delete(np.array(Taiwan_all[index-5-int((index-9)/2)-1]), [15,36])) / (np.delete(np.array(Taiwan_all[index-5-int((index-9)/2)-1]), [15,36]))*100
    residualerror_Ihep_Taiwan = np.sqrt( ((1./taiwan_ttcson)*iheperror_ttcson)**2 + (ihep_ttcson/taiwan_ttcson**2*taiwanerror_ttcson)**2 ) *100

    ax2.errorbar(draw.dates[0:15]+draw.dates[16:36]+draw.dates[37:38], residual_Ihep_Taiwan,  yerr = residualerror_Ihep_Taiwan, marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
    plt.ylabel(r'$\frac{\rm IHEP-\rm Taiwan}{\rm Taiwan}$'+" (%)",fontsize=40)
    # Fit in Residual
    #b = np.polyfit( np.delete(dates_unixsecond[0:othergroup.shape[0]],[15,36]), residual_Ihep_Taiwan, 0)
    #plt.hlines(b, dates_TaiwanRange[0], dates_TaiwanRange[-1], linestyles = "solid")
    plt.xticks( rotation=40 )
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    ax2.axes.set_ylim([ -28,  28])
    ax2.tick_params(direction='in',length=10, width=3, bottom=True, top=False, labelbottom=True)
    major_f = md.DateFormatter('%Y-%m')
    major_L = md.MonthLocator(bymonth=[1,7])
    ax2.xaxis.set_major_locator(major_L)
    ax2.xaxis.set_major_formatter(major_f)
    fig.subplots_adjust(hspace=0)

    # Save Plot
    if rigidityrange == "low":
        plt.savefig(lowworkpath +"/totalall/results/plot/" + "3BartalRotation_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + "IHEPandTaiwan" +".pdf")
    elif rigidityrange == "intermediate":
        plt.savefig(intermediateworkpath +"/total/results/plot/" + "3BartalRotation_" + str(trackerpattern) + "_" + str(richcut) + "_" + "binmerge" + str(binmerge) + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + "IHEPandTaiwan" +".pdf")

    plt.close()



