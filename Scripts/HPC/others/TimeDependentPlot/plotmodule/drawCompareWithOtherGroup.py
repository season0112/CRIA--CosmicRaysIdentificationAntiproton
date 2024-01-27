import draw
from draw import *

def drawCompareWithOtherGroup(lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, Taiwan_all, Ihep_all, Taiwan_all_6, Taiwan_error_all, Taiwan_error_all_6, Ihep_error_all, index, yaxisfactor, mode, rigidityrange):
    
    if mode == "6Bartels":
        allothergroup = [np.array(Taiwan_all_6[index-5-int((index-9)/2)-1-1])] #Only make sense from index 3:2.97-3.64.
        allothergroup_error = [np.array(Taiwan_error_all_6[index-5-int((index-9)/2)-1-1]) ]
        dates_allothergroup = [draw.dates6B_TaiwanRange]
        allgroupname = ["Taiwan"]
        allcolors = ["red"]
    elif mode == "3Bartels":
        allothergroup = [np.array(Taiwan_all[index-5-int((index-9)/2)-1]), np.array(Ihep_all[index-6-int((index-9)/2)-1])] # For Taiwan result,only make sense above 2.97GV.
        allothergroup_error = [np.array(Taiwan_error_all[index-5-int((index-9)/2)-1]), np.array(Ihep_error_all[index-6-int((index-9)/2)-1]) ]
        dates_allothergroup = [draw.dates_TaiwanRange, draw.dates_IhepRange]
        allgroupname = ["Taiwan","IHEP"]
        allcolors = ["red","magenta"]

    for othergroup, othergroup_error, dates_othergroup, groupname, color in zip (allothergroup, allothergroup_error, dates_allothergroup, allgroupname, allcolors):

    # Plot1: Plot Ratio
        fig = plt.figure(figsize=(30,18))
        grid = plt.GridSpec(4, 1, wspace=0, hspace=0)
        ax1 = fig.add_subplot(grid[0:3,0])
        if mode == "6Bartels":
            ax1.errorbar(dates_othergroup, othergroup*10000, yerr=othergroup_error*10000, marker='s', linestyle="None", markersize=12, markerfacecolor=color, ecolor=color, label=groupname)
            ax1.errorbar(draw.dates6B, draw.B6result*10000, yerr=draw.Error6b*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue",label="Aachen")
        elif mode == "3Bartels":
            ax1.errorbar(np.delete(dates_othergroup, [15,36]), np.delete(othergroup, [15,36])*10000, yerr=np.delete(othergroup_error, [15,36])*10000, marker='s', linestyle="None",markersize=12, markerfacecolor=color,ecolor=color, label=groupname)
            ax1.errorbar(draw.dates[0:15]+draw.dates[16:], np.delete(draw.B3result,15)*10000, yerr=np.delete(draw.Error3b,15)*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue",label="Aachen")
        plt.ylabel("Antiproton ratio (10$^{-4}$)",fontsize=40)
        ax1.axes.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
        ax1.tick_params(direction='in',length=10, width=3)
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=40)
        plt.legend(loc='best',fontsize=35)

    # Plot Residual of ratio
        ax2 = fig.add_subplot(grid[3,0])
        if mode == "6Bartels":
            othergroup = othergroup*10000
            othergroup_error = othergroup_error*10000
            aachen = draw.B6result[0:othergroup.shape[0]]*10000
            aachenerror = draw.Error6b[0:othergroup.shape[0]]*10000
        elif mode == "3Bartels":
            aachen = np.delete(draw.B3result[0:othergroup.shape[0]],[15,36]) *10000
            othergroup_ttcson = np.delete(othergroup, [15,36]) *10000
            othergroup_error_ttcson = np.delete(othergroup_error, [15,36]) *10000
            aachenerror = np.delete(draw.Error3b[0:othergroup.shape[0]],[15,36]) *10000
        # Choice1: if residual is (Aachen-Other)/Aachen
        #if mode == "3Bartels":
        #residual = np.append( (np.delete(B3result[0:othergroup.shape[0]],15)-np.delete(othergroup, 15))/(np.delete(B3result[0:othergroup.shape[0]],15))*100, np.ones(B3result.shape[0]-othergroup.shape[0])*999 )
        #residualerror = np.sqrt( (othergroup_ttcson/aachen**2*aachenerror)**2 + (othergroup_error_ttcson/aachen)**2 ) *100 #residualerror = othergroup_ttcson/aachen**2*aachenerror *100 ?
        #ax2.errorbar(dates[0:15]+dates[16:], residual,  yerr=np.append(residualerror, np.ones(B3result.shape[0]-othergroup.shape[0])*1), marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
        # Choice2: if residual is (Aachen-Other)/Other
        if mode == "6Bartels":
            residual = np.append( (draw.B6result[0:othergroup.shape[0]]*10000-othergroup)/othergroup*100, np.ones(draw.B6result.shape[0]-othergroup.shape[0])*99999 )
            residualerror = np.append ( np.sqrt( ((1./othergroup)*aachenerror)**2 +  (aachen/othergroup**2*othergroup_error)**2 ) *100, np.ones(draw.B6result.shape[0]-othergroup.shape[0])*0 )
            ax2.errorbar ( draw.dates6B[0:residual.shape[0]], residual, yerr=residualerror, marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
        elif mode == "3Bartels":
            residual = np.append( (np.delete(draw.B3result[0:othergroup.shape[0]],[15,36])-np.delete(othergroup, [15,36]))/(np.delete(othergroup, [15,36]))*100, np.ones(draw.B3result.shape[0]-othergroup.shape[0])*99999 )
            residualerror = np.append ( np.sqrt( ((1./othergroup_ttcson)*aachenerror)**2 + (aachen/othergroup_ttcson**2*othergroup_error_ttcson)**2 ) *100, np.ones(draw.B3result.shape[0]-othergroup.shape[0])*0 )
            ax2.errorbar (draw.dates[0:15]+draw.dates[16:36]+draw.dates[37:42], residual,  yerr=residualerror, marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
        plt.ylabel(r'$\frac{\rm Aachen-\rm Taiwan}{\rm Taiwan}$'+" (%)",fontsize=40)
        # Fit
        if mode == "6Bartels":
            print("Fit in resudual: in progress")
        elif mode == "3Bartels":
            b = np.polyfit( np.delete(draw.dates_unixsecond[0:othergroup.shape[0]],[15,36]), residual[0:(othergroup.shape[0]-1-1)], 0)
            plt.hlines(b, dates_othergroup[0], dates_othergroup[-1], linestyles = "solid")
            #print("Now is:" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]))
            #print("average difference is:" + str(b) + "%")
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=25)
        ax2.axes.set_ylim([ -40,  40])
        fig.subplots_adjust(hspace=0)

        # Save Plot
        if rigidityrange == "low":
            if mode == "6Bartels":
                plt.savefig(lowworkpath +"/totalall/results/plot/" + "6BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
            elif mode == "3Bartels":
                plt.savefig(lowworkpath +"/totalall/results/plot/" + "3BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
        elif rigidityrange == "intermediate":
            if mode == "6Bartels":
                plt.savefig(intermediateworkpath +"/total/results/plot/" + "6BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
            elif mode == "3Bartels":
                plt.savefig(intermediateworkpath +"/total/results/plot/" + "3BartalRotation_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")

        plt.close()



    # Plot2: Plot Relative Ratio
        fig_rel = plt.figure(figsize=(30,18))
        grid = plt.GridSpec(4, 1, wspace=0, hspace=0)
        ax_rel = fig_rel.add_subplot(grid[0:3,0])

        if mode == "6Bartels":
            othergroup_relative_6B = othergroup/np.mean(othergroup)
            othergroup_relativeError_6B = othergroup_error/np.mean(othergroup)
            Aachen_relative_6B = draw.B6result/np.mean(draw.B6result)
            Aachen_relativeError_6B = draw.Error6b/np.mean(draw.B6result)

            ax_rel.errorbar(dates_othergroup, othergroup_relative_6B, yerr=othergroup_relativeError_6B, marker='s', linestyle="None", markersize=12, markerfacecolor=color, ecolor=color, label=groupname)
            ax_rel.errorbar(draw.dates6B, Aachen_relative_6B, yerr=Aachen_relativeError_6B, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue",label="Aachen")

        elif mode == "3Bartels":
            othergroup_relative_3B = np.delete(othergroup, [15,36])*10000 / np.mean(np.delete(othergroup, [15,36])*10000)
            othergroup_relativeError_3B = np.delete(othergroup_error, [15,36])*10000 / np.mean(np.delete(othergroup, [15,36])*10000)
            Aachen_relative_3B = np.delete(draw.B3result, [15,36])*10000 / np.mean(np.delete(draw.B3result, [15,36])*10000)
            Aachen_relativeError_3B = np.delete(draw.Error3b, [15,36])*10000 / np.mean(np.delete(draw.B3result, [15,36])*10000)

            ax_rel.errorbar(np.delete(dates_othergroup, [15,36]), othergroup_relative_3B, yerr=othergroup_relativeError_3B, marker='s', linestyle="None",markersize=12, markerfacecolor=color,ecolor=color, label=groupname)
            ax_rel.errorbar(draw.dates[0:15]+draw.dates[16:36]+draw.dates[37:42], Aachen_relative_3B, yerr=Aachen_relativeError_3B, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue",label="Aachen")

        plt.ylabel("Antiproton relative ratio (10$^{-4}$)",fontsize=40)
        ax_rel.tick_params(direction='in',length=10, width=3)
        plt.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=True)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=40)
        plt.legend(loc='best',fontsize=35)

    # Plot Residual of relative ratio
        ax_rel2 = fig_rel.add_subplot(grid[3,0])
        if mode == "6Bartels":
            othergroup = othergroup*10000
            othergroup_error = othergroup_error*10000
            aachen = draw.B6result[0:othergroup.shape[0]]*10000
            aachenerror = draw.Error6b[0:othergroup.shape[0]]*10000
        elif mode == "3Bartels":
            aachen = np.delete(draw.B3result[0:othergroup.shape[0]],[15,36]) *10000
            othergroup_ttcson = np.delete(othergroup, [15,36]) *10000
            othergroup_error_ttcson = np.delete(othergroup_error, [15,36]) *10000
            aachenerror = np.delete(draw.Error3b[0:othergroup.shape[0]],[15,36]) *10000

        ''' 
        othergroup_relative_6B = othergroup*10000/np.mean(othergroup*10000)
        othergroup_relativeError_6B = othergroup_error*10000/np.mean(othergroup_error*10000)
        Aachen_relative_6B = draw.B6result*10000/np.mean(draw.B6result*10000)
        Aachen_relativeError_6B = draw.Error6b*10000/np.mean(draw.Error6b)

        othergroup_relative_3B = np.delete(othergroup, [15,36])*10000/np.mean(np.delete(othergroup, [15,36])*10000)
        othergroup_relativeError_3B = np.delete(othergroup_error, [15,36])*10000/np.mean(np.delete(othergroup_error, [15,36])*10000)
        Aachen_relative_3B = np.delete(draw.B3result,15)*10000/np.mean(np.delete(draw.B3result,15)*10000)
        Aachen_relativeError_3B = np.delete(draw.Error3b,15)*10000/np.mean(np.delete(draw.Error3b,15)*10000)
        '''

        # Choice1: if residual is (Aachen-Other)/Aachen
        # Choice2: if residual is (Aachen-Other)/Other
        if mode == "6Bartels":
            residual = np.append( (Aachen_relative_6B[0:othergroup_relative_6B.shape[0]]-othergroup_relative_6B)/othergroup_relative_6B*100, np.ones(Aachen_relative_6B.shape[0]-othergroup_relative_6B.shape[0])*99999 )
            residualerror = 0
            ax_rel2.errorbar ( draw.dates6B[0:residual.shape[0]], residual, yerr = residualerror, marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
        elif mode == "3Bartels":
            residual = np.append( (Aachen_relative_3B[0:othergroup_relative_3B.shape[0]]-othergroup_relative_3B)/othergroup_relative_3B*100, np.ones(Aachen_relative_3B.shape[0]-othergroup_relative_3B.shape[0])*99999 )
            residualerror = 0
            ax_rel2.errorbar (draw.dates[0:15]+draw.dates[16:36]+draw.dates[37:42], residual,  yerr = residualerror, marker='s',linestyle="None",markersize=16, markerfacecolor="k", ecolor="k" )
        plt.ylabel(r'$\frac{\rm Aachen-\rm Taiwan}{\rm Taiwan}$'+" (%)",fontsize=40)

        plt.xticks(fontsize=40)
        plt.yticks(fontsize=25)
        ax_rel2.axes.set_ylim([ -40,  40])
        fig_rel.subplots_adjust(hspace=0)


        # Save Plot
        if rigidityrange == "low":
            if mode == "6Bartels":
                plt.savefig(lowworkpath +"/totalall/results/plot/" + "6BartalRotation_Relative_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
            elif mode == "3Bartels":
                plt.savefig(lowworkpath +"/totalall/results/plot/" + "3BartalRotation_Relative_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
        elif rigidityrange == "intermediate":
            if mode == "6Bartels":
                plt.savefig(intermediateworkpath +"/total/results/plot/" + "6BartalRotation_Relative_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")
            elif mode == "3Bartels":
                plt.savefig(intermediateworkpath +"/total/results/plot/" + "3BartalRotation_Relative_" + trackerpattern + "_" + richcut + "_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_" + groupname +".pdf")

        plt.close()


