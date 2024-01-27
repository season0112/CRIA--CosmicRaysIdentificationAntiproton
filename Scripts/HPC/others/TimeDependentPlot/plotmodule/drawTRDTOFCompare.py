import draw
from draw import *



def drawTRDTOFCompare(intermediateworkpath, lowworkpath, yaxisfactor, binmerge, index):
    ## "4.43_5.37", "3.64_4.43" tof trd compare

    #1. Load TRD Effective acceptance
    EffectiveAcceptance_antipro = TFile(intermediateworkpath + "/total" + "/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root")
    EffectiveAcceptance_proton  = TFile(intermediateworkpath + "/total" + "/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root")
    Acceptance_antiproton = EffectiveAcceptance_antipro.Get("QualityCuts/effectiveAcceptanceAfterAllCuts") ## 62 points in total,0.8-1130, GetX()[0]=0.9, GetX()[61]=976. fitst and last points are 0, since MC generated R range.
    Acceptance_proton = EffectiveAcceptance_proton.Get("QualityCuts/effectiveAcceptanceAfterAllCuts")
    effective_acceptance_correction=np.array([])
    effective_acceptance_correction = np.append(effective_acceptance_correction, 0.0) ## first point gives 0.
    for i in range(1,Acceptance_antiproton.GetN()-1):  ## from index 1 and 60, index 1 and 61 are manually given as 0.
        effective_acceptance_correction = np.append(effective_acceptance_correction, Acceptance_proton.GetY()[i]/Acceptance_antiproton.GetY()[i])
    effective_acceptance_correction = np.append(effective_acceptance_correction, 0.0) ## last point gives 0.

    #2. Load TRD result
    with open (intermediateworkpath + "/total/" + "results/fit_results_" + "0124" + "_" + "free" + "_" + "4.43_5.37" + "_3BartalRotation.txt") as f3b:
        lines3b=f3b.readlines()
    B3result_trd = np.array(list(map(lambda s: float(s.strip()), lines3b)))
    B3result_trd = B3result_trd * (effective_acceptance_correction[14]+effective_acceptance_correction[15])/2  #14:4.43-4.88 15:4.88-5.37
    with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + "0124" + "_" + "free" + "_" + "4.43_5.37" + "_3BartalRotation.txt") as e3b:
        error3b_trd = e3b.readlines()
    Error3b_trd = np.array(list(map(lambda s: float(s.strip()), error3b_trd)))

    #3. Plot
    dates_unixsecond_trd = binning.Bartals3Unixtime[0:B3result_trd.shape[0]]
    dates_trd = [datetime.fromtimestamp(i) for i in dates_unixsecond_trd]
    fig_toftrdcom, ax_toftrdcom = plt.subplots(figsize=(30,18))
    #plt.errorbar(draw.dates[0:15]+draw.dates[16:], np.delete(draw.B3result_all[5,:],15)*10000, yerr=np.delete(draw.Error3b_all[5,:],15)*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="red",ecolor="red",label= "4.43_5.37GV (TOF Fit)" )
    plt.errorbar(dates_trd[0:15]+dates_trd[16:], np.delete(B3result_trd,15)*10000, yerr=np.delete(Error3b_trd,15)*10000, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue",label= "4.43_5.37GV (TRD Fit)" )
    plt.legend(loc='best',fontsize=35)
    ax_toftrdcom.set_ylim([ plt.ylim()[0]-((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor,  plt.ylim()[1]+((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor ])
    plt.xticks(fontsize=40)
    plt.yticks(fontsize=40)
    plt.savefig(lowworkpath +"/totalall/results/plot/" + "3BartalRotation_" + "binmerge" + binmerge + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) +"_tof_trd_compare.png")
