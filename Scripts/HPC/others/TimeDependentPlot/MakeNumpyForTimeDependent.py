import draw
from draw import *
import load

def SaveTimeDependentResult(mode, lowworkpath, intermediateworkpath, binmerge, trackerpattern, richcut, SplitTotal, mergestep, modename, effective_acceptance_correction_low, effective_acceptance_correction_intermediate):

    if mode == "6months":
        x_data = draw.datesM6
        nameindex = "_6months"
    elif mode == "3Bartels":
        x_data = draw.dates[0:15]+draw.dates[16:38]+draw.dates[39:]
        nameindex = "_3BartalRotation"
    elif mode == "6Bartels":
        x_data = draw.dates6B
        nameindex = "_6BartalRotation"

    plotrange_low          = range(1, 17, 2)
    plotrange_intermediate = range(9, 29, 2)

    PbarRatioLowAll                = []
    PbarStaErrorLowAll             = []
    PbarACCSysErrorLowAll          = []
    PbarFitSysErrorLowAll          = []
    PbarTotSysErrorLowAll          = []
    PbarTotErrorLowAll             = []

    PbarRatioIntermediateAll       = []
    PbarStaErrorIntermediateAll    = []
    PbarACCSysErrorIntermediateAll = []
    PbarFitSysErrorIntermediateAll = []
    PbarTotSysErrorIntermediateAll = []
    PbarTotErrorIntermediateAll    = []


    #### Load Pbar ratio and Statistical Error
    for index in plotrange_low:
        with open ( lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + nameindex + ".txt") as RatioLow:
            lineRatioLow = RatioLow.readlines()
        PbarRatio_Low = np.array(list(map(lambda s: float(s.strip()), lineRatioLow)))
        PbarRatioLowAll.append(PbarRatio_Low)
        with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + nameindex + ".txt") as ErrorLow:
            lineErrorLow = ErrorLow.readlines()
        PbarRatioError_Low = np.array(list(map(lambda s: float(s.strip()), lineErrorLow)))
        PbarStaErrorLowAll.append(PbarRatioError_Low)

    for index in plotrange_intermediate:
        with open ( intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + nameindex + ".txt" ) as RatioIntermediate:
            lineRatioIntermediate = RatioIntermediate.readlines()
        PbarRatio_Intermediate = np.array(list(map(lambda s: float(s.strip()), lineRatioIntermediate)))
        PbarRatioIntermediateAll.append(PbarRatio_Intermediate)
        with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + nameindex + ".txt") as ErrorIntermediate:
            lineErrorIntermediate = ErrorIntermediate.readlines()
        PbarRatioError_Intermediate = np.array(list(map(lambda s: float(s.strip()), lineErrorIntermediate)))
        PbarStaErrorIntermediateAll.append(PbarRatioError_Intermediate)


    #### Effective Acceptance Correction
    for index in plotrange_low: 
        PbarRatioLowAll[int((index-1)/2)]          = PbarRatioLowAll[int((index-1)/2)]          * (effective_acceptance_correction_low[index+1] + effective_acceptance_correction_low[index+1+(int(binmerge)-1)])/2
    for index in plotrange_intermediate:
        PbarRatioIntermediateAll[int((index-9)/2)] = PbarRatioIntermediateAll[int((index-9)/2)] * (effective_acceptance_correction_low[index+1] + effective_acceptance_correction_low[index+1+(int(binmerge)-1)])/2


    #### Load Systematic Error due to Acceptance
    SystematicRelativeErrorAccpetance = load.LoadSystematicErrorFromAcceptance(2)
    '''
    time_averaged_ratio_low         , time_averaged_ratio_with_effective_low         , time_averaged_error_low         , time_averaged_proton_low         , time_averaged_antiproton_low          = load.LoadTimeAveragedResult("low"         , lowworkpath, intermediateworkpath, "0124", "free", 2) # need to / 100000.
    time_averaged_ratio_intermediate, time_averaged_ratio_with_effective_intermediate, time_averaged_error_intermediate, time_averaged_proton_intermediate, time_averaged_antiproton_intermediate = load.LoadTimeAveragedResult("intermediate", lowworkpath, intermediateworkpath, "0124", "free", 2) # need to / 100000.    
    SysErrorAccTimeIndependentAll_low          = time_averaged_ratio_with_effective_low/100000          * SystematicRelativeErrorAccpetance[0:8]
    SysErrorAccTimeIndependentAll_intermediate = time_averaged_ratio_with_effective_intermediate/100000 * SystematicRelativeErrorAccpetance[4:14]
    '''
    '''
    print("shape:" + str(len(SystematicRelativeErrorAccpetance))) # 29 (already binmerge=2, just because include high rigidity range, so go up to 29, original is around 60 points)
    '''
    for i in range(len(PbarRatioLowAll)):          # 0-8
        PbarACCSysErrorLowAll         .append( SystematicRelativeErrorAccpetance[i] * PbarRatioLowAll[i] )
    for i in range(len(PbarRatioIntermediateAll)): # 0-10
        PbarACCSysErrorIntermediateAll.append( SystematicRelativeErrorAccpetance[i+4] * PbarRatioIntermediateAll[i] )


    #### Load Systematic Error due to Fit Range
    for timeindex in range(0, SplitTotal, mergestep):
        PbarFitSysErrorLow          = []
        PbarFitSysErrorIntermediate = []
        pointnumber_low          = 8
        pointnumber_intermediate = 10

        for rigidityindex in range(pointnumber_low):
            if timeindex == 123 and mode == "3Bartels":
                PbarFitSysErrorLow.append(10000000)
            else:
                rootfile = uproot.open(lowworkpath + "/totalall/SysError_Eff/SysErr_SigEff_" + modename + '_' + str(timeindex) + ".root")
                PbarFitSysErrorLow.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
        for rigidityindex in range(pointnumber_intermediate):
            if timeindex == 123 and mode == "3Bartels":
                PbarFitSysErrorIntermediate.append(10000000)
            else:
                rootfile = uproot.open(intermediateworkpath + "/total/SysError_Eff/SysErr_SigEff_" + modename + '_' + str(timeindex) + ".root")
                PbarFitSysErrorIntermediate.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])

        PbarFitSysErrorLowAll.append(PbarFitSysErrorLow)
        PbarFitSysErrorIntermediateAll.append(PbarFitSysErrorIntermediate)

    #### Convert from list to numpy
    PbarFitSysErrorLowAll          = np.array(PbarFitSysErrorLowAll)/100000
    PbarFitSysErrorIntermediateAll = np.array(PbarFitSysErrorIntermediateAll)/100000
    PbarFitSysErrorLowAll          = np.transpose(PbarFitSysErrorLowAll)
    PbarFitSysErrorIntermediateAll = np.transpose(PbarFitSysErrorIntermediateAll)

    PbarACCSysErrorLowAll = np.array(PbarACCSysErrorLowAll)
    PbarACCSysErrorIntermediateAll = np.array(PbarACCSysErrorIntermediateAll)

    PbarStaErrorLowAll          = np.array(PbarStaErrorLowAll)
    PbarStaErrorIntermediateAll = np.array(PbarStaErrorIntermediateAll)

    #### Calculat total error
    ## Time-dependent error
    PbarTotErrorTimeDependentLowAll          = np.sqrt( PbarStaErrorLowAll**2          + PbarFitSysErrorLowAll**2 )
    PbarTotErrorTimeDependentIntermediateAll = np.sqrt( PbarStaErrorIntermediateAll**2 + PbarFitSysErrorIntermediateAll**2)
    ## Time-dependent and Time-independent error
    PbarTotSysErrorLowAll          = np.sqrt( PbarACCSysErrorLowAll**2          + PbarFitSysErrorLowAll**2          )
    PbarTotSysErrorIntermediateAll = np.sqrt( PbarACCSysErrorIntermediateAll**2 + PbarFitSysErrorIntermediateAll**2 ) 
    PbarTotErrorLowAll             = np.sqrt( PbarStaErrorLowAll**2          + PbarTotSysErrorLowAll**2 )
    PbarTotErrorIntermediateAll    = np.sqrt( PbarStaErrorIntermediateAll**2 + PbarTotSysErrorIntermediateAll**2)

    #### Save in Numpy 
    np.savez(lowworkpath + "/totalall/results/" + "TimeDependentResult_" + mode + ".npz", x_data=x_data, 
        PbarRatioLowAll       = PbarRatioLowAll      ,   PbarRatioIntermediateAll       = PbarRatioIntermediateAll      ,

        PbarStaErrorLowAll    = PbarStaErrorLowAll   ,   PbarStaErrorIntermediateAll    = PbarStaErrorIntermediateAll   ,

        PbarACCSysErrorLowAll = PbarACCSysErrorLowAll,   PbarACCSysErrorIntermediateAll = PbarACCSysErrorIntermediateAll,
        PbarFitSysErrorLowAll = PbarFitSysErrorLowAll,   PbarFitSysErrorIntermediateAll = PbarFitSysErrorIntermediateAll, 
        PbarTotSysErrorLowAll = PbarTotSysErrorLowAll,   PbarTotSysErrorIntermediateAll = PbarTotSysErrorIntermediateAll,

        PbarTotErrorTimeDependentLowAll = PbarTotErrorTimeDependentLowAll, PbarTotErrorTimeDependentIntermediateAll = PbarTotErrorTimeDependentIntermediateAll,
        PbarTotErrorLowAll              = PbarTotErrorLowAll             , PbarTotErrorIntermediateAll              = PbarTotErrorIntermediateAll   ,

        RelativeStatisticalErrorLowAll = PbarStaErrorLowAll/PbarRatioLowAll   , RelativeStatisticalErrorIntermediateAll = PbarStaErrorIntermediateAll/PbarRatioIntermediateAll,
        RelativeSystematicErrorTimeDependentLowAll   = PbarFitSysErrorLowAll/PbarRatioLowAll, RelativeSystematicErrorTimeDependentIntermediateAll   = PbarFitSysErrorIntermediateAll/PbarRatioIntermediateAll,
        RelativeSystematicErrorTimeInDependentLowAll = PbarACCSysErrorLowAll/PbarRatioLowAll, RelativeSystematicErrorTimeInDependentIntermediateAll = PbarACCSysErrorIntermediateAll/PbarRatioIntermediateAll, 
        RelativeSystematicErrorLowAll  = PbarTotSysErrorLowAll/PbarRatioLowAll, RelativeSystematicErrorIntermediateAll  = PbarTotSysErrorIntermediateAll/PbarRatioIntermediateAll, 
    )






    #if mode == "6Bartels":
        #print("PbarRatioLowAll:"       + str(PbarRatioLowAll      ) )
        #print("PbarStaErrorLowAll:"    + str(PbarStaErrorLowAll   ) )
        #print("PbarTotSysErrorLowAll:" + str(PbarTotSysErrorLowAll) )
        #print('\n')
        #print("PbarRatioIntermediateAll:"       + str(PbarRatioIntermediateAll      ) )
        #print("PbarStaErrorIntermediateAll:"    + str(PbarStaErrorIntermediateAll   ) )
        #print("PbarTotSysErrorIntermediateAll:" + str(PbarTotSysErrorIntermediateAll) )      







  


