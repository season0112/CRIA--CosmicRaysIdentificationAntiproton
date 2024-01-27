
## Section 4.1: DataSelection
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/others/B1042_antipr.pl1.1800_7.6_all_mc_iss_ratio.pdf /Users/sichenli/Desktop/PhDThesis/Figures/chapter4/DataSelection/


## Section 4.2: CC Estimator
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/CCEstimatorplots/CorrectAndConfused_1D_Logy_259_330_binnumber_30_CCCut_0.20_CCN_20_TRDN_20_InROOT.pdf chapter4/ChargeConfusion/CCEstimator/
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/CCEstimatorplots/RejectionPower259_330_binnumber_30_CCCut_0.00_CCN_20_TRDN_20_InROOT.pdf chapter4/ChargeConfusion/RejectionPower/


## Section 4.3: TRD Likelihood
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/HighRange/TRD_ChargeCorrectProtonTemplate_MC_B1042_pr.pl1.1800_7.6_all_Tree_positive_rigidity_14.1_15.3_pattern_0_bin_100_LogY.pdf chapter4/TRDLikelihood/


## Section 4.4: Template Fit
<<'COMMENT'
# Time Averaged
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/NumberPlot_AllPatternsOverRigidityBinWidthpass78.pdf  chapter4/TemplateFit/TimeAveragedPlot/high/NumberPlot_AllPatternsOverRigidityBinWidthpass78.pdf
# Time Dependent
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/TemplateFit_TRD_6.47_7.76_78.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TRD/TemplateFit_TRD_6.47_7.76_78.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Antiproton_TRD_6.47_7.76_miu.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TRD/Antiproton/Antiproton_TRD_6.47_7.76_miu.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Antiproton_TRD6.47_7.76_sigma.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TRD/Antiproton/Antiproton_TRD6.47_7.76_sigma.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Electron_TRD_6.47_7.76_miu.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TRD/Electron/Electron_TRD_6.47_7.76_miu.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Electron_TRD6.47_7.76_sigma.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TRD/Electron/Electron_TRD6.47_7.76_sigma.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Antiproton_TOF3.64_4.43_sigma.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Antiproton/Antiproton_TOF3.64_4.43_sigma.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Antiproton_TOF_3.64_4.43_miu.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Antiproton/Antiproton_TOF_3.64_4.43_miu.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Electron_TOF3.64_4.43_sigma.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Electron/Electron_TOF3.64_4.43_sigma.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Electron_TOF_3.64_4.43_miu.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Electron/Electron_TOF_3.64_4.43_miu.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Pion_TOF3.64_4.43_sigma.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Pion/Pion_TOF3.64_4.43_sigma.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/IntermediateRange/Pion_TOF_3.64_4.43_miu.pdf chapter4/TemplateFit/TimeDependentPlot/TimeDependentTemplate/TOF/Pion/Pion_TOF_3.64_4.43_miu.pdf
#scp copy:/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v6.0/total/results/plot/AntiprotonNumber_6BartalRotation_0124_free_binmerge2_6.47_7.76.pdf chapter4/TemplateFit/TimeDependentPlot/intermediate/FitResult/AntiprotonNumber_6BartalRotation_0124_free_binmerge2_6.47_7.76.pdf

#scp copy:/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v6.0/total/results/plot/Chi2dof_6BartalRotation_0124_free_binmerge2_6.47_7.76.pdf chapter4/TemplateFit/TimeDependentPlot/intermediate/FitResult/Chi2dof_6BartalRotation_0124_free_binmerge2_6.47_7.76.pdf
COMMENT


## Section 4.5: Acceptance
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/others/Acceptance_Antiproton.pdf chapter4/Acceptance/EffectiveAcceptance/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/others/Acceptance_Proton.pdf chapter4/Acceptance/EffectiveAcceptance/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/AcceptancePlot.pdf chapter4/Acceptance/EffectiveAcceptanceRatio/


## Section 4.6 MeasuringTime 
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/others/TimeDependentMeasuringTime/TimeDependentMeasuringTime.pdf chapter4/MeasuringTime/TimeDependent/


## Section 4.7: Unfolding
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/NumberPlot_pass78.pdf chapter4/Unfolding/RawUnfoldCompare/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/ProtonNumberPlot_pass78.pdf chapter4/Unfolding/RawUnfoldCompare/
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/unfolding/plot_ratio/MM_Pattern_0.pdf chapter4/Unfolding/MM/MM_Pattern_0.pdf


## Section 4.8: AntiprotonToProtonFluxRatioCalculation
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/fullratio_StaErrOnly_NotOverlapped_ThisAnallysisOnly_pass78.pdf chapter4/AntiprotonToProtonFluxRatioCalculation/TimeAveraged/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/fullratio_StaErrOnly_ThisAnallysisOnly_pass78.pdf chapter4/AntiprotonToProtonFluxRatioCalculation/TimeAveraged/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/StaRelErr_pass78.pdf chapter4/AntiprotonToProtonFluxRatioCalculation/TimeAveraged/
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/6Bartels_binmerge2_1.92_2.4_StaErrOnly.pdf chapter4/AntiprotonToProtonFluxRatioCalculation/TimeDependent/6Bartels_binmerge2_1.92_2.4_StaErrOnly.pdf
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/StatisticalAbsoluteError_6Bartels_binmerge_2_1.92_2.4.pdf chapter4/AntiprotonToProtonFluxRatioCalculation/TimeDependent/StatisticalAbsoluteError_6Bartels_binmerge_2_1.92_2.4.pdf


## Section 4.9: Systematic Error
#Systematic Error: Acceptance
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysErr_Acceptance_LogY.pdf chapter4/SystematicError/Acceptance_SystemError/SysErr_Acceptance_LogY.pdf
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_ACC_* chapter4/SystematicError/Acceptance_SystemError/

# Systematic Error: CC Level
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/templatefit/negative/FixedCC/Plot_UncertaintyBand/CCLevel_Pattern_0_VGG16NN_525version_LogY.pdf chapter4/SystematicError/CC_SystemError/CCLevel_Pattern_0_VGG16NN_525version_LogY.pdf
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_CCLevel_LinearFit_Pattern_0_VGG16NN_525version_LinearX.pdf chapter4/SystematicError/CC_SystemError/UncertaintyBand_CCLevel_LinearFit_Pattern_0_VGG16NN_525version_LinearX.pdf
#scp copy:/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_CCLevel_Option3CurveFit* chapter4/SystematicError/CC_SystemError/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_CC_pass78_LogX.pdf chapter4/SystematicError/CC_SystemError/SysRelErr_CC_pass78_LogX.pdf
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysErr_CC_LogY_pass78.pdf chapter4/SystematicError/CC_SystemError/SysErr_CC_LogY_pass78.pdf

# Systematic Error: Template Fit Range
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_FitRange_* chapter4/SystematicError/FitRange_SystemError/

# Systematic Error: Total 
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_pass78.pdf chapter4/SystematicError/TotalSystemError/SysRelErr_pass78.pdf

# Systematic Error: TimeDependent
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/RelativeBreakDownOfTotalSysError_6Bartels_binmerge_2_1.92_2.4.pdf chapter4/SystematicError/TimeDependent/RelativeBreakDownOfTotalSysError_6Bartels_binmerge_2_1.92_2.4.pdf
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysErr_Tot_LogY_pass78.pdf chapter4/SystematicError/TotalSystemError/

#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/RelativeBreakDownOfTotalSysError_6Bartels_binmerge_2_1.92_2.4.pdf chapter4/SystematicError/TimeDependent/


## Section 5.1 Time averaged result
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/StaSysRelErrCompare_pass78.pdf chapter5/timeaveraged/StaSysRelErrCompare_pass78.pdf
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/fullratio_NOT_overlapped_PhysicsReport2021_PhyRep2021.pdf chapter5/timeaveraged/fullratio_NOT_overlapped_PhysicsReport2021_PhyRep2021.pdf   #( This needs change issversion to Physics Report)
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/fullratio_ThisAnallysisOnly_pass78*.pdf chapter5/timeaveraged/
#scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/TotalRelErr_BreakDown_pass78.pdf chapter5/timeaveraged/

## Section 5.2 Time dependent result
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/PbarRatioInAllRigidityBins6Bartels.pdf chapter5/timedependent/low
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/CompareWithPbarOverProAslamModel1.92_2.4.pdf chapter5/timedependent/model
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/RelativeErrorvsRigidity_6Bartels.pdf chapter5/timedependent/low
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/RelativeBreakDownOfTotalError_6Bartels_binmerge_2_1.92_2.4.pdf chapter5/timedependent/low
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/Antiproton_Electron_Aachen_1.92_2.4_6Bartels.pdf chapter5/timedependent/low
#scp copy:/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/results/plot/Antiproton_Electron_Aachen_2.4_2.97_6Bartels.pdf chapter5/timedependent/low
#scp copy:/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v6.0/total/results/plot/Antiproton_Electron_Aachen_4.43_5.37_6Bartels.pdf chapter5/timedependent/intermediate
#scp copy:/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v6.0/total/results/plot/Antiproton_Electron_Aachen_6.47_7.76_6Bartels.pdf chapter5/timedependent/intermediate

# test 
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_CC_pass78.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/SysRelErr_CC_pass78.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysRelErr_CC_pass78_LogX.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/SysRelErr_CC_pass78_LogX.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/SysErr_CC_LogY_pass78.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/SysErr_CC_LogY_pass78.pdf

scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/StaSysRelErrCompare_pass78.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/StaSysRelErrCompare_pass78.pdf
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/fullratio_ThisAnallysisOnly_pass78*.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/
scp copy:/home/bo791269/Software/AntiprotonAnalysis/Macros/TotalRelErr_BreakDown_pass78.pdf chapter4/SystematicError/CC_SystemError/CC_CCcut0.20_OneSigmaWithBump_Intersect1.58/FinalResult/


