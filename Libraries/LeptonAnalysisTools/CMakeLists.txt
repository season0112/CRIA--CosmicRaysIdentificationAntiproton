set(source_files
  AnalysisSettings.C
  EcalBDTEfficiencyParametrization.C
  EcalChiSquareLateralCalibration.C
  EfficiencyTemplateFits.C
  ElectronModelTwoPowerLawsWithSolarModulation.C
  ElectronAndPositronFluxModel.C
  ElectronPositronModelWithSmoothDoubleBreak.C
  ElectronPositronDoublePowerLawApproximation.C
  FluxModel.C
  GraphTools.C
  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff.C
  SolarFlares.C
  SimplePowerLawModel.C
  SpectralIndex.C
  SystematicUncertainties.C
  MeasuringTimeParametrization.C
  ModulatedPowerLawModel.C
  TemplateFitTools.C
  TemplateFit2DFitFunction.C
  TrdLikelihoodRatioCalibration.C
  MultiTrackChargeConfusionStudyTree.C
  RooFitTools.C
)

set(header_files
  AveragingTools.hh
  AnalysisSettings.hh
  CanvasTools.hh
  EcalBDTEfficiencyParametrization.hh
  EcalChiSquareLateralCalibration.hh
  ChargeConfusionTools.hh
  ElectronModelTwoPowerLawsWithSolarModulation.hh
  ElectronAndPositronFluxModel.hh
  ElectronPositronModelWithSmoothDoubleBreak.hh
  ElectronPositronDoublePowerLawApproximation.hh
  FluxModel.hh
  MultiTrackChargeConfusionStudyTree.hh
  SolarFlares.hh
  SimplePowerLawModel.hh
  SpectralIndex.hh
  SystematicUncertainties.hh
  MeasuringTimeParametrization.hh
  ModulatedPowerLawModel.hh
  PositronRatioTools.hh
  PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff.hh
  RooFitTools.hh
  SplineFitTools.hh
  TemplateFitTools.hh
  TemplateFit2DFitFunction.hh
  TimeTools.hh
  TrdLikelihoodRatioCalibration.hh
  TriggerEfficiency.hh
)

set(header_files_with_classdef
  EfficiencyTemplateFits.hh
)

# Work-around for gbatch support, see ACANA-24.
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
endif()

# These additional ROOT libraries are required.
set(ADDITIONAL_ROOT_LIBS TreePlayer RooFitCore)

add_project_library(LeptonAnalysisTools SOURCES ${source_files} HEADERS ${header_files} HEADERS_WITH_CLASSDEF ${header_files_with_classdef} LINKDEF Linkdef.h)
target_link_libraries(LeptonAnalysisTools ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} ${ADDITIONAL_ROOT_LIBS})
