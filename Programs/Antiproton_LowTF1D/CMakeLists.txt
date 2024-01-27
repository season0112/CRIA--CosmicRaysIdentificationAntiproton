set(source_files MakeHistogram_TimeAveraged.C)
add_project_executable(Antiproton_LowMakeTemplate SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonLowEnergyTree AntiprotonLowEnergyTemplatefitter AntiprotonBinning AntiprotonAnalysisTools)

set(source_files MakeHistogram_TimeDependent.C)
add_project_executable(Antiproton_LowMakeTimeDependentTemplate SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonLowEnergyTree AntiprotonLowEnergyTemplatefitter AntiprotonBinning AntiprotonAnalysisTools)

set(source_files 2DTemplateFit_TimeAveraged.C)
add_project_executable(Antiproton_LowTF1D SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonLowEnergyTree AntiprotonLowEnergyTemplatefitter AntiprotonBinning AntiprotonAnalysisTools)

set(source_files 2DTemplateFit_TimeDependent.C)
add_project_executable(Antiproton_LowTF1D_TimeStamp_3and6_BartalRotations SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonLowEnergyTree AntiprotonLowEnergyTemplatefitter AntiprotonBinning AntiprotonAnalysisTools)


