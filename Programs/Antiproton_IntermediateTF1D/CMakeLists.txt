set(source_files MakeHistogram_TimeAveraged.C)
add_project_executable(Antiproton_IntermediateMakeTemplate SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonIntermediateEnergyTree AntiprotonIntermediateEnergyTemplatefitter AntiprotonBinning)

set(source_files MakeHistogram_TimeDependent.C)
add_project_executable(Antiproton_IntermediateMakeTemplate_TimeDependent SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonIntermediateEnergyTree AntiprotonIntermediateEnergyTemplatefitter AntiprotonBinning)

set(source_files 1DTemplateFit_TimeAveraged.C)
add_project_executable(Antiproton_IntermediateTF1D SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonIntermediateEnergyTree AntiprotonIntermediateEnergyTemplatefitter AntiprotonBinning)

set(source_files 1DTemplateFit_TimeDependent.C)
add_project_executable(Antiproton_IntermediateTF1D_TimeStamp_3and6_BartalRotations SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonIntermediateEnergyTree AntiprotonIntermediateEnergyTemplatefitter AntiprotonBinning)

set(source_files Timestamp_edge.C)
add_project_executable(Timestamp_edge SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonBinning)


