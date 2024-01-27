set(source_files main.C)
add_project_executable(Unfolding_Bayes SOURCES ${source_files} HEADERS ${header_files} LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ExampleAnalysisTree LeptonUnfoldingTools RooUnfold LeptonAnalysisTools LeptonRooPdfs AntiprotonBinning)

