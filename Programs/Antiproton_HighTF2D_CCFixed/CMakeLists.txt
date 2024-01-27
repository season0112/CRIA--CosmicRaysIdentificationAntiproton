set(source_files main.C)
add_project_executable(Antiproton_HighTF2D_Unbinned SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ExampleAnalysisTree Templatefitter AntiprotonAnalysisTools)

