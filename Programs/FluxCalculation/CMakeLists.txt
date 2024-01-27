set(source_files FluxCalculation.C)
add_project_executable(FluxCalculation SOURCES ${source_files} ${source_files_exe} HEADERS ${header_files} LIBS ExampleAnalysisTree AntiprotonBinning)

