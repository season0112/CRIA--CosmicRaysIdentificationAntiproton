set(source_files
  AntiprotonAnalysisTools.C
)

set(header_files
  AntiprotonAnalysisTools.hh
)

add_project_library(AntiprotonAnalysisTools SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonAnalysisTools ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY})
