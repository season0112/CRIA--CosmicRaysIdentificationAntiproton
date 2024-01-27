set(source_files AntiprotonBinning.C)
set(header_files AntiprotonBinning.hh)
add_project_library(AntiprotonBinning SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonBinning ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} TreePlayer)
