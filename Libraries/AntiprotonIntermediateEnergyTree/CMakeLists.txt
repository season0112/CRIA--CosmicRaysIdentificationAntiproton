set(source_files AntiprotonIntermediateEnergyTree.C)
set(header_files AntiprotonIntermediateEnergyTree.hh)
add_project_library(AntiprotonIntermediateEnergyTree SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonIntermediateEnergyTree ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} TreePlayer)
