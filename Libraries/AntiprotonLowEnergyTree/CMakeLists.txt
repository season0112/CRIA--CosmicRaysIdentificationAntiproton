set(source_files AntiprotonLowEnergyTree.C)
set(header_files AntiprotonLowEnergyTree.hh)
add_project_library(AntiprotonLowEnergyTree SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonLowEnergyTree ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} TreePlayer)
