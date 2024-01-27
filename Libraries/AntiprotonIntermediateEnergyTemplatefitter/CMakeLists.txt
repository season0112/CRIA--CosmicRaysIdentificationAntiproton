set(source_files TemplateFitterforIntermediateEnergy.C Statistics.C Quantity.C)
set(header_files TemplateFitterforIntermediateEnergy.hh Statistics.hh Quantity.hh)
add_project_library(AntiprotonIntermediateEnergyTemplatefitter SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonIntermediateEnergyTemplatefitter ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY})

