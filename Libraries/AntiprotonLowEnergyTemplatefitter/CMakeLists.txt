set(source_files TemplateFitterforLowEnergy.C Statistics.C Quantity.C TemplateFitterforLowEnergy2D.C)
set(header_files TemplateFitterforLowEnergy.hh Statistics.hh Quantity.hh TemplateFitterforLowEnergy2D.hh)
add_project_library(AntiprotonLowEnergyTemplatefitter SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(AntiprotonLowEnergyTemplatefitter ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY})

