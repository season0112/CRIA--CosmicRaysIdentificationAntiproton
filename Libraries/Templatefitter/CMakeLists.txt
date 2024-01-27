set(source_files TemplateFitter2DforProtonCC.C TemplateFitterforProtonCC.C Statistics.C Quantity.C TemplateFitter2DforProtonCC_CCfixed_dev.C)
set(header_files TemplateFitter2DforProtonCC.hh TemplateFitterforProtonCC.hh Statistics.hh Quantity.hh TemplateFitter2DforProtonCC_CCfixed_dev.hh)

add_project_library(Templatefitter SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(Templatefitter ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY})

