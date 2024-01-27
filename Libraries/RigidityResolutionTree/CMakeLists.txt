set(source_files RigidityResolutionTree.C)
set(header_files RigidityResolutionTree.hh)
add_project_library(RigidityResolutionTree SOURCES ${source_files} HEADERS ${header_files})
target_link_libraries(RigidityResolutionTree ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} ${ACSOFT_GBATCH_LIBRARY} TreePlayer)
