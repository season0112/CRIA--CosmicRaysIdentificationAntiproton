set(source_files main.C)
add_project_executable(BartelRotation SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES})

set(source_files Timestamp_edge.C)
add_project_executable(Timestamp_edge SOURCES ${source_files} HEADERS ${header_files} /usr/include  LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES})

