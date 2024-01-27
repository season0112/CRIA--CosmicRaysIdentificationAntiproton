set(source_files main.C)
add_project_executable(TriggerEfficiency SOURCES ${source_files} HEADERS ${header_files} /usr/include LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} )


