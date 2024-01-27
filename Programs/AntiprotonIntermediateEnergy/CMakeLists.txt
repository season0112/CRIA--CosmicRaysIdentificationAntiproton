set(source_files main.C)
add_project_executable(AntiprotonIntermediateEnergy SOURCES ${source_files} HEADERS ${header_files} /usr/include LIBS ${ACSOFT_LIBRARIES} ${ACSOFT_AUX_LIBRARIES} AntiprotonIntermediateEnergyTree AntiprotonBinning)



#set(source_files_lib AntiprotonIntermediateEnergy.C)
#set(header_files_lib_with_classdef AntiprotonIntermediateEnergy.hh)
#add_project_library(AntiprotonIntermediateEnergyLibrary SOURCES ${source_files_lib} HEADERS_WITH_CLASSDEF ${header_files_lib_with_classdef} LINKDEF Linkdef.h LIBS )


