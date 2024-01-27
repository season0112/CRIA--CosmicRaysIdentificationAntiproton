cmake_minimum_required(VERSION 2.8)

# find acsoft first, to override compiler
if (ACsoft_DIR)
  find_package(ACsoft REQUIRED NO_DEFAULT_PATH)
else()
  find_package(ACsoft REQUIRED)
endif()

PROJECT(ExampleAnalysis)
complete_ACSoft_setup()

# set ExampleAnalysis includes
include_directories("${PROJECT_SOURCE_DIR}/Libraries/ExampleAnalysisTree")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/LeptonUnfoldingTools")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/RooUnfold")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/LeptonRooPdfs")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/LeptonAnalysisTools")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/Templatefitter")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonIntermediateEnergyTree")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/RigidityResolutionTree")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonIntermediateEnergyTemplatefitter")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonLowEnergyTemplatefitter")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonBinning")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonAnalysisTools")
include_directories("${PROJECT_SOURCE_DIR}/Libraries/AntiprotonLowEnergyTree")

# doxygen target
include(Auxiliary/doxygen.cmake)

set(FOLDER_NAME "Libraries")
add_subdirectory(Libraries)

set(FOLDER_NAME "Programs")
add_subdirectory(Programs)

set(FOLDER_NAME "ROOT Macros")
add_subdirectory(Macros)

