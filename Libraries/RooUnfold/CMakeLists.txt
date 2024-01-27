set(source_files
  RooUnfold.cxx
  RooUnfoldBayes.cxx
  RooUnfoldBinByBin.cxx
  RooUnfoldErrors.cxx
  RooUnfoldIds.cxx
  RooUnfoldInvert.cxx
  RooUnfoldParms.cxx
  RooUnfoldResponse.cxx
  RooUnfoldSvd.cxx
  RooUnfoldTUnfold.cxx
  TSVDUnfold.cxx
)

# RooUnfoldDagostini.cxx excluded, as it needs fortran, and we don't use it.

set(header_files_with_classdef
  RooUnfold.h
  RooUnfoldBayes.h
  RooUnfoldBinByBin.h
  RooUnfoldDagostini.h
  RooUnfoldErrors.h
  RooUnfoldIds.h
  RooUnfoldInvert.h
  RooUnfoldParms.h
  RooUnfoldResponse.h
  RooUnfoldSvd.h
  RooUnfoldTUnfold.h
  TSVDUnfold_local.h
)

set(LINKED_LIBRARIES ${ROOT_LIBRARIES})
string(REGEX REPLACE "/" "." _ROOT_VERSION ${ROOT_VERSION})

set(ADDITIONAL_CFLAGS )

# Only provide RooUnfoldTUnfold if TUnfold is available.
# TUnfold is included in ROOT 5.22 and later.
if(${_ROOT_VERSION} VERSION_GREATER 5.22.0)
  set(ADDITIONAL_CFLAGS ${ADDITIONAL_CFLAGS} -DHAVE_TUNFOLD=1)
endif()

# Only use TUnfoldSys if it is available (ROOT 5.23.02 and later).
if(${_ROOT_VERSION} VERSION_GREATER 5.23.02)
  set(ADDITIONAL_CFLAGS ${ADDITIONAL_CFLAGS} -DHAVE_TUNFOLDSYS=1)
endif()

# In ROOT version 6.09.02 and above TUnfold was split from libHist to libUnfold
if(${_ROOT_VERSION} VERSION_GREATER 6.09.01)
  set(LINKED_LIBRARIES ${LINKED_LIBRARIES} Unfold)
endif()

add_project_library(RooUnfold SOURCES ${source_files} HEADERS_WITH_CLASSDEF ${header_files_with_classdef} LINKDEF Linkdef.h)
target_compile_definitions(RooUnfold PUBLIC ${ADDITIONAL_CFLAGS})
target_link_libraries(RooUnfold ${LINKED_LIBRARIES})
