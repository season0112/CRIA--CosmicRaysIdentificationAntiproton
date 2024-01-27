set(source_files
  RooAsymmLandauLikePdf.C
  RooNovosibirskPdf.C
  RooLangausPdf.C
  RooGaussianPdf.C
  RooGaussExpPdf.C
  RooExpGaussExpPdf.C
)

set(header_files_with_classdef
  RooAsymmLandauLikePdf.hh
  RooNovosibirskPdf.hh
  RooLangausPdf.hh
  RooGaussianPdf.hh
  RooGaussExpPdf.hh
  RooExpGaussExpPdf.hh
)

add_project_library(LeptonRooPdfs SOURCES ${source_files} HEADERS_WITH_CLASSDEF ${header_files_with_classdef} LINKDEF Linkdef.h)
target_link_libraries(LeptonRooPdfs ${ROOT_LIBRARIES} RooFitCore)
