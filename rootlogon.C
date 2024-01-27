{
  cout << " ++++ ExampleAnalysis rootlogon ++++ " << endl;
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/ExampleAnalysisTree");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/LeptonUnfoldingTools");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/RooUnfold");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/LeptonAnalysisTools");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/LeptonRooPdfs");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/Templatefitter");
  gInterpreter->AddIncludePath("$MY_ANALYSIS/Libraries/AntiprotonIntermediateEnergyTree");
  gSystem->Load("$MY_ANALYSIS/lib/libExampleAnalysisLibrary");

}


