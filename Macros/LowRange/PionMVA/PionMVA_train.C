void PionMVA_train(){

//// TCut 
// Tracker
std::string trackerpattern = "0";
TCut PatternCut = (string("Pattern==") + trackerpattern).c_str();
TCut ExtrapolatedPhotoElectronsFromTracker ="NPhotoElectrons-ExtrapolatedRichExpectedPhotoElectronsProton<20 && NPhotoElectrons>-999 && ExtrapolatedRichTileIndex == TileIndex || NPhotoElectrons==-999 || ExtrapolatedRichTileIndex != TileIndex";
TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";
// Trigger
TCut TriggerPhysics = "TriggerFlags != 1 && TriggerFlags != 64 && TriggerFlags != 65";
// TRD
TCut TrdLikelihoodCut = "TrdLogLikelihoodRatioElectronProtonTracker > 0.7 || TrdLogLikelihoodRatioElectronProtonTracker==-1.5";    // Antiproton:1to1.2, Electron:0.5
TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
// ECAL
TCut EcalBDT_EnergyDCut = "EcalBDT_EnergyD < -0.9";
// RICH
TCut RichBetaCut = "RichBeta==0 || RichIsNaF==1";
TCut RichNaFCut = "RichIsNaF==1";
TCut RichAglCut = "RichIsNaF==0";
TCut RichPhotoElectron = "NPhotoElectrons-NExpectedPhotoElectrons<10  || NPhotoElectrons==-999";
TCut RichChargeCut = "RichCharge<2 || RichCharge==0";
// TOF
//TCut BetaConverted = "BetaConverted<0.9";
//TCut TofMassonecharge = "TofMassonecharge>0.5";
TCut TofMassonechargeCut = "TofMassonecharge>0";
//TCut TOFBETALikelihood = "TOFBETALikelihood<0.7";
// All Cuts
TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonechargeCut && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;

//// Load training data
string bin[17] = {"0.8_1.0", "1.0_1.16", "1.16_1.33", "1.33_1.51", "1.51_1.71", "1.71_1.92", "1.92_2.15", "2.15_2.4", "2.4_2.67", "2.67_2.97",  "2.97_3.29","3.29_3.64","3.64_4.02","4.02_4.43","4.43_4.88","4.88_5.37","5.37_5.9"  };
//string bin[1] = {"1.51_1.71"};

chdir( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/test2/")).c_str());

for (int i=0; i<17; i++){
    TChain *fb1 = new TChain ("AntiprotonLowEnergyTree");
    TChain *fb2 = new TChain ("AntiprotonLowEnergyTree");
    TChain *fs = new TChain ("AntiprotonLowEnergyTree");
    fb1->AddFile( (string("../B1091_el.pl1.0_25200_7.6_all_Tree_negative_")  + bin[i] + string(".root")).c_str() );
    fb2->AddFile( (string("../B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + bin[i] + string(".root")).c_str() );
    fs->AddFile(  (string("../B1042_antipr.pl1.1800_7.6_all_Tree_") + bin[i] + string(".root")).c_str() );
//    fs->AddFile(  (string("../B1130_pass7_7.8_all_Tree_positive_May2015_") + bin[i] + string(".root")).c_str() ); 

    cout<<  fb1->GetEntries() <<endl;
    cout<<  fb2->GetEntries() <<endl;
    cout<<  fs->GetEntries() <<endl;

    TTree *signalTree = fs->GetTree();
    TTree *backgroundA = fb1->GetTree();
    TTree *backgroundB = fb2->GetTree();

    TString outfileName( (string("TMVA_") + bin[i] + string(".root")).c_str() );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("PionMVA_dataset_" + bin[i]);
    dataloader->AddVariable( "RichCharge",                                                  "Variable 1", "units", 'F' );
    dataloader->AddVariable( "NCollectedPhotoElectrons",                                    "Variable 2", "units", 'F' );
    dataloader->AddVariable( "RICHNumberOfHits",                                            "Variable 3", "units", 'F' );
    dataloader->AddVariable( "TrdLogLikelihoodRatioProtonHeliumTracker",                    "Variable 4", "units", 'F' );
    dataloader->AddVariable( "UpperTofCharge",                                              "Variable 5", "units", 'F' );
    dataloader->AddVariable( "LowerTofCharge",                                              "Variable 6", "units", 'F' );
    dataloader->AddVariable( "CCBDTLapp",                                                   "Variable 7", "units", 'F' );
    dataloader->AddVariable( "(UpperTofBeta-LowerTofBeta)/(UpperTofBeta+LowerTofBeta)",     "Variable 8", "units", 'F' );
//    dataloader->AddVariable( "TofMassonecharge",      "Variable 3", "units", 'F' );
    dataloader->AddSpectator( "Rigidity",  "Spectator 1", "units", 'F' );

    Double_t signalWeight     = 1.0/signalTree->GetEntries()/2;
    Double_t backgroundAWeight = 1.0/backgroundA->GetEntries()/4;
    Double_t backgroundBWeight = 1.0/backgroundB->GetEntries()/4;

    dataloader->AddSignalTree    ( signalTree,     signalWeight );
    dataloader->AddBackgroundTree( backgroundA, backgroundAWeight );
    dataloader->AddBackgroundTree( backgroundB, backgroundBWeight );

    TCut mycutsignal = NegativeCut;
    TCut mycutbackground = NegativeCut;

    dataloader->PrepareTrainingAndTestTree( mycutsignal, mycutbackground,
                                        "nTrain_Signal=0.8:nTrain_Background=0.8:nTest_Signal=0.2:nTest_Background=0.2:SplitMode=Random:NormMode=EqualNumEvents:!V" );

    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:GradBaggingFraction=0.5:nCuts=50:MaxDepth=3:IgnoreNegWeightsInTraining:VarTransform=Norm");

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    TCanvas *canvas = (TCanvas*)factory->GetROCCurve(dataloader);
    canvas->Draw();
    canvas->Print( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/test2/result_bdt_") + bin[i] + string(".png")).c_str() );

    fb1->Reset();
    fb2->Reset();
    fs->Reset();

    }
}





