void PionMVA_App(){

//TCut
// Allcut:     TRD_electron_proton: all, TRD_He_Proton: Antiproton Template, positivedata,
// Tracker
std::string trackerpattern = "0";
TCut PatternCut = (string("Pattern==") + trackerpattern).c_str();
TCut ExtrapolatedPhotoElectronsFromTracker ="NPhotoElectrons-ExtrapolatedRichExpectedPhotoElectronsProton<20 && NPhotoElectrons>-999 && ExtrapolatedRichTileIndex == TileIndex || NPhotoElectrons==-999 || ExtrapolatedRichTileIndex != TileIndex";
TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";
//Trigger
TCut TriggerPhysics = "TriggerFlags != 1 && TriggerFlags != 64 && TriggerFlags != 65";
//TRD
TCut TrdLikelihoodCut = "TrdLogLikelihoodRatioElectronProtonTracker > 0.7 || TrdLogLikelihoodRatioElectronProtonTracker==-1.5";
TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
//ECAL
TCut EcalBDT_EnergyDCut = "EcalBDT_EnergyD < -0.9";
//RICH
TCut RichBetaCut = "RichBeta==0 || RichIsNaF==1";
TCut RichNaFCut = "RichIsNaF==1";
TCut RichAglCut = "RichIsNaF==0";
TCut RichPhotoElectron = "NPhotoElectrons-NExpectedPhotoElectrons<10  || NPhotoElectrons==-999";
TCut RichChargeCut = "RichCharge<2 || RichCharge==0";
//TOF
//TCut BetaConverted = "BetaConverted<0.9";
//TCut TofMassonechargecut = "TofMassonecharge>0.5";
TCut TofMassonechargecut = "TofMassonecharge>0";
//TCut TOFBETALikelihood = "TOFBETALikelihood<0.7";

TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonechargecut && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
//TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;

//string bin[14] = {"1.33_1.51","1.51_1.71", "1.71_1.92", "1.92_2.15", "2.15_2.4", "2.4_2.67", "2.67_2.97","2.97_3.29","3.29_3.64","3.64_4.02","4.02_4.43","4.43_4.88","4.88_5.37","5.37_5.9"};
string bin[17] = {"0.8_1.0", "1.0_1.16", "1.16_1.33", "1.33_1.51","1.51_1.71", "1.71_1.92", "1.92_2.15", "2.15_2.4", "2.4_2.67", "2.67_2.97",  "2.97_3.29","3.29_3.64","3.64_4.02","4.02_4.43","4.43_4.88","4.88_5.37","5.37_5.9"  };

for (int i=0; i<17; i++){

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    Float_t RichCharge, TofBeta, BetaConverted, TOFBETALikelihood, UpperTofCharge, LowerTofCharge, Rigidity, CCBDTLapp, NCollectedPhotoElectrons, TofMassonecharge, RICHNumberOfHits, TrdLogLikelihoodRatioProtonHeliumTracker, UpperTofBeta, LowerTofBeta, TofBetaSymmetry;

    reader->AddVariable( "RichCharge",                                  &RichCharge );
    reader->AddVariable( "NCollectedPhotoElectrons",                    &NCollectedPhotoElectrons );
    reader->AddVariable( "RICHNumberOfHits",                            &RICHNumberOfHits );
    reader->AddVariable( "TrdLogLikelihoodRatioProtonHeliumTracker", &TrdLogLikelihoodRatioProtonHeliumTracker);
    reader->AddVariable( "UpperTofCharge", &UpperTofCharge);
    reader->AddVariable( "LowerTofCharge", &LowerTofCharge);
    reader->AddVariable( "CCBDTLapp", &CCBDTLapp);
    reader->AddVariable( "(UpperTofBeta-LowerTofBeta)/(UpperTofBeta+LowerTofBeta)", &TofBetaSymmetry);
    //reader->AddVariable( "TofMassonecharge",   &TofMassonecharge );
    reader->AddSpectator( "Rigidity",   &Rigidity );
    reader->BookMVA("BDTG", (string("/rwthfs/rz/cluster/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/PionMVA_dataset_") + bin[i] + string("/weights/TMVAClassification_BDTG.weights.xml")).c_str() );


    //TFile *input = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1130_pass7_7.8_all_Tree_positive_May2015_1.71_1.92.root");
    //TFile *input = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_1.71_1.92.root");
    //TFile *input = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1091_el.pl1.0_25200_7.6_all_Tree_1.71_1.92.root");
    //TFile *input = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1042_antipr.pl1.1800_7.6_all_Tree_1.71_1.92.root");

//    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1130_pass7_7.8_all_Tree_negative_May2015_") + bin[i] + string(".root")).c_str() );
    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1130_pass7_7.8_all_Tree_negative_") + bin[i] + string(".root")).c_str() );
//    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1042_antipr.pl1.1800_7.6_all_Tree_") + bin[i] + string(".root")).c_str() );
//    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1130_pass7_7.8_all_Tree_positive_May2015_") + bin[i] + string(".root")).c_str() );
//    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + bin[i] + string(".root")).c_str() );
//    TFile *input = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + bin[i] + string(".root")).c_str() );

    TTree* OriginalTree = (TTree*)input->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *theTree = OriginalTree->CopyTree(NegativeCut);
    input->Close();

    Float_t PionBDT;
    auto PionMVA = theTree->Branch("PionBDT", &PionBDT, "PionBDT/F");

    theTree->SetBranchAddress( "RichCharge", &RichCharge );
    theTree->SetBranchAddress( "NCollectedPhotoElectrons", &NCollectedPhotoElectrons );
    theTree->SetBranchAddress( "RICHNumberOfHits", &RICHNumberOfHits);
    theTree->SetBranchAddress( "TrdLogLikelihoodRatioProtonHeliumTracker", &TrdLogLikelihoodRatioProtonHeliumTracker );
    theTree->SetBranchAddress( "UpperTofCharge", &UpperTofCharge );
    theTree->SetBranchAddress( "LowerTofCharge", &LowerTofCharge );
    theTree->SetBranchAddress( "CCBDTLapp", &CCBDTLapp );
    theTree->SetBranchAddress( "UpperTofBeta", &UpperTofBeta);
    theTree->SetBranchAddress( "LowerTofBeta", &LowerTofBeta);
    //theTree->SetBranchAddress( "(UpperTofBeta-LowerTofBeta)/(UpperTofBeta+LowerTofBeta)", &TofBetaSymmetry );
    //theTree->SetBranchAddress( "TofMassonecharge", &TofMassonecharge );
    theTree->SetBranchAddress( "Rigidity", &Rigidity );

    std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();
    for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
      if (ievt%100000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
      theTree->GetEntry(ievt);
      TofBetaSymmetry = (UpperTofBeta-LowerTofBeta)/(UpperTofBeta+LowerTofBeta);
      PionBDT = reader->EvaluateMVA( "BDTG");
      PionMVA->Fill();
      if (ievt>20000) break;
    }
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

//    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1130_pass7_7.8_all_Tree_positive_May2015_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
//    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1130_pass7_7.8_all_Tree_negative_May2015_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1130_pass7_7.8_all_Tree_negative_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
//    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1042_antipr.pl1.1800_7.6_all_Tree_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
//    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
//    TFile *target  = new TFile( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/test2/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + bin[i] + string("_withMVA.root")).c_str(),"RECREATE" );
    theTree->Write();
    delete reader;
}

}
