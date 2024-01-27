//root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x TrdHit_ForPositiveISSDataStudy.C
#include "AntiprotonAnalysisTools.hh"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void TrdHit_ForPositiveISSDataStudy(){

vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());

vector<double> v_TRDHit = {1,2,3,4,5,6};

for (int index = 1; index < 2; index=index+1){
    cout<< "test: " << v_TRDHit.at(index) <<endl;

    //// Load TrdlikelihoodCut
    vector<double> v_TRDcut;
    ifstream TRDfile;
    TRDfile.open( string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/TRDLogLikelihood_CutValue_eff_0.9") + string("_1") + string(".txt") );
    string trdcutvalue;
    assert(v_TRDcut.empty());
    while (getline(TRDfile, trdcutvalue)) {
    v_TRDcut.push_back(stod(trdcutvalue));
    }
    TRDfile.close();

    double RECITOFBETALOW = -0.7;
    double RECITOFBETAHIGH = 0.3;
    int RECITOFBETANUMBER = 40;
    double TrdLOW = ( v_TRDcut.at( (index-1) ));
    double TrdHIGH = 1.7;
    int TrdBINNUMBER = 20;

    TCut TrdLikelihoodCut         = (string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW)).c_str(); // Antiproton:1to1.2, Electron:0.5
    TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
    TCut TRDVTracksSizeCut        = "TRDVTracksSize==1";
    TCut TrdNumberOfHitsCut_P     = "TrdNumberOfHits<40";
    TCut TrdNumberOfHitsCut_test  = "TrdNumberOfHits<41";
    TCut PositiveCut      = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut && TrdNumberOfHitsCut_P;
    TCut PositiveCut_test = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut && TrdNumberOfHitsCut_test;

    TCut TrdCut = ( string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW) + string("&&") + string("TrdLogLikelihoodRatioElectronProtonTracker<") + doubleToString(TrdHIGH) ).c_str();
    TCut TOFBETACut =  (string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ) ") + string(">") + to_string(RECITOFBETALOW) + string("&&") + string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ) ") + string("<") + to_string(RECITOFBETAHIGH) ).c_str();

    std::string left = doubleToString(binning[index]);
    std::string right = doubleToString(binning[index+1]);
    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
        left = to_string_with_precision(binning[index], 1);}
    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
        right = to_string_with_precision(binning[index+1], 1);}
    cout << "Now is " <<left  << "_" << right << " GV: " << endl;

    TChain *fpass7_positive = new TChain("AntiprotonLowEnergyTree");
    fpass7_positive->AddFile((string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/B1130_pass7_7.8_all_Tree_positive_") + left + string("_") + right + string("_test.root")).c_str());

    //fpass7_positive->Draw( (string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_data(") + to_string(TrdBINNUMBER) + string(", ") + to_string(TrdLOW) + string(",") + to_string(TrdHIGH) + string(", ") + to_string(RECITOFBETANUMBER) + string(", ") + to_string(RECITOFBETALOW) + string(", ") + to_string(RECITOFBETAHIGH) + string(")")).c_str(), TOFBETACut && TrdCut && PositiveCut);

    fpass7_positive->Draw("TofBeta>>ori(100,-100,100)" , TOFBETACut && TrdCut && PositiveCut);
    fpass7_positive->Draw("TofBeta>>test(100,-100,100)", TOFBETACut && TrdCut && PositiveCut_test);
    TH1F *h_ori  = (TH1F*)gDirectory->Get("ori");
    TH1F *h_test = (TH1F*)gDirectory->Get("test");

    cout<< (h_test->GetEntries() - h_ori->GetEntries()) / h_ori->GetEntries() << endl;

    fpass7_positive->Reset();

    }

}


