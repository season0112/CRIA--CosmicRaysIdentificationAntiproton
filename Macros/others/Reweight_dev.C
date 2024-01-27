void Reweight(){

TCut Pattern0Cut = "Pattern==0";

//// Positive
TFile* outputFile_positive = TFile::Open("B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_250_330.root", "open" );
TTree *outputTree_positive = (TTree*)outputFile_positive->Get("ExampleAnalysisTree");

Float_t Rigidity_p;
Double_t Weight_p;
Short_t Pattern_p;

outputTree_positive->SetBranchAddress("Weight"  , &Weight_p);
outputTree_positive->SetBranchAddress("Rigidity", &Rigidity_p);
outputTree_positive->SetBranchAddress("Pattern" , &Pattern_p);

TH1D *B1042positive_RereweightedRigidity_p = new TH1D("", "", 10, 200, 350);

double totalweight_p = 0;
int nentries_p = outputTree_positive->GetEntries(Pattern0Cut);
for (int i=0; i<nentries_p; i++) {
    outputTree_positive->GetEntry(i);
    //cout<< Rigidity_p <<endl;
    //cout<< Weight_p <<endl;
    if (Pattern_p==0){
        B1042positive_RereweightedRigidity_p->Fill(Rigidity_p, Weight_p); 
        totalweight_p = totalweight_p + Weight_p;}
}
cout<< "totalweight_p is " << totalweight_p <<endl; 
cout<< "totalevent_p is" << nentries_p <<endl;

cout<< "Part2 :" << endl;
double totalcontent_p = 0;
for (int p=1; p<B1042positive_RereweightedRigidity_p->GetNbinsX(); p++) {
    totalcontent_p = totalcontent_p + B1042positive_RereweightedRigidity_p->GetBinContent(p);
}
cout<< "totalcontent_p is " << totalcontent_p <<endl;


//// Negative
TFile* outputFile_negative = TFile::Open("B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_250_330.root", "open" );
TTree *outputTree_negative = (TTree*)outputFile_negative->Get("ExampleAnalysisTree");

Float_t Rigidity_n;
Double_t Weight_n;
Short_t Pattern_n;

outputTree_negative->SetBranchAddress("Weight"  , &Weight_n);
outputTree_negative->SetBranchAddress("Rigidity", &Rigidity_n);
outputTree_negative->SetBranchAddress("Pattern" , &Pattern_n);

TH1D *B1042positive_RereweightedRigidity_n = new TH1D("", "", 10, -350, -200);

double totalweight_n = 0;
int nentries_n = outputTree_negative->GetEntries(Pattern0Cut);
for (int i=0; i<nentries_n; i++) {
    outputTree_negative->GetEntry(i);
    //cout<< Rigidity_p <<endl;
    //cout<< Weight_p <<endl;
    if (Pattern_n==0){
        B1042positive_RereweightedRigidity_n->Fill(Rigidity_n, Weight_n);
        totalweight_n = totalweight_n + Weight_n;}
}
cout<< "totalweight_n is " << totalweight_n <<endl;
cout<< "totalevent_n is " << nentries_n <<endl;

cout<< "Part2 :" << endl;
double totalcontent_n = 0;
for (int p=1; p<B1042positive_RereweightedRigidity_n->GetNbinsX(); p++) {
    totalcontent_n = totalcontent_n + B1042positive_RereweightedRigidity_n->GetBinContent(p);
}
cout<< "totalcontent_n is" << totalcontent_n <<endl;








/*
outputTree_positive->Draw("Rigidity_p>>B1042positive_Rigidity_p");
//outputFile_positive->Draw("Rigidity_p>>B1042positive_Rigidity_p", Pattern0Cut);
//outputFile_positive->Draw("Weight_p>>B1042positive_Weight_p"    , Pattern0Cut);
TH1D *B1042positive_Rigidity_p = (TH1D*)gDirectory->Get("B1042positive_Rigidity_p");
//TH1D *B1042positive_Weight_p   = (TH1D*)gDirectory->Get("B1042positive_Weight_p");

//TH1D *B1042positive_RereweightedRigidity_p = new TH1D(*B1042positive_Rigidity_p);
//TH1D *B1042positive_RereweightedRigidity_p = new TH1D();


//for (int i=1; i<B1042positive_Rigidity_p->GetEntries(); i++){
//for (int i=1; i<100; i++){
//    cout<< B1042positive_Rigidity_p->GetBinContent(i) <<endl;
//    B1042positive_RereweightedRigidity_p->Fill(B1042positive_Rigidity_p->GetBinContent(i), B1042positive_Weight_p->GetBinContent(i));

//    B1042positive_RereweightedRigidity_p->Fill(outputFile_positive
//}
*/


/*
//// Plot
B1042positive_Rigidity_p->SetLineColor(2);
B1042positive_RereweightedRigidity_p->SetLineColor(4);

//TCanvas * c112 = new TCanvas;
//B1042positive_Rigidity_p->Draw("");
//B1042positive_RereweightedRigidity_p->Draw("same");
//c112->SaveAs( string("Ztest.pdf").c_str());

TCanvas * c113 = new TCanvas;
B1042positive_RereweightedRigidity_p->Draw("HIST");
B1042positive_Rigidity_p->Draw("same");
c113->SaveAs( string("Ztest2.pdf").c_str());
*/

}


