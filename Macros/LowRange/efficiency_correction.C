

TFile *f0 = new TFile("AntiprotonLowEnergy_Auxiliary_00000_00015.root");

Cuts::Selector *Preselection = (Cuts::Selector*)f0->Get("Preselection");
Cuts::Selector *QualityCuts = (Cuts::Selector*)f0->Get("QualityCuts");

Preselection->GetCut(0)->FindAttachment("Cuts::EfficiencyHistograms");


TEfficiency *QualityCuts1MC = static_cast<const Cuts::EfficiencyHistograms*>(QualityCuts->GetCut(1)->FindAttachment("Cuts::EfficiencyHistograms"))->ProduceTagAndProbeEfficiencyMc("Quality");

TH1 *eff1 = QualityCuts1MC->GetCopyPassedHisto();
eff1->Divide(QualityCuts1MC->GetCopyPassedHisto(), QualityCuts1MC->GetTotalHistogram(),1,1,"B"); // checked: Error is same with QualityCuts1MC.




