
void Plot_TimeDependentMeasuringTime() {

    const std::vector<string> TimeEdge{
    "1305417600", "1319414400", "1333411200", "1347408000", "1361404800", "1375401600",
    "1389398400", "1403395200", "1417392000", "1431388800", "1445385600", "1459382400",
    "1473379200", "1487376000", "1501372800", "1515369600", "1529366400", "1543363200",
    "1557360000", "1571356800", "1585353600", "1599350400", "1613347200", "1627344000",
    };

    for (int j = 0; j < 23; j++) {
        cout<< j <<endl;
        //// Open File and get histograms
        TFile *t1 = new TFile( (string("/hpcwork/jara0052/sichen/Measuringtime/TimeDependent/MeasuringTime_pass7.8_11_2021_GEOMETRIC25_1.2_") + TimeEdge[j] + string("_") + TimeEdge[j+1] + string(".root")).c_str() );
        TH1D *h1 = (TH1D*)t1->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff"); 
        h1->Rebin(2);
        cout<< h1->GetBinContent(10) << endl; // h1->GetBinLowEdge(10) = 6.47 GV; h1->GetBinLowEdge(25) = 69.7 GV
    }

    return 0;    
}

