
// Rich
TCut HasRichBetaCut = "0<RichBeta";


// TRD
TCut TrdLikelihoodHeProtonCut_proton = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.2";
TCut TrdLikelihoodHeProtonCut_electron = "TrdLogLikelihoodRatioProtonHeliumTracker > 0.1";
TCut TrdSegmentsXZNumberCut = "TrdSegmentsXZNumber==1";
TCut TrdSegmentsYZNumberCut = "TrdSegmentsYZNumber==1";
TCut TRDVTracksSizeCut = "TRDVTracksSize==1";
TCut TrdNumberOfHitsCut = "TrdNumberOfHits<35";


// ECal
TCut EcalBDT_EnergyDCut_proton = "EcalBDT_EnergyD < -0.5";                             //-1:proton like, 1:electron like. not applicable:-999.
TCut EcalBDT_EnergyDCut_electron = "EcalBDT_EnergyD > 0";                              //-1:proton like, 1:electron like. not appplicable:-999.
TCut EcalBDT_EnergyDCut_electron_v2 = "EcalBDT_EnergyD > 0 || EcalBDT_EnergyD < -100"; //-1:proton like, 1:electron like. not appplicable:-999.


// Tracker
std::tuple<TCut> patterncut(std::string trackerpattern){
    TCut PatternCut;
    if (trackerpattern==std::string("01")){
        std::cout<< "Current tracker pattern is 0 and 1." <<std::endl;
        PatternCut = ("Pattern==0 || Pattern==1");
    }
    else if (trackerpattern==std::string("0124")){
        std::cout<< "Current tracker pattern is 0,1,2,4" <<std::endl;
        PatternCut = ("Pattern==0 || Pattern==1 || Pattern==2 || Pattern==4");
    }
    else{
        std::cout<< "Current tracker pattern is "<< trackerpattern <<std::endl;
        PatternCut = (std::string("Pattern==") + trackerpattern).c_str();
    }
    return {PatternCut};
}
TCut Pattern0and1 = ("Pattern==0 || Pattern==1");
TCut Pattern0 = ("Pattern==0");
TCut Pattern1 = ("Pattern==1");
TCut Pattern2 = ("Pattern==2");
TCut Pattern4 = ("Pattern==4");


// ACC
TCut ACCHitsCut = "ACCHits==0";


// Others
TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";



////
// In progress
////

//// Templates selectron
//TCut ElectronCut       = RichBetaCut && PatternCut && ProtonCCMVABDCut && TrdLikelihoodHeProtonCut_electron && EcalBDT_EnergyDCut_electron;
//TCut ElectronCut_v2    = RichBetaCut && PatternCut && ProtonCCMVABDCut && TrdLikelihoodHeProtonCut_electron && EcalBDT_EnergyDCut_electron_v2;
//TCut ElectronCut_p0    = RichBetaCut && Pattern0 && ProtonCCMVABDCut && TrdLikelihoodHeProtonCut_electron && EcalBDT_EnergyDCut_electron;
//TCut ElectronCut_p0_v2 = RichBetaCut && Pattern0 && ProtonCCMVABDCut && TrdLikelihoodHeProtonCut_electron && EcalBDT_EnergyDCut_electron_v2;
//TCut AllCut     = RichBetaCut && PatternCut && ProtonCCMVABDCut && TrdLikelihoodHeProtonCut_proton && EcalBDT_EnergyDCut_proton;

//// Pion and electron template from Data  (Mass_electron=0.000511)
// Set1:
//TCut PionTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72";
//TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002 && TrdLogLikelihoodRatioElectronProtonTracker<0.55";
//TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002";
//TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002";
//TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker<0.55";

// Set2:
//TCut ElectronTemplateDataCut = "EcalBDT_EnergyD>-999 && EcalBDT_EnergyD>0.5";
//TCut PionTemplateDataCut = "EcalBDT_EnergyD>-999 && EcalBDT_EnergyD<-0.5";

// Set3:
//TCut ElectronTemplateDataCut = PatternCut && ProtonCCMVABDCut && "abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002" && "RichBeta>0" && "EcalBDT_EnergyD>-999" && "EcalBDT_EnergyD>0.5" && "TrdLogLikelihoodRatioProtonHeliumTracker>-1" && "TrdLogLikelihoodRatioProtonHeliumTracker>0.1";
//TCut PionTemplateDataCut = PatternCut && ProtonCCMVABDCut && "abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002" && "RichBeta>0" && "EcalBDT_EnergyD>-999" && "EcalBDT_EnergyD<-0.5" && "TrdLogLikelihoodRatioProtonHeliumTracker>-1" && "TrdLogLikelihoodRatioProtonHeliumTracker<0.1";

// Set4:
//TCut ElectronTemplateDataCut = PatternCut && ProtonCCMVABDCut && "EcalBDT_EnergyD>-999" && "EcalBDT_EnergyD>0.5" && "TrdLogLikelihoodRatioProtonHeliumTracker<0.2" && "abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002";
//TCut PionTemplateDataCut = PatternCut && ProtonCCMVABDCut && "abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002" && "RichBeta>0" && "EcalBDT_EnergyD>-999" && "EcalBDT_EnergyD<-0.2" && "TrdLogLikelihoodRatioProtonHeliumTracker>-1" && "TrdLogLikelihoodRatioProtonHeliumTracker<0.07";

// Set5: (same as LowenErgyRange)
//TCut ElectronTemplateDataCut = PatternCut && "EcalBDT_EnergyD>0.5" && "abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002" && "TrdSegmentsXZNumber==1" && "TrdSegmentsYZNumber==1" && "TrdNumberOfHits<35" && "RichBeta>0"; //"RichBeta>0"?
//TCut PionTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5 && TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>70 && TrdLogLikelihoodRatioProtonHeliumTracker<0.02 && TrdLogLikelihoodRatioProtonHeliumTracker>-1";

// setting two,four
//TCut ElectronTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5" && "EcalBDT_EnergyD>0.5" && "TrdSegmentsXZNumber==1" && "TrdSegmentsYZNumber==1" && "TrdNumberOfHits<35";
//TCut PionTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5" && "EcalBDT_EnergyD>-1" && "EcalBDT_EnergyD<0.5" && "TrdSegmentsXZNumber>1" && "TrdSegmentsYZNumber>1" && "TrdNumberOfHits>40";








