import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import matplotlib.dates as md
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
#import arrow
import os
import argparse
import binning
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import numpy.polynomial.polynomial as poly
from scipy import optimize



p0_2016  = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_0_0.997_2016paperbinmerge1.root")
p1_2016  = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_1_0.997_2016paperbinmerge1.root")
p2_2016  = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_2_0.997_2016paperbinmerge1.root")
p4_2016  = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_4_0.997_2016paperbinmerge1.root")
#p01_2016  = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_01_0.997_2016paperbinmerge1.root")

p0ratio_2016 = p0_2016.Get("g_ratio_with_effective_acceptance")
publishedratio = p0_2016.Get("publisehd_ratio")
p1ratio_2016 = p1_2016.Get("g_ratio_with_effective_acceptance")
p2ratio_2016 = p2_2016.Get("g_ratio_with_effective_acceptance")
p4ratio_2016 = p4_2016.Get("g_ratio_with_effective_acceptance")
#p01ratio_2016 = p01_2016.Get("g_ratio_with_effective_acceptance")

p0_full = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_0_0.997_pass7.8binmerge1.root")
p1_full = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_1_0.997_pass7.8binmerge1.root")
p2_full = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_2_0.997_pass7.8binmerge1.root")
p4_full = TFile("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/intermediate_4_0.997_pass7.8binmerge1.root")
#p01_full = TFile("intermediate_01_0.997_pass7.8binmerge1.root")

p0ratio_full = p0_full.Get("g_ratio_with_effective_acceptance")
p1ratio_full = p1_full.Get("g_ratio_with_effective_acceptance")
p2ratio_full = p2_full.Get("g_ratio_with_effective_acceptance")
p4ratio_full = p4_full.Get("g_ratio_with_effective_acceptance")
#p01ratio_full = p01_full.Get("g_ratio_with_effective_acceptance")

'''
c1 = TCanvas("c1","c1",1000,500)
p0ratio_2016.Draw('AP')
p1ratio_2016.Draw("P same")
publishedratio.Draw("P same")
p0ratio_2016.SetMarkerColor(2);
p1ratio_2016.SetMarkerColor(4);
publishedratio.SetMarkerColor(1);
p0ratio_2016.SetMarkerStyle(8);
p1ratio_2016.SetMarkerStyle(8);
publishedratio.SetMarkerStyle(8);
p0ratio_2016.GetYaxis().SetRangeUser(0.00005,0.00022);
p0ratio_2016.GetXaxis().SetTitle("Rigidity (GV)");
p0ratio_2016.GetYaxis().SetTitle("Antiproton ratio");
p0ratio_2016.SetTitle("")
leg =TLegend(0.5,0.2,0.9,0.4)
leg.SetFillColor(0)
leg.AddEntry(p0ratio_2016, "Pattern 0 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(p1ratio_2016, "Pattern 1 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(publishedratio, "PRL paper in 2016","lp")
gStyle.SetLegendTextSize(0.03);
leg.Draw()
c1.Update()
c1.SaveAs("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/Ratio_2016_p0_p1.pdf")

c2 = TCanvas("c2","c2",1000,500)
p0ratio_full.Draw('AP')
p1ratio_full.Draw("P same")
publishedratio.Draw("P same")
p0ratio_full.SetMarkerColor(2);
p1ratio_full.SetMarkerColor(4);
publishedratio.SetMarkerColor(1);
p0ratio_full.SetMarkerStyle(8);
p1ratio_full.SetMarkerStyle(8);
publishedratio.SetMarkerStyle(8);
p0ratio_full.GetYaxis().SetRangeUser(0.00005,0.00022);
p0ratio_full.GetXaxis().SetTitle("Rigidity (GV)");
p0ratio_full.GetYaxis().SetTitle("Antiproton ratio");
p0ratio_full.SetTitle("")
leg =TLegend(0.5,0.2,0.9,0.4)
leg.SetFillColor(0)
leg.AddEntry(p0ratio_full, "Pattern 0 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(p1ratio_full, "Pattern 1 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(publishedratio, "PRL paper in 2016","lp")
gStyle.SetLegendTextSize(0.03);
leg.Draw()
c2.Update()
c2.SaveAs("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/Ratio_full_p0_p1.pdf")



c3 = TCanvas("c3","c3",1000,500)
p2ratio_full.Draw('AP')
p4ratio_full.Draw("P same")
publishedratio.Draw("P same")
p2ratio_full.SetMarkerColor(2);
p4ratio_full.SetMarkerColor(4);
publishedratio.SetMarkerColor(1);
p2ratio_full.SetMarkerStyle(8);
p4ratio_full.SetMarkerStyle(8);
publishedratio.SetMarkerStyle(8);
p2ratio_full.GetYaxis().SetRangeUser(0.00005,0.00022);
p2ratio_full.GetXaxis().SetTitle("Rigidity (GV)");
p2ratio_full.GetYaxis().SetTitle("Antiproton ratio");
p2ratio_full.SetTitle("")
leg =TLegend(0.5,0.2,0.9,0.4)
leg.SetFillColor(0)
leg.AddEntry(p2ratio_full, "Pattern 2 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(p4ratio_full, "Pattern 4 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(publishedratio, "PRL paper in 2016","lp")
gStyle.SetLegendTextSize(0.03);
leg.Draw()
c3.Update()
c3.SaveAs("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/Ratio_full_p2_p4.pdf")
'''

c4 = TCanvas("c4","c4",1000,500)
p2ratio_full.Draw('AP')
p4ratio_full.Draw("P same")
p0ratio_full.Draw("P same")
p1ratio_full.Draw("P same")
publishedratio.Draw("P same")
p2ratio_full.SetMarkerColor(2);
p4ratio_full.SetMarkerColor(4);
publishedratio.SetMarkerColor(1);
p0ratio_full.SetMarkerColor(3);
p1ratio_full.SetMarkerColor(6);
p2ratio_full.SetMarkerStyle(8);
p4ratio_full.SetMarkerStyle(8);
publishedratio.SetMarkerStyle(8);
p0ratio_full.SetMarkerStyle(8);
p1ratio_full.SetMarkerStyle(8);
p2ratio_full.GetYaxis().SetRangeUser(0.00005,0.00022);
p2ratio_full.GetXaxis().SetTitle("Rigidity (GV)");
p2ratio_full.GetYaxis().SetTitle("Antiproton ratio");
p2ratio_full.SetTitle("")
leg =TLegend(0.5,0.2,0.9,0.4)
leg.SetFillColor(0)
leg.AddEntry(p2ratio_full, "Pattern 2 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(p4ratio_full, "Pattern 4 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(p0ratio_full, "Pattern 0 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(p1ratio_full, "Pattern 1 (Full Range) (Statistic Error only) ","lp")
leg.AddEntry(publishedratio, "PRL paper in 2016","lp")
gStyle.SetLegendTextSize(0.03);
leg.Draw()
c4.Update()
c4.SaveAs("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/Ratio_full_all.pdf")





c5 = TCanvas("c5","c5",1000,500)
p2ratio_2016.Draw('AP')
p4ratio_2016.Draw("P same")
p0ratio_2016.Draw("P same")
p1ratio_2016.Draw("P same")
publishedratio.Draw("P same")
p2ratio_2016.SetMarkerColor(2);
p4ratio_2016.SetMarkerColor(4);
publishedratio.SetMarkerColor(1);
p0ratio_2016.SetMarkerColor(3);
p1ratio_2016.SetMarkerColor(6);
p2ratio_2016.SetMarkerStyle(8);
p4ratio_2016.SetMarkerStyle(8);
publishedratio.SetMarkerStyle(8);
p0ratio_2016.SetMarkerStyle(8);
p1ratio_2016.SetMarkerStyle(8);
p2ratio_2016.GetYaxis().SetRangeUser(0.00005,0.00022);
p2ratio_2016.GetXaxis().SetTitle("Rigidity (GV)");
p2ratio_2016.GetYaxis().SetTitle("Antiproton ratio");
p2ratio_2016.SetTitle("")
leg =TLegend(0.5,0.2,0.9,0.4)
leg.SetFillColor(0)
leg.AddEntry(p2ratio_2016, "Pattern 2 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(p4ratio_2016, "Pattern 4 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(p0ratio_2016, "Pattern 0 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(p1ratio_2016, "Pattern 1 (Published Range) (Statistic Error only) ","lp")
leg.AddEntry(publishedratio, "PRL paper in 2016","lp")
gStyle.SetLegendTextSize(0.03);
leg.Draw()
c5.Update()
c5.SaveAs("$HPCINTERMEDIATEDIR/total/Time_Averaged_ratio/Ratio_2016_all.pdf")






