#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad
import argparse
import os


def main():

    for binmerge in binmergeAll:

        f_2016paper = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_"+ arguments.trackerpattern +"_" + arguments.richcut + "_" + "2016paper" + "binmerge" + binmerge + ".root")
        f_pass78 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_"+ arguments.trackerpattern +"_" + arguments.richcut + "_" + "pass7.8" + "binmerge" + binmerge + ".root")

        g_ratio_with_effective_acceptance_pass78 = f_pass78.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_2016paper = f_2016paper.Get("g_ratio_with_effective_acceptance")
        h_error_relative_pass78 = f_pass78.Get("h_error_relative")
        h_error_relative_2016paper = f_2016paper.Get("h_error_relative")

        '''
        f_2016paper_0 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_0_" + arguments.richcut + "_" + "2016paper" + ".root")
        f_2016paper_1 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_1_" + arguments.richcut + "_" + "2016paper" + ".root")
        f_2016paper_2 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_2_" + arguments.richcut + "_" + "2016paper" + ".root")
        f_2016paper_4 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_4_" + arguments.richcut + "_" + "2016paper" + ".root")
        f_2016paper_01 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_01_" + arguments.richcut + "_" + "2016paper" + ".root")
        f_pass78_0 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_0_" + arguments.richcut + "_" + "pass7.8" + ".root")
        f_pass78_1 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_1_" + arguments.richcut + "_" + "pass7.8" + ".root")
        f_pass78_2 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_2_" + arguments.richcut + "_" + "pass7.8" + ".root")
        f_pass78_4 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_4_" + arguments.richcut + "_" + "pass7.8" + ".root")
        f_pass78_01 = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/intermediate_01_" + arguments.richcut + "_" + "pass7.8" + ".root")
        '''

        '''
        g_ratio_with_effective_acceptance_2016paper_0 = f_2016paper_0.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_2016paper_1 = f_2016paper_1.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_2016paper_2 = f_2016paper_2.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_2016paper_4 = f_2016paper_4.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_2016paper_01 = f_2016paper_01.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_pass78_0 = f_pass78_0.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_pass78_1 = f_pass78_1.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_pass78_2 = f_pass78_2.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_pass78_4 = f_pass78_4.Get("g_ratio_with_effective_acceptance")
        g_ratio_with_effective_acceptance_pass78_01 = f_pass78_01.Get("g_ratio_with_effective_acceptance")
        '''

        publisehd_ratio = f_2016paper.Get("publisehd_ratio")
        h_Published_Statistic_Error_Relative = f_2016paper.Get("h_Published_Statistic_Error_Relative")
        h_Published_Systematic_Error_Relative = f_2016paper.Get("h_Published_Systematic_Error_Relative")

        c1 = TCanvas("c1","c1",1000,500)
        publisehd_ratio.Draw('AP')
        g_ratio_with_effective_acceptance_2016paper.Draw("P same")
        g_ratio_with_effective_acceptance_pass78.Draw("P same")
        publisehd_ratio.GetXaxis().SetTitle("Rigidity (GV)");
        publisehd_ratio.GetYaxis().SetTitle("Antiproton ratio");
        publisehd_ratio.SetTitle("")
        publisehd_ratio.GetXaxis().SetLabelSize(0.04);
        publisehd_ratio.GetXaxis().SetTitleSize(0.04);
        publisehd_ratio.GetYaxis().SetLabelSize(0.04);
        publisehd_ratio.GetYaxis().SetTitleSize(0.04);
        g_ratio_with_effective_acceptance_2016paper.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_2016paper.SetMarkerColor(2);
        g_ratio_with_effective_acceptance_pass78.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_pass78.SetMarkerColor(4);
        gPad.SetLogx()
        leg =TLegend(0.15,0.7,0.6,0.9)
        leg.SetFillColor(0)
        leg.AddEntry(g_ratio_with_effective_acceptance_2016paper, "Pattern" + arguments.trackerpattern +" (Publisehd Range) (Statistic Error only)","lp")
        leg.AddEntry(g_ratio_with_effective_acceptance_pass78, "Pattern" + arguments.trackerpattern + " (Full Range) (Statistic Error only)","lp")
        leg.AddEntry(publisehd_ratio, "AMS PRL paper in 2016","lp")
        gStyle.SetLegendTextSize(0.03);
        leg.Draw()
        c1.Update()
        c1.SaveAs(intermediateworkpath + "/total/Time_Averaged_ratio/plot/pattern" + arguments.trackerpattern + "_" + arguments.richcut + "_" + "binmerge" + binmerge + ".pdf")


        c1e = TCanvas("c1e","c1e",1000,500)
        h_Published_Systematic_Error_Relative.Draw('')
        h_Published_Statistic_Error_Relative.Draw('same')
        h_error_relative_2016paper.Draw('same')
        h_error_relative_pass78.Draw('same')
        h_Published_Systematic_Error_Relative.GetYaxis().SetLimits(0,5.6);
        h_Published_Systematic_Error_Relative.GetYaxis().SetRangeUser(0,5.6);
        h_Published_Systematic_Error_Relative.GetXaxis().SetTitle("Rigidity (GV)");
        h_Published_Systematic_Error_Relative.GetYaxis().SetTitle("Relative Error (%)");
        h_Published_Systematic_Error_Relative.SetTitle("")
        h_Published_Systematic_Error_Relative.GetXaxis().SetLabelSize(0.04);
        h_Published_Systematic_Error_Relative.GetXaxis().SetTitleSize(0.04);
        h_Published_Systematic_Error_Relative.GetYaxis().SetLabelSize(0.04);
        h_Published_Systematic_Error_Relative.GetYaxis().SetTitleSize(0.04);
        h_Published_Systematic_Error_Relative.SetFillColor(0);
        h_Published_Systematic_Error_Relative.SetLineColor(1);
        h_Published_Statistic_Error_Relative.SetLineColor(2);
        h_error_relative_2016paper.SetLineColor(3);
        h_error_relative_pass78.SetLineColor(4);
        gPad.SetLogx()
        h_Published_Systematic_Error_Relative.GetXaxis().SetMoreLogLabels();
        gStyle.SetOptStat("00000000")
        leg =TLegend(0.15,0.7,0.6,0.9)
        leg.SetFillColor(0)
        leg.AddEntry(h_Published_Statistic_Error_Relative, "Statistic_Error_Relative (AMS PRL paper in 2016)","lp")
        leg.AddEntry(h_Published_Systematic_Error_Relative, "Systematic_Error_Relative (AMS PRL paper in 2016)","lp")
        leg.AddEntry(h_error_relative_2016paper, "Statistic_Error_Relative (Pattern" + arguments.trackerpattern + ") (Publisehd Range)","lp")
        leg.AddEntry(h_error_relative_pass78, "Statistic_Error_Relative (Pattern" + arguments.trackerpattern + ") (Full Range)","lp")
        gStyle.SetLegendTextSize(0.03);
        leg.Draw()
        c1e.Update()
        c1e.SaveAs(intermediateworkpath + "/total/Time_Averaged_ratio/plot/error_pattern" + arguments.trackerpattern + "_" + arguments.richcut + "_" + "binmerge" + binmerge +".pdf")



        '''
        c2 = TCanvas("c2","c2",1000,500)
        publisehd_ratio.Draw('AP')
        g_ratio_with_effective_acceptance_2016paper_2.Draw("P same")
        g_ratio_with_effective_acceptance_2016paper_4.Draw("P same")
        g_ratio_with_effective_acceptance_pass78_2.Draw("P same")
        g_ratio_with_effective_acceptance_pass78_4.Draw("P same")
        publisehd_ratio.GetXaxis().SetTitle("Rigidity (GV)");
        publisehd_ratio.GetYaxis().SetTitle("Antiproton ratio");
        publisehd_ratio.SetTitle("")
        publisehd_ratio.GetXaxis().SetLabelSize(0.04);
        publisehd_ratio.GetXaxis().SetTitleSize(0.04);
        publisehd_ratio.GetYaxis().SetLabelSize(0.04);
        publisehd_ratio.GetYaxis().SetTitleSize(0.04);
        g_ratio_with_effective_acceptance_2016paper_2.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_2016paper_2.SetMarkerColor(2);
        g_ratio_with_effective_acceptance_2016paper_4.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_2016paper_4.SetMarkerColor(3);
        g_ratio_with_effective_acceptance_pass78_2.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_pass78_2.SetMarkerColor(4);
        g_ratio_with_effective_acceptance_pass78_4.SetMarkerStyle(8);
        g_ratio_with_effective_acceptance_pass78_4.SetMarkerColor(6);
        gPad.SetLogx()
        leg =TLegend(0.15,0.7,0.6,0.9)
        leg.SetFillColor(0)
        leg.AddEntry(g_ratio_with_effective_acceptance_2016paper_2, "Pattern2 (Publisehd Range) (Statistic Error only)","lp")
        leg.AddEntry(g_ratio_with_effective_acceptance_2016paper_4, "Pattern4 (Publisehd Range) (Statistic Error only)","lp")
        leg.AddEntry(g_ratio_with_effective_acceptance_pass78_2, "Pattern2 (Full Range) (Statistic Error only)","lp")
        leg.AddEntry(g_ratio_with_effective_acceptance_pass78_4, "Pattern4 (Full Range) (Statistic Error only)","lp")
        leg.AddEntry(publisehd_ratio, "AMS PRL paper in 2016","lp")
        gStyle.SetLegendTextSize(0.03);
        leg.Draw()
        c2.Update()
        c2.SaveAs(intermediateworkpath + "/total/Time_Averaged_ratio/plot/pattern" + "24" + "_" + arguments.richcut + ".pdf")
        '''


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-C','--cluster', help='which cluster you choose')
    parser.add_argument('-P','--trackerpattern', help='which trackerpattern you choose.')
    parser.add_argument('-R','--richcut', help='what richcut you choose.')
    arguments = parser.parse_args()

    binmergeAll = [1, 2]

    if arguments.cluster == "JUAMS":
        intermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
    elif arguments.cluster == "HPC":
        intermediateworkpath = os.getenv('HPCINTERMEDIATEDIR')

    main()


