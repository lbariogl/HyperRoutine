import numpy as np
import uproot
import ROOT
import sys
sys.path.append('../utils')
import utils as utils

ROOT.gROOT.SetBatch(True)

df1 = uproot.open("/data/shared/hyp_run_3/mc/AO2D_MC.root")["O2mchypcands"].arrays(library="pd")
df2 = df1.copy(deep=True)
df1.query("fIsReco > 0", inplace=True)
df2.query("fIsReco > 0", inplace=True)

utils.correct_and_convert_df(df1)

h2MomResoVsPtHe3, h_pt_shift = utils.create_pt_shift_histo(df2)
utils.correct_and_convert_df(df2, h_pt_shift)

h_df1_mass_pt = ROOT.TH2D("h_df1_mass_pt", "; M (^{3}He + #pi) (GeV/#it{c}^{2}); #it{p}_{T}^{gen}", 10, 1, 5, 30, 2.96, 3.04)
h_df2_mass_pt = ROOT.TH2D("h_df2_mass_pt", "; M (^{3}He + #pi) (GeV/#it{c}^{2}); #it{p}_{T}^{gen}", 10, 1, 5, 30, 2.96, 3.04)

utils.fill_th2_hist(h_df1_mass_pt, df1, "fGenPt", "fMassH3L")
utils.fill_th2_hist(h_df2_mass_pt, df2, "fGenPt", "fMassH3L")
h_df2_pt_reso = ROOT.TH2D("h_df2_pt_reso", "h2_pt_reso", 50, 1, 5, 50, -0.4, 0.4)
df2["PtReso"] = (df2["fPtHe3"] - df2["fGenPtHe3"]) / df2["fGenPtHe3"]
utils.fill_th2_hist(h_df2_pt_reso, df2, "fGenPtHe3", "PtReso")


outf = ROOT.TFile("../utils/he3_pt_calibration.root", "RECREATE")
h2MomResoVsPtHe3.Write()
h_pt_shift.Write()
h_df2_pt_reso.Write()
h_df2_mass_pt.Write()

## fit the invariant mass slices with a gaussian, save the sigma
h_df1_mass_pt.FitSlicesY()
h_df1_mass_pt_sigma = ROOT.gDirectory.Get("h_df1_mass_pt_2")
h_df2_mass_pt.FitSlicesY()
h_df2_mass_pt_sigma = ROOT.gDirectory.Get("h_df2_mass_pt_2")
## Transform the sigma in MeV
for i in range(1, h_df1_mass_pt_sigma.GetNbinsX()+1):
    h_df1_mass_pt_sigma.SetBinContent(i, h_df1_mass_pt_sigma.GetBinContent(i)*1000)
    h_df2_mass_pt_sigma.SetBinContent(i, h_df2_mass_pt_sigma.GetBinContent(i)*1000)
    h_df1_mass_pt_sigma.SetBinError(i, h_df1_mass_pt_sigma.GetBinError(i)*1000)
    h_df2_mass_pt_sigma.SetBinError(i, h_df2_mass_pt_sigma.GetBinError(i)*1000)
## plot in the same canvas the sigma vs pt
c_sigma = ROOT.TCanvas("c_sigma", "", 800, 600)
h_df2_mass_pt_sigma.SetLineColor(ROOT.kBlue)
h_df1_mass_pt_sigma.SetLineColor(ROOT.kRed)
h_df2_mass_pt_sigma.SetMarkerColor(ROOT.kBlue)
h_df1_mass_pt_sigma.SetMarkerColor(ROOT.kRed)
h_df2_mass_pt_sigma.SetMarkerStyle(20)
h_df1_mass_pt_sigma.SetMarkerStyle(20)
h_df2_mass_pt_sigma.SetMarkerSize(1.2)
h_df1_mass_pt_sigma.SetMarkerSize(1.2)
h_df2_mass_pt_sigma.SetStats(0)
h_df1_mass_pt_sigma.SetStats(0)
h_df1_mass_pt_sigma.GetXaxis().SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})")
h_df1_mass_pt_sigma.GetYaxis().SetTitle("#sigma (MeV/#it{c}^{2})")
h_df1_mass_pt_sigma.GetXaxis().SetTitleSize(0.05)
h_df1_mass_pt_sigma.GetYaxis().SetTitleSize(0.05)
h_df1_mass_pt_sigma.SetTitle("")
h_df1_mass_pt_sigma.Draw("pe")
h_df2_mass_pt_sigma.Draw("pe same")
leg_sigma = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
leg_sigma.AddEntry(h_df1_mass_pt_sigma, "Default", "pe")
leg_sigma.AddEntry(h_df2_mass_pt_sigma, "W/ ^{3}He #it{p}_{T} calibration", "pe")
leg_sigma.Draw("same")
c_sigma.Write()



### select invariant masses with genPt < 2 GeV/c and fit them with a gaussian
df1_lowpt = df1.query("abs(fGenPt) < 2")
df2_lowpt = df2.query("abs(fGenPt) < 2")
h_df1_lowpt_mass = ROOT.TH1D("h_df1_lowpt_mass", ";M (^{3}He + #pi) (GeV/#it{c}^{2}); Counts", 40, 2.97, 3.02)
h_df2_lowpt_mass = ROOT.TH1D("h_df2_lowpt_mass", ";M (^{3}He + #pi) (GeV/#it{c}^{2}); Counts", 40, 2.97, 3.02)
h_df1_lowpt_mass.SetLineColor(ROOT.kRed)
utils.fill_th1_hist(h_df1_lowpt_mass, df1_lowpt, "fMassH3L")
utils.fill_th1_hist(h_df2_lowpt_mass, df2_lowpt, "fMassH3L")
h_df1_lowpt_mass.Fit("gaus", "S", "", 2.9775, 3.002)
h_df2_lowpt_mass.Fit("gaus", "S", "", 2.983, 2.9948)
##change color of the fit functions
h_df1_lowpt_mass.GetFunction("gaus").SetLineColor(ROOT.kRed)
h_df2_lowpt_mass.GetFunction("gaus").SetLineColor(ROOT.kBlue)
h_df1_lowpt_mass.SetStats(0)
h_df2_lowpt_mass.SetStats(0)
h_df1_lowpt_mass.SetMarkerStyle(20)
h_df2_lowpt_mass.SetMarkerStyle(20)
h_df1_lowpt_mass.SetMarkerSize(1.2)
h_df2_lowpt_mass.SetMarkerSize(1.2)
h_df1_lowpt_mass.SetMarkerColor(ROOT.kRed)
h_df2_lowpt_mass.SetMarkerColor(ROOT.kBlue)
## add legend
leg = ROOT.TLegend(0.6, 0.6, 0.89, 0.89)
leg.AddEntry(h_df1_lowpt_mass, "Default", "pe")
leg.AddEntry(h_df2_lowpt_mass, "W/ ^{3}He #it{p}_{T} calibration", "pe")
h_df2_lowpt_mass.GetXaxis().SetTitle("M (^{3}He + #pi) (GeV/#it{c}^{2})")
h_df2_lowpt_mass.GetYaxis().SetTitle("Entries")
h_df2_lowpt_mass.GetXaxis().SetTitleSize(0.05)
h_df2_lowpt_mass.GetYaxis().SetTitleSize(0.05)
## Add TPaveText
pt = ROOT.TPaveText(0.6, 0.4, 0.89, 0.6, "NDC")
pt.AddText("MC Simulation")
pt.AddText("#it{p}_{T}^{gen} < 2 GeV/#it{c}")
pt.SetBorderSize(0)
pt.SetFillStyle(0)
pt.SetTextFont(42)
pt.SetTextAlign(12)
pt_width = ROOT.TPaveText(0.6, 0.2, 0.89, 0.4, "NDC")
##add widths of the gaussian fits
sigma1 = h_df1_lowpt_mass.GetFunction("gaus").GetParameter(2)*1000
sigma2 = h_df2_lowpt_mass.GetFunction("gaus").GetParameter(2)*1000
t1 = pt_width.AddText(f"#sigma = {sigma1:.1f}" + " MeV/#it{c^{2}}")
t1.SetTextColor(ROOT.kRed)
t2 = pt_width.AddText("#sigma = " + f"{sigma2:.1f}" + " MeV/#it{c^{2}}")
t2.SetTextColor(ROOT.kBlue)
pt_width.SetBorderSize(0)
pt_width.SetFillStyle(0)
pt_width.SetTextFont(42)
pt_width.SetTextAlign(12)
## plot in the same canvas the two mass distributions
c_slice = ROOT.TCanvas("c_slice", "", 800, 600)
h_df2_lowpt_mass.Draw("pe")
h_df1_lowpt_mass.Draw("pe same")
leg.Draw("same")
pt.Draw("same")
pt_width.Draw("same")
c_slice.Write()


outf.Close()
