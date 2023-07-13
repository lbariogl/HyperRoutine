import numpy as np
import uproot
import ROOT
import sys
sys.path.append('../utils')
import utils as utils

ROOT.gROOT.SetBatch(True)


def create_pt_shift_histo(df):
    h2MomResoVsPtHe3 = ROOT.TH2F("h2MomResoVsPtHe3", ";^{3}He #it{p}_{T} (GeV/#it{c});  ^{3}He #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c})", 50, 1., 5, 50, -0.4, 0.4)
    df.eval("PtResHe3 = (fPtHe3- fGenPtHe3)", inplace=True)
    utils.fill_th2_hist(h2MomResoVsPtHe3, df, "fPtHe3", "PtResHe3")
    h2MomResoVsPtHe3.FitSlicesY()
    h2MomResoVsPtHe3_mean = ROOT.gDirectory.Get("h2MomResoVsPtHe3_1")
    return h2MomResoVsPtHe3, h2MomResoVsPtHe3_mean

def create_derived_columns(df):
    df["fPxHe3"] = df["fPtHe3"] * np.cos(df["fPhiHe3"])
    df["fPyHe3"] = df["fPtHe3"] * np.sin(df["fPhiHe3"])
    df["fPzHe3"] = df["fPtHe3"] * np.sinh(df["fEtaHe3"])
    df["fPxPi"] = df["fPtPi"] * np.cos(df["fPhiPi"])
    df["fPyPi"] = df["fPtPi"] * np.sin(df["fPhiPi"])
    df["fPzPi"] = df["fPtPi"] * np.sinh(df["fEtaPi"])
    df["fPt"] = np.sqrt((df["fPxHe3"] + df["fPxPi"])**2 + (df["fPyHe3"] + df["fPyPi"])**2)
    df["fPx"] = df["fPxHe3"] + df["fPxPi"]
    df["fPy"] = df["fPyHe3"] + df["fPyPi"]
    df["fPz"] = df["fPzHe3"] + df["fPzPi"]
    df["fP"] = np.sqrt(df["fPx"]**2 + df["fPy"]**2 + df["fPz"]**2)
    df["fDecayL"] = np.sqrt(df["fXDecVtx"]**2 + df["fYDecVtx"]**2 + df["fZDecVtx"]**2)
    df["fCosPA"] = (df["fPx"]*df["fXDecVtx"] + df["fPy"]*df["fYDecVtx"] + df["fPz"]*df["fZDecVtx"]) / (df["fP"]*df["fDecayL"])
    df["fEnHe3"] = np.sqrt(df["fPxHe3"]**2 + df["fPyHe3"]**2 + df["fPzHe3"]**2 + 2.8083916**2)
    df["fEnPi"] = np.sqrt(df["fPxPi"]**2 + df["fPyPi"]**2 + df["fPzPi"]**2 + 0.139570**2)
    df["fMassH3L"] = np.sqrt((df["fEnHe3"] + df["fEnPi"])**2 - (df["fPxHe3"] + df["fPxPi"])**2 - (df["fPyHe3"] + df["fPyPi"])**2 - (df["fPzHe3"] + df["fPzPi"])**2)


df1 = uproot.open("/data/shared/hyp_run_3/mc/AO2D_new_task.root")["O2mchypcands"].arrays(library="pd")
df2 = df1.copy(deep=True)
df1.query("fIsReco > 0", inplace=True)
df2.query("fIsReco > 0", inplace=True)

create_derived_columns(df1)

### loop over df2 and shift the pT according to the resolution
h2MomResoVsPtHe3, h_pt_shift = create_pt_shift_histo(df2)
cloned_pt_arr = np.array(df2["fPtHe3"])
for i in range(len(cloned_pt_arr)):
    pt_shift = h_pt_shift.GetBinContent(h_pt_shift.FindBin(cloned_pt_arr[i]))
    cloned_pt_arr[i] -= pt_shift
df2["fPtHe3"] = cloned_pt_arr
create_derived_columns(df2)


h_df1_mass_pt = ROOT.TH2D("h_df1_mass_pt", "; M (^{3}He + #pi) (GeV/#it{c}^{2}); #it{p}_{T}^{gen}", 10, 1, 5, 30, 2.96, 3.04)
h_df2_mass_pt = ROOT.TH2D("h_df2_mass_pt", "; M (^{3}He + #pi) (GeV/#it{c}^{2}); #it{p}_{T}^{gen}", 10, 1, 5, 30, 2.96, 3.04)

utils.fill_th2_hist(h_df1_mass_pt, df1, "fGenPt", "fMassH3L")
utils.fill_th2_hist(h_df2_mass_pt, df2, "fGenPt", "fMassH3L")
h_df2_pt_reso = ROOT.TH2D("h_df2_pt_reso", "h2_pt_reso", 50, 1, 5, 50, -0.4, 0.4)
df2["PtReso"] = (df2["fPtHe3"] - df2["fGenPtHe3"]) / df2["fGenPtHe3"]
utils.fill_th2_hist(h_df2_pt_reso, df2, "fGenPtHe3", "PtReso")


outf = ROOT.TFile("../../results/pt_mass_resolutions.root", "RECREATE")
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
