

import ROOT


f_1x = ROOT.TFile("../../results/absorption/absorption_histos_x1.root")
f_15x = ROOT.TFile("../../results/absorption/absorption_histos_x1.5.root")
f_2x = ROOT.TFile("../../results/absorption/absorption_histos_x2.root")
f_5x = ROOT.TFile("../../results/absorption/absorption_histos_x5.root")

print(f_1x, f_15x, f_2x, f_5x)
h_1x_mat = f_1x.Get("h_abso_frac_pt_mat")
h_1x_mat.SetDirectory(0)
h_1x_mat.SetName("h_abso_frac_pt_mat_1x")

h_15x_mat = f_15x.Get("h_abso_frac_pt_mat")
h_15x_mat.SetDirectory(0)
h_15x_mat.SetName("h_abso_frac_pt_mat_15x")

h_2x_mat = f_2x.Get("h_abso_frac_pt_mat")
h_2x_mat.SetDirectory(0)
h_2x_mat.SetName("h_abso_frac_pt_mat_2x")

h_5x_mat = f_5x.Get("h_abso_frac_pt_mat")
h_5x_mat.SetDirectory(0)
h_5x_mat.SetName("h_abso_frac_pt_mat_5x")


## repeat for antimatter
h_1x_antimat = f_1x.Get("h_abso_frac_pt_antimat")
h_1x_antimat.SetDirectory(0)
h_1x_antimat.SetName("h_abso_frac_pt_mat_1x_antimat")

h_15x_antimat = f_15x.Get("h_abso_frac_pt_antimat")
h_15x_antimat.SetDirectory(0)
h_15x_antimat.SetName("h_abso_frac_pt_mat_15x_antimat")

h_2x_antimat = f_2x.Get("h_abso_frac_pt_antimat")
h_2x_antimat.SetDirectory(0)
h_2x_antimat.SetName("h_abso_frac_pt_mat_2x_antimat")

h_5x_antimat = f_5x.Get("h_abso_frac_pt_antimat")
h_5x_antimat.SetDirectory(0)
h_5x_antimat.SetName("h_abso_frac_pt_mat_5x_antimat")


## canvas 15x for mat and antimat corrections
c_15x = ROOT.TCanvas("c_15x", "c_15x", 800, 600)

h_15x_mat.SetLineColor(ROOT.kBlue)
h_15x_mat.SetMarkerColor(ROOT.kBlue)
h_15x_mat.SetMaximum(1.0)
h_15x_mat.SetMinimum(0.95)
h_15x_antimat.SetLineColor(ROOT.kRed)
h_15x_antimat.SetMarkerColor(ROOT.kRed)

h_15x_mat.Draw("pe")
h_15x_antimat.Draw("pe same")

leg = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
leg.AddEntry(h_15x_antimat, "1.5 #times #sigma(^{3}#bar{He})", "l")
leg.AddEntry(h_15x_mat, "1.5 #times #sigma(^{3}He)", "l")
leg.SetBorderSize(0)
leg.Draw()



## canvas 1,1.5, 2, 5 x for antimat corrections
c_all = ROOT.TCanvas("c_all", "c_all", 800, 600)

h_1x_antimat.SetLineColor(ROOT.kBlack)
h_1x_antimat.SetMarkerColor(ROOT.kBlack)
h_1x_antimat.SetMaximum(1.0)
h_1x_antimat.SetMinimum(0.8)
h_15x_antimat.SetLineColor(ROOT.kRed)
h_15x_antimat.SetMarkerColor(ROOT.kRed)
h_2x_antimat.SetLineColor(ROOT.kBlue)
h_2x_antimat.SetMarkerColor(ROOT.kBlue)
h_5x_antimat.SetLineColor(ROOT.kGreen)
h_5x_antimat.SetMarkerColor(ROOT.kGreen)

h_1x_antimat.Draw("pe")
h_15x_antimat.Draw("pe same")
h_2x_antimat.Draw("pe same")
h_5x_antimat.Draw("pe same")

leg_all = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
leg_all.AddEntry(h_1x_antimat, "1 #times #sigma(^{3}#bar{He})", "l")
leg_all.AddEntry(h_15x_antimat, "1.5 #times #sigma(^{3}#bar{He})", "l")
leg_all.AddEntry(h_2x_antimat, "2 #times #sigma(^{3}#bar{He})", "l")
leg_all.AddEntry(h_5x_antimat, "5 #times #sigma(^{3}#bar{He})", "l")
leg_all.SetBorderSize(0)
leg_all.Draw()






outfile = ROOT.TFile("absorption_ratios.root", "RECREATE")
c_15x.Write()
c_all.Write()
