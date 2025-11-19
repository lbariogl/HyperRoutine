import ROOT
import numpy as np

## val, stat, syst

#new_yield = [2.1e-08, 0.3e-08, 0.1e-08]
#trd_yield = [2.1e-08, 0.6e-08, 0.4e-08]

#define kBlueC TColor::GetColor("#2077b4")
#define kAzureC TColor::GetColor("#18becf")
#define kGreenC TColor::GetColor("#2ba02b")
#kOrangeC  = ROOT.TColor.GetColor("#ff7f0f");

kOrangeC = ROOT.TColor.GetColor("#ff7f0f")
kGreenC = ROOT.TColor.GetColor('#2ca02c')
kBlueC = ROOT.TColor.GetColor("#2077b4")
kAzureC = ROOT.TColor.GetColor("#18becf")


lambda_yield = 0.2502 / 2
lambda_yield_hm = 1.147 / 2

x_trd = np.array([6.89], dtype=np.float64)
x_hm = np.array([30.81], dtype=np.float64)

x_new = np.array([7.12], dtype=np.float64)
x_err_stat = np.array([0.], dtype=np.float64)
x_err = np.array([0.1], dtype=np.float64)
x_err_hm = np.array([0.44], dtype=np.float64)


yield_new = np.array([1.4e-08 / lambda_yield], dtype=np.float64)
yield_trd = np.array([2.1e-08 / lambda_yield], dtype=np.float64)
yield_hm = np.array([2.4e-07 / lambda_yield_hm], dtype=np.float64)

print(yield_new, yield_trd, yield_hm)

yield_new_stat = np.array([0.1e-08 / lambda_yield], dtype=np.float64)
yield_trd_stat = np.array([0.7e-08 / lambda_yield], dtype=np.float64)
yield_hm_stat = np.array([0.5e-07 / lambda_yield_hm], dtype=np.float64)

yield_new_syst = np.array([0.3e-08 / lambda_yield], dtype=np.float64)
yield_trd_syst = np.array([0.4e-08 / lambda_yield], dtype=np.float64)
yield_hm_syst = np.array([0.3e-07 / lambda_yield_hm], dtype=np.float64)


## add values in tgrapherrors and plot

gr_stat_new = ROOT.TGraphErrors(len(x_new), x_new, yield_new, x_err_stat, yield_new_stat)
gr_stat_trd = ROOT.TGraphErrors(len(x_trd), x_trd, yield_trd, x_err_stat, yield_trd_stat)

gr_syst_new = ROOT.TGraphErrors(len(x_new), x_new, yield_new, x_err, yield_new_syst)
gr_syst_trd = ROOT.TGraphErrors(len(x_trd), x_trd, yield_trd, x_err, yield_trd_syst)

gr_stat_hm = ROOT.TGraphErrors(len(x_hm), x_hm, yield_hm, x_err_hm, yield_hm_stat)
gr_syst_hm = ROOT.TGraphErrors(len(x_hm), x_hm, yield_hm, x_err_hm, yield_hm_syst)

gr_stat_new.SetMarkerColor(ROOT.kRed)
gr_stat_new.SetLineColor(ROOT.kRed)
gr_stat_new.SetMarkerSize(2)
gr_stat_new.SetMarkerStyle(ROOT.kFullDiamond)
gr_stat_new.SetLineWidth(1)

gr_syst_new.SetMarkerColor(ROOT.kRed)
gr_syst_new.SetLineColor(ROOT.kRed)
gr_syst_new.SetFillStyle(0)
gr_syst_new.SetMarkerSize(2)
gr_syst_new.SetMarkerStyle(ROOT.kFullDiamond)
gr_syst_new.SetLineWidth(1)

gr_stat_trd.SetMarkerColor(kOrangeC)
gr_stat_trd.SetLineColor(kOrangeC)
gr_stat_trd.SetMarkerSize(2)
gr_stat_trd.SetMarkerStyle(ROOT.kStar)
gr_stat_trd.SetLineWidth(1)


gr_syst_trd.SetMarkerColor(kOrangeC)
gr_syst_trd.SetLineColor(kOrangeC)
gr_syst_trd.SetMarkerSize(2)
gr_syst_trd.SetMarkerStyle(ROOT.kStar)
gr_syst_trd.SetLineWidth(1)
gr_syst_trd.SetFillStyle(0)

gr_stat_hm.SetMarkerColor(kOrangeC)
gr_stat_hm.SetLineColor(kOrangeC)
gr_stat_hm.SetMarkerSize(2)
gr_stat_hm.SetMarkerStyle(ROOT.kStar)
gr_stat_hm.SetLineWidth(1)

gr_syst_hm.SetMarkerColor(kOrangeC)
gr_syst_hm.SetLineColor(kOrangeC)
gr_syst_hm.SetMarkerSize(2)
gr_syst_hm.SetMarkerStyle(ROOT.kStar)
gr_syst_hm.SetLineWidth(1)
gr_syst_hm.SetFillStyle(0)


## adding SHM models
hp_ratio_csm_1 = ROOT.TGraphErrors("../utils/models/VanillaCSM.Yields.Vc.eq.1.6dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")
hp_ratio_csm_3 = ROOT.TGraphErrors("../utils/models/VanillaCSM.Yields.Vc.eq.3dVdy.dat","%*s %*s %*s %lg %*s %*s %*s %*s %lg %*s")

hp_ratio_csm_1.SetLineColor(922)
hp_ratio_csm_1.SetLineWidth(2)
hp_ratio_csm_1.SetTitle("SHM, #it{Vc} = 1.6 d#it{V}/d#it{y}")
hp_ratio_csm_1.SetMarkerSize(0)

hp_ratio_csm_3.SetLineColor(922)
hp_ratio_csm_3.SetLineWidth(2)
hp_ratio_csm_3.SetLineStyle(2)
hp_ratio_csm_3.SetMarkerSize(0)
hp_ratio_csm_3.SetTitle("SHM, #it{Vc} = 3d#it{V}/d#it{y}")


## adding coalescence models
hp_2body = ROOT.TGraphErrors("../utils/models/coal2b.csv","%lg %lg %lg")
# hp_2body.SetFillStyle(3014)
hp_2body.SetMarkerSize(0)
hp_2body.SetLineWidth(2)
hp_2body.SetLineColor(kAzureC)
hp_2body.SetMarkerColor(kAzureC)
hp_2body.SetFillColorAlpha(kAzureC, 0.3)

hp_3body = ROOT.TGraphErrors("../utils/models/coal3b.csv","%lg %lg %lg")
# hp_3body.SetFillStyle(3014)
hp_3body.SetMarkerSize(0)
hp_3body.SetLineWidth(2)
hp_3body.SetLineColor(kAzureC)
hp_3body.SetMarkerColor(kAzureC)
hp_3body.SetFillColorAlpha(kAzureC, 0.3)




cv_out = ROOT.TCanvas("cv", "cv")
frame = cv_out.DrawFrame(3, 2e-08, 80, 0.6e-06)
cv_out.SetLogx()
# cv_out.SetLogy()
frame.GetXaxis().SetTitleOffset(0.8)
frame.GetYaxis().SetTitleOffset(0.9)
frame.GetYaxis().SetTitleSize(0.06)
frame.GetXaxis().SetTitleSize(0.06)
frame.GetYaxis().SetLabelSize(0.04)
frame.GetXaxis().SetLabelSize(0.04)

frame.GetXaxis().SetTitle("#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}")
frame.GetYaxis().SetTitle("{}_{#Lambda}^{3}H/#Lambda")
hp_2body.Draw('3 same')

gr_stat_new.Draw('Pz')
gr_syst_new.Draw('P2')

# gr_stat_trd.Draw('Pz')
# gr_syst_trd.Draw('P2')

# gr_stat_hm.Draw('Pz')
# gr_syst_hm.Draw('P2')

hp_ratio_csm_1.Draw('L same')

leg = ROOT.TLegend(0.57, 0.3, 0.94, 0.47)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetHeader("ALICE Preliminary")
## reduce space between symbols and text
leg.SetMargin(0.1)
leg.AddEntry(gr_stat_new, "Run 3, pp #sqrt{#it{s}} = 13.6 TeV", "P")
# leg.AddEntry(gr_stat_trd, "Run 2, pp #sqrt{#it{s}} = 13 TeV", "P")

leg_theory = ROOT.TLegend(0.57, 0.15, 0.94, 0.26)
leg_theory.AddEntry(hp_ratio_csm_1, "SHM, #it{Vc} = 1.6 d#it{V}/d#it{y}", "L")
leg_theory.AddEntry(hp_2body, "Coalescence, 2-body", "F")
leg.Draw()
leg_theory.Draw()




outfile = ROOT.TFile("yield.root", "RECREATE")
cv_out.Write()
outfile.Close()

cv_out.SaveAs("yield.pdf")
cv_out.SaveAs("yield.png")

