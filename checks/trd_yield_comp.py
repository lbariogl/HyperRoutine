import ROOT
import numpy as np

## val, stat, syst

#new_yield = [2.1e-08, 0.3e-08, 0.1e-08]
#trd_yield = [2.1e-08, 0.6e-08, 0.4e-08]

x_trd = np.array([0.5], dtype=np.float64)
x_new = np.array([1.5], dtype=np.float64)
x_err_stat = np.array([0.], dtype=np.float64)
x_err = np.array([0.1], dtype=np.float64)


yield_new = np.array([1.45e-08], dtype=np.float64)
yield_trd = np.array([2.1e-08], dtype=np.float64)

yield_new_stat = np.array([0.1e-08], dtype=np.float64)
yield_trd_stat = np.array([0.7e-08], dtype=np.float64)

yield_new_syst = np.array([0.3e-08], dtype=np.float64)
yield_trd_syst = np.array([0.4e-08], dtype=np.float64)


## add values in tgrapherrors and plot

gr_stat_new = ROOT.TGraphErrors(len(x_new), x_new, yield_new, x_err_stat, yield_new_stat)
gr_stat_trd = ROOT.TGraphErrors(len(x_trd), x_trd, yield_trd, x_err_stat, yield_trd_stat)

gr_syst_new = ROOT.TGraphErrors(len(x_new), x_new, yield_new, x_err, yield_new_syst)
gr_syst_trd = ROOT.TGraphErrors(len(x_trd), x_trd, yield_trd, x_err, yield_trd_syst)

gr_stat_new.SetMarkerColor(ROOT.kAzure + 2)
gr_stat_new.SetLineColor(ROOT.kAzure + 2)
gr_stat_new.SetMarkerSize(2)
gr_stat_new.SetMarkerStyle(ROOT.kOpenDiamond)
gr_stat_new.SetLineWidth(1)

gr_syst_new.SetMarkerColor(ROOT.kAzure + 2)
gr_syst_new.SetLineColor(ROOT.kAzure + 2)
gr_syst_new.SetFillStyle(0)
gr_syst_new.SetMarkerSize(2)
gr_syst_new.SetMarkerStyle(ROOT.kOpenDiamond)
gr_syst_new.SetLineWidth(1)

gr_stat_trd.SetMarkerColor(ROOT.kRed)
gr_stat_trd.SetLineColor(ROOT.kRed)
gr_stat_trd.SetMarkerSize(2)
gr_stat_trd.SetMarkerStyle(ROOT.kFullDiamond)
gr_stat_trd.SetLineWidth(1)


gr_syst_trd.SetMarkerColor(ROOT.kRed)
gr_syst_trd.SetLineColor(ROOT.kRed)
gr_syst_trd.SetMarkerSize(2)
gr_syst_trd.SetMarkerStyle(ROOT.kFullDiamond)
gr_syst_trd.SetLineWidth(1)
gr_syst_trd.SetFillStyle(0)

gr_stat_new.GetYaxis().SetTitle("#frac{1}{N_{ev}}#frac{#it{d}N}{#it{d}y} (GeV/#it{c})^{-1}")



cv_out = ROOT.TCanvas("cv", "cv", 700,700)
frame = cv_out.DrawFrame(0, 0.5e-08, 2, 2.9e-8)
frame.GetXaxis().SetLabelSize(0.)
frame.GetYaxis().SetTitle("#frac{1}{N_{ev}}#frac{#it{d}N}{#it{d}y} (GeV/#it{c})^{-1}")

gr_stat_new.Draw('Pz')
gr_syst_new.Draw('P2')

gr_stat_trd.Draw('Pz')
gr_syst_trd.Draw('P2')

leg = ROOT.TLegend(0.6, 0.6, 0.9, 0.9)
leg.AddEntry(gr_stat_new, "Run 3, MB", "P")
leg.AddEntry(gr_stat_trd, "Run 2, TRD trigger", "P")
leg.Draw()




outfile = ROOT.TFile("yield.root", "RECREATE")
cv_out.Write()