import ROOT
import numpy as np


kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor('#ff7f00')

hyp_file = ROOT.TFile("../results/2024_pt_analysis/pt_analysis_antimat_merged_ivan.root")
hyp_hist_stat = hyp_file.Get("std/hStat")
hyp_hist_syst = hyp_file.Get("std/hSyst")

he3_file = ROOT.TFile("../../non-prompt-nuclei/utils/he3_spectrum.root")
he3_hist_stat = he3_file.Get("fStatTPCA")
he3_hist_syst = he3_file.Get("fSystTPCA")

## first do hyp/he ratio
ratio_stat = hyp_hist_stat.Clone("ratio_stat")
## remove attached fit function
ratio_stat.GetListOfFunctions().Clear()

ratio_syst = hyp_hist_syst.Clone("ratio_syst")
ratio_stat.SetMarkerColor(ROOT.kAzure + 2)
ratio_stat.SetLineColor(ROOT.kAzure + 2)
ratio_syst.SetMarkerColor(ROOT.kAzure + 2)
ratio_syst.SetLineColor(ROOT.kAzure + 2)

ratio_stat.GetYaxis().SetTitle("{}^{3}_{#bar{#Lambda}}#bar{H} / ^{3}#bar{He}")

## sum in quadrature stat and syst errors

for i in range(1, ratio_stat.GetNbinsX()+1):
    hyp_val = hyp_hist_stat.GetBinContent(i)
    he3_val = he3_hist_stat.GetBinContent(i + 1)  # he3 histogram is shifted by 1 bin

    hyp_stat_err = hyp_hist_stat.GetBinError(i)
    he3_stat_err = he3_hist_stat.GetBinError(i + 1)

    hyp_syst_err = hyp_hist_syst.GetBinError(i)
    he3_syst_err = he3_hist_syst.GetBinError(i + 1)

    ratio_val = hyp_val / he3_val
    ratio_stat_err = ratio_val * np.sqrt((hyp_stat_err / hyp_val)**2 + (he3_stat_err / he3_val)**2)
    ratio_syst_err = ratio_val * np.sqrt((hyp_syst_err / hyp_val)**2 + (he3_syst_err / he3_val)**2)
    ratio_err = np.sqrt(ratio_stat_err**2 + ratio_syst_err**2)
    ratio_stat.SetBinContent(i, ratio_val)
    ratio_stat.SetBinError(i, ratio_err)


## enable stats
# drop fit function
ratio_stat.SetStats(1)

ratio_stat.Fit('pol0', 'RMI+', '', 1.5, 5)
fit_func = ratio_stat.GetFunction('pol0')
fit_func.SetLineColor(kOrangeC)

c = ROOT.TCanvas()
ratio_stat.Draw('PEX0')
ratio_syst.Draw('E2 same')

outfile = ROOT.TFile("ratio.root", "RECREATE")
ratio_stat.Write()
c.Write()
outfile.Close()



