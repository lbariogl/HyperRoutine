import ROOT
from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import yaml
import argparse
import numpy as np
from itertools import product
import os

import sys

sys.path.append("../utils")
import utils as utils

p_bins = np.linspace(0.5, 5, 25)
n_p_bins = len(p_bins) - 1
p_bins_array = np.array(p_bins, dtype=np.float64)


func_string = '([1] - TMath::Power((TMath::Abs(2*x/2.80839160743) / TMath::Sqrt(1 + 2*2*x*x/2.80839160743/2.80839160743)),[3]) - TMath::Log([2] + TMath::Power(1/TMath::Abs(2*x/2.80839160743),[4]))) * [0] / TMath::Power((TMath::Abs(2*x/2.80839160743) / TMath::Sqrt(1 + 2*2*x*x/2.80839160743/2.80839160743)),[3])'

tree_list = ["/data3/fmazzasc/hyp_run_3/pbpb/pass4/AO2D.root"]    
hdl = TreeHandler(tree_list, 'O2hypcands', folder_name='DF*')
utils.correct_and_convert_df(hdl, True, False, True)

hdl.apply_preselections('fAvgClusterSizeHe > 6 and fNSigmaHe > -3 and fNTPCclusHe > 90 and fIsMatter==False')
df = hdl.get_data_frame()

hTPCdEdXvsP = ROOT.TH2D(f"hTPCdEdXvsP", r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)", 25, 1.0, 6, 175, 0, 1400,)
hTPCdEdXvsP_toFit = ROOT.TH1D(f"hTPCdEdXvsP_toFit", r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)", n_p_bins, p_bins_array)


outfile = ROOT.TFile("output.root", "RECREATE")
dataset_dir = outfile.mkdir("TPCdEdXvsP")

for i_p in range(0, n_p_bins):

    p_bin = [p_bins[i_p], p_bins[i_p + 1]]
    p_sel = f"abs(fTPCmomHe) > {p_bin[0]} and abs(fTPCmomHe) < {p_bin[1]}"
    print(f"psel: {p_sel}")
    bin_df = df.query(p_sel, inplace=False)

    p_label = (f"{p_bin[0]:.2f} "+ r"#leq #it{p}/z < " + f"{p_bin[1]:.2f}"+ r" GeV/#it{c}")

    hTPCdEdX_pbin = ROOT.TH1D(f"hTPCdEdX_p{i_p}", p_label + r";d#it{E}/d#it{X} (a. u.); counts", 175, 0, 1400,)

    for dEdX, p in zip(bin_df["fTPCsignalHe"], bin_df[f"fTPCmomHe"]):
        hTPCdEdXvsP.Fill(p, dEdX)
        hTPCdEdX_pbin.Fill(dEdX)

    mean = hTPCdEdX_pbin.GetMean()
    rms = hTPCdEdX_pbin.GetRMS()
    hTPCdEdX_pbin.Fit("gaus", "MQRL+", "", mean - 3 * rms, mean + 3 * rms)
    fitFunc = hTPCdEdX_pbin.GetFunction("gaus")
    hTPCdEdXvsP_toFit.SetBinContent(i_p + 1, fitFunc.GetParameter(1))
    hTPCdEdXvsP_toFit.SetBinError(i_p + 1, fitFunc.GetParameter(2))
    hTPCdEdX_pbin.Write()

# Drawing default BB function
func_BB_default = ROOT.TF1("func_BB_default", func_string, 0.5, 6, 5)
func_BB_default.SetParameters(-321.34, 0.6539, 1.591, 0.8225, 2.363)
func_BB_default.SetLineColor(ROOT.kRed)

# Defining BB function for fit
func_BB_fit = ROOT.TF1("func_BB_fit", func_string, 0.5, 6, 5)
func_BB_fit.SetLineColor(ROOT.kBlue)
func_BB_fit.SetParameters(-321.34, 0.6539, 1.591, 0.8225, 2.363)
hTPCdEdXvsP_toFit.Fit(func_BB_fit)

outfile.cd()

cTPCdEdXvsP = ROOT.TCanvas(f"cTPCdEdXvsP", f"cTPCdEdXvsP", 800, 600)
cTPCdEdXvsP.DrawFrame(0.8, 0.0, 4.0, 1400.0, r";#it{p}/z (GeV/#it{c}); d#it{E}/d#it{X} (a. u.)")

hTPCdEdXvsP.Draw("colz same")
func_BB_default.Draw("L same")
func_BB_fit.Draw("L same")

hTPCdEdXvsP.Write()
cTPCdEdXvsP.Write()
hTPCdEdXvsP_toFit.Write()