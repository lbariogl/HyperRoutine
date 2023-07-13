import ROOT
import uproot
from hipe4ml.tree_handler import TreeHandler
import numpy as np

import argparse
import yaml

import sys
sys.path.append('../utils')
import utils as utils


hdlMC = TreeHandler(['/data/shared/hyp_run_3/mc/AO2D_MC_BDT.root'], 'O2mchypcands')
hdlMC.print_summary()

spectra_file = ROOT.TFile.Open('../utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
he3_spectrum.SetRange(1, 5)

hdlMC.eval_data_frame("fAbsGenPt = abs(fGenPt)")

utils.reweight_pt_spectrum(hdlMC, 'fAbsGenPt', he3_spectrum)

hdlMC_rew = hdlMC.apply_preselections('rej == True', inplace=False)
## plot normalised pt spectra before and after reweighting

hPtShapeBefore = ROOT.TH1F('hPtShapeBefore', 'hPtShapeBefore; #it{p}_{T}^{gen}; Counts', 500, 1, 5)
hPtShapeAfter = ROOT.TH1F('hPtShapeAfter', 'hPtShapeAfter; #it{p}_{T}^{gen} (GeV/#it{c}); Counts', 500, 1, 5)
hRecoSpectrumAfter = ROOT.TH1F('hRecoSpectrumAfter', 'hRecoSpectrumAfter; #it{p}_{T}^{rec} (GeV/#it{c}); Counts', 100, 1, 5)
hCosPABefore = ROOT.TH1F('hCosPABefore', 'hCosPABefore; cos(#theta_{PA}); Counts', 500, 0.985, 1)
hCosPAAfter = ROOT.TH1F('hCosPAAfter', 'hCosPAAfter; cos(#theta_{PA}); Counts', 500, 0.985, 1)

utils.fill_th1_hist(hPtShapeBefore, hdlMC, 'fAbsGenPt')
utils.fill_th1_hist(hPtShapeAfter, hdlMC_rew, 'fAbsGenPt')
utils.fill_th1_hist(hCosPABefore, hdlMC, 'fCosPA')
utils.fill_th1_hist(hCosPAAfter, hdlMC_rew.apply_preselections('fPt<2', inplace=False), 'fCosPA')
utils.fill_th1_hist(hRecoSpectrumAfter, hdlMC_rew, 'fPt')


cv6 = ROOT.TCanvas("cv6", "cv6", 800, 600)
he3_spectrum.SetRange(1, 5)
he3_spectrum.SetNpx(500)
he3_spectrum.SetNormalized(True)
hPtShapeAfter.DrawNormalized()
histo_from_tf1 = he3_spectrum.CreateHistogram()
histo_from_tf1.SetLineColor(2)
histo_from_tf1.SetLineWidth(2)
histo_from_tf1.SetLineColor(ROOT.kRed)
histo_from_tf1.DrawNormalized("same")
legend = ROOT.TLegend(0.6, 0.6, 0.78, 0.8)
legend.AddEntry(hPtShapeAfter, "Reweighted sample", "l")
legend.AddEntry(histo_from_tf1, "^{3} He Levy-Tsallis", "l")
legend.Draw()





outfile = ROOT.TFile("../../results/test_rej.root", "recreate")
outfile.cd()
hPtShapeBefore.Write()
hPtShapeAfter.Write()
hRecoSpectrumAfter.Write()

hCosPABefore.Write()
hCosPAAfter.Write()

cv6.Write()

outfile.Close()

