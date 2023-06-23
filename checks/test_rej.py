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
hdlMC.eval_data_frame("fAbsGenPt = abs(fGenPt)")

utils.reweight_pt_spectrum(hdlMC, 'fAbsGenPt', he3_spectrum)

hdlMC_rew = hdlMC.apply_preselections('rej == True', inplace=False)
## plot normalised pt spectra before and after reweighting

hPtShapeBefore = ROOT.TH1F('hPtShapeBefore', 'hPtShapeBefore; #it{p}_{T}^{gen}; Counts', 100, 0, 5)
hPtShapeAfter = ROOT.TH1F('hPtShapeAfter', 'hPtShapeAfter; #it{p}_{T}^{gen}; Counts', 100, 0, 5)
hCosPABefore = ROOT.TH1F('hCosPABefore', 'hCosPABefore; cos(#theta_{PA}); Counts', 100, 0.985, 1)
hCosPAAfter = ROOT.TH1F('hCosPAAfter', 'hCosPAAfter; cos(#theta_{PA}); Counts', 100, 0.985, 1)

utils.fill_th1_hist(hPtShapeBefore, hdlMC, 'fAbsGenPt')
utils.fill_th1_hist(hPtShapeAfter, hdlMC_rew, 'fAbsGenPt')
utils.fill_th1_hist(hCosPABefore, hdlMC, 'fCosPA')
utils.fill_th1_hist(hCosPAAfter, hdlMC_rew, 'fCosPA')

outfile = ROOT.TFile("../../results/test_rej.root", "recreate")
outfile.cd()
hPtShapeBefore.Write()
hPtShapeAfter.Write()

hCosPABefore.Write()
hCosPAAfter.Write()
outfile.Close()

