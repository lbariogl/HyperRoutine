import ROOT
import uproot
import pandas as pd
import numpy as np
import argparse
import os
import sys
sys.path.append("../")
import utils

ROOT.gROOT.SetBatch(True)

def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(var1_val, var2_val)


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--dau_pdg', dest='daughter_pdg', help="daughter pdg code.")
parser.add_argument('--suffix', dest='suffix', help="suffix of the input/output files.")


args = parser.parse_args()
daughter_pdg = args.daughter_pdg
suffix = args.suffix
if suffix!=None:
    suffix = "_" + suffix

path_dir = "../../match_res"
tree = uproot.open(path_dir + "/" + f"DauTreeMC{suffix}.root")["DauTreeMC"].arrays(library="pd")
outfile_name = path_dir + "/" + "ab_efficiency_" + str(daughter_pdg) + suffix  + ".root"
##
tree.query(f"abs(pdg) == {daughter_pdg}", inplace=True)

h_gen_radius_hist   = ROOT.TH1F("h_gen_radius", ";Radius (cm)", 50, 0, 40)
h_gen_pt_hist = ROOT.TH1F("h_gen_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)

h_avail_tpc_radius_hist = ROOT.TH1F("h_avail_tpc_radius", "has TPCtrack & !hasITStrack;Radius (cm)", 50, 0, 40)
h_avail_ab_radius_hist = ROOT.TH1F("h_avail_ab_radius", "has L5, L6, TPCtrack & !hasITStrack;Radius (cm)", 50, 0, 40)
h_avail_tpc_pt_hist = ROOT.TH1F("h_avail_tpc_pt", "has TPCtrack & !hasITStrack;#it{p}_T (GeV/#it{c})", 50, 0, 10)
h_avail_ab_pt_hist = ROOT.TH1F("h_avail_ab_pt", "has L5, L6, TPCtrack & !hasITStrack;#it{p}_T (GeV/#it{c})", 50, 0, 10)

h_reco_ab_radius_hist = ROOT.TH1F("h_reco_ab_radius", ";Radius (cm)", 50, 0, 40)
h_reco_ab_pt_hist = ROOT.TH1F("h_reco_ab_pt", ";#it{p}_T (GeV/#it{c})", 50, 0, 10)
h_chi2_match_hist = ROOT.TH1F("h_chi2_match", ";#chi^{2}", 100, 0, 100)


fill_th1_hist(h_gen_radius_hist, tree, "genRad")
fill_th1_hist(h_gen_pt_hist, tree, "genPt")

## select only candidates that could be found by the TPC 
tree.query("tpcRef!=-1 and itsRef==-1", inplace=True)
fill_th1_hist(h_avail_tpc_radius_hist, tree, "genRad")
fill_th1_hist(h_avail_tpc_pt_hist, tree, "genPt")

## select only candidates that could be found by the AB
tree.query("clRefL5!=-1 and clRefL6!=-1", inplace=True)
tree.query("not(isAB==False and chi2Match!=-1)", inplace=True)
print("**** Findable AB tree **** \n", tree.query('tfNum==3 and isAB==False')[['tfNum', 'itsRef','tpcRef',  'clRefL5', 'clRefL6', 'clL5tracked', 'clL6tracked', 'isAB', 'isITSTPCfake', 'nClus','chi2Match']])
fill_th1_hist(h_avail_ab_radius_hist, tree, "genRad")
fill_th1_hist(h_avail_ab_pt_hist, tree, "genPt")


## select candidates that reco by the AB
tree.query("isAB==True", inplace=True)
fill_th1_hist(h_reco_ab_radius_hist, tree, "genRad")
fill_th1_hist(h_reco_ab_pt_hist, tree, "genPt")
fill_th1_hist(h_chi2_match_hist, tree, "chi2Match")

h_tpc_eff_rad = utils.computeEfficiency(h_avail_tpc_radius_hist, h_reco_ab_radius_hist, 'hEffTPCRadius')
h_eff_rad = utils.computeEfficiency(h_avail_ab_radius_hist, h_reco_ab_radius_hist, 'hEffRadius')

h_tpc_eff_pt = utils.computeEfficiency(h_avail_tpc_pt_hist, h_reco_ab_pt_hist, 'hEffTPCPt')
h_eff_pt = utils.computeEfficiency(h_avail_ab_pt_hist, h_reco_ab_pt_hist, 'hEffPt')





## dump the histograms to a file
outfile = ROOT.TFile(outfile_name, "RECREATE")
outfile.cd()
h_gen_radius_hist.Write()
h_gen_pt_hist.Write()
h_avail_ab_radius_hist.Write()
h_reco_ab_radius_hist.Write()
h_chi2_match_hist.Write()

## efficiency
h_tpc_eff_rad.Write()
h_tpc_eff_pt.Write()
h_eff_rad.Write()
h_eff_pt.Write()

outfile.Close()