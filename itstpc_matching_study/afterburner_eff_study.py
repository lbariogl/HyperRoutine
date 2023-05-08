import ROOT
import uproot
import pandas as pd
import numpy as np
import argparse

ROOT.gROOT.SetBatch(True)

def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(var1_val, var2_val)


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--dau_pdg', dest='daughter_pdg', help="daughter pdg code.")
parser.add_argument('--debug', dest='debug', help="fill debug histograms.", action='store_true')

args = parser.parse_args()
daughter_pdg = args.daughter_pdg
debug = args.debug
path_dir = "../../match_res"
tree = uproot.open(path_dir + "/" + "DauTreeMC.root")["DauTreeMC"].arrays(library="pd")
outfile_name = path_dir + "/" + "ab_efficiency_" + str(daughter_pdg) + ".root"
##
tree.query(f"abs(pdg) == {daughter_pdg}", inplace=True)

h_gen_radius_hist   = ROOT.TH1F("h_gen_radius", ";Radius (cm)", 50, 0, 40)
h_gen_pt_hist = ROOT.TH1F("h_gen_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)
h_avail_ab_radius_hist = ROOT.TH1F("h_avail_ab_radius", "has L5, L6, TPCtrack & !hasITStrack;Radius (cm)", 50, 0, 40)
h_reco_ab_radius_hist = ROOT.TH1F("h_reco_ab_radius", ";Radius (cm)", 50, 0, 40)


fill_th1_hist(h_gen_radius_hist, tree, "genRad")
fill_th1_hist(h_gen_pt_hist, tree, "genPt")

## select only candidates that could be found by the AB
tree.query("clRefL5!=-1 and clRefL6!=-1 and tpcRef!=-1 and itsRef==-1 ", inplace=True)

print("**** Findable AB tree **** \n", tree[['tfNum', 'itsRef','tpcRef',  'clRefL5', 'clRefL6', 'clL5tracked', 'clL6tracked', 'isAB']])

fill_th1_hist(h_avail_ab_radius_hist, tree, "genRad")

## select candidates that were found by the AB
tree.query("isAB==True", inplace=True)

fill_th1_hist(h_reco_ab_radius_hist, tree, "genRad")



## dump the histograms to a file
outfile = ROOT.TFile(outfile_name, "RECREATE")
outfile.cd()
h_gen_radius_hist.Write()
h_gen_pt_hist.Write()
h_avail_ab_radius_hist.Write()
h_reco_ab_radius_hist.Write()
outfile.Close()



