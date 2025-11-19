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


h2MomResoVsPtHe3 = ROOT.TH2F("h2MomResoVsPtHe3", ";#it{p}_{T} (GeV/#it{c});#sigma(#it{p}_{T})/#it{p}_{T}", 50, 1.5, 5, 50, -0.4, 0.4)
h2MomResoVsPtPi = ROOT.TH2F("h2MomResoVsPtPi", ";#it{p}_{T} (GeV/#it{c});#sigma(#it{p}_{T})/#it{p}_{T}", 50, 0, 2, 50, -0.4, 0.4)

tree = uproot.open("../../match_res/DauTreeMC.root")["DauTreeMC"].arrays(library="pd")
print(tree.columns)

tree.query("itsTPCPt>0", inplace=True)
 


tree_he3 = tree.query(f"abs(pdg) == 1000020030")
tree_he3.loc[:,'itsTPCPt'] *= 2
tree_pi = tree.query(f"abs(pdg) == 211")

tree_he3.eval("PtRes = (itsTPCPt- genPt) / genPt", inplace=True)
tree_pi.eval("PtRes = (itsTPCPt - genPt) / genPt", inplace=True)


fill_th2_hist(h2MomResoVsPtHe3, tree_he3, "itsTPCPt", "PtRes")
fill_th2_hist(h2MomResoVsPtPi, tree_pi, "itsTPCPt", "PtRes")


### fit slices and save mean and sigma
h2MomResoVsPtHe3.FitSlicesY()
h2MomResoVsPtPi.FitSlicesY()

h2MomResoVsPtHe3_mean = ROOT.gDirectory.Get("h2MomResoVsPtHe3_1")
h2MomResoVsPtHe3_sigma = ROOT.gDirectory.Get("h2MomResoVsPtHe3_2")
h2MomResoVsPtPi_mean = ROOT.gDirectory.Get("h2MomResoVsPtPi_1")
h2MomResoVsPtPi_sigma = ROOT.gDirectory.Get("h2MomResoVsPtPi_2")

outfile = ROOT.TFile("../../match_res/dau_mom_reso.root", "recreate")
outfile.cd()

h2MomResoVsPtHe3.Write()
h2MomResoVsPtPi.Write()
h2MomResoVsPtHe3_mean.Write()
h2MomResoVsPtHe3_sigma.Write()
h2MomResoVsPtPi_mean.Write()
h2MomResoVsPtPi_sigma.Write()

outfile.Close()
