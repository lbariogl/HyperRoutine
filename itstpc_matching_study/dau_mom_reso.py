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


h2MomResoVsPtHe3 = ROOT.TH2F("h2MomResoVsPtHe3", ";#it{p}_{T} (GeV/#it{c});#sigma(#it{p}_{T})/#it{p}_{T}", 50, 0, 5, 50, -0.2, 0.2)
h2MomResoVsPtPi = ROOT.TH2F("h2MomResoVsPtPi", ";#it{p}_{T} (GeV/#it{c});#sigma(#it{p}_{T})/#it{p}_{T}", 50, 0, 2, 50, -0.2, 0.2)

tree = uproot.open("../../match_res/DauTreeMC.root")["DauTreeMC"].arrays(library="pd")
print(tree.columns)

tree.query("itsTPCPt>0", inplace=True)
 


tree_he3 = tree.query(f"abs(pdg) == 1000020030")
tree_he3.loc[:,'itsTPCPt'] *= 2
tree_pi = tree.query(f"abs(pdg) == 211")

tree_he3.eval("PtRes = (itsTPCPt- genPt) / genPt", inplace=True)
tree_pi.eval("PtRes = (itsTPCPt - genPt) / genPt", inplace=True)



fill_th2_hist(h2MomResoVsPtHe3, tree_he3, "genPt", "PtRes")
fill_th2_hist(h2MomResoVsPtPi, tree_pi, "genPt", "PtRes")

outfile = ROOT.TFile("../../match_res/dau_mom_reso.root", "recreate")
outfile.cd()
h2MomResoVsPtHe3.Write()
h2MomResoVsPtPi.Write()
outfile.Close()
