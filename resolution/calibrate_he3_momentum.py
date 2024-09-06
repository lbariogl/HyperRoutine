import numpy as np
import uproot
import ROOT
import sys
sys.path.append('../utils')
import utils as utils

ROOT.gROOT.SetBatch(True)

df = uproot.open("/home/fmazzasc/AO2D_h3l.root")["DF_2316988992350369/O2mchypcands;1"].arrays(library="pd")
df['fHePIDHypo'] = np.right_shift(df['fFlags'], 4)
df.query("fIsReco > 0 and fHePIDHypo==6", inplace=True)
utils.correct_and_convert_df(df)

h2MomResoVsPtHe3, h_pt_shift = utils.create_pt_shift_histo(df)
h_pt_shift.Fit("pol2", "Q")

h2MomResoUpdated = ROOT.TH2F('updated', ';{}^{3}He #it{p}_{T} (GeV/#it{c});{}^{3}He #it{p}_{T}^{gen} - #it{p}_{T}^{reco} (GeV/#it{c})', 30, 1.3, 5, 50, -0.4, 0.4)

df["fPtHe3"] += -0.1286 - 0.1269 * df["fPtHe3"] + 0.06 * df["fPtHe3"]**2
print(df["fPtHe3"])
df.eval('PtResHe3 = (fGenPtHe3 - fPtHe3)', inplace=True)
print(df["PtResHe3"])
utils.fill_th2_hist(h2MomResoUpdated, df, 'fPtHe3', 'PtResHe3')


outfilename = "he3_momentum_calibration.root"
outfile = ROOT.TFile(outfilename, "RECREATE")
h2MomResoVsPtHe3.Write()
h2MomResoUpdated.Write()
h_pt_shift.Write()
outfile.Close()

