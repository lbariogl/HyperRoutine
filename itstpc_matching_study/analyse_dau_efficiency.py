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

tree = uproot.open("DauTreeMC.root")["DauTreeMC"].arrays(library="pd")

outfile_name = "dau_efficiency_" + str(daughter_pdg) + ".root"
##
tree.query(f"abs(pdg) == {daughter_pdg}", inplace=True)

h_gen_radius_hist   = ROOT.TH1F("h_gen_radius", ";Radius (cm)", 50, 0, 40)
h_gen_pt_hist = ROOT.TH1F("h_gen_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)
h_rec_radius_hist_list = []
h_rec_pt_hist_list = []

h_detector_list = ["ITS", "TPC", "ITS-TPC", "TPC-TOF", "TPC-TRD", "ITS-TPC-TOF", "ITS-TPC-TRD","TPC-TRD-TOF","ITS-TPC-TRD-TOF"]

for idet, det in enumerate(h_detector_list):
    h_rec_radius_hist_list.append(ROOT.TH1F("h_rec_radius_" + det, ";Radius (cm)", 50, 0, 40))
    h_rec_pt_hist_list.append(ROOT.TH1F("h_rec_pt_" + det, ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10))

fill_th1_hist(h_gen_radius_hist, tree, "genRad")
fill_th1_hist(h_gen_pt_hist, tree, "genPt")
h_gen_radius_hist.Sumw2()
h_gen_pt_hist.Sumw2()

for idet, det in enumerate(h_detector_list):
    ## bitwise operation to check if the detector is present
    mask = 1 << idet
    bool_mask = tree['detectorBMap'] & mask
    ## move into bool mask
    bool_mask = bool_mask.astype(bool)
    df_det = tree[bool_mask]
    hist_rad = h_rec_radius_hist_list[idet]
    hist_pt = h_rec_pt_hist_list[idet]
    fill_th1_hist(hist_rad, df_det, "genRad")
    fill_th1_hist(hist_pt, df_det, "genPt")

    hist_rad.Sumw2()
    hist_pt.Sumw2()

    hist_rad.Divide(h_gen_radius_hist)
    hist_pt.Divide(h_gen_pt_hist)



### filter ITS + TPC 
if debug:
    h_rej_flag_hist = ROOT.TH1F("h_rej_flag", ";Rejection flag", 20, -1.5, 18.5)
    h2_rej_flag_pt = ROOT.TH2F("h2_rej_flag_pt", ";#it{p}_{T} (GeV/#it{c}); Rejection flag", 50, 0, 10, 20, -1.5, 18.5)
    h_rec_its_tpc_rad = ROOT.TH1F("h_rec_its_tpc_rad", ";Radius (cm)", 50, 0, 40)
    h_rec_its_tpc_pt = ROOT.TH1F("h_rec_its_tpc_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)
    h_rec_its_tpc_eff_rad = ROOT.TH1F("its_tpc_algo_eff_rad", "; Radius (cm); Algorithm Efficiency", 50, 0, 40)
    h_rec_its_tpc_eff_pt = ROOT.TH1F("its_tpc_algo_eff_pt", "; #it{p}_{T} (GeV/#it{c}); Algorithm Efficiency", 50, 0, 10)
    h_rec_match_chi2 = ROOT.TH1F("h_rec_match_chi2", ";Matching #chi^{2}", 50, 0, 100)
    h_fake_itstpc_rad = ROOT.TH1F("h_fake_itstpc_rad", ";Radius (cm)", 50, 0, 40)
    h_fake_itstpc_pt = ROOT.TH1F("h_fake_itstpc_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)
    h_true_itstpc_rad = ROOT.TH1F("h_true_itstpc_rad", ";Radius (cm)", 50, 0, 40)
    h_true_itstpc_pt = ROOT.TH1F("h_true_itstpc_pt", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)


    maskITS = 1 << 0
    maskTPC = 1 << 1
    maskITSTPC = 1 << 2

    ### select only ITS + TPC
    bool_mask_its = (tree['detectorBMap'] & maskITS).astype(bool)
    tree = tree[bool_mask_its]
    bool_mask_tpc = (tree['detectorBMap'] & maskTPC).astype(bool)
    tree = tree[bool_mask_tpc]
    fill_th1_hist(h_rec_its_tpc_rad, tree, "genRad")
    fill_th1_hist(h_rec_its_tpc_pt, tree, "genPt")
    fill_th1_hist(h_rec_match_chi2, tree, "chi2Match")

    bool_mask_its_tpc = (tree['detectorBMap'] & maskITSTPC).astype(bool)
    fill_th1_hist(h_rec_its_tpc_eff_rad, tree[bool_mask_its_tpc], "genRad")
    fill_th1_hist(h_rec_its_tpc_eff_pt, tree[bool_mask_its_tpc], "genPt")
    fill_th1_hist(h_fake_itstpc_rad, tree[bool_mask_its_tpc].query('isITSTPCfake==True'), "genRad")
    fill_th1_hist(h_fake_itstpc_pt, tree[bool_mask_its_tpc].query('isITSTPCfake==True'), "genPt")
    fill_th1_hist(h_true_itstpc_rad, tree[bool_mask_its_tpc].query('isITSTPCfake==False'), "genRad")
    fill_th1_hist(h_true_itstpc_pt, tree[bool_mask_its_tpc].query('isITSTPCfake==False'), "genPt")
    h_rec_its_tpc_eff_rad.Sumw2()
    h_rec_its_tpc_eff_pt.Sumw2()
    h_rec_its_tpc_eff_rad.Divide(h_rec_its_tpc_rad)
    h_rec_its_tpc_eff_pt.Divide(h_rec_its_tpc_pt)

    bool_mask_not_its_tpc = np.logical_not(bool_mask_its_tpc)
    fill_th1_hist(h_rej_flag_hist, tree[bool_mask_not_its_tpc], "rejFlag")
    fill_th2_hist(h2_rej_flag_pt, tree[bool_mask_not_its_tpc], "genPt", "rejFlag")

# print(tree[bool_mask_its_tpc].query("rejFlag==0")[['tfNum', 'itsRef','tpcRef','rejFlag', 'chi2Match', 'nRefs', 'detectorBMap', 'isAB', 'isITSTPCfake', 'isITSfake', 'isTPCfake']])
# print(tree[bool_mask_not_its_tpc].query("tfNum==1")[['tfNum', 'itsRef','tpcRef','rejFlag', 'genPt', 'genEta']].head(20))



n_det = 3
## eff vs radius
c_rad = ROOT.TCanvas("c_rad", "c_rad", 800, 800)
leg = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for i in range(0, n_det):
    h_rec_radius_hist_list[i].SetLineColor(i+1)
    h_rec_radius_hist_list[i].SetMarkerColor(i+1)
    h_rec_radius_hist_list[i].SetMarkerStyle(20)
    h_rec_radius_hist_list[i].SetMarkerSize(1)
    h_rec_radius_hist_list[i].SetMinimum(0)
    h_rec_radius_hist_list[i].SetMaximum(1.2)
    h_rec_radius_hist_list[i].GetXaxis().SetTitle("Production radius (cm)")
    h_rec_radius_hist_list[i].GetYaxis().SetTitle("Efficiency")
    h_rec_radius_hist_list[i].SetTitle(f"Secondary {daughter_pdg} efficiency")
    h_rec_radius_hist_list[i].SetStats(0)
    h_rec_radius_hist_list[i].Draw("same")
    sa_string = ''
    if i<2:
        sa_string = '_SA'
    leg.AddEntry(h_rec_radius_hist_list[i], f'{h_detector_list[i]}{sa_string}', "lep")
leg.Draw()

## eff vs pt
c_pt = ROOT.TCanvas("c_pt", "c_pt", 800, 800)
leg2 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for i in range(0, n_det):
    h_rec_pt_hist_list[i].SetLineColor(i+1)
    h_rec_pt_hist_list[i].SetMarkerColor(i+1)
    h_rec_pt_hist_list[i].SetMarkerStyle(20)
    h_rec_pt_hist_list[i].SetMarkerSize(1)
    h_rec_pt_hist_list[i].SetMinimum(0)
    h_rec_pt_hist_list[i].SetMaximum(1.2)
    h_rec_pt_hist_list[i].GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
    h_rec_pt_hist_list[i].GetYaxis().SetTitle("Efficiency")
    h_rec_pt_hist_list[i].SetTitle(f"Secondary {daughter_pdg} efficiency")
    h_rec_pt_hist_list[i].SetStats(0)
    h_rec_pt_hist_list[i].Draw("same")
    sa_string = ''
    if i<2:
        sa_string = '_SA'
    leg2.AddEntry(h_rec_pt_hist_list[i], f'{h_detector_list[i]}{sa_string}', "lep")
leg2.Draw()


### itstpc efficiency for fully matched and AB tracks
h_radius_rec_ab = ROOT.TH1F("h_radius_rec_ab", ";Radius (cm)", 50, 0, 40)
h_radius_rec_full = ROOT.TH1F("h_radius_rec_full", ";Radius (cm)", 50, 0, 40)

mask = 1 << 2
bool_mask = tree['detectorBMap'] & mask
## move into bool mask
bool_mask = bool_mask.astype(bool)
df_det = tree[bool_mask]
fill_th1_hist(h_radius_rec_ab, df_det.query("isAB==True"), "genRad")
fill_th1_hist(h_radius_rec_full, df_det.query("isAB==False"), "genRad")

h_radius_rec_ab.Sumw2()
h_radius_rec_full.Sumw2()

h_radius_rec_ab.Divide(h_gen_radius_hist)
h_radius_rec_full.Divide(h_gen_radius_hist)

c_rad_ab = ROOT.TCanvas("c_rad_ab", "c_rad_ab", 800, 800)
h_radius_rec_ab.SetLineColor(1)
h_radius_rec_ab.SetMarkerColor(1)
h_radius_rec_ab.SetMarkerStyle(20)
h_radius_rec_ab.SetMarkerSize(1)
h_radius_rec_ab.SetMinimum(0)
h_radius_rec_ab.SetMaximum(1.2)
h_radius_rec_ab.GetXaxis().SetTitle("Production radius (cm)")
h_radius_rec_ab.GetYaxis().SetTitle("Efficiency")
h_radius_rec_ab.SetStats(0)
h_radius_rec_ab.Draw()
h_radius_rec_full.SetLineColor(2)
h_radius_rec_full.SetMarkerColor(2)
h_radius_rec_full.SetMarkerStyle(20)
h_radius_rec_full.SetMarkerSize(1)
h_radius_rec_full.Draw("same")

h_rec_radius_hist_list[2].SetLineColor(3)
h_rec_radius_hist_list[2].SetMarkerColor(3)
h_rec_radius_hist_list[2].SetMarkerStyle(20)
h_rec_radius_hist_list[2].SetMarkerSize(1)
h_rec_radius_hist_list[2].Draw("same")


leg3 = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
leg3.AddEntry(h_radius_rec_ab, "AB", "lep")
leg3.AddEntry(h_radius_rec_full, "Fully matched", "lep")
leg3.AddEntry(h_rec_radius_hist_list[2], "Total ITS-TPC", "lep")
leg3.Draw()



############ write to file
outfile = ROOT.TFile(outfile_name, "RECREATE")
for h in h_rec_radius_hist_list:
    h.Write()

c_rad.Write()
c_pt.Write()
c_rad_ab.Write()

if debug:
    h_rec_its_tpc_eff_pt.Write()
    h_rec_its_tpc_eff_rad.Write()
    h_rej_flag_hist.Write()
    h2_rej_flag_pt.Write()
    h_rec_match_chi2.Write()
    h_fake_itstpc_pt.Write()
    h_fake_itstpc_rad.Write()
    h_true_itstpc_pt.Write()
    h_true_itstpc_rad.Write()

outfile.Close()

