import ROOT
import uproot
import pandas as pd

def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for i in range(df.shape[0]):
        h.Fill(df[var1][i], df[var2][i])


## create histograms
h_px = ROOT.TH1F("h_px", ";Resolution Px", 50, -0.2, 0.2)
h_py = ROOT.TH1F("h_py", ";Resolution Py", 50, -0.2, 0.2)
h_pz = ROOT.TH1F("h_pz", ";Resolution Pz", 50, -0.2, 0.2)

h_x = ROOT.TH1F("h_x", "; Resolution Dec X", 50, -0.2, 0.2)
h_y = ROOT.TH1F("h_y", "; Resolution Dec Y", 50, -0.2, 0.2)
h_z = ROOT.TH1F("h_z", "; Resolution Dec Z", 50, -0.2, 0.2)

h_cosPA = ROOT.TH1F("h_cosPA", "; cosPA", 50, 0.95, 1)
h_nTPCClus = ROOT.TH1F("h_nTPCClus", "; nTPCClus", 50, 60, 200)
h_mass_3lh = ROOT.TH1F("h_3lh_mass", "; 3LH mass", 50, 2.96, 3.04)
h_mass_4lh = ROOT.TH1F("h_4lh_mass", "; 4LH mass", 50, 3.96, 4.04)

h_pt_gen = ROOT.TH1F("h_pt_gen", "; Pt gen", 50, 0, 10)
h_pt_rec = ROOT.TH1F("h_pt_rec", "; Pt rec", 50, 0, 10)

file = uproot.open("AnalysisResults_trees.root")
dfs = file.keys()

mc = True
tree_name = "O2datahypcands" if not mc else "O2mchypcands"

pd_arr = []

for df in dfs:
    if tree_name in df:
        pd_arr.append(file[df].arrays(library='pd'))

df = pd.concat(pd_arr, ignore_index=True)
df.eval('pt = sqrt(fPx**2 + fPy**2)', inplace=True)
if mc:
    df.eval('pt_gen = sqrt(fGenPx**2 + fGenPy**2)', inplace=True)
    fill_th1_hist(h_pt_gen, df, 'pt_gen')
    df.query('fIsReco==True', inplace=True)
    df.eval('resPx = (fPx - fGenPx)/fGenPx', inplace=True)
    df.eval('resPy = (fPy - fGenPy)/fGenPy', inplace=True)
    df.eval('resPz = (fPz - fGenPz)/fGenPz', inplace=True)
    df.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)
    ## fill histograms
    fill_th1_hist(h_px, df, 'resPx')
    fill_th1_hist(h_py, df, 'resPy')
    fill_th1_hist(h_pz, df, 'resPz')
    fill_th1_hist(h_x, df, 'ResDecX')
    fill_th1_hist(h_y, df, 'ResDecY')
    fill_th1_hist(h_z, df, 'ResDecZ')



fill_th1_hist(h_pt_rec, df, 'pt')
fill_th1_hist(h_cosPA, df, 'fCosPA')
fill_th1_hist(h_nTPCClus, df, 'fNTPCclusHe')
fill_th1_hist(h_mass_3lh, df, 'fMassH3L')
fill_th1_hist(h_mass_4lh, df, 'fMassH4L')


## save to file root
f = ROOT.TFile("histos.root", "RECREATE")

if mc:
    h_px.Write()
    h_py.Write()
    h_pz.Write()
    h_x.Write()
    h_y.Write()
    h_z.Write()
    h_pt_gen.Write()
    h_eff_pt = h_pt_rec.Clone("h_eff_pt")
    h_eff_pt.Divide(h_pt_gen)
    h_eff_pt.SetTitle("; #it{p}_T (GeV/#it{c}); Efficiency")
    h_eff_pt.Write()

h_cosPA.Write()
h_nTPCClus.Write()
h_pt_rec.Write()
h_mass_3lh.Write()
h_mass_4lh.Write()



