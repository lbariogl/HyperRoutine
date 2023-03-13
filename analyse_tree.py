import ROOT
import uproot
import pandas as pd

def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for i in range(df.shape[0]):
        h.Fill(df[var1][i], df[var2][i])

input_file_name = '../data/AO2D_merged.root'
output_dir_name = '../results/'
tree_name_data = "O2datahypcands"
tree_name_mc = "O2mchypcands"
mc = False

## create histograms
hCosPA = ROOT.TH1F("hCosPA", ";cos(#theta_{PA})", 50, 0.95, 1)
hNTPCclus = ROOT.TH1F("hNTPCclus", ";n TPC clusters", 50, 60, 200)
hMass3LH = ROOT.TH1F("h_3lh_mass", "; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 2.96, 3.04)
hMass4LH = ROOT.TH1F("h_4lh_mass", "; 4LH mass", 50, 3.96, 4.04)
hPtRec = ROOT.TH1F("hPtRec", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)

## for MC only
hPtGen = ROOT.TH1F("hPtGen", "; Pt gen", 50, 0, 10)
hResolutionPx = ROOT.TH1F("hResolutionPx", ";Resolution #it{p_{x}}", 50, -0.2, 0.2)
hResolutionPy = ROOT.TH1F("hResolutionPy", ";Resolution #it{p_{y}}", 50, -0.2, 0.2)
hResolutionPz = ROOT.TH1F("hResolutionPz", ";Resolution #it{p_{z}}", 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F("hResolutionDecVtxX", "; Resolution Dec X", 50, -0.2, 0.2)
hResolutionDecVtxY = ROOT.TH1F("hResolutionDecVtxY", "; Resolution Dec Y", 50, -0.2, 0.2)
hResolutionDecVtxZ = ROOT.TH1F("hResolutionDecVtxZ", "; Resolution Dec Z", 50, -0.2, 0.2)

# creating the dataframe
pd_arr = []
data_events = uproot.open(f'{input_file_name}:{tree_name_data}')
pd_arr.append(data_events.arrays(library='pd'))
mc_events = None
if mc:
    pd_arr.append(mc_events.arrays(library='pd'))

df = pd.concat(pd_arr, ignore_index=True)
df.eval('pt = sqrt(fPx**2 + fPy**2)', inplace=True)
if mc:
    df.eval('pt_gen = sqrt(fGenPx**2 + fGenPy**2)', inplace=True)
    df.query('fIsReco==True', inplace=True)
    df.eval('resPx = (fPx - fGenPx)/fGenPx', inplace=True)
    df.eval('resPy = (fPy - fGenPy)/fGenPy', inplace=True)
    df.eval('resPz = (fPz - fGenPz)/fGenPz', inplace=True)
    df.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)

# filtering
df_filtered = df.query("fCosPA > 0.998 & fNTPCclusHe > 110 & abs(fDcaHe) < 0.1")

## fill histograms
fill_th1_hist(hPtRec, df_filtered, 'pt')
fill_th1_hist(hCosPA, df_filtered, 'fCosPA')
fill_th1_hist(hNTPCclus, df_filtered, 'fNTPCclusHe')
fill_th1_hist(hMass3LH, df_filtered, 'fMassH3L')
fill_th1_hist(hMass4LH, df_filtered, 'fMassH4L')

if mc:
    fill_th1_hist(hPtGen, df_filtered, 'pt_gen')
    fill_th1_hist(hResolutionPx, df_filtered, 'resPx')
    fill_th1_hist(hResolutionPy, df_filtered, 'resPy')
    fill_th1_hist(hResolutionPz, df_filtered, 'resPz')
    fill_th1_hist(hResolutionDecVtxX, df_filtered, 'ResDecX')
    fill_th1_hist(hResolutionDecVtxY, df_filtered, 'ResDecY')
    fill_th1_hist(hResolutionDecVtxZ, df_filtered, 'ResDecZ')

## save to file root
f = ROOT.TFile(f"{output_dir_name}/HypertritonResults.root", "RECREATE")

hCosPA.Write()
hNTPCclus.Write()
hPtRec.Write()
hMass3LH.Write()
hMass4LH.Write()

if mc:
    hResolutionPx.Write()
    hResolutionPy.Write()
    hResolutionPz.Write()
    hResolutionDecVtxX.Write()
    hResolutionDecVtxY.Write()
    hResolutionDecVtxZ.Write()
    hPtGen.Write()
    # compute efficiencies
    hEffPt = hPtRec.Clone("hEffPt")
    hEffPt.Divide(hPtGen)
    hEffPt.SetTitle("; #it{p}_T (GeV/#it{c}); Efficiency")
    hEffPt.Write()


