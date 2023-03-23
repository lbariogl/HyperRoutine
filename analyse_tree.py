import ROOT
import uproot
import pandas as pd
import argparse
import yaml
from hipe4ml.tree_handler import TreeHandler

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--mc', dest='mc', action='store_true', help="if True MC information is stored.")
parser.set_defaults(mc=False)
parser.add_argument('--input-files', dest='input_files', help="path to the input files.")
parser.set_defaults(input_file='../data/AO2D_merged.root')
parser.add_argument('--output-dir', dest='output_dir', help="path to the directory in which the output is stored.")
parser.set_defaults(output_dir='../results/')
parser.add_argument('--output-file', dest='output_file', help="name of the output file.")
parser.set_defaults(output_file='HypertritonResults.root')
parser.add_argument('--selection', dest='selection', help="selections to be bassed as query.")
parser.set_defaults(selection='fCosPA > 0.998 & fNTPCclusHe > 110 & abs(fDcaHe) > 0.1')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.")
parser.set_defaults(config_file='')
args = parser.parse_args()

# initialise parameters from parser (can be overwritten by external yaml file)
mc = args.mc
input_file_name = args.input_files
output_dir_name = args.output_dir
output_file_name = args.output_file
selections = args.selection

if  args.config_file != "":
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    mc = config['mc']
    input_files_name = config['input_files']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']
    selections = config['selection']

# utils
def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for i in range(df.shape[0]):
        h.Fill(df[var1][i], df[var2][i])

# tree names
tree_name_data = 'O2datahypcands'
tree_name_mc = 'O2mchypcands'

## create histograms
hCosPA = ROOT.TH1F("hCosPA", ";cos(#theta_{PA})", 50, 0.95, 1)
hNTPCclus = ROOT.TH1F("hNTPCclus", ";n TPC clusters", 50, 60, 200)
hMass3LH = ROOT.TH1F("h_3lh_mass", "; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 2.96, 3.04)
hMass4LH = ROOT.TH1F("h_4lh_mass", "; m({}^{4}_{#Lambda}H) (GeV/#it{c})", 50, 3.96, 4.04)
hPtRec = ROOT.TH1F("hPtRec", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 10)
hRadius = ROOT.TH1F("hRadius", ";Radius (cm)", 100, 0, 40)
hNSigHe = ROOT.TH1F("hNSigmaHe", ";n #sigma He3", 50, -3, 3)

## for MC only
hPtGen = ROOT.TH1F("hPtGen", "; Pt gen", 50, 0, 10)
hResolutionPx = ROOT.TH1F("hResolutionPx", ";Resolution #it{p_{x}}", 50, -0.2, 0.2)
hResolutionPy = ROOT.TH1F("hResolutionPy", ";Resolution #it{p_{y}}", 50, -0.2, 0.2)
hResolutionPz = ROOT.TH1F("hResolutionPz", ";Resolution #it{p_{z}}", 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F("hResolutionDecVtxX", "; Resolution Dec X", 50, -0.2, 0.2)
hResolutionDecVtxY = ROOT.TH1F("hResolutionDecVtxY", "; Resolution Dec Y", 50, -0.2, 0.2)
hResolutionDecVtxZ = ROOT.TH1F("hResolutionDecVtxZ", "; Resolution Dec Z", 50, -0.2, 0.2)

# creating the dataframe
tree_hdl = TreeHandler(input_files_name, tree_name_data)
df = tree_hdl.get_data_frame()


df.eval('pt = sqrt(fPx**2 + fPy**2)', inplace=True)
df.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)
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
df_filtered = df.query(selections)


## fill histograms
fill_th1_hist(hPtRec, df_filtered, 'pt')
fill_th1_hist(hCosPA, df_filtered, 'fCosPA')
fill_th1_hist(hRadius, df_filtered, 'fDecRad')
fill_th1_hist(hNTPCclus, df_filtered, 'fNTPCclusHe')
fill_th1_hist(hNSigHe, df_filtered, 'fNSigmaHe')
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
f = ROOT.TFile(f"{output_dir_name}/{output_file_name}", "RECREATE")

hCosPA.Write()
hRadius.Write()
hNTPCclus.Write()
hNSigHe.Write()
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


