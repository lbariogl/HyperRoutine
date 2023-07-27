import ROOT
import uproot
import pandas as pd
import argparse
import yaml
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils


parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--mc', dest='mc', action='store_true',
                    help="if True MC information is stored.", default=False)
parser.add_argument('--input-files', dest='input_files',
                    help="path to the input files.", default='../data/AO2D_merged.root')
parser.add_argument('--output-dir', dest='output_dir',
                    help="path to the directory in which the output is stored.", default='../results/')
parser.add_argument('--output-file', dest='output_file',
                    help="name of the output file.", default='HypertritonResults.root')
parser.add_argument('--selection', dest='selection', help="selections to be bassed as query.",
                    default='fCosPA > 0.998 & fNTPCclusHe > 110 & abs(fDcaHe) > 0.1')
parser.add_argument('--is-matter', dest='is_matter',
                    help="path to the YAML file with configuration.", default='matter')

parser.add_argument('--skip-out-tree', dest='skip_out_tree', action='store_true', help="if True do not save output tree.")

parser.add_argument('--correction-file', dest='correction_file',
                    help="path to the file use for 3He pt correction.", default=None)

parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.", default='')
args = parser.parse_args()

# initialise parameters from parser (can be overwritten by external yaml file)
mc = args.mc
skip_out_tree = args.skip_out_tree
input_file_name = args.input_files
output_dir_name = args.output_dir
output_file_name = args.output_file
selections = args.selection
is_matter = args.is_matter
correction_file = args.correction_file

if args.config_file != "":
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    mc = config['mc']
    input_files_name = config['input_files']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']
    selections = config['selection']
    is_matter = config['is_matter']
    correction_file = config['correction_file']

matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')


print('**********************************')
print('    Running analyse_tree.py')
print('**********************************')

# create histograms
hCosPA = ROOT.TH1F("hCosPA", ";cos(#theta_{PA})", 50, 0.95, 1)
hNTPCclus = ROOT.TH1F("hNTPCclus", ";n TPC clusters", 50, 60, 200)
hMass3LH = ROOT.TH1F(
    "h_3lh_mass", "; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 40, 2.96, 3.04)
hMass4LH = ROOT.TH1F(
    "h_4lh_mass", ";  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})", 30, 3.89, 3.97)
hPtRec = ROOT.TH1F("hPtRec", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 5)
hRadius = ROOT.TH1F("hRadius", ";Radius (cm)", 100, 0, 40)
hDecLen = ROOT.TH1F("hDecLen", ";Decay length (cm)", 100, 0, 40)
hNSigHe = ROOT.TH1F("hNSigmaHe", ";n_{#sigma}^{TPC}({}^{3}He)", 50, -3, 3)

h2MassCosPA = ROOT.TH2F(
    "h2MassCosPA", ";cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0.99, 1, 50, 2.96, 3.04)
h2MassDecLen = ROOT.TH2F(
    "h2MassDecLen", ";Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 40, 50, 2.96, 3.04)
h2MassDCADaughters = ROOT.TH2F(
    "h2MassDCADaughters", ";DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 200, 0, 0.3, 50, 2.96, 3.04)
h2MassDCAHePv = ROOT.TH2F(
    "h2MassDCAHe", ";DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 2, 50, 2.96, 3.04)
h2MassPt = ROOT.TH2F(
    "h2MassPt", ";#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 0, 7, 50, 2.96, 3.04)
h2Mass4LHnSigmaHe = ROOT.TH2F(
    "h2Mass4LHnSigmaHe", ";n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})", 50, -4, 4, 30, 3.89, 3.97)
 

# for MC only
hPtGen = ROOT.TH1F("hPtGen", "; Pt gen", 50, 0, 5)
hResolutionPt = ROOT.TH1F(
    "hResolutionPt", ";(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, -0.2, 0.2)
hResolutionPtvsPt = ROOT.TH2F(
    "hResolutionPtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionHe3PtvsPt = ROOT.TH2F(
    "hResolutionHe3PtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionPiPtvsPt = ROOT.TH2F(
    "hResolutionPiPtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 1, 50, -0.2, 0.2)
hResolutionP = ROOT.TH1F(
    "hResolutionP", ";(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}", 50, -0.2, 0.2)
hResolutionPvsP = ROOT.TH2F(
    "hResolutionPvsP", ";#it{p}^{gen} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F(
    "hResolutionDecVtxX", "; Resolution Dec X", 50, -0.2, 0.2)
hResolutionDecVtxY = ROOT.TH1F(
    "hResolutionDecVtxY", "; Resolution Dec Y", 50, -0.2, 0.2)
hResolutionDecVtxZ = ROOT.TH1F(
    "hResolutionDecVtxZ", "; Resolution Dec Z", 50, -0.2, 0.2)

# creating the dataframe
tree_name = 'O2datahypcands' if not mc else 'O2mchypcands'
tree_hdl = TreeHandler(input_files_name, tree_name)
df = tree_hdl.get_data_frame()

# import correction file
correction_hist = None
if correction_file:
    corr_file = ROOT.TFile(correction_file)
    correction_hist = corr_file.Get('hShiftVsPtHe3')
    correction_hist.SetDirectory(0)

# try to convert
utils.correct_and_convert_df(df, correction_hist)

# add new columns
df.eval('fP = fPt * cosh(fEta)', inplace=True)
df.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)

# for MC only
if mc:
    df.eval('fGenP = fGenPt * cosh(fGenEta)', inplace=True)

    ##apply pT rejection
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    df.eval("fAbsGenPt = abs(fGenPt)", inplace=True)
    utils.reweight_pt_spectrum(df, 'fAbsGenPt', he3_spectrum)
    df.query('rej==True', inplace=True)
    n_gen = len(df)



    if is_matter == 'matter':
        utils.fill_th1_hist_abs(hPtGen, df.query(
            'fGenPt>0', inplace=False), 'fGenPt')
    elif is_matter == 'antimatter':
        utils.fill_th1_hist_abs(hPtGen, df.query(
            'fGenPt<0', inplace=False), 'fGenPt')
    else:
        utils.fill_th1_hist_abs(hPtGen, df, 'fGenPt')
    # select only reconstructed candidates
    df.query('fIsReco==True', inplace=True)
    n_reco = len(df)
    print("MC, integrated efficiency: ", n_reco/n_gen)

if selections == '':
    if is_matter == 'matter':
        selections = 'fIsMatter == True'
    elif is_matter == 'antimatter':
        selections = 'fIsMatter == False'
else:
    if is_matter == 'matter':
        selections = selections + ' & fIsMatter == True'
    elif is_matter == 'antimatter':
        selections = selections + ' & fIsMatter == False'

# filtering
if selections != '':
    df_filtered = df.query(selections)
else:
    df_filtered = df

# fill histograms
utils.fill_th1_hist(hPtRec, df_filtered, 'fPt')
utils.fill_th1_hist(hCosPA, df_filtered, 'fCosPA')
utils.fill_th1_hist(hRadius, df_filtered, 'fDecRad')
utils.fill_th1_hist(hDecLen, df_filtered, 'fDecLen')
utils.fill_th1_hist(hNTPCclus, df_filtered, 'fNTPCclusHe')
utils.fill_th1_hist(hNSigHe, df_filtered, 'fNSigmaHe')
utils.fill_th1_hist(hMass3LH, df_filtered, 'fMassH3L')
utils.fill_th1_hist(hMass4LH, df_filtered, 'fMassH4L')

utils.fill_th2_hist(h2MassCosPA, df_filtered, 'fCosPA', 'fMassH3L')
utils.fill_th2_hist(h2MassDecLen, df_filtered, 'fDecLen', 'fMassH3L')
utils.fill_th2_hist(h2MassDCADaughters, df_filtered, 'fDcaV0Daug', 'fMassH3L')
utils.fill_th2_hist(h2MassDCAHePv, df_filtered, 'fDcaHe', 'fMassH3L')
utils.fill_th2_hist(h2MassPt, df_filtered, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHnSigmaHe, df_filtered, 'fNSigmaHe', 'fMassH4L')

# for MC only
if mc:
    df_filtered.eval('resPt = (fPt - abs(fGenPt))/abs(fGenPt)', inplace=True)
    df_filtered.eval('resP = (fP - abs(fGenP))/abs(fGenP)', inplace=True)
    df_filtered.eval(
        'ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df_filtered.eval(
        'ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df_filtered.eval(
        'ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)
    utils.fill_th1_hist(hResolutionPt, df_filtered, 'resPt')
    utils.fill_th1_hist(hResolutionP, df_filtered, 'resP')
    utils.fill_th1_hist(hResolutionDecVtxX, df_filtered, 'ResDecX')
    utils.fill_th1_hist(hResolutionDecVtxY, df_filtered, 'ResDecY')
    utils.fill_th1_hist(hResolutionDecVtxZ, df_filtered, 'ResDecZ')
    utils.fill_th2_hist_abs(hResolutionPtvsPt, df_filtered, 'fGenPt', 'resPt')
    utils.fill_th2_hist_abs(hResolutionPvsP, df_filtered, 'fGenP', 'resP')

# save to file root
f = ROOT.TFile(f"{output_dir_name}/{output_file_name}", "RECREATE")

hCosPA.Write()
hRadius.Write()
hDecLen.Write()
hNTPCclus.Write()
hNSigHe.Write()
hPtRec.Write()
hMass3LH.Write()
hMass4LH.Write()

h2MassCosPA.Write()
h2MassDecLen.Write()
h2MassDCADaughters.Write()
h2MassDCAHePv.Write()
h2Mass4LHnSigmaHe.Write()
h2MassPt.Write()

if not skip_out_tree:
    df_filtered.to_parquet(f"{output_dir_name}/{output_file_name}.parquet")

if mc:
    f.mkdir("MC")
    f.cd("MC")
    hResolutionPt.Write()
    hResolutionPtvsPt.Write()
    hResolutionP.Write()
    hResolutionPvsP.Write()
    hResolutionDecVtxX.Write()
    hResolutionDecVtxY.Write()
    hResolutionDecVtxZ.Write()
    hPtGen.Write()
