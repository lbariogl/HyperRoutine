import ROOT
import uproot
import pandas as pd
import argparse
import yaml
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--mc', dest='mc', action='store_true', help="if True MC information is stored.", default=False)
parser.add_argument('--input-files', dest='input_files', help="path to the input files.", default='../data/AO2D_merged.root')
parser.add_argument('--output-dir', dest='output_dir', help="path to the directory in which the output is stored.", default='../results/')
parser.add_argument('--output-file', dest='output_file', help="name of the output file without extension", default='HypertritonResults')
parser.add_argument('--selection', dest='selection', help="selections to be bassed as query.", default='fCosPA > 0.998 & fNTPCclusHe > 110 & abs(fDcaHe) > 0.1')
parser.add_argument('--is-matter', dest='is_matter', help="path to the YAML file with configuration.", default='matter')
parser.add_argument('--skip-out-tree', dest='skip_out_tree', action='store_true', help="if True do not save output tree.")
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()

# initialise parameters from parser (can be overwritten by external yaml file)
mc = args.mc
skip_out_tree = args.skip_out_tree
input_file_name = args.input_files
output_dir_name = args.output_dir
output_file_name = args.output_file
selections_string = args.selection
is_matter = args.is_matter

if args.config_file != "":
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    mc = config['mc']
    input_files_name = config['input_files']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']
    selections = config['selection']
    selections_string = utils.convert_sel_to_string(selections)
    is_matter = config['is_matter']


matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')


print('**********************************')
print('    Running analyse_tree.py')
print('**********************************')

############# Create histograms #############
hCosPA = ROOT.TH1F("hCosPA", ";cos(#theta_{PA})", 50, 0.95, 1)
hNTPCclus = ROOT.TH1F("hNTPCclus", ";n TPC clusters", 50, 60, 200)
hMass3LH = ROOT.TH1F("h_3lh_mass", "; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 40, 2.96, 3.04)
hMass4LH = ROOT.TH1F("h_4lh_mass", ";  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})", 30, 3.89, 3.97)
hPtRec = ROOT.TH1F("hPtRec", ";#it{p}_{T} (GeV/#it{c})", 50, 0, 5)
hCtRec = ROOT.TH1F("hCtRec", ";#it{c#tau} (cm)", 50, 0, 40)
hRadius = ROOT.TH1F("hRadius", ";Radius (cm)", 100, 0, 40)
hDecLen = ROOT.TH1F("hDecLen", ";Decay length (cm)", 100, 0, 40)
hNSigHe = ROOT.TH1F("hNSigmaHe", ";n_{#sigma}^{TPC}({}^{3}He)", 50, -3, 3)
h2MassCosPA = ROOT.TH2F("h2MassCosPA", ";cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0.99, 1, 50, 2.96, 3.04)
h2MassDecLen = ROOT.TH2F("h2MassDecLen", ";Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 40, 50, 2.96, 3.04)
h2MassDCADaughters = ROOT.TH2F("h2MassDCADaughters", ";DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 200, 0, 0.3, 50, 2.96, 3.04)
h2MassDCAHePv = ROOT.TH2F("h2MassDCAHe", ";DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 2, 50, 2.96, 3.04)
h2MassPt = ROOT.TH2F("h2MassPt", ";#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 0, 7, 50, 2.96, 3.04)
h2Mass4LHnSigmaHe = ROOT.TH2F("h2Mass4LHnSigmaHe", ";n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})", 50, -4, 4, 30, 3.89, 3.97)
# for MC only
hPtGen = ROOT.TH1F("hPtGen", "; Pt gen", 50, 0, 5)
hCtGen = ROOT.TH1F("hCtGen", "; Ct gen", 50, 0, 40)
hResolutionPt = ROOT.TH1F("hResolutionPt", ";(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, -0.2, 0.2)
hResolutionPtvsPt = ROOT.TH2F("hResolutionPtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionHe3PtvsPt = ROOT.TH2F("hResolutionHe3PtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionPiPtvsPt = ROOT.TH2F("hResolutionPiPtvsPt", ";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 1, 50, -0.2, 0.2)
hResolutionP = ROOT.TH1F("hResolutionP", ";(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}", 50, -0.2, 0.2)
hResolutionPvsP = ROOT.TH2F("hResolutionPvsP", ";#it{p}^{gen} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F("hResolutionDecVtxX", "; Resolution Dec X", 50, -0.2, 0.2)
hResolutionDecVtxY = ROOT.TH1F("hResolutionDecVtxY", "; Resolution Dec Y", 50, -0.2, 0.2)
hResolutionDecVtxZ = ROOT.TH1F("hResolutionDecVtxZ", "; Resolution Dec Z", 50, -0.2, 0.2)

############# Read trees #############
tree_name = 'O2datahypcands' if not mc else 'O2mchypcands'
tree_hdl = TreeHandler(input_files_name, tree_name)
df = tree_hdl.get_data_frame()
# correct and convert dataframe
utils.correct_and_convert_df(df, True, mc)


############# Apply pre-selections to MC #############
if mc:
    mc_pre_sels = ""
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    utils.reweight_pt_spectrum(df, 'fAbsGenPt', he3_spectrum)
    mc_pre_sels += 'rej==True'
    if is_matter == 'matter':
        mc_pre_sels += 'and fGenPt>0'
    elif is_matter == 'antimatter':
        mc_pre_sels += 'and fGenPt<0'
    

    df.query(mc_pre_sels, inplace=True)
    df.query("fGenCt<28.5 or fGenCt>28.6", inplace=True) ### Needed to remove the peak at 28.5 cm in the anchored MC
    ## fill histograms to be put at denominator of efficiency
    utils.fill_th1_hist(hPtGen, df, 'fAbsGenPt')
    utils.fill_th1_hist(hCtGen, df, 'fGenCt')
    ## now we select only the reconstructed particles
    df.query('fIsReco==True', inplace=True)

############# Apply pre-selections to data #############        
else:
    data_pre_sels = ""
    if is_matter == 'matter':
        data_pre_sels += 'fIsMatter == True'
    elif is_matter == 'antimatter':
        data_pre_sels += 'fIsMatter == False'
    if data_pre_sels != '':
        df.query(data_pre_sels, inplace=True)


############# Common filtering #############
if selections_string != '':
    df.query(selections_string, inplace=True)


############# Fill output histograms #############
utils.fill_th1_hist(hPtRec, df, 'fPt')
utils.fill_th1_hist(hCtRec, df, 'fCt')
utils.fill_th1_hist(hCosPA, df, 'fCosPA')
utils.fill_th1_hist(hRadius, df, 'fDecRad')
utils.fill_th1_hist(hDecLen, df, 'fDecLen')
utils.fill_th1_hist(hNTPCclus, df, 'fNTPCclusHe')
utils.fill_th1_hist(hNSigHe, df, 'fNSigmaHe')
utils.fill_th1_hist(hMass3LH, df, 'fMassH3L')
utils.fill_th1_hist(hMass4LH, df, 'fMassH4L')
utils.fill_th2_hist(h2MassCosPA, df, 'fCosPA', 'fMassH3L')
utils.fill_th2_hist(h2MassDecLen, df, 'fDecLen', 'fMassH3L')
utils.fill_th2_hist(h2MassDCADaughters, df, 'fDcaV0Daug', 'fMassH3L')
utils.fill_th2_hist(h2MassDCAHePv, df, 'fDcaHe', 'fMassH3L')
utils.fill_th2_hist(h2MassPt, df, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHnSigmaHe, df, 'fNSigmaHe', 'fMassH4L')

# for MC only
if mc:
    df.eval('resPt = (fPt - fAbsGenPt)/fAbsGenPt', inplace=True)
    df.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)
    utils.fill_th1_hist(hResolutionPt, df, 'resPt')
    utils.fill_th1_hist(hResolutionDecVtxX, df, 'ResDecX')
    utils.fill_th1_hist(hResolutionDecVtxY, df, 'ResDecY')
    utils.fill_th1_hist(hResolutionDecVtxZ, df, 'ResDecZ')
    utils.fill_th2_hist_abs(hResolutionPtvsPt, df, 'fAbsGenPt', 'resPt')


# save to file root
f = ROOT.TFile(f"{output_dir_name}/{output_file_name}.root", "RECREATE")

hPtRec.Write()
hCtRec.Write()
hCosPA.Write()
hRadius.Write()
hDecLen.Write()
hNTPCclus.Write()
hNSigHe.Write()
hMass3LH.Write()
hMass4LH.Write()

h2MassCosPA.Write()
h2MassDecLen.Write()
h2MassDCADaughters.Write()
h2MassDCAHePv.Write()
h2Mass4LHnSigmaHe.Write()
h2MassPt.Write()

if mc:
    f.mkdir("MC")
    f.cd("MC")
    hResolutionPt.Write()
    hResolutionPtvsPt.Write()
    hResolutionDecVtxX.Write()
    hResolutionDecVtxY.Write()
    hResolutionDecVtxZ.Write()
    hPtGen.Write()
    hCtGen.Write()

    h_eff = hPtRec.Clone("hEfficiencyPt")
    h_eff.SetTitle(";#it{p}_{T} (GeV/#it{c}); Efficiency")
    h_eff.Divide(hPtGen)
    h_eff.Write()

    h_eff_ct = hCtRec.Clone("hEfficiencyCt")
    h_eff_ct.SetTitle(";#it{c#tau} (cm); Efficiency")
    h_eff_ct.Divide(hCtGen)
    h_eff_ct.Write()
    


if not skip_out_tree:
    df.to_parquet(f"{output_dir_name}/{output_file_name}.parquet")

