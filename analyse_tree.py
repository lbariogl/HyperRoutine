import ROOT
import numpy as np
import uproot
import pandas as pd
import argparse
import os
import yaml
from hipe4ml.tree_handler import TreeHandler
from signal_extraction import SignalExtraction


import sys
sys.path.append('utils')
import utils as utils

ROOT.gStyle.SetOptStat(0)


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help='path to the YAML file with configuration.', default='')
args = parser.parse_args()

# initialise parameters from parser (can be overwritten by external yaml file)

if args.config_file == '':
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
mc = config['mc']
input_files_name = config['input_files']

if 'input_analysis_results_file' not in config:
    input_analysis_results_file = ''
else:
    input_analysis_results_file = config['input_analysis_results_file']

if 'is_trigger' not in config:
    is_trigger = False
else:
    is_trigger = config['is_trigger']

output_dir_name = config['output_dir']
output_file_name = config['output_file']
selections = config['selection']
selections_string = utils.convert_sel_to_string(selections)
is_matter = config['is_matter']
is_h4l = config['is_h4l']
skip_out_tree = config['skip_out_tree']
calibrate_he_momentum = config['calibrate_he_momentum']
do_signal_extraction = config['do_signal_extraction']


matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')


print('**********************************')
print('    Running analyse_tree.py')
print('**********************************')

mass_low_limit = 3.87 if is_h4l else 2.96
mass_high_limit = 3.98 if is_h4l else 3.04

############# Create histograms #############
hCosPA = ROOT.TH1F('hCosPA', r';cos(#theta_{PA})', 500, 0.95, 1)
hNTPCclus = ROOT.TH1F('hNTPCclus', r';n TPC clusters', 80, 79.5, 159.5)
h2NTPCclusPt = ROOT.TH2F('h2NTPCclusPt', r';#it{p}_{T} (GeV/#it{c}); n TPC clusters', 50, 0, 10, 80, 79.5, 159.5)

hMass3LH = ROOT.TH1F('h_3lh_mass', r'; m({}^{3}_{#Lambda}H) (GeV/#it{c})', 40, 2.96, 3.04)
hMass4LH = ROOT.TH1F('h_4lh_mass', r';  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})', 32, 3.87, 3.98)
h2Mass3LHVvsMass4LH = ROOT.TH2F('h2Mass3LHVvsMass4LH', r'; m({}^{3}_{#Lambda}H) (GeV/#it{c}); m({}^{4}_{#Lambda}H) (GeV/#it{c})', 40, 2.96, 3.04, 32, 3.87, 3.98)
hPtRec = ROOT.TH1F('hPtRec', r';#it{p}_{T} (GeV/#it{c})', 8, np.array([1.4, 1.7, 2., 2.3, 2.6, 3, 3.5, 4., 5], dtype=np.float64))
hCtRec = ROOT.TH1F('hCtRec', r';#it{c#tau} (cm)', 50, 0, 40)
hRadius = ROOT.TH1F('hRadius', r';Radius (cm)', 50, 0, 40)
hDecLen = ROOT.TH1F('hDecLen', r';Decay length (cm)', 100, 0, 40)
hPVz = ROOT.TH1F('hPVz', r';Primary Vertex z (cm)', 100, -20, 20)
hNSigHe = ROOT.TH1F('hNSigmaHe', r';n_{#sigma}^{TPC}({}^{3}He)', 50, -3, 3)
h2NSigHe3VsMom = ROOT.TH2F('h2NSigHe3VsMom', r';{}^{3}He #it{p}_{T} (GeV/#it{c});n_{#sigma}^{TPC}({}^{3}He)', 50, -10, 10, 50, -3, 3)
h2NSigHe4VsMom = ROOT.TH2F('h2NSigHe4VsMom', r';{}^{4}He #it{p}_{T} (GeV/#it{c});n_{#sigma}^{TPC}({}^{4}He)', 50, -10, 10, 50, -3, 3)
h2Mass3HLvsPt = ROOT.TH2F('h2Mass3HLvsPt', r';#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0, 10, 50, 2.96, 3.04)
h2Mass4LHvsPt = ROOT.TH2F('h2Mass4LHvsPt', r';#it{p}_{T} (GeV/#it{c}); m({}^{4}_{#Lambda}H) (GeV/#it{c})', 100, 0, 10, 50, 3.87, 3.98)
h2TPCSigHe3VsMom = ROOT.TH2F('h2TPCSigHe3VsMom', r';{}^{3}He #it{p}_{T} (GeV/#it{c});TPC signal', 50, 0, 7, 100, 0.5, 1000)
h2TPCChi2VsPt = ROOT.TH2F('h2TPCChi2VsPt', r';#it{p}_{T} (GeV/#it{c});TPC #chi^{2}', 50, 0, 7, 100, 0, 2)
hClusterSizeHe = ROOT.TH1F('hClusterSizeHe', r';#LT Cluster size #GT', 15, 0.5, 15.5)
hTrackedClSize = ROOT.TH1F('hTrackedClSize', r';#LT Cluster size #GT', 15, 0.5, 15.5)
hClusterSizeHeCosLam = ROOT.TH1F('hClusterSizeHeCosLam', r';#LT Cluster size #GT x cos(#lambda)', 15, 0.5, 15.5)
h2ClusSizeVsCosLam = ROOT.TH2F('h2ClusSizeVsCosLam', r'; Cos(#lambda); #LT Cluster size #GT', 100, 0.95, 1, 15, 0.5, 15.5)
h2NSigClusSizeHe = ROOT.TH2F('h2NSigClusSizeHe', r';n_{#sigma}^{TPC}({}^{3}He);<Cluster size>', 50, -3, 3, 15, 0.5, 15.5)
h2TPCSigClusSize = ROOT.TH2F('h2TPCSigClusSize', r';<Cluster size>; TPC signal', 50, 0.5, 15.5, 100, 0.5, 1000)
hClusterSizePi = ROOT.TH1F('hClusterSizePi', r';#LT Cluster size #GT', 15, 0.5, 15.5)
h2NSigClusSizePi = ROOT.TH2F('h2NSigClusSizePi', r';n_{#sigma}^{TPC}(#pi); #LT Cluster size #GT', 50, -3, 3, 15, 0.5, 15.5)
hHeMomTPCMinusMomGlo = ROOT.TH2F('hHeMomTPCMinusMomGlo', r';#it{p}^{glo}/z (GeV/#it{c});(#it{p}^{TPC} - #it{p}^{Glo}) / z (GeV/#it{c})', 50, -5, 5, 50, -2, 2)
hHeMomTPCMinusMomGloTritHyp = ROOT.TH2F('hHeMomTPCMinusMomGloTritHyp', r';#it{p}^{glo}/z (GeV/#it{c});(#it{p}^{TPC} - #it{p}^{Glo}) / z (GeV/#it{c})', 50, -5, 5, 50, -2, 2)
hHeMomTPCMinusMomGloHeHyp = ROOT.TH2F('hHeMomTPCMinusMomGloHeHyp', r';#it{p}^{glo}/z (GeV/#it{c});(#it{p}^{TPC} - #it{p}^{Glo}) / z (GeV/#it{c})', 50, -5, 5, 50, -2, 2)

h2MassV2 = ROOT.TH2F('h2MassV2', r';m({}^{3}_{#Lambda}H) (GeV/#it{c}); v2', 30, mass_low_limit, mass_high_limit, 500, -1, 1)
hMeanV2VsMass = ROOT.TH1F('hMeanV2VsMass', r';m({}^{3}_{#Lambda}H) (GeV/#it{c}); #LT v2 #GT', 30, mass_low_limit, mass_high_limit)
h2MassRadius = ROOT.TH2F('h2MassRadius', r';Radius (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0, 40, 50, mass_low_limit, mass_high_limit)
h2MassCosPA = ROOT.TH2F('h2MassCosPA', r';cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0.99, 1, 50, mass_low_limit, mass_high_limit)
h2PtCosPA = ROOT.TH2F('h2PtCosPA', r';#it{p}_{T} (GeV/#it{c}); cos(#theta_{PA})', 10, 1, 6, 100, 0.98, 1)
h2MassDecLen = ROOT.TH2F('h2MassDecLen', r';Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0, 40, 50, mass_low_limit, mass_high_limit)
h2MassDCADaughters = ROOT.TH2F('h2MassDCADaughters', r';DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 200, 0, 0.3, 50, mass_low_limit, mass_high_limit)
h2MassDCAHePv = ROOT.TH2F('h2MassDCAHe', r';DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 400, -2,2, 50, mass_low_limit, mass_high_limit)
h2MassDCAPiPv = ROOT.TH2F('h2MassDCAPi', r';DCA #pi PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 400, -20,20, 50, mass_low_limit, mass_high_limit)
h2MassPt = ROOT.TH2F('h2MassPt', r';#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 50, 0, 7, 50, mass_low_limit, mass_high_limit)
h2MassPIDHypo = ROOT.TH2F('h2MassPIDHypo', r';Hypothesis; m({}^{3}_{#Lambda}H) (GeV/#it{c})', 16, 0.5, 16.5, 50, mass_low_limit, mass_high_limit)
h2Mass4LHnSigmaHe = ROOT.TH2F('h2Mass4LHnSigmaHe', r';n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})', 50, -4, 4, 30, 3.89, 3.97)
# for MC only
hPtGen = ROOT.TH1F('hPtGen', r';#it{p}_{T}^{gen} (GeV/#it{c})', 8, np.array([1.4, 1.7, 2., 2.3, 2.6, 3, 3.5, 4., 5], dtype=np.float64))
hPtGenFineBins = ROOT.TH1F('hPtGenFineBins', r';#it{p}_{T}^{gen} (GeV/#it{c})', 100, 0, 6)
hPtGenFull = ROOT.TH1F('hPtGenFull', r';#it{p}_{T}^{gen} (GeV/#it{c})', 8, np.array([1.4, 1.7, 2., 2.3, 2.6, 3, 3.5, 4., 5], dtype=np.float64))

hRadGen = ROOT.TH1F('hRadGen', r';Radius^{gen} (cm)', 50, 0, 40)
hCtGen = ROOT.TH1F('hCtGen', r';#it{c}#tau (cm)', 50, 0, 40)
hResolutionPt = ROOT.TH1F('hResolutionPt', r';(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, -0.2, 0.2)
hResolutionPtvsPt = ROOT.TH2F('hResolutionPtvsPt', r';#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, 0, 5, 50, -0.2, 0.2)
hResolutionRadvsRad = ROOT.TH2F('hResolutionRadvsRad', r';Radius^{gen} (cm);(Radius^{rec} - Radius^{gen}) / Radius^{gen}', 50, 0, 40, 50, -0.2, 0.2)
hResolutionHe3PtvsPt = ROOT.TH2F('hResolutionHe3PtvsPt', r';#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, 0, 5, 50, -0.2, 0.2)
hResolutionPiPtvsPt = ROOT.TH2F('hResolutionPiPtvsPt', r';#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, 0, 1, 50, -0.2, 0.2)
hResolutionP = ROOT.TH1F('hResolutionP', r';(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}', 50, -0.2, 0.2)
hResolutionPvsP = ROOT.TH2F('hResolutionPvsP', r';#it{p}^{gen} (GeV/#it{c});(#it{p}^{rec} - #it{p}^{gen}) / #it{p}^{gen}', 50, 0, 5, 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F('hResolutionDecVtxX', r'; Resolution Dec X', 50, -0.2, 0.2)
hResolutionDecVtxY = ROOT.TH1F('hResolutionDecVtxY', r'; Resolution Dec Y', 50, -0.2, 0.2)
hResolutionDecVtxZ = ROOT.TH1F('hResolutionDecVtxZ', r'; Resolution Dec Z', 50, -0.2, 0.2)
hHeliumPIDHypo = ROOT.TH1F('hHeliumPIDHypo', r';Hypothesis', 16, 0.5, 16.5)
hPiPIDHypo = ROOT.TH1F('hPiPIDHypo', r';Hypothesis', 16, 0.5, 16.5)


############# Read trees #############
tree_names = ['O2datahypcands','O2hypcands', 'O2hypcandsflow', 'O2mchypcands']
tree_keys = uproot.open(input_files_name[0]).keys()
for tree in tree_names:
    for key in tree_keys:
        if tree in key:
            tree_name = tree
            break
print(f'Tree name: {tree_name}')
tree_hdl = TreeHandler(input_files_name, tree_name, folder_name='DF*')
# tree_hdl = TreeHandler(input_files_name, tree_name)

df = tree_hdl.get_data_frame()
utils.correct_and_convert_df(df, calibrate_he_momentum, mc, is_h4l)


############# Apply pre-selections to MC #############
if mc:
    mc_pre_sels = 'fIsSurvEvSel==True'
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    # utils.reweight_pt_spectrum(df, 'fAbsGenPt', he3_spectrum)
    # mc_pre_sels += ' and rej==True'
    if is_matter == 'matter':
        mc_pre_sels += ' and fGenPt>0'
    elif is_matter == 'antimatter':
        mc_pre_sels += ' and fGenPt<0'
    
    utils.fill_th1_hist(hPtGenFull, df, 'fAbsGenPt')
    df.query(mc_pre_sels, inplace=True)
    ## fill histograms to be put at denominator of efficiency
    utils.fill_th1_hist(hPtGen, df, 'fAbsGenPt')
    utils.fill_th1_hist(hPtGenFineBins, df, 'fAbsGenPt')
    utils.fill_th1_hist(hCtGen, df, 'fGenCt')
    utils.fill_th1_hist(hRadGen, df, 'fGenDecRad')
    ## now we select only the reconstructed particles
    df.query('fIsReco==True', inplace=True)

############# Apply pre-selections to data #############
else:
    data_pre_sels = ''
    if is_matter == 'matter':
        data_pre_sels += 'fIsMatter == True'
    elif is_matter == 'antimatter':
        data_pre_sels += 'fIsMatter == False'
    if data_pre_sels != '':
        df.query(data_pre_sels, inplace=True)


############# Common filtering #############
if selections_string != '':
    df.query(selections_string, inplace=True)



# df.query('fAvgClusterSizeHe>4', inplace=True)
mass_string = 'fMassH3L' if not is_h4l else 'fMassH4L'

############# Fill output histograms #############
utils.fill_th1_hist(hPtRec, df, 'fPt')
utils.fill_th1_hist(hCtRec, df, 'fCt')
utils.fill_th1_hist(hCosPA, df, 'fCosPA')
utils.fill_th1_hist(hRadius, df, 'fDecRad')
utils.fill_th1_hist(hDecLen, df, 'fDecLen')
utils.fill_th1_hist(hPVz, df, 'fZPrimVtx')
utils.fill_th1_hist(hNTPCclus, df, 'fNTPCclusHe')
utils.fill_th2_hist(h2NTPCclusPt, df, 'fPt', 'fNTPCclusHe')
utils.fill_th2_hist(h2TPCChi2VsPt, df, 'fPt', 'fTPCChi2He')
utils.fill_th1_hist(hNSigHe, df, 'fNSigmaHe3')
utils.fill_th1_hist(hMass3LH, df, 'fMassH3L')
utils.fill_th1_hist(hMass4LH, df, 'fMassH4L')
utils.fill_th2_hist(h2Mass3LHVvsMass4LH, df, 'fMassH3L', 'fMassH4L')
utils.fill_th2_hist(h2MassCosPA, df, 'fCosPA', mass_string)
utils.fill_th2_hist(h2PtCosPA, df, 'fPt', 'fCosPA')
utils.fill_th2_hist(h2MassRadius, df, 'fDecRad', mass_string)
utils.fill_th2_hist(h2MassDecLen, df, 'fDecLen', mass_string)
utils.fill_th2_hist(h2MassDCADaughters, df, 'fDcaV0Daug', mass_string)
utils.fill_th2_hist(h2MassDCAHePv, df, 'fDcaHe', mass_string)
utils.fill_th2_hist(h2MassDCAPiPv, df, 'fDcaPi', mass_string)
utils.fill_th2_hist(h2MassPt, df, 'fPt', mass_string)
utils.fill_th2_hist(h2Mass4LHnSigmaHe, df, 'fNSigmaHe3', mass_string)
utils.fill_th2_hist(h2NSigHe3VsMom, df, 'fTPCSignMomHe3', 'fNSigmaHe3')
utils.fill_th2_hist(h2NSigHe4VsMom, df, 'fTPCSignMomHe3', 'fNSigmaHe4')
utils.fill_th2_hist(h2TPCSigHe3VsMom, df, 'fTPCmomHe', 'fTPCsignalHe')
utils.fill_th2_hist(h2Mass3HLvsPt, df, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHvsPt, df, 'fPt', 'fMassH4L')

df.eval('MomDiffHe3 = fTPCmomHe - fPHe3/2', inplace=True)
utils.fill_th2_hist(hHeMomTPCMinusMomGlo, df, 'fGloSignMomHe3', 'MomDiffHe3')

if 'fITSclusterSizesHe' in df.columns:
    utils.fill_th1_hist(hClusterSizeHe, df, 'fAvgClusterSizeHe')
    utils.fill_th1_hist(hClusterSizeHeCosLam, df, 'fAvgClSizeCosLambda')
    utils.fill_th1_hist(hClusterSizePi, df, 'fAvgClusterSizePi')
    utils.fill_th2_hist(h2NSigClusSizeHe, df, 'fNSigmaHe', 'fAvgClusterSizeHe')
    utils.fill_th2_hist(h2TPCSigClusSize, df, 'fAvgClusterSizeHe', 'fTPCsignalHe')
    utils.fill_th2_hist(h2TPCSigClusSize, df, 'fAvgClusterSizePi', 'fTPCsignalPi')
    utils.fill_th2_hist(h2ClusSizeVsCosLam, df, 'fCosLambdaHe', 'fAvgClusterSizeHe')

if 'fTrackedClSize' in df.columns:
    utils.fill_th1_hist(hTrackedClSize, df, 'fTrackedClSize')

if 'fFlags' in df.columns:
    utils.fill_th1_hist(hHeliumPIDHypo, df, 'fHePIDHypo')
    utils.fill_th1_hist(hPiPIDHypo, df, 'fPiPIDHypo')
    df_He_PID = df.query('fHePIDHypo==7')
    df_Trit_PID = df.query('fHePIDHypo==6')
    utils.fill_th2_hist(hHeMomTPCMinusMomGloTritHyp, df_Trit_PID, 'fGloSignMomHe3', 'MomDiffHe3')
    utils.fill_th2_hist(hHeMomTPCMinusMomGloHeHyp, df_He_PID, 'fGloSignMomHe3', 'MomDiffHe3')
    utils.fill_th2_hist(h2MassPIDHypo, df, 'fHePIDHypo', mass_string)

if 'fV2' in df.columns:
    utils.fill_th2_hist(h2MassV2, df, mass_string, 'fV2')
    ## fill the mean v2 vs mass starting from the 2D hist
    for i in range(1, h2MassV2.GetNbinsX()+1):
        bin_entries = []
        v2 = []
        for j in range(1, h2MassV2.GetNbinsY()+1):
            bin_entries.append(h2MassV2.GetBinContent(i, j))
            v2.append(h2MassV2.GetYaxis().GetBinCenter(j))
        if len(bin_entries) > 0 and np.sum(bin_entries) > 0:
            mean = np.average(v2, weights=bin_entries)
            std = np.sqrt(np.average((v2 - mean)**2, weights=bin_entries))
            hMeanV2VsMass.SetBinContent(i, mean)
            hMeanV2VsMass.SetBinError(i, std/np.sqrt(np.sum(bin_entries)))
            # print(f'bin {i}: {mean}, {std}, {np.sum(bin_entries)}, {std/np.sqrt(np.sum(bin_entries))} ')

        else:
            hMeanV2VsMass.SetBinContent(i, 0)
            hMeanV2VsMass.SetBinError(i, 0)

# for MC only
if mc:
    df.eval('resPt = (fPt - fAbsGenPt)/fAbsGenPt', inplace=True)
    df.eval('resRad = (fDecRad - fGenDecRad)/fGenDecRad', inplace=True)
    df.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)
    utils.fill_th1_hist(hResolutionPt, df, 'resPt')
    utils.fill_th1_hist(hResolutionDecVtxX, df, 'ResDecX')
    utils.fill_th1_hist(hResolutionDecVtxY, df, 'ResDecY')
    utils.fill_th1_hist(hResolutionDecVtxZ, df, 'ResDecZ')
    utils.fill_th2_hist_abs(hResolutionPtvsPt, df, 'fAbsGenPt', 'resPt')
    utils.fill_th2_hist_abs(hResolutionRadvsRad, df, 'fGenDecRad', 'resRad')


# save to file root
f = ROOT.TFile(f'{output_dir_name}/{output_file_name}.root', 'RECREATE')


hPtRec.Write()
hCtRec.Write()
hCosPA.Write()
hRadius.Write()
hDecLen.Write()
hPVz.Write()
hNTPCclus.Write()
h2NTPCclusPt.Write()
h2TPCChi2VsPt.Write()
hNSigHe.Write()
hMass3LH.Write()
hMass4LH.Write()
h2MassV2.Write()
hMeanV2VsMass.Write()
h2MassCosPA.Write()
h2PtCosPA.Write()
h2MassRadius.Write()
h2MassDecLen.Write()
h2MassDCADaughters.Write()
h2MassDCAHePv.Write()
h2MassDCAPiPv.Write()
h2Mass4LHnSigmaHe.Write()
h2MassPt.Write()
h2NSigClusSizePi.Write()
h2TPCSigClusSize.Write()
h2NSigHe3VsMom.Write()
h2TPCSigHe3VsMom.Write()
h2NSigHe4VsMom.Write()
hHeMomTPCMinusMomGlo.Write()
h2Mass3LHVvsMass4LH.Write()
h2Mass3HLvsPt.Write()
h2Mass4LHvsPt.Write()

if 'fFlags' in df.columns:
    hHeliumPIDHypo.Write()
    hPiPIDHypo.Write()
    hHeMomTPCMinusMomGloTritHyp.Write()
    hHeMomTPCMinusMomGloHeHyp.Write()
    h2MassPIDHypo.Write()

hClusterSizeHe.Write()
hTrackedClSize.Write()
hClusterSizeHeCosLam.Write()
hClusterSizePi.Write()
h2NSigClusSizeHe.Write()
h2ClusSizeVsCosLam.Write()

if mc:
    f.mkdir('MC')
    f.cd('MC')
    hResolutionPt.Write()
    hResolutionPtvsPt.Write()
    hResolutionRadvsRad.Write()
    hResolutionDecVtxX.Write()
    hResolutionDecVtxY.Write()
    hResolutionDecVtxZ.Write()
    hPtGen.Write()
    hPtGenFineBins.Write()
    hCtGen.Write()
    hRadGen.Write()

    h_eff = utils.computeEfficiency(hPtGen, hPtRec, 'hEffPt')
    h_eff.SetTitle(';#it{p}_{T} (GeV/#it{c}); Efficiency')
    h_eff.Write()

    h_eff_ct = hCtRec.Clone('hEfficiencyCt')
    h_eff_ct.SetTitle(';#it{c#tau} (cm); Efficiency')
    h_eff_ct.Divide(hCtGen)
    h_eff_ct.Write()

    h_eff_rad = hRadius.Clone('hEfficiencyRad')
    h_eff_rad.SetTitle(';Radius (cm); Efficiency')
    h_eff_rad.Divide(hRadGen)
    h_eff_rad.Write()


    h_sigloss = utils.computeEfficiency(hPtGenFull, hPtGen, 'hSigLoss')
    h_sigloss.SetTitle(';#it{p}_{T} (GeV/#it{c}); Signal loss')
    h_sigloss.Write()


    ### check if the gCt values are repeated



if not skip_out_tree:
    df.to_parquet(f'{output_dir_name}/{output_file_name}.parquet')


if do_signal_extraction:
    sign_extr_dir = f.mkdir('SignalExtraction')
    f.cd('SignalExtraction')
    signal_extraction = SignalExtraction(df)
    signal_extraction.bkg_fit_func = "pol2"
    signal_extraction.signal_fit_func = "gaus"
    signal_extraction.n_bins_data = 25
    signal_extraction.n_evts = utils.getNEvents(input_analysis_results_file, is_trigger) if input_analysis_results_file != '' else 0
    signal_extraction.matter_type = is_matter
    signal_extraction.performance = True
    signal_extraction.is_3lh = not is_h4l
    signal_extraction.out_file =  sign_extr_dir
    signal_extraction.process_fit()
f.Close()

if input_analysis_results_file != '':
    print("Number of events analysed: ", utils.getNEvents(input_analysis_results_file, is_trigger))
