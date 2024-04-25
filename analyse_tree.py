import ROOT
import numpy as np
import uproot
import pandas as pd
import argparse
import os
import yaml
from hipe4ml.tree_handler import TreeHandler

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
output_dir_name = config['output_dir']
output_file_name = config['output_file']
selections = config['selection']
selections_string = utils.convert_sel_to_string(selections)
is_matter = config['is_matter']
is_h4l = config['is_h4l']
skip_out_tree = config['skip_out_tree']
calibrate_he_momentum = config['calibrate_he_momentum']


matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')


print('**********************************')
print('    Running analyse_tree.py')
print('**********************************')

############# Create histograms #############
hCosPA = ROOT.TH1F('hCosPA', r';cos(#theta_{PA})', 50, 0.95, 1)
hNTPCclus = ROOT.TH1F('hNTPCclus', r';n TPC clusters', 80, 79.5, 159.5)
h2NTPCclusPt = ROOT.TH2F('h2NTPCclusPt', r';#it{p}_{T} (GeV/#it{c}); n TPC clusters', 50, 0, 5, 80, 79.5, 159.5)

hMass3LH = ROOT.TH1F('h_3lh_mass', r'; m({}^{3}_{#Lambda}H) (GeV/#it{c})', 60, 2.96, 3.04)
hMass4LH = ROOT.TH1F('h_4lh_mass', r';  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})', 30, 3.89, 3.97)
hPtRec = ROOT.TH1F('hPtRec', r';#it{p}_{T} (GeV/#it{c})', 50, 0, 5)
hCtRec = ROOT.TH1F('hCtRec', r';#it{c#tau} (cm)', 50, 0, 40)
hRadius = ROOT.TH1F('hRadius', r';Radius (cm)', 100, 0, 40)
hDecLen = ROOT.TH1F('hDecLen', r';Decay length (cm)', 100, 0, 40)
hNSigHe = ROOT.TH1F('hNSigmaHe', r';n_{#sigma}^{TPC}({}^{3}He)', 50, -3, 3)
h2NSigHe3VsMom = ROOT.TH2F('h2NSigHe3VsMom', r';{}^{3}He #it{p}_{T} (GeV/#it{c});n_{#sigma}^{TPC}({}^{3}He)', 50, -5, 5, 50, -3, 3)
hClusterSizeHe = ROOT.TH1F('hClusterSizeHe', r';#LT Cluster size #GT', 15, 0.5, 15.5)
hClusterSizeHeCosLam = ROOT.TH1F('hClusterSizeHeCosLam', r';#LT Cluster size #GT x cos(#lambda)', 15, 0.5, 15.5)
h2NSigClusSizeHe = ROOT.TH2F('h2NSigClusSizeHe', r';n_{#sigma}^{TPC}({}^{3}He);<Cluster size>', 50, -3, 3, 15, 0.5, 15.5)
h2TPCSigClusSize = ROOT.TH2F('h2TPCSigClusSize', r';<Cluster size>; TPC signal', 50, 0.5, 15.5, 100, 0.5, 1000)
hClusterSizePi = ROOT.TH1F('hClusterSizePi', r';#LT Cluster size #GT', 15, 0.5, 15.5)
h2NSigClusSizePi = ROOT.TH2F('h2NSigClusSizePi', r';n_{#sigma}^{TPC}(#pi); #LT Cluster size #GT', 50, -3, 3, 15, 0.5, 15.5)
hHeMomTPCMinusMomGlo = ROOT.TH2F('hHeMomTPCMinusMomGlo', r';#it{p}^{glo}/z (GeV/#it{c});(#it{p}^{TPC} - #it{p}^{Glo}) / z (GeV/#it{c})', 50, -5, 5, 50, -2, 2)
h2MassV2 = ROOT.TH2F('h2MassV2', r';m({}^{3}_{#Lambda}H) (GeV/#it{c}); v2', 30, 2.96, 3.04, 500, -1, 1)
hMeanV2VsMass = ROOT.TH1F('hMeanV2VsMass', r';m({}^{3}_{#Lambda}H) (GeV/#it{c}); #LT v2 #GT', 30, 2.96, 3.04)
h2MassCosPA = ROOT.TH2F('h2MassCosPA', r';cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0.99, 1, 50, 2.96, 3.04)
h2MassDecLen = ROOT.TH2F('h2MassDecLen', r';Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0, 40, 50, 2.96, 3.04)
h2MassDCADaughters = ROOT.TH2F('h2MassDCADaughters', r';DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 200, 0, 0.3, 50, 2.96, 3.04)
h2MassDCAHePv = ROOT.TH2F('h2MassDCAHe', r';DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 100, 0, 2, 50, 2.96, 3.04)
h2MassPt = ROOT.TH2F('h2MassPt', r';#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})', 50, 0, 7, 50, 2.96, 3.04)
h2Mass4LHnSigmaHe = ROOT.TH2F('h2Mass4LHnSigmaHe', r';n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})', 50, -4, 4, 30, 3.89, 3.97)
# for MC only
hPtGen = ROOT.TH1F('hPtGen', r';#it{p}_{T}^{gen} (GeV/#it{c})', 50, 0, 5)
hCtGen = ROOT.TH1F('hCtGen', r';#it{c}#tau (cm)', 50, 0, 40)
hResolutionPt = ROOT.TH1F('hResolutionPt', r';(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, -0.2, 0.2)
hResolutionPtvsPt = ROOT.TH2F('hResolutionPtvsPt', r';#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}', 50, 0, 5, 50, -0.2, 0.2)
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
df = tree_hdl.get_data_frame()
print('Tree columns:', df.columns)
# correct and convert dataframe
utils.correct_and_convert_df(df, calibrate_he_momentum, mc, is_h4l)


############# Apply pre-selections to MC #############
if mc:
    mc_pre_sels = ''
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    utils.reweight_pt_spectrum(df, 'fAbsGenPt', he3_spectrum)
    mc_pre_sels += 'rej==True'
    if is_matter == 'matter':
        mc_pre_sels += 'and fGenPt>0'
    elif is_matter == 'antimatter':
        mc_pre_sels += 'and fGenPt<0'

    ## fill histograms to be put at denominator of efficiency
    utils.fill_th1_hist(hPtGen, df, 'fAbsGenPt')
    utils.fill_th1_hist(hCtGen, df, 'fGenCt')
    ## now we select only the reconstructed particles
    print(len(np.unique(df.round(decimals=5)['fGenCt'])))
    print(len(df))
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

############# Fill output histograms #############
utils.fill_th1_hist(hPtRec, df, 'fPt')
utils.fill_th1_hist(hCtRec, df, 'fCt')
utils.fill_th1_hist(hCosPA, df, 'fCosPA')
utils.fill_th1_hist(hRadius, df, 'fDecRad')
utils.fill_th1_hist(hDecLen, df, 'fDecLen')
utils.fill_th1_hist(hNTPCclus, df, 'fNTPCclusHe')
utils.fill_th2_hist(h2NTPCclusPt, df, 'fPt', 'fNTPCclusHe')
utils.fill_th1_hist(hNSigHe, df, 'fNSigmaHe')
utils.fill_th1_hist(hMass3LH, df, 'fMassH3L')
utils.fill_th1_hist(hMass4LH, df, 'fMassH4L')
utils.fill_th2_hist(h2MassCosPA, df, 'fCosPA', 'fMassH3L')
utils.fill_th2_hist(h2MassDecLen, df, 'fDecLen', 'fMassH3L')
utils.fill_th2_hist(h2MassDCADaughters, df, 'fDcaV0Daug', 'fMassH3L')
utils.fill_th2_hist(h2MassDCAHePv, df, 'fDcaHe', 'fMassH3L')
utils.fill_th2_hist(h2MassPt, df, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHnSigmaHe, df, 'fNSigmaHe', 'fMassH4L')
utils.fill_th2_hist(h2NSigHe3VsMom, df, 'fTPCSignMomHe3', 'fNSigmaHe')

df.eval('MomDiffHe3 = fTPCmomHe - fPHe3/2', inplace=True)
utils.fill_th2_hist(hHeMomTPCMinusMomGlo, df, 'fGloSignMomHe3', 'MomDiffHe3')

if 'fITSclusterSizesHe' in df.columns:
    utils.fill_th1_hist(hClusterSizeHe, df, 'fAvgClusterSizeHe')
    utils.fill_th1_hist(hClusterSizeHeCosLam, df, 'fAvgClSizeCosLambda')
    utils.fill_th1_hist(hClusterSizePi, df, 'fAvgClusterSizePi')
    utils.fill_th2_hist(h2NSigClusSizeHe, df, 'fNSigmaHe', 'fAvgClusterSizeHe')
    utils.fill_th2_hist(h2TPCSigClusSize, df, 'fAvgClusterSizeHe', 'fTPCsignalHe')
    utils.fill_th2_hist(h2TPCSigClusSize, df, 'fAvgClusterSizePi', 'fTPCsignalPi')


if 'fFlags' in df.columns:
    df['fHePIDHypo'] = np.right_shift(df['fFlags'], 4)
    df['fPiPIDHypo'] = np.bitwise_and(df['fFlags'], 0b1111)
    utils.fill_th1_hist(hHeliumPIDHypo, df, 'fHePIDHypo')
    utils.fill_th1_hist(hPiPIDHypo, df, 'fPiPIDHypo')

if 'fV2' in df.columns:
    utils.fill_th2_hist(h2MassV2, df, 'fMassH3L', 'fV2')
    ## fill the mean v2 vs mass starting from the 2D hist
    for i in range(1, h2MassV2.GetNbinsX()+1):
        bin_entries = []
        v2 = []
        for j in range(1, h2MassV2.GetNbinsY()+1):
            bin_entries.append(h2MassV2.GetBinContent(i, j))
            v2.append(h2MassV2.GetYaxis().GetBinCenter(j))
        if len(bin_entries) > 0:
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
    df.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
    df.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
    df.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)
    utils.fill_th1_hist(hResolutionPt, df, 'resPt')
    utils.fill_th1_hist(hResolutionDecVtxX, df, 'ResDecX')
    utils.fill_th1_hist(hResolutionDecVtxY, df, 'ResDecY')
    utils.fill_th1_hist(hResolutionDecVtxZ, df, 'ResDecZ')
    utils.fill_th2_hist_abs(hResolutionPtvsPt, df, 'fAbsGenPt', 'resPt')


# save to file root
f = ROOT.TFile(f'{output_dir_name}/{output_file_name}.root', 'RECREATE')


hPtRec.Write()
# utils.saveCanvasAsPDF(hPtRec, f'{output_dir_name}/qa_plots')
hCtRec.Write()
# utils.saveCanvasAsPDF(hCtRec, f'{output_dir_name}/qa_plots')
hCosPA.Write()
# utils.saveCanvasAsPDF(hCosPA, f'{output_dir_name}/qa_plots')
hRadius.Write()
# utils.saveCanvasAsPDF(hRadius, f'{output_dir_name}/qa_plots')
hDecLen.Write()
# utils.saveCanvasAsPDF(hDecLen, f'{output_dir_name}/qa_plots')
hNTPCclus.Write()
# utils.saveCanvasAsPDF(hNTPCclus, f'{output_dir_name}/qa_plots')
h2NTPCclusPt.Write()
hNSigHe.Write()
# utils.saveCanvasAsPDF(hNSigHe, f'{output_dir_name}/qa_plots')
hMass3LH.Write()
# utils.saveCanvasAsPDF(hMass3LH, f'{output_dir_name}/qa_plots')
hMass4LH.Write()
# utils.saveCanvasAsPDF(hMass4LH, f'{output_dir_name}/qa_plots')
h2MassV2.Write()
# utils.saveCanvasAsPDF(hPtRec, f'{output_dir_name}/qa_plots', is2D=True)
hMeanV2VsMass.Write()
# utils.saveCanvasAsPDF(hMeanV2VsMass, f'{output_dir_name}/qa_plots', is2D=True)

h2MassCosPA.Write()
# utils.saveCanvasAsPDF(h2MassCosPA, f'{output_dir_name}/qa_plots', is2D=True)
h2MassDecLen.Write()
# utils.saveCanvasAsPDF(h2MassDecLen, f'{output_dir_name}/qa_plots', is2D=True)
h2MassDCADaughters.Write()
# utils.saveCanvasAsPDF(h2MassDCADaughters, f'{output_dir_name}/qa_plots', is2D=True)
h2MassDCAHePv.Write()
# utils.saveCanvasAsPDF(h2MassDCAHePv, f'{output_dir_name}/qa_plots', is2D=True)
h2Mass4LHnSigmaHe.Write()
# utils.saveCanvasAsPDF(h2Mass4LHnSigmaHe, f'{output_dir_name}/qa_plots', is2D=True)
h2MassPt.Write()
# utils.saveCanvasAsPDF(h2NSigClusSizeHe, f'{output_dir_name}/qa_plots', is2D=True)
h2NSigClusSizeHe.Write()
# utils.saveCanvasAsPDF(h2NSigClusSizeHe, f'{output_dir_name}/qa_plots', is2D=True)
h2NSigClusSizePi.Write()
# utils.saveCanvasAsPDF(h2NSigClusSizePi, f'{output_dir_name}/qa_plots', is2D=True)
h2TPCSigClusSize.Write()
# utils.saveCanvasAsPDF(h2TPCSigClusSize, f'{output_dir_name}/qa_plots', is2D=True)
h2NSigHe3VsMom.Write()
# utils.saveCanvasAsPDF(h2NSigHe3VsMom, f'{output_dir_name}/qa_plots', is2D=True)
hHeMomTPCMinusMomGlo.Write()
# utils.saveCanvasAsPDF(hHeMomTPCMinusMomGlo, f'{output_dir_name}/qa_plots', is2D=True)

if 'fFlags' in df.columns:
    hHeliumPIDHypo.Write()
    # utils.saveCanvasAsPDF(hHeliumPIDHypo, f'{output_dir_name}/qa_plots')
    hPiPIDHypo.Write()
    # utils.saveCanvasAsPDF(hPiPIDHypo, f'{output_dir_name}/qa_plots')

hClusterSizeHe.Write()
# utils.saveCanvasAsPDF(hClusterSizeHe, f'{output_dir_name}/qa_plots')
hClusterSizeHeCosLam.Write()
# utils.saveCanvasAsPDF(hClusterSizeHeCosLam, f'{output_dir_name}/qa_plots')
hClusterSizePi.Write()
# utils.saveCanvasAsPDF(hClusterSizePi, f'{output_dir_name}/qa_plots')

h2NSigClusSizeHe.Write()

if mc:
    f.mkdir('MC')
    f.cd('MC')
    hResolutionPt.Write()
    # utils.saveCanvasAsPDF(hResolutionPt, f'{output_dir_name}/qa_plots')
    hResolutionPtvsPt.Write()
    # utils.saveCanvasAsPDF(hResolutionPtvsPt, f'{output_dir_name}/qa_plots', is2D=True)
    hResolutionDecVtxX.Write()
    # utils.saveCanvasAsPDF(hResolutionDecVtxX, f'{output_dir_name}/qa_plots')
    hResolutionDecVtxY.Write()
    # utils.saveCanvasAsPDF(hResolutionDecVtxY, f'{output_dir_name}/qa_plots')
    hResolutionDecVtxZ.Write()
    # utils.saveCanvasAsPDF(hResolutionDecVtxZ, f'{output_dir_name}/qa_plots')
    hPtGen.Write()
    # utils.saveCanvasAsPDF(hPtGen, f'{output_dir_name}/qa_plots')
    hCtGen.Write()
    # utils.saveCanvasAsPDF(hCtGen, f'{output_dir_name}/qa_plots')

    h_eff = hPtRec.Clone('hEfficiencyPt')
    h_eff.SetTitle(';#it{p}_{T} (GeV/#it{c}); Efficiency')
    h_eff.Divide(hPtGen)
    h_eff.Write()

    h_eff_ct = hCtRec.Clone('hEfficiencyCt')
    h_eff_ct.SetTitle(';#it{c#tau} (cm); Efficiency')
    h_eff_ct.Divide(hCtGen)
    h_eff_ct.Write()


    ### check if the gCt values are repeated



if not skip_out_tree:
    df.to_parquet(f'{output_dir_name}/{output_file_name}.parquet')

