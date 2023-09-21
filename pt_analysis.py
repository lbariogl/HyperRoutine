import ROOT
import uproot
import argparse
import yaml
import numpy as np
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils

from spectra import SpectraMaker

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='config_file',
                        help="path to the YAML file with configuration.", default='')

    args = parser.parse_args()
    if args.config_file == "":
        print('** No config file provided. Exiting. **')
        exit()

    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    input_file_name_data = config['input_files_data']
    input_file_name_mc = config['input_files_mc']
    input_analysis_results_file = config['input_analysis_results_file']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']
    pt_bins = config['pt_bins']
    selections = config['selection']
    is_matter = config['is_matter']
    pt_correction_file = config['pt_correction_file']

    matter_options = ['matter', 'antimatter', 'both']
    if is_matter not in matter_options:
        raise ValueError(
            f'Invalid is-matter option. Expected one of: {matter_options}')

    print('**********************************')
    print('    Running pt_analysis.py')
    print('**********************************\n')

    print("** Loading data and apply preselections **")
    # load data and MC
    data_hdl = TreeHandler(input_file_name_data, 'O2datahypcands')
    mc_hdl = TreeHandler(input_file_name_mc, 'O2mchypcands')

    # import correction file
    correction_hist = None
    if pt_correction_file:
        corr_file = ROOT.TFile(pt_correction_file)
        correction_hist = corr_file.Get('hShiftVsPtHe3')
        correction_hist.SetDirectory(0)

    # declare output file
    output_file = ROOT.TFile.Open(
        f'{output_dir_name}/{output_file_name}', 'recreate')
    output_dir_std = output_file.mkdir('std')

    # Add columns to the handlers
    utils.correct_and_convert_df(data_hdl, correction_hist)
    utils.correct_and_convert_df(mc_hdl, correction_hist, isMC=True)

    # apply preselections
    matter_sel = ''
    mc_matter_sel = ''
    if is_matter == 'matter':
        matter_sel = 'fIsMatter == True'
        mc_matter_sel = 'fGenPt > 0'

    elif is_matter == 'antimatter':
        matter_sel = 'fIsMatter == False'
        mc_matter_sel = 'fGenPt < 0'

    if matter_sel != '':
        data_hdl.apply_preselections(matter_sel)
        mc_hdl.apply_preselections(mc_matter_sel)

    # reweight MC pT spectrum
    mc_hdl.eval_data_frame("fGenAbsPt = abs(fGenPt)")
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    utils.reweight_pt_spectrum(mc_hdl, 'fGenAbsPt', he3_spectrum)

    mc_hdl.apply_preselections('rej==True')
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)

    print("** Data loaded. ** \n")
    print("** Starting pt analysis **")

    # Spectra routine

    spectra_maker = SpectraMaker()

    spectra_maker.data_hdl = data_hdl
    spectra_maker.mc_hdl = mc_hdl
    spectra_maker.mc_reco_hdl = mc_reco_hdl

    spectra_maker.n_ev = uproot.open(input_analysis_results_file)[
        'hyper-reco-task']['hZvtx'].values().sum()
    spectra_maker.branching_ratio = 0.25
    spectra_maker.delta_rap = 2.0

    spectra_maker.var = 'fPt'
    spectra_maker.bins = pt_bins
    spectra_maker.selections = selections
    spectra_maker.is_matter = is_matter

    spectra_maker.output_dir = output_dir_std

    # create raw spectra
    spectra_maker.make_spectra()

    # create corrected spectra
    spectra_maker.make_histos()

    # define fit function

    he3_spectrum.SetParameter(0, he3_spectrum.GetParameter(0))
    he3_spectrum.FixParameter(1, he3_spectrum.GetParameter(1))
    he3_spectrum.FixParameter(2, he3_spectrum.GetParameter(2))
    he3_spectrum.FixParameter(3, 2.99131)
    he3_spectrum.SetLineColor(ROOT.kRed)

    spectra_maker.fit_func = he3_spectrum
    spectra_maker.fit()

    output_file.Close()
