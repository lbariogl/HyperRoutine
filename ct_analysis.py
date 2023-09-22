import ROOT
import uproot
import argparse
import yaml
import copy
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils

from spectra import SpectraMaker

if __name__ == '__main__':

    ROOT.gROOT.SetBatch(True)
    # silent mode for fits
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

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
    ct_bins = config['ct_bins']
    selections = config['selection']
    is_matter = config['is_matter']
    pt_correction_file = config['pt_correction_file']

    matter_options = ['matter', 'antimatter', 'both']
    if is_matter not in matter_options:
        raise ValueError(
            f'Invalid is-matter option. Expected one of: {matter_options}')

    print('**********************************')
    print('    Running ct_analysis.py')
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

    #########################
    #     standard cuts
    #########################

    output_dir_std = output_file.mkdir('std')

    # Spectra routine

    spectra_maker = SpectraMaker()

    spectra_maker.data_hdl = copy.deepcopy(data_hdl)
    spectra_maker.mc_hdl = copy.deepcopy(mc_hdl)
    spectra_maker.mc_reco_hdl = copy.deepcopy(mc_reco_hdl)

    spectra_maker.n_ev = uproot.open(input_analysis_results_file)[
        'hyper-reco-task']['hZvtx'].values().sum()
    spectra_maker.branching_ratio = 0.25
    spectra_maker.delta_rap = 2.0

    spectra_maker.var = 'fCt'
    spectra_maker.bins = ct_bins
    spectra_maker.selections = selections
    spectra_maker.is_matter = is_matter

    spectra_maker.output_dir = output_dir_std

    fit_range = [ct_bins[0], ct_bins[-1]]
    spectra_maker.fit_range = fit_range

    # create raw spectra
    spectra_maker.make_spectra()

    # create corrected spectra
    spectra_maker.make_histos()

    # define fit function
    expo = ROOT.TF1(
        'myexpo', '[0]*exp(-x/([1]*0.029979245800))/((exp(-[2]/([1]*0.029979245800)) - exp(-[3]/([1]*0.029979245800))) * [1]*0.029979245800)', fit_range[0], fit_range[1])
    expo.SetParLimits(1, 230, 500)
    start_bin = spectra_maker.h_corrected_counts.FindBin(fit_range[0])
    end_bin = spectra_maker.h_corrected_counts.FindBin(fit_range[1])
    expo.FixParameter(0, spectra_maker.h_corrected_counts.Integral(
        start_bin, end_bin, "width"))
    expo.FixParameter(2, fit_range[0])
    expo.FixParameter(3, fit_range[1])
    expo.SetLineColor(ROOT.kRed)

    spectra_maker.fit_func = expo
    spectra_maker.fit_options = 'MI+'
    spectra_maker.fit()

    del spectra_maker

    #########################
    #     varied cuts
    #########################

    cut_dict = {
        'fCosPA' : [0.996, 0.997, 0.998, 0.999],
        'fNSigmaHe' : [-1.99, -1.98, 1.97]
    }

    sel_kind_dict = {
        'fCosPA' : r'fCosPA > {}',
        'fNSigmaHe' : r'fNSigmaHe > {}'
    }

    print("** Starting systematic variations **")

    for var, cuts in cut_dict.items():

        for i, val in enumerate(cuts):

            print(f'{var}: {i} / {len(cuts)}')

            output_dir_varied = output_file.mkdir(f'{var}_{val}')

            spectra_maker = SpectraMaker()

            spectra_maker.data_hdl = copy.deepcopy(data_hdl)
            spectra_maker.mc_hdl = copy.deepcopy(mc_hdl)
            spectra_maker.mc_reco_hdl = copy.deepcopy(mc_reco_hdl)

            spectra_maker.n_ev = uproot.open(input_analysis_results_file)[
                'hyper-reco-task']['hZvtx'].values().sum()
            spectra_maker.branching_ratio = 0.25
            spectra_maker.delta_rap = 2.0

            spectra_maker.var = 'fCt'
            spectra_maker.bins = ct_bins
            local_sel = sel_kind_dict[var].format(val)
            if type(selections) == str:
                selections + ' and ' + local_sel
            else:
                for sel in selections:
                    sel = sel + ' and ' + local_sel
            spectra_maker.selections = selections
            spectra_maker.is_matter = is_matter

            spectra_maker.output_dir = output_dir_varied

            fit_range = [ct_bins[0], ct_bins[-1]]
            spectra_maker.fit_range = fit_range

            # create raw spectra
            spectra_maker.make_spectra()

            # create corrected spectra
            spectra_maker.make_histos()

            expo = ROOT.TF1(
                'myexpo', '[0]*exp(-x/([1]*0.029979245800))/((exp(-[2]/([1]*0.029979245800)) - exp(-[3]/([1]*0.029979245800))) * [1]*0.029979245800)', fit_range[0], fit_range[1])
            expo.SetParLimits(1, 230, 500)
            start_bin = spectra_maker.h_corrected_counts.FindBin(fit_range[0])
            end_bin = spectra_maker.h_corrected_counts.FindBin(fit_range[1])
            expo.FixParameter(0, spectra_maker.h_corrected_counts.Integral(
                start_bin, end_bin, "width"))
            expo.FixParameter(2, fit_range[0])
            expo.FixParameter(3, fit_range[1])
            expo.SetLineColor(ROOT.kRed)

            spectra_maker.fit_func = copy.deepcopy(expo)
            spectra_maker.fit_options = 'MI+'
            spectra_maker.fit()

            del spectra_maker

    output_file.Close()
