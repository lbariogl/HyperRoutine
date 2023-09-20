from signal_extraction import SignalExtraction
import ROOT
import uproot
import argparse
import yaml
import numpy as np
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils


def ct_bin_analysis(data_hdl, mc_hdl, mc_reco_hdl, ct_bins, selections, output_dir):
    raw_counts = []
    raw_counts_err = []
    efficiency = []

    for ibin in range(0, len(ct_bins) - 1):
        ct_bin = [ct_bins[ibin], ct_bins[ibin + 1]]
        ct_sel = f'fCt > {ct_bin[0]} & fCt < {ct_bin[1]}'

        # count generated per ct bin
        ct_mc_hdl = mc_hdl.apply_preselections(ct_sel, inplace=False)

        if type(selections) == str:
            ct_sel = f'{ct_sel} & {selections}'
        else:
            ct_sel = f'{ct_sel} & {selections[ibin]}'

        # select reconstructed in data and mc
        ct_data_hdl = data_hdl.apply_preselections(ct_sel, inplace=False)
        ct_mc_reco_hdl = mc_reco_hdl.apply_preselections(ct_sel, inplace=False)

        # compute efficiency
        eff = len(ct_mc_reco_hdl) / len(ct_mc_hdl)
        efficiency.append(eff)

        output_dir.mkdir(f'ct_{ct_bin[0]}_{ct_bin[1]}')
        output_dir.cd(f'ct_{ct_bin[0]}_{ct_bin[1]}')

        signal_extraction = SignalExtraction(ct_data_hdl, ct_mc_hdl)
        signal_extraction.n_bins = 30
        signal_extraction.n_evts = n_ev_plot
        signal_extraction.matter_type = is_matter
        signal_extraction.performance = False
        signal_extraction.is_3lh = True
        fit_stats = signal_extraction.process_fit()
        signal_extraction.data_frame_fit.Write()
        signal_extraction.mc_frame_fit.Write()

        raw_counts.append(fit_stats['signal'][0])
        raw_counts_err.append(fit_stats['signal'][1])

    return efficiency, raw_counts, raw_counts_err


def make_output_objects(ct_bins, efficiency, raw_counts, raw_counts_err, output_dir):

    corrected_counts = []
    corrected_counts_err = []

    h_raw_counts = ROOT.TH1D('h_raw_counts', 'h_raw_counts', len(
        ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
    h_efficiency = ROOT.TH1D('h_efficiency', 'h_efficiency', len(
        ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
    h_corrected_counts = ROOT.TH1D('h_corrected_counts', 'h_corrected_counts;#it{ct} (cm); #frac{dN}{d(#it{ct})}', len(
        ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
    h_corrected_counts.GetXaxis().SetTitleSize(0.05)
    h_corrected_counts.GetYaxis().SetTitleSize(0.05)

    for ibin in range(0, len(ct_bins) - 1):
        print(f'ibin: {ibin}, raw_count: {raw_counts[ibin]}, raw_count_err: {raw_counts[ibin]}')
        bin_width = ct_bins[ibin + 1] - ct_bins[ibin]
        h_raw_counts.SetBinContent(ibin + 1, raw_counts[ibin]/bin_width)
        h_raw_counts.SetBinError(ibin + 1, raw_counts_err[ibin]/bin_width)
        h_efficiency.SetBinContent(ibin + 1, efficiency[ibin])

        local_corrected_counts = raw_counts[ibin] / \
            efficiency[ibin] / bin_width
        local_corrected_counts_err = raw_counts_err[ibin] / \
            efficiency[ibin] / bin_width

        h_corrected_counts.SetBinContent(
            ibin + 1, local_corrected_counts)
        h_corrected_counts.SetBinError(
            ibin + 1, local_corrected_counts_err)

        corrected_counts.append(local_corrected_counts)
        corrected_counts_err.append(local_corrected_counts_err)

    fit_range = [ct_bins[0], ct_bins[-1]]
    expo = ROOT.TF1(
        'myexpo', '[0]*exp(-x/([1]*0.029979245800))/((exp(-[2]/([1]*0.029979245800)) - exp(-[3]/([1]*0.029979245800))) * [1]*0.029979245800)', fit_range[0], fit_range[1])
    expo.SetParLimits(1, 230, 500)
    start_bin = h_corrected_counts.FindBin(fit_range[0])
    end_bin = h_corrected_counts.FindBin(fit_range[1])
    expo.FixParameter(0, h_corrected_counts.Integral(
        start_bin, end_bin, "width"))
    expo.FixParameter(2, fit_range[0])
    expo.FixParameter(3, fit_range[1])

    h_corrected_counts.Fit(expo, 'MI+', '', fit_range[0], fit_range[1])
    expo.SetLineColor(ROOT.kRed)

    output_dir.cd()
    h_raw_counts.Write()
    h_efficiency.Write()
    h_corrected_counts.Write()

    return corrected_counts, corrected_counts_err


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
    output_dir_std = output_file.mkdir('std')

    # Add columns to the handlers
    utils.correct_and_convert_df(data_hdl, correction_hist)
    utils.correct_and_convert_df(mc_hdl, correction_hist, isMC=True)

    # get number of events, branching ratio and delta rapidity
    n_ev = uproot.open(input_analysis_results_file)[
        'hyper-reco-task']['hZvtx'].values().sum()
    n_ev_plot = n_ev / 1e9
    n_ev_plot = round(n_ev, 0)
    branching_ratio = 0.25
    delta_rap = 2.0

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

    ct_bin_selections = True
    if type(selections) == str:
        data_hdl.apply_preselections(selections)
        mc_hdl.apply_preselections(selections)
        ct_bin_selections = False

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

    efficiency, raw_counts, raw_counts_err = ct_bin_analysis(
        data_hdl, mc_hdl, mc_reco_hdl, ct_bins, selections, output_dir_std)

    make_output_objects(ct_bins, efficiency, raw_counts,
                        raw_counts_err, output_dir_std)

    output_file.Close()
