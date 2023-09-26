import ROOT
import uproot
import argparse
import yaml
import copy
from itertools import product
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
    parser.add_argument('--separated', dest='separated', action='store_true',
                    help="if True one variable a time is varied.", default=False)

    args = parser.parse_args()
    if args.config_file == "":
        print('** No config file provided. Exiting. **')
        exit()

    separated = args.separated

    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    input_file_name_data = config['input_files_data']
    input_file_name_mc = config['input_files_mc']
    input_analysis_results_file = config['input_analysis_results_file']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']
    pt_bins = config['pt_bins']
    selections_std = config['selection']
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

    spectra_maker.var = 'fPt'
    spectra_maker.bins = pt_bins
    spectra_maker.selections = selections_std
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

    corrected_counts_std = copy.deepcopy(spectra_maker.corrected_counts)
    corrected_counts_err_std = copy.deepcopy(spectra_maker.corrected_counts_err)

    del spectra_maker

    yield_histos = []
    bin_labels = []

    for ibin in range(0, len(pt_bins) - 1):
        bin = [pt_bins[ibin], pt_bins[ibin + 1]]
        bin_label = f'fPt_{bin[0]}_{bin[1]}'
        histo_title = str(bin[0]) + r' < #it{{ct}} < ' + str(bin[1]) + r'GeV/#it{c}; #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}'
        histo = ROOT.TH1D(f'hYield_{bin_label}', histo_title, 50, 0.5 * corrected_counts_std[ibin], 1.5 * corrected_counts_std[ibin])
        yield_histos.append(histo)
        bin_labels.append(bin_label)

    #########################
    #     varied cuts
    #########################

    cut_fCosPA = [f'fCosPA > {0.98 + i * 0.001}' for i in range(0, 20)]
    cut_fNSigmaHe = [f'fNSigmaHe > {-3.5 + i * 0.25}' for i in range(0, 9)]

    cut_dict = {'fCosPA': cut_fCosPA,
                'fNSigmaHe': cut_fNSigmaHe}

    print("** Starting systematic variations **")

    if separated:

        print("  ** separated cuts **")

        for var, cuts in cut_dict.items():

            for i_cut, cut in enumerate(cuts):

                print(f'{var}: {i_cut} / {len(cuts)}')

                output_dir_varied = output_file.mkdir(f'{var}_{i_cut}')

                spectra_maker = SpectraMaker()

                spectra_maker.data_hdl = copy.deepcopy(data_hdl)
                spectra_maker.mc_hdl = copy.deepcopy(mc_hdl)
                spectra_maker.mc_reco_hdl = copy.deepcopy(mc_reco_hdl)

                spectra_maker.n_ev = uproot.open(input_analysis_results_file)[
                    'hyper-reco-task']['hZvtx'].values().sum()
                spectra_maker.branching_ratio = 0.25
                spectra_maker.delta_rap = 2.0

                spectra_maker.var = 'fPt'
                spectra_maker.bins = pt_bins
                spectra_maker.selections = copy.deepcopy(selections_std)
                spectra_maker.vary_selection(var, cut)
                spectra_maker.is_matter = is_matter

                spectra_maker.output_dir = output_dir_varied

                # create raw spectra
                spectra_maker.make_spectra()

                # create corrected spectra
                spectra_maker.make_histos()

                for ibin in range(0, len(pt_bins) - 1):
                    yield_histos[ibin].Fill(spectra_maker.corrected_counts[ibin])

                he3_spectrum.SetParameter(0, he3_spectrum.GetParameter(0))
                he3_spectrum.FixParameter(1, he3_spectrum.GetParameter(1))
                he3_spectrum.FixParameter(2, he3_spectrum.GetParameter(2))
                he3_spectrum.FixParameter(3, 2.99131)
                he3_spectrum.SetLineColor(ROOT.kRed)

                spectra_maker.fit_func = copy.deepcopy(he3_spectrum)
                spectra_maker.fit()

                del spectra_maker

    else:

        print("  ** mixed cuts **")

        # create all possible combinations within cut_dict
        combos = list(product(*list(cut_dict.values())))

        for i_combo, combo in enumerate(combos):

            print(f'mixed: {i_combo} / {len(combos)}')

            output_dir_varied = output_file.mkdir(f'mixed_{i_combo}')

            spectra_maker = SpectraMaker()

            spectra_maker.data_hdl = copy.deepcopy(data_hdl)
            spectra_maker.mc_hdl = copy.deepcopy(mc_hdl)
            spectra_maker.mc_reco_hdl = copy.deepcopy(mc_reco_hdl)

            spectra_maker.n_ev = uproot.open(input_analysis_results_file)[
                'hyper-reco-task']['hZvtx'].values().sum()
            spectra_maker.branching_ratio = 0.25
            spectra_maker.delta_rap = 2.0

            spectra_maker.var = 'fPt'
            spectra_maker.bins = pt_bins
            spectra_maker.selections = copy.deepcopy(selections_std)
            for i_var, var in enumerate(cut_dict.keys()):
                spectra_maker.vary_selection(var, combo[i_var])
            spectra_maker.is_matter = is_matter

            spectra_maker.output_dir = output_dir_varied

            # create raw spectra
            spectra_maker.make_spectra()

            # create corrected spectra
            spectra_maker.make_histos()

            for ibin in range(0, len(pt_bins) - 1):
                yield_histos[ibin].Fill(spectra_maker.corrected_counts[ibin])

            he3_spectrum.SetParameter(0, he3_spectrum.GetParameter(0))
            he3_spectrum.FixParameter(1, he3_spectrum.GetParameter(1))
            he3_spectrum.FixParameter(2, he3_spectrum.GetParameter(2))
            he3_spectrum.FixParameter(3, 2.99131)
            he3_spectrum.SetLineColor(ROOT.kRed)

            spectra_maker.fit_func = copy.deepcopy(he3_spectrum)
            spectra_maker.fit()

            del spectra_maker

    for ibin in range(0, len(pt_bins) - 1):
        output_dir_std.cd(bin_labels[ibin])
        yield_histos[ibin].Write()

    output_file.Close()
