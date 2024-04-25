from spectra import SpectraMaker
from hipe4ml.tree_handler import TreeHandler
import yaml
import argparse
import uproot
import numpy as np
import copy
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

import sys
sys.path.append('utils')
import utils as utils

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
    output_file_name = config['output_file'] + '_separated.root'

    if 'ct_bins' in config:
        analysis_var = 'fCt'
    else:
        analysis_var = 'fPt'

    if analysis_var == 'fCt':
        analysis_bins = config['ct_bins']
    else:
        analysis_bins = config['pt_bins']

    selections_std = config['selection']
    is_matter = config['is_matter']
    cut_dict_syst = config['cut_dict_syst']

    signal_fit_func = config['signal_fit_func']
    bkg_fit_func = config['bkg_fit_func']
    do_syst = config['do_syst']
    n_trials = config['n_trials']
    n_bins_mass_data = config['n_bins_mass_data']
    n_bins_mass_mc = config['n_bins_mass_mc']

    matter_options = ['matter', 'antimatter', 'both']
    if is_matter not in matter_options:
        raise ValueError(
            f'Invalid is-matter option. Expected one of: {matter_options}')

    print('**********************************')
    if analysis_var == 'fCt':
        print('    Running ct_analysis.py')
    else:
        print('    Running pt_analysis.py')
    print('**********************************\n')
    print("----------------------------------")
    print("** Loading data and apply preselections **")

    data_hdl = TreeHandler(input_file_name_data, 'O2hypcands', folder_name='DF*')
    mc_hdl = TreeHandler(input_file_name_mc, 'O2mchypcands', folder_name='DF*')

    # declare output file
    output_file = ROOT.TFile.Open(f'{output_dir_name}/{output_file_name}', 'recreate')

    # Add columns to the handlers
    print("Data summary:", data_hdl.print_summary())
    utils.correct_and_convert_df(data_hdl, calibrate_he3_pt=True)
    utils.correct_and_convert_df(mc_hdl, calibrate_he3_pt=True, isMC=True)

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

    # get Standard Spectrum
    standard_file = ROOT.TFile(
        f"{output_dir_name}/{config['output_file']}.root")
    std_spectrum = standard_file.Get('std/h_corrected_counts')
    std_spectrum.SetDirectory(0)
    utils.setHistStyle(std_spectrum, ROOT.kRed)
    std_corrected_counts = []
    std_corrected_counts_err = []
    for i_bin in range(1, std_spectrum.GetNbinsX()+1):
        std_corrected_counts.append(std_spectrum.GetBinContent(i_bin))
        std_corrected_counts_err.append(std_spectrum.GetBinError(i_bin))

    # reweight MC pT spectrum
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
    spectra_file.Close()
    utils.reweight_pt_spectrum(mc_hdl, 'fAbsGenPt', he3_spectrum)

    mc_hdl.apply_preselections('rej==True')
    # Needed to remove the peak at 28.5 cm in the anchored MC
    mc_hdl.apply_preselections('fGenCt < 28.5 or fGenCt > 28.6')
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)

    print("** Data loaded. ** \n")
    print("----------------------------------")

    if analysis_var == 'fCt':
        print("** Starting ct analysis **")
    else:
        print("** Starting pt analysis **")

    # get number of events
    n_ev = uproot.open(input_analysis_results_file)[
        'hyper-reco-task']['hZvtx'].values().sum()

    #########################
    #     varied cuts
    #########################

    print("** Starting systematic variations **")

    # create a dictionary with all the possible selections for a specific variable
    cut_string_dict = {}
    for var in cut_dict_syst:
        var_dict = cut_dict_syst[var]
        cut_greater = var_dict['cut_greater']
        cut_greater_string = " > " if cut_greater else " < "

        cut_list = var_dict['cut_list']
        cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
        cut_string_dict[var] = []
        for cut in cut_arr:
            cut_string_dict[var].append(var + cut_greater_string + str(cut))

    print("  ** separated cuts **")

    spectra_dict = {}
    canvas_dict = {}
    legend_dict = {}

    chi2_selection_dict = {}
    relative_error_selection_dict = {}
    outlier_selection_dict = {}

    for var, cuts in cut_string_dict.items():
        var_dir = output_file.mkdir(f'{var}')

        spectra_dict[var] = []
        chi2_selection_dict[var] = []
        relative_error_selection_dict[var] = []
        outlier_selection_dict[var] = []
        canvas_dict[var] = ROOT.TCanvas(f'c{var}', f'c{var}', 800, 600)
        legend_dict[var] = ROOT.TLegend(0.45, 0.52, 0.92, 0.86, '', 'brNDC')

        for i_cut, cut in enumerate(cuts):

            print(f'{var}: {i_cut} / {len(cuts)} ==> {cut}')

            output_dir_varied = var_dir.mkdir(f'{i_cut}')

            spectra_maker = SpectraMaker()

            spectra_maker.data_hdl = data_hdl
            spectra_maker.mc_hdl = mc_hdl
            spectra_maker.mc_reco_hdl = mc_reco_hdl

            spectra_maker.n_ev = n_ev
            spectra_maker.branching_ratio = 0.25
            spectra_maker.delta_rap = 2.0

            spectra_maker.var = analysis_var
            spectra_maker.bins = analysis_bins
            # varying the standard selections with the cut of interest
            selections_new = copy.deepcopy(selections_std)
            for element in selections_new:
                element[var] = cut
            sel_string_list = [utils.convert_sel_to_string(
                sel) for sel in selections_new]
            spectra_maker.selection_string = sel_string_list
            spectra_maker.is_matter = is_matter
            spectra_maker.n_bins_mass_data = n_bins_mass_data
            spectra_maker.n_bins_mass_mc = n_bins_mass_mc

            spectra_maker.output_dir = output_dir_varied

            fit_range = [analysis_bins[0], analysis_bins[-1]]
            spectra_maker.fit_range = fit_range

            # create raw spectra
            spectra_maker.make_spectra()
            chi2_check = spectra_maker.chi2_selection()
            chi2_selection_dict[var].append(chi2_check)

            if not chi2_check:
                print('   Rejeted for chi2')

            # draw plot for signal extraction in each bin
            data_output_dir_varied = output_dir_varied.mkdir('data')
            mc_output_dir_varied = output_dir_varied.mkdir('mc')

            data_output_dir_varied.cd()
            for i, frame in enumerate(spectra_maker.h_signal_extractions_data):
                frame.Write(f'fInvariantMass_{i}')

            mc_output_dir_varied.cd()
            for i, frame in enumerate(spectra_maker.h_signal_extractions_mc):
                frame.Write(f'fInvariantMass_{i}')

            # create corrected spectra
            spectra_maker.make_histos()
            histo = copy.deepcopy(spectra_maker.h_corrected_counts)

            relative_error_check = spectra_maker.relative_error_selection()
            relative_error_selection_dict[var].append(relative_error_check)
            if not relative_error_check:
                print('   Rejeted for large relative error of corrected counts')

            if analysis_var == 'fCt':
                histo.SetName(f'hCt{var}_{i_cut}')
            else:
                histo.SetName(f'hPt{var}_{i_cut}')
            spectra_dict[var].append(histo)
            data_output_dir_varied.cd()
            histo.Write()

            outlier_check = spectra_maker.outlier_selection(
                std_corrected_counts, std_corrected_counts_err)
            outlier_selection_dict[var].append(outlier_check)
            if not outlier_check:
                print('   Rejeted for outlier')

            del spectra_maker

    # get color paletter
    cols = ROOT.TColor.GetPalette()

    output_file.cd()
    output_file.mkdir('std')
    # std_spectrum.Write()

    for var, histos in spectra_dict.items():
        output_file.cd(f'{var}')
        canvas_dict[var].cd()
        if analysis_var == 'fCt':
            canvas_dict[var].DrawFrame(
                0., 0., 20., 3000., r';#it{ct} (cm);#frac{d#it{N}}{d(#it{ct})} (cm^{-1})')
        else:
            canvas_dict[var].DrawFrame(
                1., 0., 5., 1.5e-8, r';#it{p}_{T} (GeV/#it{c});#frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}')
        for i_histo, histo in enumerate(histos):

            if not chi2_selection_dict[var][i_histo]:
                continue

            if not relative_error_selection_dict[var][i_histo]:
                continue

            if not outlier_selection_dict[var][i_histo]:
                continue

            utils.setHistStyle(histo, cols.At(i_histo*4))
            legend_dict[var].AddEntry(
                histo, f'{cut_string_dict[var][i_histo]}', 'PE')
            histo.Draw('PE SAME')
        legend_dict[var].AddEntry(
            std_spectrum, 'std', 'PE')
        std_spectrum.Draw('PE SAME')
        legend_dict[var].Draw()
        legend_dict[var].SetNColumns(5)
        canvas_dict[var].Write()
