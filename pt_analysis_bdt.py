from hipe4ml.tree_handler import TreeHandler
from hipe4ml.model_handler import ModelHandler
import hipe4ml.analysis_utils as au
import hipe4ml.plot_utils as pu
import xgboost as xgb
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('PDF')

from signal_extraction import SignalExtraction



from sklearn.utils import shuffle
from itertools import product
import copy
import yaml
import argparse
import uproot
import pandas as pd
import numpy as np
import os
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
kOrangeC  = ROOT.TColor.GetColor('#ff7f00')

import sys
sys.path.append('utils')
import utils as utils


## levy-tsallis is defined in the file AdditionalFunctions.h
ROOT.gROOT.ProcessLine('.L utils/AdditionalFunctions.h++')
from ROOT import LevyTsallis

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)

## input and output files
input_file_name_data = config['input_files_data']
input_file_name_mc = config['input_files_mc']
additional_training_datasets_mc = config['additional_training_datasets_mc']
input_analysis_results_file = config['input_analysis_results_file']
is_trigger = config['is_trigger']
output_dir_name = config['output_dir']
output_file_suffix = config['output_file_suffix']

debug = config['debug']
## actions to perform
do_training = config['do_training']
do_application = config['do_application']
do_yield = config['do_yield']

## training settings
preselections = config['preselections']
training_variables = config['training_variables']
train_bkg_fraction = config['train_bkg_fraction']
random_state = config['random_state']
n_cores = config['n_cores']
models_dir = config['models_dir']
ml_plots_dir = config['ml_plots_dir']

## signal extraction settings
signal_fit_func = config['signal_fit_func']
bkg_fit_func = config['bkg_fit_func']
n_bins_mass_data = config['n_bins_mass_data']
n_bins_mass_mc = config['n_bins_mass_mc']
sigma_range_mc_to_data = config['sigma_range_mc_to_data']

## general settings
calibrate_he_momentum = config['calibrate_he_momentum']
coal_based_mc = config['coal_based_mc']
pt_bins = config['pt_bins']
bin_structure = config['bin_structure']
is_matter = config['is_matter']
matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(f'Invalid is-matter option. Expected one of: {matter_options}')

## yield settings
absorption_histo_file = config['absorption_histo_file']
event_loss = config['event_loss']
signal_loss = config['signal_loss']
syst_multi_trials = config['syst_multi_trials']
absorption_syst = config['absorption_syst']
br_syst = config['br_syst']

print('**********************************')
print('    Running pt_analysis_bdt.py')
print('**********************************\n')
print("----------------------------------")
print("** Loading data and apply preselections **")

tree_names = ['O2datahypcands','O2hypcands', 'O2hypcandsflow']
tree_keys = uproot.open(input_file_name_data[0]).keys()
for tree in tree_names:
    for key in tree_keys:
        if tree in key:
            tree_name = tree
            break
print(f"Data tree found: {tree_name}")
data_hdl = TreeHandler(input_file_name_data, tree_name, folder_name='DF*')
mc_hdl = TreeHandler(input_file_name_mc, 'O2mchypcands', folder_name='DF*')

# Add columns to the handlers
utils.correct_and_convert_df(data_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=False)
utils.correct_and_convert_df(mc_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=True)
if additional_training_datasets_mc != []:
    mc_add_hdl = TreeHandler(additional_training_datasets_mc, 'O2mchypcands', folder_name='DF*')
    utils.correct_and_convert_df(mc_add_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=True)

# apply preselections + get absorption histo
matter_sel = ''
mc_matter_sel = ''
absorption_histo = None

if absorption_histo_file != '':
    absorption_file = ROOT.TFile.Open(absorption_histo_file)
    absorption_histo_mat = absorption_file.Get('x1.5/h_abso_frac_pt_mat')
    absorption_histo_anti = absorption_file.Get('x1.5/h_abso_frac_pt_antimat')
    absorption_histo_mat.SetDirectory(0)
    absorption_histo_anti.SetDirectory(0)

if is_matter == 'matter':
    matter_sel = 'fIsMatter == True'
    mc_matter_sel = 'fGenPt > 0'
    if absorption_histo_file != '':
        absorption_histo = absorption_histo_mat

elif is_matter == 'antimatter':
    matter_sel = 'fIsMatter == False'
    mc_matter_sel = 'fGenPt < 0'
    if absorption_histo_file != '':
        absorption_histo = absorption_histo_anti

if matter_sel != '':
    data_hdl.apply_preselections(matter_sel)
    mc_hdl.apply_preselections(mc_matter_sel)
    if additional_training_datasets_mc != []:
        mc_add_hdl.apply_preselections(mc_matter_sel)
    if absorption_histo_file != '':      ## get average between matter and antimatter absorption
        absorption_histo = absorption_histo_mat.Clone('h_abso_frac_pt')
        absorption_histo.Add(absorption_histo_anti)
        absorption_histo.Scale(0.5)

# reweight MC pT spectrum
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()
utils.reweight_pt_spectrum(mc_hdl, 'fAbsGenPt', he3_spectrum)
if additional_training_datasets_mc != []:
    utils.reweight_pt_spectrum(mc_add_hdl, 'fAbsGenPt', he3_spectrum)
if not coal_based_mc:
    mc_hdl.apply_preselections('rej==True')
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)
    mc_hdl_evsel = mc_hdl.apply_preselections('fIsSurvEvSel==True', inplace=False)
else:
    signal_loss = 1
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)
    mc_hdl_evsel = mc_hdl

preselections = utils.convert_sel_to_string(preselections)
data_hdl.apply_preselections(preselections)
mc_reco_hdl.apply_preselections(preselections)
if additional_training_datasets_mc != []:
    mc_add_hdl.apply_preselections(preselections)

print("** Data loaded. ** \n")
print("----------------------------------")

if do_training:
    print("** Starting training for BDT **")
    if not os.path.exists(models_dir):
        os.makedirs(models_dir)
    if not os.path.exists(ml_plots_dir):
        os.makedirs(ml_plots_dir)
    bkg_hdl = data_hdl.apply_preselections('(fMassH3L<2.95 or fMassH3L>3.02)', inplace=False)
    mc_training_hdl = copy.deepcopy(mc_reco_hdl)
    if additional_training_datasets_mc != []:
        mc_training_hdl._full_data_frame = pd.concat([mc_reco_hdl.get_data_frame(), mc_add_hdl.get_data_frame()], ignore_index=True)
    for ptbin in range(len(pt_bins) - 1):
        pt_range = f'(fPt >= {pt_bins[ptbin]}) & (fPt < {pt_bins[ptbin +1]})'
        bkg_hdl_pt = bkg_hdl.apply_preselections(pt_range, inplace=False)
        signal_hdl_pt = mc_training_hdl.apply_preselections(pt_range, inplace=False)
        if train_bkg_fraction!=None and train_bkg_fraction*len(signal_hdl_pt)<len(bkg_hdl_pt):
            bkg_hdl_pt.shuffle_data_frame(size=int(train_bkg_fraction*len(signal_hdl_pt)), inplace=True, random_state=random_state)
        print(f"Training for pT bin {ptbin}: {pt_bins[ptbin]} - {pt_bins[ptbin+1]} GeV/c")
        print("Signal candidates: ", len(signal_hdl_pt))
        print("Background candidates: ", len(bkg_hdl_pt))
        train_test_data = au.train_test_generator([signal_hdl_pt, bkg_hdl_pt], [1,0], test_size=0.5, random_state=random_state)
        distr = pu.plot_distr([bkg_hdl_pt, signal_hdl_pt], training_variables, bins=63, labels=['Signal',"Background"],colors=["blue","red"], log=True, density=True, figsize=(18, 13), alpha=0.3, grid=False)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(ml_plots_dir + f"/features_distributions_ptbin_{ptbin}.png", bbox_inches='tight')
        corr = pu.plot_corr([signal_hdl_pt,bkg_hdl_pt], training_variables + ["fMassH3L"], ['Signal',"Background"])
        corr[0].savefig(ml_plots_dir + f"/correlations_ptbin_{ptbin}.png",bbox_inches='tight')
        model_hdl = ModelHandler(xgb.XGBClassifier(n_jobs=n_cores, random_state=random_state, max_depth=3, n_estimators=100), training_variables)
        y_pred_test = model_hdl.train_test_model(train_test_data, True, True)
        print("Model trained and tested. Saving results ...")
        bdt_out_plot = pu.plot_output_train_test(model_hdl, train_test_data, 100, True, ["Signal", "Background"], True, density=True)
        bdt_out_plot.savefig(ml_plots_dir + f"/bdt_output_ptbin_{ptbin}.png")
        feature_importance_plot = pu.plot_feature_imp(train_test_data[2], train_test_data[3], model_hdl)
        feature_importance_plot[0].savefig( ml_plots_dir + f"/feature_importance_ptbin_{ptbin}.png")
        ## dump model handler and efficiencies vs score
        model_hdl.dump_model_handler(models_dir + f"/model_handler_ptbin_{ptbin}.pkl")
        eff_arr = np.round(np.arange(0.4,0.99,0.01),2)
        score_eff_arr = au.score_from_efficiency_array(train_test_data[3], y_pred_test, eff_arr)
        np.save(models_dir + f"/efficiency_arr_ptbin_{ptbin}.npy", eff_arr)
        np.save(models_dir + f"/score_efficiency_arr_ptbin_{ptbin}.npy",score_eff_arr)
        true_eff, true_score = au.bdt_efficiency_array(train_test_data[3], y_pred_test, n_points=300)
        ## plot real and interpolated efficiency vs score in the same figure
        fig = plt.figure(figsize=(8,6))
        plt.plot(true_score, true_eff, 'r-o', label='Real efficiency', markersize=2)
        plt.plot(score_eff_arr, eff_arr, 'b-', label='Interpolated efficiency')
        plt.xlabel('BDT score')
        plt.ylabel('Efficiency')
        fig.savefig(ml_plots_dir + f"/efficiency_interpolation_ptbin_{ptbin}.png")

        print("---------------------------------------------")
        del signal_hdl_pt, bkg_hdl_pt, train_test_data
    print("** Training finished ** \n")
    print("----------------------------------")

if do_application:
    output_file = ROOT.TFile.Open(output_dir_name + f"/signal_extraction_{output_file_suffix}.root", "RECREATE")
    print("** Proceeding to application **")
    for ibin in range(len(pt_bins) - 1):
        pt_dir = output_file.mkdir(f'pt_bin_{ibin}')
        ptbin = [pt_bins[ibin], pt_bins[ibin + 1]]
        for bin_struct in bin_structure:
            sel_bins = bin_struct['pt_bins']
            if ptbin in sel_bins:
                break
        output_file.cd(f'pt_bin_{ibin}')
        bdt_eff_array = np.load(models_dir + f"/efficiency_arr_ptbin_{ibin}.npy")
        bdt_score_eff_array = np.load(models_dir + f"/score_efficiency_arr_ptbin_{ibin}.npy")
        pt_range = f'(fPt >= {ptbin[0]}) & (fPt < {ptbin[1]})'
        data_hdl_pt = data_hdl.apply_preselections(pt_range, inplace=False)
        mc_hdl_pt = mc_reco_hdl.apply_preselections(pt_range, inplace=False)
        model_hdl = ModelHandler()
        model_hdl.load_model_handler(models_dir + f"/model_handler_ptbin_{ibin}.pkl")
        data_hdl_pt.apply_model_handler(model_hdl, column_name="BDT_score")
        mc_hdl_pt.apply_model_handler(model_hdl, column_name="BDT_score")

        sig_functions = bin_struct['signal_fit_func']
        bkg_functions = bin_struct['bkg_fit_func']
        for sig_func, bkg_func in product(sig_functions, bkg_functions):
            raw_counts_hist = ROOT.TH1D(f'raw_counts_ptbin_{ibin}_sig_{sig_func}_bkg_{bkg_func}', f';BDT Efficiency; Raw counts / BDT efficiency', len(bdt_eff_array), bdt_eff_array[0] - 0.005, bdt_eff_array[-1] + 0.005)
            outdir_fits = pt_dir.mkdir(f'fits_sig_{sig_func}_bkg_{bkg_func}')
            outdir_fits.cd()
            for eff, score in zip(bdt_eff_array, bdt_score_eff_array):
                data_hdl_eff = data_hdl_pt.apply_preselections(f'BDT_score >= {score}', inplace=False)
                mc_hdl_eff = mc_hdl_pt.apply_preselections(f'BDT_score >= {score}', inplace=False)
                signal_extraction = SignalExtraction(data_hdl_eff, mc_hdl_eff)
                signal_extraction.bkg_fit_func = bkg_func
                signal_extraction.signal_fit_func = sig_func
                signal_extraction.sigma_range_mc_to_data = sigma_range_mc_to_data
                signal_extraction.n_bins_data = 40
                signal_extraction.matter_type = is_matter
                signal_extraction.performance = False
                signal_extraction.out_file = outdir_fits
                signal_extraction.data_frame_fit_name = f'data_frame_fit_ptbin_eff_{int(eff*100)}'
                signal_extraction.mc_frame_fit_name = f'mc_frame_fit_ptbin_eff_{int(eff*100)}'
                fit_stats = signal_extraction.process_fit()
                if fit_stats['chi2'] < 1.2:
                    raw_counts_hist.SetBinContent(raw_counts_hist.FindBin(eff), fit_stats['signal'][0] / eff)
                    raw_counts_hist.SetBinError(raw_counts_hist.FindBin(eff), fit_stats['signal'][1] / eff)
                else:
                    raw_counts_hist.SetBinContent(raw_counts_hist.FindBin(eff), 0)
                    raw_counts_hist.SetBinError(raw_counts_hist.FindBin(eff), 0)
            pt_dir.cd()
            raw_counts_hist.Write()
    output_file.Close()
    print("** Application finished ** \n")


if do_yield:
    print("** Starting yield calculation **")
    output_yield_file = ROOT.TFile.Open(output_dir_name + f"/final_yield_{output_file_suffix}.root", "RECREATE")
    signal_extraction_file = ROOT.TFile.Open(output_dir_name + f"/signal_extraction_{output_file_suffix}.root", "READ")
    ## first we compute preselection efficiency from MC
    h_gen_pt = ROOT.TH1D('h_gen_pt', ';p_{T} (GeV/c);Counts', len(pt_bins) - 1, np.array(pt_bins))
    h_reco_pt = ROOT.TH1D('h_reco_pt', ';p_{T} (GeV/c);Counts', len(pt_bins) - 1, np.array(pt_bins))
    utils.fill_th1_hist(h_gen_pt, mc_hdl_evsel, 'fAbsGenPt')
    utils.fill_th1_hist(h_reco_pt, mc_reco_hdl, 'fPt')
    eff_hist = utils.computeEfficiency(h_gen_pt, h_reco_pt, 'PreselEff')
    output_yield_file.cd()
    eff_hist.Write()
    corr_counts_dir = output_yield_file.mkdir('corrected_counts_vs_bdt_eff')
    syst_variations_dir = output_yield_file.mkdir('syst_variations')
    n_events = utils.getNEvents(input_analysis_results_file, is_trigger)
    br = 0.25
    delta_rap = 2.0
    h_default_spectrum_stat = ROOT.TH1D('h_default_spectrum_stat', ';p_{T} (GeV/c);Counts', len(pt_bins) - 1, np.array(pt_bins))
    h_default_spectrum_syst = ROOT.TH1D('h_default_spectrum_syst', ';p_{T} (GeV/c);Counts', len(pt_bins) - 1, np.array(pt_bins))

    trials = []

    for ibin in range(len(pt_bins) - 1):
        ptbin = [pt_bins[ibin], pt_bins[ibin + 1]]
        for bin_struct in bin_structure:
            sel_bins = bin_struct['pt_bins']
            if ptbin in sel_bins:
                break
        bin_width = pt_bins[ibin + 1] - pt_bins[ibin]
        eff = eff_hist.GetBinContent(eff_hist.FindBin((pt_bins[ibin] + pt_bins[ibin + 1]) / 2))
        eff_error = eff_hist.GetBinError(eff_hist.FindBin((pt_bins[ibin] + pt_bins[ibin + 1]) / 2))
        absorption = absorption_histo.GetBinContent(absorption_histo.FindBin((pt_bins[ibin] + pt_bins[ibin + 1]) / 2)) if absorption_histo else 1

        sig_functions = bin_struct['signal_fit_func']
        bkg_functions = bin_struct['bkg_fit_func']
        bdt_eff_cut = bin_struct['bdt_eff_cut']
        bdt_syst_range = [bdt_eff_cut - bin_struct['bdt_syst_range'], bdt_eff_cut + bin_struct['bdt_syst_range']]

        syst_bin_list = []
        syst_bin_list_err = []

        for sig_func, bkg_func in product(sig_functions, bkg_functions):
            signal_extraction_file.cd()
            print(f'pt_bin_{ibin}/raw_counts_ptbin_{ibin}_sig_{sig_func}_bkg_{bkg_func}')
            raw_counts = signal_extraction_file.Get(f'pt_bin_{ibin}/raw_counts_ptbin_{ibin}_sig_{sig_func}_bkg_{bkg_func}')
            h_corrected_counts = raw_counts.Clone(f'corrected_counts_ptbin_{ibin}_sig_{sig_func}_bkg_{bkg_func}')
            for ibin_eff in range(1, h_corrected_counts.GetNbinsX() + 1):
                raw = h_corrected_counts.GetBinContent(ibin_eff)
                raw_err = h_corrected_counts.GetBinError(ibin_eff)
                if debug:
                    print('----------------------------------')
                    print(f'Raw counts at eff {h_corrected_counts.GetBinLowEdge(ibin_eff)}: {raw} +/- {raw_err}')
                    print(f'Using eff: {eff}, absorption: {absorption}, event_loss: {event_loss}, signal_loss: {signal_loss}, br: {br}, delta_rap: {delta_rap}, n_events: {n_events}, bin_width: {bin_width}')
                corrected_yield = raw * event_loss / (bin_width * eff * absorption * signal_loss * br * delta_rap * n_events) 
                corrected_yield_relative_err = np.sqrt( (raw_err / raw)**2 + (eff_error / eff)**2 ) if raw !=0 else 0
                corrected_yield_err = corrected_yield * corrected_yield_relative_err

                h_corrected_counts.SetBinContent(ibin_eff, corrected_yield)
                h_corrected_counts.SetBinError(ibin_eff, corrected_yield_err)

                if h_corrected_counts.GetBinLowEdge(ibin_eff) >= bdt_syst_range[0] and h_corrected_counts.GetBinLowEdge(ibin_eff) < bdt_syst_range[1]:
                    syst_bin_list.append(h_corrected_counts.GetBinContent(ibin_eff))
                    syst_bin_list_err.append(h_corrected_counts.GetBinError(ibin_eff))

            if sig_func == "dscb" and bkg_func == "pol2":
                if h_corrected_counts.GetBinContent(h_corrected_counts.FindBin(bdt_eff_cut)) != 0:
                    h_default_spectrum_stat.SetBinContent(ibin + 1, h_corrected_counts.GetBinContent(h_corrected_counts.FindBin(bdt_eff_cut)))
                    h_default_spectrum_stat.SetBinError(ibin + 1, h_corrected_counts.GetBinError(h_corrected_counts.FindBin(bdt_eff_cut)))
                else:
                    h_default_spectrum_stat.SetBinContent(ibin + 1, h_corrected_counts.GetBinContent(h_corrected_counts.FindBin(bdt_eff_cut - 0.01)))
                    h_default_spectrum_stat.SetBinError(ibin + 1, h_corrected_counts.GetBinError(h_corrected_counts.FindBin(bdt_eff_cut - 0.01)))

            corr_counts_dir.cd()
            h_corrected_counts.Write()
        
        ## fill histo with syst variations
        ## remove zeros
        syst_bin_list = [x for x in syst_bin_list if x > 0]
        syst_bin_list_err = [x for x in syst_bin_list_err if x > 0]
        trials.append([x for x in zip(syst_bin_list, syst_bin_list_err)])

        bin_low = np.mean(syst_bin_list) - 3 * np.std(syst_bin_list)
        bin_high = np.mean(syst_bin_list) + 3 * np.std(syst_bin_list)
        h_syst_variations = ROOT.TH1D(f'h_syst_variations_ptbin_{ibin}', ';Yield variations;Counts', 40, bin_low, bin_high)
        for val in syst_bin_list:
            h_syst_variations.Fill(val)

        syst_variations_dir.cd()
        h_syst_variations.Write()
        ## compute syst error from variations
        syst_error = h_syst_variations.GetRMS()
        h_default_spectrum_syst.SetBinContent(ibin + 1, h_default_spectrum_stat.GetBinContent(ibin + 1))
        total_syst = np.sqrt(syst_error**2 + (absorption_syst * h_default_spectrum_stat.GetBinContent(ibin + 1))**2 + (br_syst * h_default_spectrum_stat.GetBinContent(ibin + 1))**2)
        h_default_spectrum_syst.SetBinError(ibin + 1, total_syst)

    
    ## loop over trials, choose a random combination and fill the final syst histo
    n_comb = 500 
    h3l_spectrum = LevyTsallis('levy', 2.99131)
    h3l_spectrum.SetParLimits(1, 4, 30)
    h3l_spectrum.SetParLimits(3, 1e-08, 2.5e-08)

    list_yields = []

    rng = np.random.default_rng()
    for i in range(n_comb):
        combo_vals = []
        combo_errs = []
        for trial in trials:
            choice = rng.choice(trial)
            combo_vals.append(choice[0])
            combo_errs.append(choice[1])
        
        h_temp_spectrum = ROOT.TH1D('h_temp_spectrum', ';p_{T} (GeV/c);Counts', len(pt_bins) - 1, np.array(pt_bins))
        for ibin in range(len(pt_bins) - 1):
            h_temp_spectrum.SetBinContent(ibin + 1, combo_vals[ibin])
            h_temp_spectrum.SetBinError(ibin + 1, combo_errs[ibin])
        ## fit the spectrum with levy function
        fit_result = h_temp_spectrum.Fit(h3l_spectrum, 'SRQ')
        integral = h3l_spectrum.Integral(0, 10)
        if fit_result.Prob() >= 0.05 and integral < 2.4e-08:
            list_yields.append(integral)
    
    h_yields = ROOT.TH1D('h_yields', ';Yield;Counts', 50, np.min(list_yields)*0.8, np.max(list_yields)*1.2)
    for yield_val in list_yields:
        h_yields.Fill(yield_val)

    h_default_spectrum_stat.Fit(h3l_spectrum, 'R')


    output_yield_file.cd()
    h_default_spectrum_stat.Write()
    h_default_spectrum_syst.Write()
    h_yields.Write()
    absorption_histo.Write()
    output_yield_file.Close()
    print("** Yield calculation finished ** \n")
    print("----------------------------------")
            
        


   