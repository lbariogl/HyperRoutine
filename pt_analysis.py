from spectra import SpectraMaker
from hipe4ml.tree_handler import TreeHandler
from sklearn.utils import shuffle
from itertools import product
import copy
import yaml
import argparse
import uproot
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
input_file_name_data = config['input_files_data']
input_file_name_mc = config['input_files_mc']
input_analysis_results_file = config['input_analysis_results_file']
is_trigger = config['is_trigger']

output_dir_name = config['output_dir']
output_file_name = config['output_file']

pt_bins = config['pt_bins']
bin_structure = config['bin_structure']

is_matter = config['is_matter']
calibrate_he_momentum = config['calibrate_he_momentum']

signal_fit_func = config['signal_fit_func']
bkg_fit_func = config['bkg_fit_func']
n_bins_mass_data = config['n_bins_mass_data']
n_bins_mass_mc = config['n_bins_mass_mc']
sigma_range_mc_to_data = config['sigma_range_mc_to_data']

coal_based_mc = config['coal_based_mc']
absorption_histo_file = config['absorption_histo_file']
event_loss = config['event_loss']
signal_loss = config['signal_loss']

do_syst = config['do_syst']
n_trials = config['n_trials']
absorption_syst = config['absorption_syst']
br_syst = config['br_syst']


matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(f'Invalid is-matter option. Expected one of: {matter_options}')

print('**********************************')
print('    Running pt_analysis.py')
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

# declare output file
output_file = ROOT.TFile.Open(f'{output_dir_name}/{output_file_name}.root', 'recreate')

# Add columns to the handlers
utils.correct_and_convert_df(data_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=False)
utils.correct_and_convert_df(mc_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=True)

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
    if absorption_histo_file != '':      ## get average between matter and antimatter absorption
        absorption_histo = absorption_histo_mat.Clone('h_abso_frac_pt')
        absorption_histo.Add(absorption_histo_anti)
        absorption_histo.Scale(0.5)


# reweight MC pT spectrum
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()
utils.reweight_pt_spectrum(mc_hdl, 'fAbsGenPt', he3_spectrum)
if not coal_based_mc:
    mc_hdl.apply_preselections('rej==True')
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)
    mc_hdl_evsel = mc_hdl.apply_preselections('fIsSurvEvSel==True', inplace=False)
else:
    signal_loss = 1
    mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)
    mc_hdl_evsel = mc_hdl

print("** Data loaded. ** \n")
print("----------------------------------")
print("** Starting pt analysis **")

output_dir_std = output_file.mkdir('std')
spectra_maker = SpectraMaker()

spectra_maker.data_hdl = data_hdl
spectra_maker.mc_hdl = mc_hdl_evsel
spectra_maker.mc_hdl_sign_extr = mc_hdl
spectra_maker.mc_reco_hdl = mc_reco_hdl

spectra_maker.n_ev = utils.getNEvents(input_analysis_results_file, is_trigger)
spectra_maker.branching_ratio = 0.25
spectra_maker.delta_rap = 2.0
spectra_maker.h_absorption = absorption_histo
spectra_maker.event_loss = event_loss
spectra_maker.signal_loss = signal_loss

spectra_maker.var = 'fPt'
spectra_maker.bins = pt_bins

## create an empty list of size len(pt_bins) - 1
selection_array = [None] * (len(pt_bins) - 1)
selection_string_array = [None] * (len(pt_bins) - 1)

systematic_cut_array = [None] * (len(pt_bins) - 1)
systematic_signal_shape_array = [None] * (len(pt_bins) - 1)
systematic_bkg_shape_array = [None] * (len(pt_bins) - 1)


jagged_bins = [[pt_bins[i], pt_bins[i + 1]] for i in range(len(pt_bins) - 1)]
for bin_struct in bin_structure:
    sel_bins = bin_struct['pt_bins']
    for sel_bin in sel_bins:
        for i, jagged_bin in enumerate(jagged_bins):
            if jagged_bin == sel_bin:
                selection_array[i] = bin_struct['selections']
                selection_string_array[i] = utils.convert_sel_to_string(bin_struct['selections'])
                systematic_cut_array[i] = bin_struct['systematic_cuts']
                systematic_signal_shape_array[i] = bin_struct['systematic_fit_func']['signal_fit_func']
                systematic_bkg_shape_array[i] = bin_struct['systematic_fit_func']['bkg_fit_func']


spectra_maker.selection_string = selection_string_array
spectra_maker.is_matter = is_matter
spectra_maker.inv_mass_signal_func = signal_fit_func
spectra_maker.inv_mass_bkg_func = bkg_fit_func
spectra_maker.sigma_range_mc_to_data = sigma_range_mc_to_data
spectra_maker.output_dir = output_dir_std
fit_range = [pt_bins[0], pt_bins[-1]]
spectra_maker.fit_range = fit_range
# create raw spectra
spectra_maker.make_spectra()
# create corrected spectra
spectra_maker.make_histos()
h3l_spectrum = LevyTsallis('levy', 2.99131)
# h3l_spectrum.SetParLimits(1, 10, 30)
h3l_spectrum.SetParLimits(3, 1e-08, 2.5e-08)
h3l_spectrum.SetLineColor(kOrangeC)
spectra_maker.fit_func = h3l_spectrum
spectra_maker.fit_options = 'MIQ+'
spectra_maker.fit()
spectra_maker.dump_to_output_dir()
std_corrected_counts = copy.deepcopy(spectra_maker.corrected_counts)
std_corrected_counts_err = copy.deepcopy(spectra_maker.corrected_counts_err)
final_stat = copy.deepcopy(spectra_maker.h_corrected_counts)
final_stat.SetName('hStat')
utils.setHistStyle(final_stat, ROOT.kAzure + 2)
final_syst = final_stat.Clone('hSyst')
final_syst_rms = final_stat.Clone('hSystRMS')
final_syst_rms.SetLineColor(ROOT.kAzure + 2)
final_syst_rms.SetMarkerColor(ROOT.kAzure + 2)


std_yield = spectra_maker.fit_func.Integral(0, 10)
std_yield_err = spectra_maker.fit_func.IntegralError(0, 10)
extrapolated_fraction = spectra_maker.fit_func.Integral(0, pt_bins[0]) / std_yield
std_yield_chi2 = spectra_maker.fit_func.GetChisquare() / spectra_maker.fit_func.GetNDF()
std_yield_prob = spectra_maker.fit_func.GetProb()


fit_fun_stat = copy.deepcopy(spectra_maker.fit_func)
fit_fun_stat.SetLineColor(kOrangeC)
fit_fun_stat.FixParameter(0, spectra_maker.fit_func.GetParameter(0))
fit_fun_stat.FixParameter(1, spectra_maker.fit_func.GetParameter(1))
fit_fun_stat.FixParameter(2, spectra_maker.fit_func.GetParameter(2))

yield_dist = ROOT.TH1D('hYieldSyst', ';dN/dy ;Counts', 40, 1.e-08, 3.e-08)
yield_prob = ROOT.TH1D('hYieldProb', ';prob. ;Counts', 100, 0, 1)


h_pt_syst = []
for i_bin in range(0, len(spectra_maker.bins) - 1):

    bin_label = f'{spectra_maker.bins[i_bin]}' + r' #leq #it{p}_{T} < ' f'{spectra_maker.bins[i_bin + 1]}' + r' GeV/#it{c}'
    histo = ROOT.TH1D(f'hPtSyst_{i_bin}', f'{bin_label}' + r';d#it{N} / d#it{p}_{T} (GeV/#it{c})^{-1};', 30, 0.5 * std_corrected_counts[i_bin], 2 * std_corrected_counts[i_bin])
    h_pt_syst.append(histo)

spectra_maker.del_dyn_members()

print("** pt analysis done. ** \n")


if do_syst:
    print("** Starting systematic variations **")
    n_trials = config['n_trials']
    output_dir_syst = output_file.mkdir('trials')
    # list of trial strings to be printed to a text file
    trial_strings = []
    print("----------------------------------")
    print("** Starting systematics analysis **")
    print(f'** {n_trials} trials will be tested **')
    print("----------------------------------")


    # create a dictionary with the same keys
    
    syst_combos = [None] * (len(pt_bins) - 1)
    for ipt in range(len(pt_bins) - 1):
        cut_string_dict = {}  ## sono arrivato qui
        cut_dict_syst = systematic_cut_array[ipt]
        for var in cut_dict_syst:
            var_dict = cut_dict_syst[var]
            cut_greater = var_dict['cut_greater']
            cut_greater_string = " > " if cut_greater else " < "

            cut_list = var_dict['cut_list']
            cut_arr = np.linspace(cut_list[0], cut_list[1], cut_list[2])
            cut_string_dict[var] = []
            for cut in cut_arr:
                if var_dict['cut_abs']:
                    cut_string_dict[var].append("abs(" + var + ")" + cut_greater_string + str(cut))
                else:
                    cut_string_dict[var].append(var + cut_greater_string + str(cut))
        
        cut_string_dict['signal_fit_func'] = systematic_signal_shape_array[ipt]
        cut_string_dict['bkg_fit_func'] = systematic_bkg_shape_array[ipt]

        combos = list(product(*list(cut_string_dict.values())))
        syst_combos[ipt] = combos
    

    ### check that all the combos have same length, if not take the minimum
    min_n_combos = min([len(syst_combos[ipt]) for ipt in range(len(pt_bins) - 1)])

    for ipt in range(len(pt_bins) - 1):    
        ## shuffle the combos
        if n_trials < min_n_combos:
            syst_combos[ipt] = shuffle(syst_combos[ipt], random_state=42, n_samples=n_trials)
        else:
            print(f'Warning: number of trials is smaller than the number of systematic variations for pt bin {ipt}.')
            syst_combos[ipt] = shuffle(syst_combos[ipt], random_state=42, n_samples=min_n_combos)
    

    print("Sampling ", len(syst_combos[0]), " trials, out of ", min_n_combos, " systematic variations.")
    for i_combo, _ in enumerate(syst_combos[0]):
        trial_strings.append("----------------------------------")
        trial_num_string = f'Trial: {i_combo} / {len(syst_combos[0])}'
        trial_strings.append(trial_num_string)
        print(trial_num_string)
        print("----------------------------------")

        cut_selection_list = []
        bkg_fit_func_list = []
        signal_fit_func_list = []

        for ipt in range(len(pt_bins) - 1):
            combo = syst_combos[ipt][i_combo]
            signal_fit_func = combo[-2]
            bkg_fit_func = combo[-1]
            sel_string = " & ".join(combo[: -2])

            ## add to the sel string all the variables that are not included in the syst variations but are in the std selections
            for var, sel in selection_array[ipt].items():
                if var not in sel_string:
                    sel_string += " & " + sel

            cut_selection_list.append(sel_string)
            bkg_fit_func_list.append(bkg_fit_func)
            signal_fit_func_list.append(signal_fit_func)


        trial_strings.append(str(cut_selection_list))
        trial_strings.append(str(bkg_fit_func_list))
        trial_strings.append(str(signal_fit_func_list))

        # make_spectra
        spectra_maker.selection_string = cut_selection_list
        spectra_maker.inv_mass_signal_func = signal_fit_func_list
        spectra_maker.inv_mass_bkg_func = bkg_fit_func_list
        spectra_maker.n_bins_mass_data = n_bins_mass_data
        spectra_maker.n_bins_mass_mc = n_bins_mass_mc
        spectra_maker.sigma_range_mc_to_data = sigma_range_mc_to_data
        trial_dir = output_dir_syst.mkdir(f'trial_{i_combo}')
        spectra_maker.output_dir = trial_dir
        spectra_maker.make_spectra()
        spectra_maker.make_histos()
        spectra_maker.fit()

        res_string = "Integral: " + str(spectra_maker.fit_func.Integral(0, 10)) + " Prob: " + str(spectra_maker.fit_func.GetProb())
        trial_strings.append(res_string)

        for i_bin in range(0, len(spectra_maker.bins) - 1):
            h_pt_syst[i_bin].Fill(spectra_maker.corrected_counts[i_bin])

        spectra_maker.dump_to_output_dir()
        # if spectra_maker.fit_func.GetProb() > 0.05 and spectra_maker.chi2_selection():
        yield_dist.Fill(spectra_maker.fit_func.Integral(0, 10))
        yield_prob.Fill(spectra_maker.fit_func.GetProb())

        spectra_maker.del_dyn_members()

output_dir_std.cd()

# systematic uncetrainty fo each pt bin
for i_bin in range(0, len(spectra_maker.bins) - 1):

    canvas = ROOT.TCanvas(f'cYield_{i_bin}', f'cYield_{i_bin}', 800, 600)
    canvas.SetTopMargin(0.1)
    canvas.SetBottomMargin(0.15)
    canvas.SetLeftMargin(0.08)
    canvas.SetRightMargin(0.08)
    canvas.DrawFrame(0, 0, 2 * std_corrected_counts[i_bin], 1.1 * h_pt_syst[i_bin].GetMaximum(), r';d#it{N} / d#it{p}_{T} (GeV/#it{c})^{-1};')
    # create a line for the standard value of lifetime
    std_line = ROOT.TLine(std_corrected_counts[i_bin], 0, std_corrected_counts[i_bin], 1.05 * h_pt_syst[i_bin].GetMaximum())
    std_line.SetLineColor(kOrangeC)
    std_line.SetLineWidth(2)
    # create box for statistical uncertainty
    std_errorbox = ROOT.TBox(std_corrected_counts[i_bin] - std_corrected_counts_err[i_bin], 0, std_corrected_counts[i_bin] + std_corrected_counts_err[i_bin], 1.05 * h_pt_syst[i_bin].GetMaximum())
    std_errorbox.SetFillColorAlpha(kOrangeC, 0.5)
    std_errorbox.SetLineWidth(0)
    # fitting histogram with systematic variations
    fit_func = ROOT.TF1(f'fit_func_{i_bin}', 'gaus', 0.5 *std_corrected_counts[i_bin], 1.5 * std_corrected_counts[i_bin])
    fit_func.SetLineColor(ROOT.kGreen+3)
    h_pt_syst[i_bin].Fit(fit_func, 'Q')
    syst_mu = fit_func.GetParameter(1)
    syst_mu_err = fit_func.GetParError(1)
    syst_sigma = fit_func.GetParameter(2)
    syst_rms = h_pt_syst[i_bin].GetRMS()

    if absorption_syst is not None and br_syst is not None:
        syst_sigma = np.sqrt(syst_sigma**2 + (std_corrected_counts[i_bin] * absorption_syst)**2 + (std_corrected_counts[i_bin] * br_syst)**2)
        syst_rms = np.sqrt(syst_rms**2 + (std_corrected_counts[i_bin] * absorption_syst)**2 + (std_corrected_counts[i_bin] * br_syst)**2)

    if do_syst:
        final_syst.SetBinError(i_bin+1, syst_sigma)
        final_syst_rms.SetBinError(i_bin+1, syst_rms)
    else:
        final_syst.SetBinError(i_bin+1, 0)
        final_syst_rms.SetBinError(i_bin+1, 0)


    syst_sigma_err = fit_func.GetParError(2)
    fit_param = ROOT.TPaveText(0.7, 0.6, 0.9, 0.82, 'NDC')
    fit_param.SetBorderSize(0)
    fit_param.SetFillStyle(0)
    fit_param.SetTextAlign(12)
    fit_param.SetTextFont(42)
    fit_param.AddText('#mu = ' + f'{syst_mu:.2e} #pm {syst_mu_err:.2e}' + r' (GeV/#it{c})^{-1}')
    fit_param.AddText('#sigma = ' + f'{syst_sigma:.2e} #pm {syst_sigma_err:.2e}' + r' (GeV/#it{c})^{-1}')
    fit_param.AddText('RMS = ' + f'{syst_rms:.2e} #pm {h_pt_syst[i_bin].GetRMSError():.2e}' + r' (GeV/#it{c})^{-1}')
    fit_param.AddText('standard value = ' + f'{std_corrected_counts[i_bin]:.2e} #pm {std_corrected_counts_err[i_bin]:.2e}' + r' (GeV/#it{c})^{-1}')
    # draw histogram with systematic variations
    canvas.cd()
    h_pt_syst[i_bin].Draw('HISTO SAME')
    fit_func.Draw('SAME')
    std_errorbox.Draw()
    std_line.Draw()
    fit_param.Draw()
    canvas.Write()
    canvas.SaveAs(f'{output_dir_name}/cYield_{i_bin}.pdf')

cFinalSpectrum = ROOT.TCanvas('cFinalSpectrum', 'cFinalSpectrum', 800, 600)
# define canvas between 0 and 10
cFinalSpectrum.DrawFrame(0.1e-09, 0.1e-09, 6, 1.5 * final_stat.GetMaximum(), r';#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{#it{d}N}{#it{d}y#it{d}#it{p}_{T}} (GeV/#it{c})^{-1}')
cFinalSpectrum.SetLogy()

## remove fit function from the list of functions
final_syst_rms.GetListOfFunctions().Remove(spectra_maker.fit_func)
final_stat.GetListOfFunctions().Remove(spectra_maker.fit_func)
final_stat.Draw('PEX0 SAME')
final_syst_rms.Draw('PE2 SAME')
fit_fun_stat.Draw('SAME')
cFinalSpectrum.Write()
cFinalSpectrum.SaveAs(f'{output_dir_name}/cFinalSpectrum.pdf')


final_stat.Write()
final_syst.Write()
final_syst_rms.Write()
fit_fun_stat.Write()
yield_dist.Write()
yield_prob.Write()
output_file.Close()

print("----------------------------------")
if absorption_syst is not None:
    print(f'NB: additional systematic uncertainty from absorption added: {absorption_syst}')
print("** Multi trial analysis done ** \n")

print(f'Number of events analysed: {spectra_maker.n_ev}')
print("Yield for the std selections: ", std_yield, " +- ", std_yield_err)
print("Extrapolated fraction: ", extrapolated_fraction)
print("Chi2/NDF: ", std_yield_chi2)
print("Prob: ", std_yield_prob)
print("Final fit parameters: ")
for i in range(fit_fun_stat.GetNpar()):
    print(f'Parameter {i}: {fit_fun_stat.GetParameter(i)} +- {fit_fun_stat.GetParError(i)}')


if do_syst:
    # write trial strings to a text file
    if os.path.exists(f'{output_dir_name}/{output_file_name}.txt'):
        os.remove(f'{output_dir_name}/{output_file_name}.txt')
    with open(f'{output_dir_name}/{output_file_name}.txt', 'w') as f:
        for trial_string in trial_strings:
            if isinstance(trial_string, list):
                for line in trial_string:
                    f.write("%s\n" % line)
            else:
                f.write("%s\n" % trial_string)
