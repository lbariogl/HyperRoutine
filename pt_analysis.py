import ROOT
import uproot
import pandas as pd
import argparse
import yaml
import numpy as np  
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils
from signal_extraction import SignalExtraction


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
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
bkg_fit_func = config['bkg_fit_func']
is_matter = config['is_matter']
pt_correction_file = config['pt_correction_file']

matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(f'Invalid is-matter option. Expected one of: {matter_options}')

print('**********************************')
print('    Running pt_analysis.py')
print('**********************************\n')

print("** Loading data and apply preselections **")
## load data and MC
dataH = TreeHandler(input_file_name_data, 'O2datahypcands')
mcH = TreeHandler(input_file_name_mc, 'O2mchypcands')

# import correction file
correction_hist = None
if pt_correction_file:
    corr_file = ROOT.TFile(pt_correction_file)
    correction_hist = corr_file.Get('hShiftVsPtHe3')
    correction_hist.SetDirectory(0)

## declare output file
output_file = ROOT.TFile.Open(f'{output_dir_name}/{output_file_name}', 'recreate')

## Add columns to the handlers
utils.correct_and_convert_df(dataH, correction_hist)
utils.correct_and_convert_df(mcH, correction_hist)


### get number of events, branching ratio and delta rapidity
n_ev = uproot.open(input_analysis_results_file)['hyper-reco-task']['hZvtx'].values().sum()
n_ev_plot = n_ev / 1e9
n_ev_plot = round(n_ev, 0)
branching_ratio = 0.25
delta_rap = 2.0

## apply preselections
matter_sel = ''
mc_matter_sel = ''
if is_matter == 'matter':
    matter_sel = 'fIsMatter == True'
    mc_matter_sel = 'fGenPt > 0'
   
elif is_matter == 'antimatter':
    matter_sel = 'fIsMatter == False'
    mc_matter_sel = 'fGenPt < 0'

if matter_sel != '':
    dataH.apply_preselections(matter_sel)
    mcH.apply_preselections(mc_matter_sel)

pt_bin_selections = True
if type(selections) == str:
    dataH.apply_preselections(selections)
    mcH.apply_preselections(selections)
    pt_bin_selections = False

### reweight MC pT spectrum
mcH.eval_data_frame("fGenAbsPt = abs(fGenPt)")
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()
utils.reweight_pt_spectrum(mcH, 'fGenAbsPt', he3_spectrum)

mcH.apply_preselections('rej==True')
mcH_reco = mcH.apply_preselections('fIsReco == 1', inplace=False)


print("** Data loaded. ** \n")
print("** Starting pt analysis **")

raw_counts = []
raw_counts_err = []
efficiency = []


for ibin in range(0, len(pt_bins) - 1):
    pt_bin = [pt_bins[ibin], pt_bins[ibin + 1]]
    pt_sel =  f'fPt > {pt_bin[0]} & fPt < {pt_bin[1]}'
    if pt_bin_selections:
        pt_sel = f'{pt_sel} & {selections[ibin]}'
    
    pt_dataH = dataH.apply_preselections(pt_sel, inplace=False)
    pt_mcH_reco = mcH_reco.apply_preselections(pt_sel, inplace=False)
    pt_mcH = mcH.apply_preselections(f' fGenAbsPt > {pt_bin[0]} & fGenAbsPt < {pt_bin[1]}', inplace=False)
    eff = len(pt_mcH_reco) / len(pt_mcH)
    efficiency.append(eff)

    output_file.mkdir(f'pt_{pt_bin[0]}_{pt_bin[1]}')
    output_file.cd(f'pt_{pt_bin[0]}_{pt_bin[1]}')

    signal_extraction = SignalExtraction(pt_dataH, pt_mcH)
    signal_extraction.n_bins = 30
    signal_extraction.n_evts = n_ev_plot
    signal_extraction.matter_type = is_matter
    signal_extraction.performance = False
    signal_extraction.is_3lh = True
    signal_extraction.additional_pave_text = f'{pt_bin[0]} #leq #it{{p}}_{{T}} < {pt_bin[1]} GeV/#it{{c}}'
    fit_stats = signal_extraction.process_fit()
    signal_extraction.data_frame_fit.Write()
    signal_extraction.mc_frame_fit.Write()

    raw_counts.append(fit_stats['signal'][0])
    raw_counts_err.append(fit_stats['signal'][1])


### create and fill yield histogram
output_file.cd()
h_raw_counts = ROOT.TH1D('h_raw_counts', 'h_raw_counts', len(pt_bins) - 1, np.array(pt_bins, dtype=np.float64))
h_efficiency = ROOT.TH1D('h_efficiency', 'h_efficiency', len(pt_bins) - 1, np.array(pt_bins, dtype=np.float64))
h_corrected_counts = ROOT.TH1D('h_corrected_counts', 'h_corrected_counts;#it{p}_{T} (GeV/#it{c}); #frac{1}{N_{ev}} #frac{dN}{dyd#it{p}_{T}}', len(pt_bins) - 1, np.array(pt_bins, dtype=np.float64))
h_corrected_counts.GetXaxis().SetTitleSize(0.05)
h_corrected_counts.GetYaxis().SetTitleSize(0.05)

for ibin in range(0, len(pt_bins) - 1):
    bin_width = pt_bins[ibin + 1] - pt_bins[ibin]
    h_raw_counts.SetBinContent(ibin + 1, raw_counts[ibin]/bin_width / n_ev / branching_ratio / delta_rap)
    h_raw_counts.SetBinError(ibin + 1, raw_counts_err[ibin]/bin_width / n_ev / branching_ratio / delta_rap)
    h_efficiency.SetBinContent(ibin + 1, efficiency[ibin])
    h_corrected_counts.SetBinContent(ibin + 1, raw_counts[ibin] / efficiency[ibin] / bin_width / n_ev / branching_ratio / delta_rap)
    h_corrected_counts.SetBinError(ibin + 1, raw_counts_err[ibin] / efficiency[ibin] / bin_width / n_ev / branching_ratio / delta_rap)

he3_spectrum.SetParameter(0, he3_spectrum.GetParameter(0))
he3_spectrum.FixParameter(1, he3_spectrum.GetParameter(1))
he3_spectrum.FixParameter(2, he3_spectrum.GetParameter(2))
he3_spectrum.FixParameter(3, 2.99131)
he3_spectrum.SetLineColor(ROOT.kRed)
h_corrected_counts.Fit(he3_spectrum, 'R')

he3_spectrum.Write('he3_spectrum')
h_raw_counts.Write()
h_efficiency.Write()
h_corrected_counts.Write()
output_file.Close()