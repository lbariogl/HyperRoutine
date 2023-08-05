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
ct_bins = config['ct_bins']
selections = config['selection']
is_matter = config['is_matter']
pt_correction_file = config['pt_correction_file']

matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(f'Invalid is-matter option. Expected one of: {matter_options}')

print('**********************************')
print('    Running ct_analysis.py')
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
utils.correct_and_convert_df(mcH, correction_hist, isMC=True)


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

ct_bin_selections = True
if type(selections) == str:
    dataH.apply_preselections(selections)
    mcH.apply_preselections(selections)
    ct_bin_selections = False

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


for ibin in range(0, len(ct_bins) - 1):
    ct_bin = [ct_bins[ibin], ct_bins[ibin + 1]]
    ct_sel =  f'fCt > {ct_bin[0]} & fCt < {ct_bin[1]}'
    if ct_bin_selections:
        ct_sel = f'{ct_sel} & {selections[ibin]}'
    
    ct_dataH = dataH.apply_preselections(ct_sel, inplace=False)
    ct_mcH_reco = mcH_reco.apply_preselections(ct_sel, inplace=False)
    ct_mcH = mcH.apply_preselections(f'fGenCt > {ct_bin[0]} and fGenCt < {ct_bin[1]}', inplace=False)
    eff = len(ct_mcH_reco) / len(ct_mcH)
    efficiency.append(eff)

    output_file.mkdir(f'ct_{ct_bin[0]}_{ct_bin[1]}')
    output_file.cd(f'ct_{ct_bin[0]}_{ct_bin[1]}')

    signal_extraction = SignalExtraction(ct_dataH, ct_mcH)
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




### create and fill yield histogram
output_file.cd()
h_raw_counts = ROOT.TH1D('h_raw_counts', 'h_raw_counts', len(ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
h_efficiency = ROOT.TH1D('h_efficiency', 'h_efficiency', len(ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
h_corrected_counts = ROOT.TH1D('h_corrected_counts', 'h_corrected_counts;#it{ct} (cm); #frac{dN}{d(#it{ct})}', len(ct_bins) - 1, np.array(ct_bins, dtype=np.float64))
h_corrected_counts.GetXaxis().SetTitleSize(0.05)
h_corrected_counts.GetYaxis().SetTitleSize(0.05)

for ibin in range(0, len(ct_bins) - 1):
    bin_width = ct_bins[ibin + 1] - ct_bins[ibin]
    h_raw_counts.SetBinContent(ibin + 1, raw_counts[ibin]/bin_width)
    h_raw_counts.SetBinError(ibin + 1, raw_counts_err[ibin]/bin_width)
    h_efficiency.SetBinContent(ibin + 1, efficiency[ibin])
    h_corrected_counts.SetBinContent(ibin + 1, raw_counts[ibin] / efficiency[ibin] / bin_width)
    h_corrected_counts.SetBinError(ibin + 1, raw_counts_err[ibin] / efficiency[ibin] / bin_width)


fit_range = [ct_bins[0], ct_bins[-1]]
expo = ROOT.TF1('myexpo', '[0]*exp(-x/([1]*0.029979245800))/((exp(-[2]/([1]*0.029979245800)) - exp(-[3]/([1]*0.029979245800))) * [1]*0.029979245800)', fit_range[0], fit_range[1])
expo.SetParLimits(1, 230, 500)
start_bin = h_corrected_counts.FindBin(fit_range[0])
end_bin = h_corrected_counts.FindBin(fit_range[1])
expo.FixParameter(0, h_corrected_counts.Integral(start_bin, end_bin, "width"))
expo.FixParameter(2, fit_range[0])
expo.FixParameter(3, fit_range[1])
fit_result = h_corrected_counts.Fit(expo, "MI+", '', fit_range[0], fit_range[1])


h_corrected_counts.Fit(expo, 'MI+', '', fit_range[0], fit_range[1])
expo.SetLineColor(ROOT.kRed)
output_file.cd()
h_raw_counts.Write()
h_efficiency.Write()
h_corrected_counts.Write()
output_file.Close()