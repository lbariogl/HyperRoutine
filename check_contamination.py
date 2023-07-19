import ROOT
import uproot
from hipe4ml.tree_handler import TreeHandler
import numpy as np

import argparse
import yaml

import sys
sys.path.append('utils')
import utils as utils

import signal_extraction

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--input-files', dest='input_files',
                    help="path to the input files.", default='../data/AO2D_merged.root')
parser.add_argument('--selection', dest='selection', help="selections to be bassed as query.",
                    default='fCosPA > 0.998 & fNTPCclusHe > 110 & abs(fDcaHe) > 0.1')
parser.add_argument('--is-matter', dest='is_matter',
                    help="path to the YAML file with configuration.", default='matter')

parser.add_argument('--skip-out-tree', dest='skip_out_tree', action='store_true', help="if True do not save output tree.")

parser.add_argument('--correction-file', dest='correction_file',
                    help="path to the file use for 3He pt correction.", default=None)

parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.", default='')
args = parser.parse_args()

# initialise parameters from parser (can be overwritten by external yaml file)
skip_out_tree = args.skip_out_tree
input_file_name = args.input_files
selections = args.selection
is_matter = args.is_matter
correction_file = args.correction_file

if args.config_file != "":
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    input_files_name = config['input_files']
    selections = config['selection']
    is_matter = config['is_matter']
    correction_file = config['correction_file']

matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')

# creating the dataframe
tree_name = 'O2datahypcands'
tree_hdl = TreeHandler(input_files_name, tree_name)
df = tree_hdl.get_data_frame()

# import correction file
correction_hist = None
# if correction_file:
corr_file = ROOT.TFile(correction_file)
correction_hist = corr_file.Get('hShiftVsPtHe3')
correction_hist.SetDirectory(0)

# try to convert
utils.correct_and_convert_df(df, correction_hist)

# add new columns
df.eval('fP = fPt * cosh(fEta)', inplace=True)
df.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)

if selections == '':
    if is_matter == 'matter':
        selections = 'fIsMatter == True'
    elif is_matter == 'antimatter':
        selections = 'fIsMatter == False'
else:
    if is_matter == 'matter':
        selections = selections + ' & fIsMatter == True'
    elif is_matter == 'antimatter':
        selections = selections + ' & fIsMatter == False'

# filtering
if selections != '':
    df_filtered = df.query(selections)
else:
    df_filtered = df

h2MassMass = ROOT.TH2F('h2MassMass',';m({}^{3}_{#Lambda}H) (GeV/#it{c}^{2});m({}^{4}_{#Lambda}H) (GeV/#it{c}^{2})', 40, 2.96, 3.04, 30, 3.89, 3.97)

h2NsigmaMass = ROOT.TH2F('h2NsigmaMass',';n_{#sigma}^{TPC}({}^{3}He);m({}^{4}_{#Lambda}H) (GeV/#it{c}^{2})', 50, 0.7, 6.7, 30, 3.89, 3.97)

utils.fill_th2_hist(h2MassMass, df_filtered, 'fMassH3L', 'fMassH4L')
utils.fill_th2_hist(h2NsigmaMass, df_filtered, 'fNSigmaHe', 'fMassH4L')

parquet_data = '../results/check_contamination.parquet'
df_filtered.to_parquet(parquet_data)

parquet_mc = '../results/HypertritonResults_antimatter_MC.root.parquet'

_, frame_fit, signal_counts, signal_counts_err = signal_extraction.getFitFrames(is_matter, parquet_data, '/data/shared/hyp_run_3/apass4/new_dataformat/AnalysisResults.root', parquet_mc, preselections='', print_info=False)

cContamina = ROOT.TCanvas('')

output_file = ROOT.TFile('../results/check_contamination.root', 'recreate')
h2MassMass.Write()
h2NsigmaMass.Write()
