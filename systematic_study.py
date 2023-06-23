import signal_extraction
import argparse
import yaml
import numpy as np
import ROOT
import uproot
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', default='configs/signal_extraction/config_signal_extraction_antimat_noCosPA.yaml',
                    help="path to the YAML file with configuration.")
args = parser.parse_args()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)

matter_type = config['matter_type']
input_parquet_data = config['input_parquet_data']
input_analysis_results = config['input_analysis_results']
input_analysis_results_mc = config['input_analysis_results_mc']
input_parquet_mc = config['input_parquet_mc']

output_file = ROOT.TFile('../results/systematic_study.root', 'recreate')

# check CosPA
cosPA_arr = np.linspace(0.985, 1., 75, dtype=np.float64)
nCosPa_bins = len(cosPA_arr) - 1

cCosPA = ROOT.TCanvas('cCosPA', 'cCosPA', 800, 600)
cCosPA.Print('../results/SystematicChecksCosPA.pdf[')

# silent mode for fits
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

# Data

hDataSigCosPA = ROOT.TH1F('hDataSigCosPA', ';cos(#theta_{PA}); signal fraction', nCosPa_bins, cosPA_arr)
utils.setHistStyle(hDataSigCosPA, ROOT.kRed+1)

for iCosPA, cosPA in enumerate(cosPA_arr[:-1]):
    presel = f'fCosPA > {cosPA:.4f}'
    label = r'cos(#theta_{PA}) > ' + f'{cosPA:.4f}'
    print(f'presel: {presel}')
    # perform fits
    _, frame_fit, signal_counts, signal_counts_err = signal_extraction.getFitFrames(matter_type, input_parquet_data, input_analysis_results,
                                                  input_parquet_mc, preselections=presel)
    cCosPA.cd()
    frame_fit.Draw()
    myLatex = ROOT.TLatex()
    myLatex.DrawLatexNDC(0.4, 0.95, label)
    cCosPA.Print('../results/SystematicChecksCosPA.pdf')
    hDataSigCosPA.SetBinContent(iCosPA+1, signal_counts)
    hDataSigCosPA.SetBinError(iCosPA+1, signal_counts_err)

hDataSigCosPA.Scale(1/hDataSigCosPA.GetBinContent(1))
cCosPA.Print('../results/SystematicChecksCosPA.pdf]')

# MC
an_vtx_z = uproot.open(input_analysis_results_mc)['hyper-reco-task']['hZvtx']
n_evts_mc = an_vtx_z.values().sum()

tree_hdl = TreeHandler(input_parquet_mc)
df = tree_hdl.get_data_frame()
hMcSigCosPA = ROOT.TH1D('hMcSigCosPA', ';cos(#theta_{PA}); signal fraction', nCosPa_bins, cosPA_arr)
utils.setHistStyle(hMcSigCosPA, ROOT.kAzure+2)

for iCosPA, cosPA in enumerate(cosPA_arr[:-1]):
    sel = f'fCosPA > {cosPA:.4f} and fMassH3L > 2.98 and fMassH3L < 3.'
    df_filtered = df.query(sel)
    signal_counts = df_filtered.shape[0]
    signal_counts_err = np.sqrt(signal_counts)
    hMcSigCosPA.SetBinContent(iCosPA+1, signal_counts)
    hMcSigCosPA.SetBinError(iCosPA+1, signal_counts_err)

hMcSigCosPA.Scale(1/hMcSigCosPA.GetBinContent(1))

output_file.cd()
hDataSigCosPA.Write()
hMcSigCosPA.Write()

cCosPAcomparison = ROOT.TCanvas('cCosPAcomparison', 'cCosPAcomparison', 800, 600)
hMcSigCosPA.Draw('PE')
hMcSigCosPA.GetYaxis().SetRangeUser(0.,1.2)
hDataSigCosPA.Draw('PE SAME')
hDataSigCosPA.GetYaxis().SetRangeUser(0.,1.2)
legend = ROOT.TLegend(0.2, 0.3, 0.4, 0.4, '', 'brNDC')
legend.SetLineWidth(0)
legend.AddEntry(hDataSigCosPA, 'Data', 'L')
legend.AddEntry(hMcSigCosPA, 'MC', 'L')
legend.Draw()
cCosPAcomparison.Write()