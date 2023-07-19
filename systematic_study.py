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

## ROOT batch mode
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', default='configs/signal_extraction/config_signal_extraction_antimat.yaml',
                    help="path to the YAML file with configuration.")
args = parser.parse_args()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)

matter_type = config['matter_type']
input_parquet_data = config['input_parquet_data']
input_analysis_results = config['input_analysis_results']
input_parquet_mc = config['input_parquet_mc']

output_file = ROOT.TFile('../results/systematic_study.root', 'recreate')

# common info for MC
tree_hdl = TreeHandler(input_parquet_mc)

##apply pT rejection
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()
tree_hdl.eval_data_frame("fAbsGenPt = abs(fGenPt)")
utils.reweight_pt_spectrum(tree_hdl, 'fAbsGenPt', he3_spectrum)
tree_hdl.apply_preselections("rej==True and fIsReco==True")

df = tree_hdl.get_data_frame()

# silent mode for fits
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

def systematic_routine(var, arr, sel_string, histo_data, histo_mc, canvas, normalise_to_first=True):
    for i, val in enumerate(arr):
        # data
        presel = sel_string.format(var, val)
        _, frame_fit, signal_counts, signal_counts_err = signal_extraction.getFitFrames(matter_type, input_parquet_data, input_analysis_results,
                                                    input_parquet_mc, preselections=presel, print_info=False)
        histo_data.SetBinContent(i+1, signal_counts)
        histo_data.SetBinError(i+1, signal_counts_err)

        # mc
        sel = presel
        df_filtered = df.query(sel)
        signal_counts = df_filtered.shape[0]
        signal_counts_err = np.sqrt(signal_counts)
        histo_mc.SetBinContent(i+1, signal_counts)
        histo_mc.SetBinError(i+1, signal_counts_err)

    if normalise_to_first:
        histo_mc.Scale(1./histo_mc.GetBinContent(1))
        histo_data.Scale(1./histo_data.GetBinContent(1))
    else:
        histo_mc.Scale(1./histo_mc.GetBinContent(histo_mc.GetNbinsX()))
        histo_data.Scale(1./histo_data.GetBinContent(histo_data.GetNbinsX()))

    output_file.cd()
    histo_data.Write()
    histo_mc.Write()

    canvas.cd()
    histo_mc.Draw('PE')
    histo_mc.GetYaxis().SetRangeUser(0.,1.2)
    histo_data.Draw('PE SAME')
    histo_data.GetYaxis().SetRangeUser(0.,1.2)
    legend = ROOT.TLegend(0.4, 0.3, 0.6, 0.4, '', 'brNDC')
    legend.SetLineWidth(0)
    legend.AddEntry(histo_data, 'Data', 'L')
    legend.AddEntry(histo_mc, 'MC', 'L')
    legend.Draw()
    canvas.Write()

##########################
##        CosPA         ##
##########################

print('Checking cosPA')

cosPA_arr = np.linspace(0.99, 1., 50, dtype=np.float64)
nCosPa_bins = len(cosPA_arr) - 1

hDataSigCosPA = ROOT.TH1F('hDataSigCosPA', ';cos(#theta_{PA}); signal fraction', nCosPa_bins, cosPA_arr)
utils.setHistStyle(hDataSigCosPA, ROOT.kRed+1)

hMcSigCosPA = ROOT.TH1D('hMcSigCosPA', ';cos(#theta_{PA}); signal fraction', nCosPa_bins, cosPA_arr)
utils.setHistStyle(hMcSigCosPA, ROOT.kAzure+2)

cCosPAcomparison = ROOT.TCanvas('cCosPAcomparison', 'cCosPAcomparison', 800, 600)

sel_string = r'{} > {}'

systematic_routine('fCosPA', cosPA_arr[:-1], sel_string, hDataSigCosPA, hMcSigCosPA, cCosPAcomparison)

##########################
##      DcaV0Daug       ##
##########################

print('Checking DcaV0Daug')

DcaV0Daug_arr = np.linspace(0., .3, 30, dtype=np.float64)
nDcaV0Daug_bins = len(DcaV0Daug_arr) - 1

hDataSigDcaV0Daug = ROOT.TH1F('hDataSigDcaV0Daug', ';DCA_{V0 daughters} (cm); signal fraction', nDcaV0Daug_bins, DcaV0Daug_arr)
utils.setHistStyle(hDataSigDcaV0Daug, ROOT.kRed+1)

hMcSigDcaV0Daug = ROOT.TH1D('hMcSigDcaV0Daug', ';DCA_{V0 daughters} (cm); signal fraction', nDcaV0Daug_bins, DcaV0Daug_arr)
utils.setHistStyle(hMcSigDcaV0Daug, ROOT.kAzure+2)

cDcaV0Daugcomparison = ROOT.TCanvas('cDcaV0Daugcomparison', 'cDcaV0Daugcomparison', 800, 600)

sel_string = r'abs({}) < {}'

systematic_routine('fDcaV0Daug', DcaV0Daug_arr[1:], sel_string, hDataSigDcaV0Daug, hMcSigDcaV0Daug, cDcaV0Daugcomparison, normalise_to_first=False)

##########################
##         DcaHe        ##
##########################

print('Checking DcaHe')

DcaHe_arr = np.linspace(0., .3, 30, dtype=np.float64)
nDcaHe_bins = len(DcaHe_arr) - 1

hDataSigDcaHe = ROOT.TH1F('hDataSigDcaHe', ';DCA({}^{3}He) (cm); signal fraction', nDcaHe_bins, DcaHe_arr)
utils.setHistStyle(hDataSigDcaHe, ROOT.kRed+1)

hMcSigDcaHe = ROOT.TH1D('hMcSigDcaHe', ';DCA({}^{3}He) (cm); signal fraction', nDcaHe_bins, DcaHe_arr)
utils.setHistStyle(hMcSigDcaHe, ROOT.kAzure+2)

cDcaHecomparison = ROOT.TCanvas('cDcaHecomparison', 'cDcaHecomparison', 800, 600)

sel_string = r'abs({}) > {}'

systematic_routine('fDcaHe', DcaHe_arr[1:], sel_string, hDataSigDcaHe, hMcSigDcaHe, cDcaHecomparison)

##########################
##         DcaPi        ##
##########################

print('CPicking DcaPi')

DcaPi_arr = np.linspace(0., 3., 30, dtype=np.float64)
nDcaPi_bins = len(DcaPi_arr) - 1

hDataSigDcaPi = ROOT.TH1F('hDataSigDcaPi', ';DCA(#pi) (cm); signal fraction', nDcaPi_bins, DcaPi_arr)
utils.setHistStyle(hDataSigDcaPi, ROOT.kRed+1)

hMcSigDcaPi = ROOT.TH1D('hMcSigDcaPi', ';DCA(#pi) (cm); signal fraction', nDcaPi_bins, DcaPi_arr)
utils.setHistStyle(hMcSigDcaPi, ROOT.kAzure+2)

cDcaPicomparison = ROOT.TCanvas('cDcaPicomparison', 'cDcaPicomparison', 800, 600)

sel_string = r'abs({}) > {}'

systematic_routine('fDcaPi', DcaPi_arr[1:], sel_string, hDataSigDcaPi, hMcSigDcaPi, cDcaPicomparison)

