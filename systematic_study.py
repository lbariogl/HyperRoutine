import signal_extraction
import argparse
import yaml
import numpy as np
import ROOT

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

output_file = '../results/systematic_study.root'

# check CosPA
cosPA_arr = np.linspace(0.998, 1., 20)

cCosPA = ROOT.TCanvas('cCosPA', 'cCosPA', 800, 600)
cCosPA.Print('../results/SystematicChecksCosPA.pdf[')

# silent mode for fits
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

for cosPA in cosPA_arr[:-1]:
    presel = f'fCosPA > {cosPA:.4f}'
    label = r'cos(#theta_{PA}) > ' + f'{cosPA:.4f}'
    print(f'presel: {presel}')
    # perform fits
    _, frame_fit = signal_extraction.getFitFrames(matter_type, input_parquet_data, input_analysis_results,
                                                  input_parquet_mc, preselections=presel)
    cCosPA.cd()
    frame_fit.Draw()
    myLatex = ROOT.TLatex()
    myLatex.DrawLatexNDC(0.4, 0.95, label)
    cCosPA.Print('../results/SystematicChecksCosPA.pdf')


cCosPA.Print('../results/SystematicChecksCosPA.pdf]')
