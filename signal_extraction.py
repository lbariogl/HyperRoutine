import ROOT
import uproot
from hipe4ml.tree_handler import TreeHandler
import numpy as np

import argparse
import yaml

import sys
sys.path.append('utils')
import utils as utils

utils.set_style()
kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor('#ff7f00')

ROOT.gROOT.LoadMacro('utils/RooCustomPdfs/RooDSCBShape.cxx++')

def configureMatterType3LH(matter_type, selections):

    if matter_type == 'matter':
        inv_mass_string = '#it{M}_{^{3}He+#pi^{-}}'
        if selections != '':
            extended_selections = selections + ' and fIsMatter==True'
        else:
            extended_selections = 'fIsMatter==True'
    elif matter_type == 'antimatter':
        inv_mass_string = '#it{M}_{^{3}#bar{He}+#pi^{+}}'
        if selections != '':
            extended_selections = selections + ' and fIsMatter==False'
        else:
            extended_selections = 'fIsMatter==False'
    else:
        '#it{M}_{^{3}He+#pi^{-}} + c.c.'
        extended_selections = selections

    return extended_selections, inv_mass_string

def getNevents(input_analysis_results, print_info=True):
    # get number of events
    an_vtx_z = uproot.open(input_analysis_results)['hyper-reco-task']['hZvtx']
    n_evts = an_vtx_z.values().sum()
    if print_info:
        print(f'Number of events: {n_evts}')
    n_evts = round(n_evts/1e9, 0)
    return n_evts

def create_handlers(input_parquet_data, input_parquet_mc):
    data_hdl = TreeHandler(input_parquet_data)
    mc_hdl =  TreeHandler(input_parquet_mc)
    return data_hdl, mc_hdl


def fit3LH(data_hdl, mc_hdl, matter_type, input_analysis_results, selections='', n_bins = 40, print_info = True):

    # adapt labels and selections to matter type
    extended_selections, inv_mass_string = configureMatterType3LH(matter_type, selections)

    # get number of events
    n_evts = getNevents(input_analysis_results, print_info)

    # define signal pdf
    mass = ROOT.RooRealVar('m', inv_mass_string, 2.96, 3.04, 'GeV/c^{2}')
    mu = ROOT.RooRealVar('mu', 'hypernucl mass', 2.98, 3.0, 'GeV/c^{2}')
    sigma = ROOT.RooRealVar('sigma', 'hypernucl width',
                            0.001, 0.004, 'GeV/c^{2}')
    a1 = ROOT.RooRealVar('a1', 'a1', 0, 5.)
    a2 = ROOT.RooRealVar('a2', 'a2', 0, 10.)
    n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
    n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
    signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)

    # define background pdf
    c0 = ROOT.RooRealVar('c0', 'constant c0', -1., 1)
    c1 = ROOT.RooRealVar('c1', 'constant c1', -1., 1)
    background = ROOT.RooChebychev(
        'bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0, c1))
    f = ROOT.RooRealVar('f', 'fraction of signal', 0.01, 0.4)

    # fix DSCB parameters to MC
    mass_roo_mc = utils.ndarray2roo(
        np.array(mc_hdl['fMassH3L'].values, dtype=np.float64), mass, 'histo_mc')
    signal.fitTo(mass_roo_mc, ROOT.RooFit.Range(2.97, 3.01))
    a1.setConstant()
    a2.setConstant()
    n1.setConstant()
    n2.setConstant()
    frame_prefit = mass.frame(80)
    mass_roo_mc.plotOn(frame_prefit)
    signal.plotOn(frame_prefit)
    fit_param = ROOT.TPaveText(0.6, 0.6, 0.9, 0.9, 'NDC')
    fit_param.SetBorderSize(0)
    fit_param.SetFillStyle(0)
    fit_param.SetTextAlign(12)
    fit_param.AddText(
        '#mu = ' + f'{mu.getVal()*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
    fit_param.AddText(
        '#sigma = ' + f'{sigma.getVal()*1e3:.2f} #pm {sigma.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
    frame_prefit.addObject(fit_param)

    # define the fit function and perform the actual fit
    fit_function = ROOT.RooAddPdf(
        'total_pdf', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(f))

    if extended_selections != '':
        data_hdl.apply_preselections(extended_selections)

    mass_array = np.array(data_hdl['fMassH3L'].values, dtype=np.float64)
    mass_roo_data = utils.ndarray2roo(mass_array, mass)
    frame_fit, signal_counts, signal_counts_err = utils.fit_and_plot(mass_roo_data, mass, fit_function, signal,
                                       background, sigma, mu, f, n_ev=n_evts, matter_type=matter_type, print_info=print_info, n_bins=n_bins)

    return frame_prefit, frame_fit, signal_counts, signal_counts_err

if __name__ == '__main__':

    # set parameters
    parser = argparse.ArgumentParser(
        description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='config_file', default='',
                        help='path to the YAML file with configuration.')
    parser.add_argument('--nbins', dest='n_bins', default=30,
                        help='number of bins in the final plot.')
    args = parser.parse_args()

    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)

    matter_type = config['matter_type']
    input_parquet_data = config['input_parquet_data']
    input_analysis_results = config['input_analysis_results']

    input_parquet_mc = config['input_parquet_mc']
    output_dir = config['output_dir']
    output_file = config['output_file']

    ml_efficiency_scan = config['ml_efficiency_scan']
    selections = config['preselections']
    input_eff_dir = config['input_eff_dir']
    n_bins = config['n_bins']

    data_hdl, mc_hdl = create_handlers(input_parquet_data, input_parquet_mc)

    # perform fits
    frame_prefit, frame_fit, signal_counts, signal_counts_err = fit3LH(data_hdl, mc_hdl, matter_type, input_analysis_results, selections, n_bins)

    # create output file and save frames
    out_file = ROOT.TFile(f'{output_dir}/{output_file}', 'recreate')
    out_file.cd()
    frame_prefit.Write('histo_mc')
    frame_fit.Write('fit')

    cSignalExtraction = ROOT.TCanvas('cSignalExtraction', 'cSignalExtraction', 800, 600)
    frame_fit.Draw()
    cSignalExtraction.SaveAs(f'{output_dir}/cSignalExtraction_{matter_type}.pdf')
