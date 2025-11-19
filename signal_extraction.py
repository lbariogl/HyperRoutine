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
kOrangeC  = ROOT.TColor.GetColor('#ff7f00')
ROOT.gROOT.SetBatch()

## create signal extraction class
class SignalExtraction:
    def __init__(self, input_data_hdl, input_mc_hdl=None): ## could be either a pandas or a tree handler

        self.data_hdl = input_data_hdl
        self.mc_hdl = input_mc_hdl
        self.is_3lh = True
        self.n_evts = 1e9
        self.is_matter = False
        self.performance = False
        self.additional_pave_text = '' ## additional text to be added to the ALICE performance pave
        self.colliding_system = 'pp'
        self.energy = 13.6

        ## fit-related variables
        self.pdf = None
        self.roo_dataset = None
        self.n_bins_data = 40
        self.n_bins_mc = 80
        self.signal_fit_func = 'dscb'
        self.sigma_range_mc_to_data = [1, 1.5]
        self.bkg_fit_func = 'pol1'

        ### frames to be saved to file
        self.out_file = None ## could also be a TDirectory
        self.data_frame_fit_name = 'data_frame_fit'
        self.mc_frame_fit_name = 'mc_frame_fit'
        
        ## output objects, broken members for ROOT > 6.22 (TO BE FIXED)
        self.mc_frame_fit = None
        self.data_frame_fit = None
        self.local_pvalue_graph = None



    def process_fit(self, extended_likelihood=True, rooworkspace_path=None):

        if self.is_3lh:
            self.inv_mass_string = '#it{M}_{^{3}He+#pi^{-}}' if self.is_matter else '#it{M}_{^{3}#bar{He}+#pi^{+}}'
            decay_string = '{}^{3}_{#Lambda}H #rightarrow ^{3}He + #pi^{-}' if self.is_matter else '{}^{3}_{#bar{#Lambda}}#bar{H} #rightarrow ^{3}#bar{He} + #pi^{+}'
            tree_var_name = 'fMassH3L'
        else:
            self.inv_mass_string = '#it{M}_{^{4}He+#pi^{-}}' if self.is_matter else '#it{M}_{^{4}#bar{He}+#pi^{+}}'
            decay_string = '{}^{4}_{#Lambda}H #rightarrow ^{4}He+#pi^{-}' if self.is_matter else '{}^{4}_{#bar{#Lambda}}#bar{H} #rightarrow ^{4}#bar{He} + #pi^{+}'
            tree_var_name = 'fMassH4L'

        # define signal and bkg variables
        if self.is_3lh:
            mass = ROOT.RooRealVar('m', self.inv_mass_string, 2.96, 3.04, 'GeV/#it{c}^{2}')
            mu = ROOT.RooRealVar('mu', 'hypernucl mass', 2.988, 2.992, 'GeV/#it{c}^{2}')
        else:
            mass = ROOT.RooRealVar('m', self.inv_mass_string, 3.89, 3.97, 'GeV/#it{c}^{2}')
            mu = ROOT.RooRealVar('mu', 'hypernucl mass', 3.9, 3.95, 'GeV/#it{c}^{2}')

        sigma = ROOT.RooRealVar('sigma', 'hypernucl width', 0.0011, 0.004, 'GeV/#it{c}^{2}')
        a1 = ROOT.RooRealVar('a1', 'a1', 0., 10.)
        a2 = ROOT.RooRealVar('a2', 'a2', 0., 10.)
        n1 = ROOT.RooRealVar('n1', 'n1', 0., 10.)
        n2 = ROOT.RooRealVar('n2', 'n2', 0., 10.)
        c0 = ROOT.RooRealVar('c0', 'constant c0', -1., 1)
        c1 = ROOT.RooRealVar('c1', 'constant c1', -1., 1)

        if self.signal_fit_func == 'dscb':
            signal = ROOT.RooCrystalBall('cb', 'cb', mass, mu, sigma, a1, n1, True)
        elif self.signal_fit_func == 'gaus':
            signal = ROOT.RooGaussian('gaus', 'gaus', mass, mu, sigma)
        else:
            raise ValueError(f'Invalid signal fit function. Expected one of: dscb, gaus')

        # define background pdf
        if self.bkg_fit_func == 'pol1':
            background = ROOT.RooChebychev('bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0))
        elif self.bkg_fit_func == 'pol2':
            background = ROOT.RooChebychev('bkg', 'pol2 bkg', mass, ROOT.RooArgList(c0, c1))
        elif self.bkg_fit_func == 'expo':
            background = ROOT.RooExponential('bkg', 'expo bkg', mass, c0)
        else:
            raise ValueError(f'Invalid background fit function. Expected one of: pol1, pol2, expo')

        if extended_likelihood:
            n_signal = ROOT.RooRealVar('n_signal', 'n_signal', 0., 1e4)
            n_background = ROOT.RooRealVar('n_background', 'n_background', 0., 1e6)
        else:
            f = ROOT.RooRealVar('f', 'fraction of signal', 0., 0.4)

        # fix DSCB parameters to MC
        if self.mc_hdl != None:
            mass_roo_mc = utils.ndarray2roo(np.array(self.mc_hdl['fMassH3L'].values, dtype=np.float64), mass, 'histo_mc')
            fit_results_mc = signal.fitTo(mass_roo_mc, ROOT.RooFit.Range(2.97, 3.01), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1))
            a1.setConstant()
            a2.setConstant()
            n1.setConstant()
            n2.setConstant()
            sigma.setRange(sigma.getVal(), self.sigma_range_mc_to_data[1]*sigma.getVal())
            self.mc_frame_fit = mass.frame(self.n_bins_mc)
            self.mc_frame_fit.SetName(self.mc_frame_fit_name)
            mass_roo_mc.plotOn(self.mc_frame_fit, ROOT.RooFit.Name('mc'), ROOT.RooFit.DrawOption('p'))
            signal.plotOn(self.mc_frame_fit, ROOT.RooFit.Name('signal'), ROOT.RooFit.DrawOption('p'))
            fit_param = ROOT.TPaveText(0.6, 0.43, 0.9, 0.85, 'NDC')
            fit_param.SetBorderSize(0)
            fit_param.SetFillStyle(0)
            fit_param.SetTextAlign(12)
            fit_param.AddText(r'#mu = ' + f'{mu.getVal()*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
            fit_param.AddText(r'#sigma = ' + f'{sigma.getVal()*1e3:.2f} #pm {sigma.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
            fit_param.AddText(r'alpha_{L} = ' + f'{a1.getVal():.2f} #pm {a1.getError():.2f}')
            fit_param.AddText(r'alpha_{R} = ' + f'{a2.getVal():.2f} #pm {a2.getError():.2f}')
            fit_param.AddText(r'n_{L} = ' + f'{n1.getVal():.2f} #pm {n1.getError():.2f}')
            fit_param.AddText(r'n_{R} = ' + f'{n2.getVal():.2f} #pm {n2.getError():.2f}')
            self.mc_frame_fit.addObject(fit_param)
            chi2_mc = self.mc_frame_fit.chiSquare('signal', 'mc')
            ndf_mc = self.n_bins_mc - fit_results_mc.floatParsFinal().getSize()
            fit_param.AddText('#chi^{2} / NDF = ' + f'{chi2_mc:.3f} (NDF: {ndf_mc})')

        # define the fit function and perform the actual fit
        if extended_likelihood:
            self.pdf = ROOT.RooAddPdf('total_pdf', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(n_signal, n_background))
        else:
            self.pdf = ROOT.RooAddPdf('total_pdf', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(f))

        mass_array = np.array(self.data_hdl[tree_var_name].values, dtype=np.float64)
        self.roo_dataset = utils.ndarray2roo(mass_array, mass)
        fit_results_data = self.pdf.fitTo(self.roo_dataset, ROOT.RooFit.Extended(extended_likelihood), ROOT.RooFit.Save(True), ROOT.RooFit.PrintLevel(-1))

        ## get fit parameters
        fit_pars = self.pdf.getParameters(self.roo_dataset)
        sigma_val = fit_pars.find('sigma').getVal()
        sigma_val_error = fit_pars.find('sigma').getError()
        mu_val = fit_pars.find('mu').getVal()
        mu_val_error = fit_pars.find('mu').getError()

        if extended_likelihood:
            signal_counts = n_signal.getVal()
            signal_counts_error = n_signal.getError()
            background_counts = n_background.getVal()
            background_counts_error = n_background.getError()

        else:
            signal_counts = (1-f.getVal())*self.roo_dataset.sumEntries()
            signal_counts_error = (1-f.getVal()) * self.roo_dataset.sumEntries()*f.getError()/f.getVal()
            background_counts = f.getVal()*self.roo_dataset.sumEntries()
            background_counts_error = f.getVal() * self.roo_dataset.sumEntries()*f.getError()/f.getVal()

        self.data_frame_fit = mass.frame(self.n_bins_data)
        self.data_frame_fit.SetName(self.data_frame_fit_name)

        self.roo_dataset.plotOn(self.data_frame_fit, ROOT.RooFit.Name('data'), ROOT.RooFit.DrawOption('p'))
        self.pdf.plotOn(self.data_frame_fit, ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(kOrangeC), ROOT.RooFit.Name('bkg'))
        self.pdf.plotOn(self.data_frame_fit, ROOT.RooFit.LineColor(ROOT.kAzure + 2 ), ROOT.RooFit.Name('fit_func'))


        chi2_data = self.data_frame_fit.chiSquare('fit_func', 'data')
        ndf_data = self.n_bins_data - fit_results_data.floatParsFinal().getSize()

        self.data_frame_fit.GetYaxis().SetTitleSize(0.06)
        self.data_frame_fit.GetYaxis().SetTitleOffset(0.9)
        self.data_frame_fit.GetYaxis().SetMaxDigits(2)
        self.data_frame_fit.GetXaxis().SetTitleOffset(1.1)

        # signal within 3 sigma
        mass.setRange('signal', mu_val-3*sigma_val, mu_val+3*sigma_val)
        signal_int = signal.createIntegral(ROOT.RooArgSet(mass), ROOT.RooArgSet(mass), 'signal')
        signal_int_val_3s = signal_int.getVal()*signal_counts
        signal_int_val_3s_error = signal_int_val_3s*signal_counts_error/signal_counts
        # background within 3 sigma
        mass.setRange('bkg', mu_val-3*sigma_val, mu_val+3*sigma_val)
        bkg_int = background.createIntegral(ROOT.RooArgSet(mass), ROOT.RooArgSet(mass), 'bkg')
        bkg_int_val_3s = bkg_int.getVal()*background_counts
        bkg_int_val_3s_error = bkg_int_val_3s*background_counts_error/background_counts
        significance = signal_int_val_3s / np.sqrt(signal_int_val_3s + bkg_int_val_3s)
        significance_err = utils.significance_error(signal_int_val_3s, bkg_int_val_3s, signal_int_val_3s_error, bkg_int_val_3s_error)
        s_b_ratio_err = np.sqrt((signal_int_val_3s_error/signal_int_val_3s)**2 + (bkg_int_val_3s_error/bkg_int_val_3s)**2)*signal_int_val_3s/bkg_int_val_3s

        # add pave for stats
        pinfo_vals = ROOT.TPaveText(0.632, 0.5, 0.932, 0.85, 'NDC')
        pinfo_vals.SetBorderSize(0)
        pinfo_vals.SetFillStyle(0)
        pinfo_vals.SetTextAlign(11)
        pinfo_vals.SetTextFont(42)
        pinfo_vals.AddText(f'Signal (S): {signal_counts:.0f} #pm {signal_counts_error:.0f}')
        pinfo_vals.AddText(f'S/B (3 #sigma): {signal_int_val_3s/bkg_int_val_3s:.1f} #pm {s_b_ratio_err:.1f}')
        pinfo_vals.AddText('S/#sqrt{S+B} (3 #sigma): ' + f'{significance:.1f} #pm {significance_err:.1f}')
        pinfo_vals.AddText('#mu = ' + f'{mu_val*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
        pinfo_vals.AddText('#sigma = ' + f'{sigma_val*1e3:.2f} #pm {sigma.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
        pinfo_vals.AddText('#chi^{2} / NDF = ' + f'{chi2_data:.3f} (NDF: {ndf_data})')

        ## add pave for ALICE performance
        if self.performance:
            pinfo_alice = ROOT.TPaveText(0.6, 0.5, 0.93, 0.85, 'NDC')
        else:
            pinfo_alice = ROOT.TPaveText(0.14, 0.6, 0.42, 0.85, 'NDC')
        pinfo_alice.SetBorderSize(0)
        pinfo_alice.SetFillStyle(0)
        pinfo_alice.SetTextAlign(11)
        pinfo_alice.SetTextFont(42)
        pinfo_alice.AddText('ALICE Performance')
        
        sqrtsnn = "#sqrt{#it{s}}"
        if self.colliding_system != 'pp':
            sqrtsnn = "#sqrt{#it{s_{NN}}}"
        pinfo_alice.AddText(f'Run 3 {self.colliding_system}, {sqrtsnn} = {self.energy} TeV')
        exponent = np.floor(np.log10(self.n_evts))
        self.n_evts = self.n_evts / 10**(exponent)
        pinfo_alice.AddText('#it{N}_{evt} = ' + f'{self.n_evts:.1f} ' + '#times 10^{' + f'{exponent:.0f}' + '}') 
        pinfo_alice.AddText(decay_string)
        pinfo_alice.AddText('0 #leq #it{p}_{T} < 9 GeV/#it{c}, |y| < 1')

        if self.additional_pave_text != '':
            pinfo_alice.AddText(self.additional_pave_text)

        legend = ROOT.TLegend(0.6, 0.2, 0.9, 0.45)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        ## signal and background components
        legend.AddEntry(self.data_frame_fit.findObject('fit_func'), 'Signal + Background', 'l')
        legend.AddEntry(self.data_frame_fit.findObject('bkg'), 'Background', 'l')
        ## set proper color to background legend entry
        legend.GetListOfPrimitives().At(1).SetLineColor(kOrangeC)
        legend.GetListOfPrimitives().At(1).SetLineStyle(2)
        legend.GetListOfPrimitives().At(1).SetLineWidth(2)
        self.data_frame_fit.addObject(pinfo_vals)
        self.data_frame_fit.addObject(pinfo_alice)
        self.data_frame_fit.addObject(legend)

        fit_stats = {'signal': [signal_counts, signal_counts_error],
                     'significance': [significance, significance_err], 's_b_ratio': [signal_int_val_3s/bkg_int_val_3s, s_b_ratio_err], 'chi2': chi2_data}
        

        if rooworkspace_path != None:
            w = ROOT.RooWorkspace('w')
            sb_model = ROOT.RooStats.ModelConfig('sb_model', w)
            sb_model.SetPdf(self.pdf)
            sb_model.SetParametersOfInterest(ROOT.RooArgSet(n_signal))
            sb_model.SetObservables(ROOT.RooArgSet(mass))
            getattr(w, 'import')(sb_model)
            getattr(w, 'import')(self.roo_dataset)
            w.writeToFile(rooworkspace_path + '/rooworkspace.root', True)
        
        if self.out_file != None:
            self.out_file.cd()
            self.data_frame_fit.Write()
            if self.mc_frame_fit != None:
                self.mc_frame_fit.Write()

        return fit_stats



    def compute_significance_asymptotic_calc(self, rooworkspace_path, do_local_p0plot=False):
        print("-----------------------------------------------")
        print("Computing significance with asymptotic calculator")
        ## get saved workspace
        workspace_file = ROOT.TFile(rooworkspace_path + '/rooworkspace.root', 'READ')
        w = workspace_file.Get('w')
        roo_abs_data = w.data('data')
        sb_model = w.obj('sb_model')
        poi = sb_model.GetParametersOfInterest().first()
        sb_model.SetSnapshot(ROOT.RooArgSet(poi))
        ## create the b-only model
        b_model = sb_model.Clone()
        b_model.SetName('b_model')
        poi.setVal(0)
        b_model.SetSnapshot(poi)
        b_model.Print()
        # w.var('sigma').setConstant(True)
        w.var('mu').setConstant(True)

        asymp_calc = ROOT.RooStats.AsymptoticCalculator(roo_abs_data, sb_model, b_model)
        asymp_calc.SetPrintLevel(0)
        asymp_calc_result = asymp_calc.GetHypoTest()
        null_p_value = asymp_calc_result.NullPValue()
        null_p_value_err = asymp_calc_result.NullPValueError()
        significance = asymp_calc_result.Significance()
        significance_err = asymp_calc_result.SignificanceError()



        if do_local_p0plot:
            ### perform a scan in mass and compute the significance
            masses = []
            p0_values = []
            p0_values_expected = []
            mass_array = np.linspace(w.var('mu').getMin(), w.var('mu').getMax(), 100)
            for mass in mass_array:

                w.var('mu').setVal(mass)
                w.var('mu').setConstant(True)
                asymp_calc_scan = ROOT.RooStats.AsymptoticCalculator(roo_abs_data, sb_model, b_model)
                asymp_calc_scan.SetOneSidedDiscovery(True)
                asym_calc_result_scan = asymp_calc_scan.GetHypoTest()
                null_p_value_scan = asym_calc_result_scan.NullPValue()
                masses.append(mass)
                p0_values.append(null_p_value_scan)

                print(f"Mass: {mass} MeV/c^2, p0: {null_p_value_scan:.10f}")

            ## create a graph with the p0 values
            self.local_pvalue_graph = ROOT.TGraph(len(masses), np.array(masses), np.array(p0_values))
            self.local_pvalue_graph.SetName('p0_values')
            self.local_pvalue_graph.GetXaxis().SetTitle(self.inv_mass_string)
            self.local_pvalue_graph.GetYaxis().SetTitle('Local p-value')
            # log Y axis
            self.local_pvalue_graph.SetMarkerStyle(20)
            self.local_pvalue_graph.SetMarkerColor(ROOT.kAzure + 2 )
            self.local_pvalue_graph.SetMarkerSize(0)
            self.local_pvalue_graph.SetLineColor(ROOT.kAzure + 2 )
            self.local_pvalue_graph.SetLineWidth(2)
            
            if self.out_file != None:
                self.out_file.cd()
                self.local_pvalue_graph.Write()

        print("****************************************************")
        print(f'p0: {null_p_value:.3E} +/- {null_p_value_err:.3E}')
        print(f'significance: {significance:.5f} +/- {significance_err:.5f}')
        print("****************************************************")







if __name__ == '__main__':

    # set parameters
    parser = argparse.ArgumentParser(
        description='Configure the parameters of the script.')
    parser.add_argument('--config-file', dest='config_file', default='',
                        help='path to the YAML file with configuration.')
    parser.add_argument('--performance', action='store_true',
                        help="True for performance plot", default=False)
    args = parser.parse_args()

    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)

    input_parquet_data = config['input_parquet_data']
    input_analysis_results = config['input_analysis_results']

    input_parquet_mc = config['input_parquet_mc']
    output_dir = config['output_dir']
    output_file = config['output_file']
    matter_type = config['matter_type']
    compute_significance = config['compute_significance']

    performance = args.performance
    data_hdl = TreeHandler(input_parquet_data)
    mc_hdl = None
    if input_parquet_mc != '':
        mc_hdl =  TreeHandler(input_parquet_mc)

    an_vtx_z = uproot.open(input_analysis_results)['hyper-reco-task']['hZvtx']
    n_evts = an_vtx_z.values().sum()

    signal_extraction = SignalExtraction(data_hdl, mc_hdl)
    signal_extraction.n_bins = config['n_bins']
    signal_extraction.n_evts = n_evts
    signal_extraction.is_matter = not matter_type == 'antimatter'
    signal_extraction.performance = performance
    signal_extraction.is_3lh = not config['is_4lh']
    signal_extraction.bkg_fit_func = 'pol2'
    signal_extraction.sigma_range_mc_to_data = [1, 1.3]

    signal_extraction.colliding_system = config['colliding_system']
    signal_extraction.energy = config['energy']
    
    out_file = ROOT.TFile(f'{output_dir}/{output_file}', 'recreate')
    signal_extraction.out_file = out_file.mkdir('signal_extraction')


    signal_extraction.process_fit(extended_likelihood=True, rooworkspace_path="../results")
    if compute_significance:
        signal_extraction.compute_significance_asymptotic_calc(rooworkspace_path="../results", do_local_p0plot=True)

    if config['is_4lh']:
        state_label = '4lh'
    else:
        state_label = '3lh'

    cSignalExtraction = ROOT.TCanvas('cSignalExtraction', 'cSignalExtraction', 800, 600)
    signal_extraction.data_frame_fit.SetTitle('')
    signal_extraction.data_frame_fit.Draw()
    cSignalExtraction.SaveAs(f'{output_dir}/cSignalExtraction_{matter_type}_{state_label}.pdf')
