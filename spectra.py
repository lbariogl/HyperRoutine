from signal_extraction import SignalExtraction
import ROOT
import numpy as np

import sys
sys.path.append('utils')
import utils as utils

class SpectraMaker:

    def __init__(self):

        # data related members
        self.data_hdl = None
        self.mc_hdl = None
        self.mc_reco_hdl = None

        self.n_ev = 0
        self.branching_ratio = 0.25
        self.delta_rap = 2.0

        # variable related members
        self.var = ''
        self.bins = []
        self.selection_string = ''  # could be either a string or a list of strings
        self.is_matter = ''

        self.n_bins_mass_data = 30
        self.n_bins_mass_mc = 80

        self.raw_counts = []
        self.raw_counts_err = []
        self.chi2 = []
        self.efficiency = []

        self.corrected_counts = []
        self.corrected_counts_err = []

        self.inv_mass_signal_func = "dscb"
        self.inv_mass_bkg_func = "pol1"  # could be either a string or a list of strings
        self.fit_func = None
        self.fit_options = None
        self.fit_range = []
        self.sigma_range_mc_to_data = [1., 1.5]

        self.output_dir = None

        self.h_signal_extractions_data = []
        self.h_signal_extractions_mc = []
        self.h_raw_counts = None
        self.h_efficiency = None
        self.h_corrected_counts = None

        # cuts for systematic uncertainties
        self.chi2_cut = 1.4
        self.relative_error_cut = 1.
        self.outlier_cut = 3

    def _check_members(self):

        var_options = ['fCt', 'fPt']
        if self.var not in var_options:
            raise ValueError(
                f'Invalid var option: {self.var}. Expected one of: {var_options}')

        if not self.data_hdl:
            raise ValueError(f'data_hdl not correctly set')

        if not self.mc_hdl:
            raise ValueError(f'mc_hdl not correctly set')

        if not self.mc_reco_hdl:
            raise ValueError(f'mc_reco_hdl not correctly set')

    def make_spectra(self):

        self._check_members()

        for ibin in range(0, len(self.bins) - 1):
            bin = [self.bins[ibin], self.bins[ibin + 1]]
            bin_sel = f'{self.var} > {bin[0]} & {self.var} < {bin[1]}'

            if self.var == 'fCt':
                mc_bin_sel = f'fGenCt > {bin[0]} & fGenCt < {bin[1]}'
            else:
                mc_bin_sel = f'abs(fGenPt) > {bin[0]} & abs(fGenPt) < {bin[1]}'

            # count generated per ct bin
            bin_mc_hdl = self.mc_hdl.apply_preselections(mc_bin_sel, inplace=False)

            if isinstance(self.selection_string, list):
                bin_sel = f'{bin_sel} and {self.selection_string[ibin]}'
                mc_bin_sel = f'{mc_bin_sel} and {self.selection_string[ibin]}'
            else:
                bin_sel = f'{bin_sel} and {self.selection_string}'
                mc_bin_sel = f'{mc_bin_sel} and {self.selection_string}'

            # select reconstructed in data and mc
            bin_data_hdl = self.data_hdl.apply_preselections(bin_sel, inplace=False)
            bin_mc_reco_hdl = self.mc_reco_hdl.apply_preselections(mc_bin_sel, inplace=False)

            # compute efficiency
            eff = len(bin_mc_reco_hdl) / len(bin_mc_hdl)
            print(mc_bin_sel)
            print("bin low", bin[0], "bin high", bin[1], "efficiency", eff)
            self.efficiency.append(eff)

            signal_extraction = SignalExtraction(bin_data_hdl, bin_mc_hdl)

            bkg_mass_fit_func = None
            if isinstance(self.inv_mass_bkg_func, list):
                bkg_mass_fit_func = self.inv_mass_bkg_func[ibin]
            else:
                bkg_mass_fit_func = self.inv_mass_bkg_func

            sgn_mass_fit_func = None
            if isinstance(self.inv_mass_signal_func, list):
                sgn_mass_fit_func = self.inv_mass_signal_func[ibin]
            else:
                sgn_mass_fit_func = self.inv_mass_signal_func

            signal_extraction.bkg_fit_func = bkg_mass_fit_func
            signal_extraction.signal_fit_func = sgn_mass_fit_func
            signal_extraction.n_bins_data = self.n_bins_mass_data
            signal_extraction.n_bins_mc = self.n_bins_mass_mc
            signal_extraction.n_evts = self.n_ev
            signal_extraction.matter_type = self.is_matter
            signal_extraction.performance = False
            signal_extraction.is_3lh = True

            if isinstance(self.sigma_range_mc_to_data[0], list):
                signal_extraction.sigma_range_mc_to_data = self.sigma_range_mc_to_data[ibin]
            else:
                signal_extraction.sigma_range_mc_to_data = self.sigma_range_mc_to_data

            fit_stats = signal_extraction.process_fit()

            if self.var == 'fPt':
                bin_label = f'{bin[0]} #leq #it{{p}}_{{T}} < {bin[1]} GeV/#it{{c}}'
            else:
                bin_label = f'{bin[0]} #leq #it{{ct}} < {bin[1]} cm'

            signal_extraction.additional_pave_text = bin_label

            self.h_signal_extractions_data.append(
                signal_extraction.data_frame_fit)
            self.h_signal_extractions_mc.append(signal_extraction.mc_frame_fit)
            self.raw_counts.append(fit_stats['signal'][0])
            self.raw_counts_err.append(fit_stats['signal'][1])
            self.chi2.append(fit_stats['chi2'])

    def make_histos(self):

        self._check_members()

        if not self.raw_counts:
            raise RuntimeError(
                'raw_counts is empty. You must run make_spectra first.')

        if self.var == 'fCt':
            x_label = r'#it{ct} (cm)'
            y_raw_label = r'#it{N}_{raw}'
            y_eff_label = r'#epsilon #times acc.'
            y_corr_label = r'#frac{d#it{N}}{d(#it{ct})} (cm^{-1})'
        else:
            x_label = r'#it{p}_{T} (GeV/#it{c})'
            y_raw_label = r'#it{N}_{raw}'
            y_eff_label = r'#epsilon #times acc.'
            y_corr_label = r'#frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}'

        self.h_raw_counts = ROOT.TH1D('h_raw_counts', f';{x_label};{y_raw_label}', len(
            self.bins) - 1, np.array(self.bins, dtype=np.float64))
        self.h_efficiency = ROOT.TH1D('h_efficiency', f';{x_label};{y_eff_label}', len(
            self.bins) - 1, np.array(self.bins, dtype=np.float64))

        self.h_corrected_counts = ROOT.TH1D('h_corrected_counts', f';{x_label};{y_corr_label}', len(
            self.bins) - 1, np.array(self.bins, dtype=np.float64))
        self.h_corrected_counts.GetXaxis().SetTitleSize(0.05)
        self.h_corrected_counts.GetYaxis().SetTitleSize(0.05)

        for ibin in range(0, len(self.bins) - 1):
            bin_width = self.bins[ibin + 1] - self.bins[ibin]

            self.h_raw_counts.SetBinContent(
                ibin + 1, self.raw_counts[ibin]/bin_width)
            self.h_raw_counts.SetBinError(
                ibin + 1, self.raw_counts_err[ibin]/bin_width)
            self.h_efficiency.SetBinContent(ibin + 1, self.efficiency[ibin])

            local_corrected_counts = self.raw_counts[ibin] / \
                self.efficiency[ibin] / bin_width
            local_corrected_counts_err = self.raw_counts_err[ibin] / \
                self.efficiency[ibin] / bin_width

            if self.var == 'fPt':
                local_corrected_counts = local_corrected_counts / \
                    self.n_ev / self.branching_ratio / self.delta_rap
                local_corrected_counts_err = local_corrected_counts_err / \
                    self.n_ev / self.branching_ratio / self.delta_rap

            self.h_corrected_counts.SetBinContent(
                ibin + 1, local_corrected_counts)
            self.h_corrected_counts.SetBinError(
                ibin + 1, local_corrected_counts_err)

            self.corrected_counts.append(local_corrected_counts)
            self.corrected_counts_err.append(local_corrected_counts_err)

    def fit(self):

        if not self.h_corrected_counts:
            raise ValueError(
                'h_corrected_counts not set. Use make_histos first.')

        if not self.fit_func:
            raise ValueError('Fit function not set.')

        if self.fit_range:
            self.h_corrected_counts.Fit(
                self.fit_func, self.fit_options, '', self.fit_range[0], self.fit_range[1])
        else:
            self.h_corrected_counts.Fit(self.fit_func, 'R')

    def del_dyn_members(self):
        self.raw_counts = []
        self.raw_counts_err = []
        self.efficiency = []
        self.corrected_counts = []
        self.corrected_counts_err = []
        self.h_signal_extractions_data = []
        self.h_signal_extractions_mc = []

        self.h_raw_counts = None
        self.h_efficiency = None
        self.h_corrected_counts = None

    def dump_to_output_dir(self):
        self.output_dir.cd()
        self.h_raw_counts.Write()
        self.h_efficiency.Write()
        self.h_corrected_counts.Write()
        for ibin in range(0, len(self.bins) - 1):
            bin = [self.bins[ibin], self.bins[ibin + 1]]
            self.output_dir.mkdir(f'{self.var}_{bin[0]}_{bin[1]}')
            self.output_dir.cd(f'{self.var}_{bin[0]}_{bin[1]}')
            h_data = self.h_signal_extractions_data[ibin]
            h_mc = self.h_signal_extractions_mc[ibin]
            bin_string_data = f'{self.var}_{self.bins[ibin]}_{self.bins[ibin + 1]}_data'
            bin_string_mc = f'{self.var}_{self.bins[ibin]}_{self.bins[ibin + 1]}_mc'
            h_data.SetName(bin_string_data)
            h_mc.SetName(bin_string_mc)
            h_data.Write()
            h_mc.Write()

    def chi2_selection(self):
        for el in self.chi2:
            if el > self.chi2_cut:
                return False
        return True

    def relative_error_selection(self):
        relative_errors = [
            err/val for val, err in zip(self.corrected_counts, self.corrected_counts_err)]
        for el in relative_errors:
            if el > self.relative_error_cut:
                return False
        return True

    def outlier_selection(self, std_corrected_counts, std_corrected_counts_err):
        for i, counts in enumerate(self.corrected_counts):
            distance = abs(std_corrected_counts[i] - counts) / np.hypot(
                self.corrected_counts_err[i], std_corrected_counts_err[i])
            if distance > self.outlier_cut:
                return False
        return True
