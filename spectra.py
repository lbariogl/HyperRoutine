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
        self.selections = ''
        self.is_matter = ''

        self.raw_counts = []
        self.raw_counts_err = []
        self.efficiency = []

        self.corrected_counts = []
        self.corrected_counts_err = []

        self.fit_func = None
        self.fit_options = None
        self.fit_range = []

        self.output_dir = None

        self.h_raw_counts = None
        self.h_efficiency = None
        self.h_corrected_counts = None

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

            # count generated per ct bin
            bin_mc_hdl = self.mc_hdl.apply_preselections(
                bin_sel, inplace=False)

            if type(self.selections) == str:
                bin_sel = f'{self.selections}'
            else:
                bin_sel = f'{bin_sel} & {self.selections[ibin]}'

            # select reconstructed in data and mc
            bin_data_hdl = self.data_hdl.apply_preselections(
                bin_sel, inplace=False)
            bin_mc_reco_hdl = self.mc_reco_hdl.apply_preselections(
                bin_sel, inplace=False)

            # compute efficiency
            eff = len(bin_mc_reco_hdl) / len(bin_mc_hdl)
            self.efficiency.append(eff)

            self.output_dir.mkdir(f'{self.var}_{bin[0]}_{bin[1]}')
            self.output_dir.cd(f'{self.var}_{bin[0]}_{bin[1]}')

            signal_extraction = SignalExtraction(bin_data_hdl, bin_mc_hdl)
            signal_extraction.n_bins = 30
            n_ev_plot = round(self.n_ev / 1e9, 0)
            signal_extraction.n_evts = n_ev_plot
            signal_extraction.matter_type = self.is_matter
            signal_extraction.performance = False
            signal_extraction.is_3lh = True
            fit_stats = signal_extraction.process_fit()

            if self.var == 'fPt':
                bin_label = f'{bin[0]} #leq #it{{p}}_{{T}} < {bin[1]} GeV/#it{{c}}'
            else:
                bin_label = f'{bin[0]} #leq #it{{ct}} < {bin[1]} cm'

            signal_extraction.additional_pave_text = bin_label

            signal_extraction.data_frame_fit.Write()
            signal_extraction.mc_frame_fit.Write()

            self.raw_counts.append(fit_stats['signal'][0])
            self.raw_counts_err.append(fit_stats['signal'][1])


    def make_histos(self):

        self._check_members()

        if not self.raw_counts:
            raise RuntimeError('raw_counts is empty. You must run make_spectra first.')

        if self.var == 'fCt':
            x_label = '#it{ct} (cm)'
            y_raw_label = '#it{N}_{raw}'
            y_eff_label = '#epsilon #times acc.'
            y_corr_label = '#frac{d#it{N}}{d(#it{ct})} (cm^{-1})'
        else:
            x_label = '#it{p}_{T} (GeV/#it{c})'
            y_raw_label = '#it{N}{_{raw}'
            y_eff_label = '#epsilon #times acc.'
            y_corr_label = '#frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}'

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

            self.h_raw_counts.SetBinContent(ibin + 1, self.raw_counts[ibin]/bin_width)
            self.h_raw_counts.SetBinError(ibin + 1, self.raw_counts_err[ibin]/bin_width)
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

        self.output_dir.cd()
        self.h_raw_counts.Write()
        self.h_efficiency.Write()

    def fit(self):

        if not self.h_corrected_counts:
            raise ValueError('h_corrected_counts not set. Use make_histos first.')

        if not self.fit_func:
            raise ValueError('Fit function not set.')

        if self.fit_range:
            self.h_corrected_counts.Fit(self.fit_func, self.fit_options, '', self.fit_range[0], self.fit_range[1])
        else:
            self.h_corrected_counts.Fit(self.fit_func, 'R')

        self.output_dir.cd()
        self.h_corrected_counts.Write()
        self.fit_func.Write()
