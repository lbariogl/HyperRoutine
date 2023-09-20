import utils as utils
import ROOT
import uproot
from hipe4ml.tree_handler import TreeHandler
import numpy as np

import sys
sys.path.append('utils')


class GeneralHandler:

    # constructors
    def __init__(self):
        self.data_hdl = None
        self.mc_hdl = None
        self.data_sel = ''
        self.mc_sel = ''
        self.n_evts = 0

    def __init__(self, data_hdl, mc_hdl, input_analysis_results):
        self.data_hdl = data_hdl
        self.mc_hdl = mc_hdl
        self.data_sel = ''
        self.mc_sel = ''
        self.n_evts = 0

        self.set_n_evts(input_analysis_results)

    # utilities
    def correct_and_convert_df(self, file_name, histo_name='hShiftVsPtHe3'):
        corr_file = ROOT.TFile(file_name)
        correction_hist = corr_file.Get(histo_name)
        correction_hist.SetDirectory(0)
        utils.correct_and_convert_df(self.data_hdl, correction_hist)
        utils.correct_and_convert_df(self.mc_hdl, correction_hist, isMC=True)

    def reweight_pt_spectrum(self, file_name='utils/heliumSpectraMB.root', histo_name='fCombineHeliumSpecLevyFit_0-100'):

        spectra_file = ROOT.TFile.Open(file_name)
        he3_spectrum = spectra_file.Get(histo_name)
        spectra_file.Close()

        self.mc_hdl.eval_data_frame("fAbsGenPt = abs(fGenPt)")
        utils.reweight_pt_spectrum(self.mc_hdl, 'fAbsGenPt', he3_spectrum)
        self.mc_hdl.apply_preselections('rej == True', inplace=True)

    def apply_selections(self):
        self.data_hdl.apply_preselections(self.data_sel)
        self.data_hdl.apply_preselections(self.mc_sel)

    # setters
    def set_data_hdl(self, data_hdl):
        self.data_hdl = data_hdl

    def set_mc_hdl(self, mc_hdl):
        self.mc_hdl = mc_hdl

    def set_n_evts(self, input_analysis_results):
        an_vtx_z = uproot.open(input_analysis_results)[
            'hyper-reco-task']['hZvtx']
        n_evts = an_vtx_z.values().sum() / 1e9
        n_evts = round(n_evts, 0)

    def set_data_sel(self, sel):
        self.data_sel = sel

    def set_mc_sel(self, sel):
        self.mc_sel = sel

    # getters
    def get_n_evts(self):
        return self.n_evts

    def get_data_sel(self):
        return self.data_sel

    def get_mc_sel(self):
        return self.mc_sel

    def get_data_hdl(self):
        return self.data_hdl

    def get_mc_hdl(self):
        return self.mc_hdl
