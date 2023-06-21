import ROOT
import yaml
from hipe4ml.tree_handler import TreeHandler

import sys
sys.path.append('utils')
import utils as utils

config_file_mc = open('configs/data_vs_mc/config_mc.yaml', 'r')
config_mc = yaml.full_load(config_file_mc)
input_file_name_mc = config_mc['input_files']

config_file_data = open('configs/data_vs_mc/config_antimatter.yaml', 'r')
config_data = yaml.full_load(config_file_data)
input_file_name_data = config_data['input_files']

# MC
tree_hdl_mc = TreeHandler(input_file_name_mc, 'O2mchypcands')
df_mc = tree_hdl_mc.get_data_frame()
df_mc.eval('fP = fPt * cosh(fEta)', inplace=True)
df_mc.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)
df_mc.eval('fDecLen = sqrt(fXDecVtx**2 + fYDecVtx**2 + fZDecVtx**2)', inplace=True)
df_mc_filtered = df_mc.query(
    'fPt>1 and fCosPA>0.998 and fNSigmaHe>-2. and fTPCsignalPi<1000 and abs(fMassH3L - 2.992) < 0.006')

# Data
tree_hdl_data = TreeHandler(input_file_name_data, 'O2datahypcands')
df_data = tree_hdl_data.get_data_frame()
df_data.eval('fP = fPt * cosh(fEta)', inplace=True)
df_data.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)
df_data.eval(
    'fDecLen = sqrt(fXDecVtx**2 + fYDecVtx**2 + fZDecVtx**2)', inplace=True)
df_data_filtered = df_data.query(
    'fPt>1 and fCosPA>0.998 and fNSigmaHe>-2. and fTPCsignalPi<1000 and abs(fMassH3L - 2.992) < 0.006')

# Data histograms
hDataMassH3L = ROOT.TH1F(
    'hDataMassH3L', '; m({}^{3}_{#Lambda}H) (GeV/#it{c})', 40, 2.96, 3.04)
hDataCosPA = ROOT.TH1F('hDataCosPA', ';cos(#theta_{PA})', 200, 0.998, 1)
hDataDCAv0Daugh = ROOT.TH1F('hDataDCAv0Daugh', ';DCA (cm)', 200, 0., 1.2)
hDataDCAhe = ROOT.TH1F('hDataDCAhe', ';DCA (cm)', 200, 0., 1.2)
hDataDCApi = ROOT.TH1F('hDataDCApi', ';DCA (cm)', 200, 0., 1.2)
hDataPtRec = ROOT.TH1F('hDataPtRec', ';#it{p}_{T} (GeV/#it{c})', 50, 0, 10)
hDataRadius = ROOT.TH1F('hDataRadius', ';Radius (cm)', 100, 0, 40)
hDataDecLen = ROOT.TH1F('hDataDecLen', ';Decay length (cm)', 100, 0, 40)

hist_data_dict = {'fMassH3L': hDataMassH3L,
                  'fCosPA': hDataCosPA,
                  'fDcaV0Daug': hDataDCAv0Daugh,
                  'fDcaHe': hDataDCAhe,
                  'fDcaPi': hDataDCApi,
                  'fPt': hDataPtRec,
                  'fDecRad': hDataRadius,
                  'fDecLen': hDataDecLen}

for var, hist in hist_data_dict.items():
    utils.fill_th1_hist(hist, df_data_filtered, var)
    hist.Scale(1./hist.Integral())
    utils.setHistStyle(hist, ROOT.kAzure+2, linewidth=2)

# MC histograms
hMcMassH3L = ROOT.TH1F(
    'hMcMassH3L', '; m({}^{3}_{#Lambda}H) (GeV/#it{c})', 40, 2.96, 3.04)
hMcCosPA = ROOT.TH1F('hMcCosPA', ';cos(#theta_{PA})', 200, 0.998, 1)
hMcDCAv0Daugh = ROOT.TH1F('hMcDCAv0Daugh', ';DCA (cm)', 200, 0., 1.2)
hMcDCAhe = ROOT.TH1F('hMcDCAhe', ';DCA (cm)', 200, 0., 1.2)
hMcDCApi = ROOT.TH1F('hMcDCApi', ';DCA (cm)', 200, 0., 1.2)
hMcPtRec = ROOT.TH1F('hMcPtRec', ';#it{p}_{T} (GeV/#it{c})', 50, 0, 10)
hMcRadius = ROOT.TH1F('hMcRadius', ';Radius (cm)', 100, 0, 40)
hMcDecLen = ROOT.TH1F('hMcDecLen', ';Decay length (cm)', 100, 0, 40)

hist_mc_dict = {'fMassH3L': hMcMassH3L,
                'fCosPA': hMcCosPA,
                'fDcaV0Daug': hMcDCAv0Daugh,
                'fDcaHe': hMcDCAhe,
                'fDcaPi': hMcDCApi,
                'fPt': hMcPtRec,
                'fDecRad': hMcRadius,
                'fDecLen': hMcDecLen}

for var, hist in hist_mc_dict.items():
    utils.fill_th1_hist(hist, df_mc_filtered, var)
    hist.Scale(1./hist.Integral())
    utils.setHistStyle(hist, ROOT.kRed+1, linewidth=2)

# Create canvases
cMassH3L = ROOT.TCanvas('cMassH3L', 'cMassH3L', 800, 600)
cCosPA = ROOT.TCanvas('cCosPA', 'cCosPA', 800, 600)
cDCAv0Daugh = ROOT.TCanvas('cDCAv0Daugh', 'cDCAv0Daugh', 800, 600)
cDCAhe = ROOT.TCanvas('cDCAhe', 'cDCAhe', 800, 600)
cDCApi = ROOT.TCanvas('cDCApi', 'cDCAhe', 800, 600)
cPtRec = ROOT.TCanvas('cPtRec', 'cPtRec', 800, 600)
cRadius = ROOT.TCanvas('cRadius', 'cRadius', 800, 600)
cDecLen = ROOT.TCanvas('cDecLen', 'cDecLen', 800, 600)

canvas_dict = {'fMassH3L': cMassH3L,
               'fCosPA': cCosPA,
               'fDcaV0Daug': cDCAv0Daugh,
               'fDcaHe': cDCAhe,
               'fDcaPi': cDCApi,
               'fPt': cPtRec,
               'fDecRad': cRadius,
               'fDecLen': cDecLen}

y_limit_dict = {'fMassH3L': 0.5,
                'fCosPA': 0.4,
                'fDcaV0Daug': 0.22,
                'fDcaHe': 0.06,
                'fDcaPi': 0.02,
                'fPt': 0.5,
                'fDecRad': 0.07,
                'fDecLen': 0.05}

output_file = ROOT.TFile('../results/data_vs_mc.root', 'recreate')

for var, canvas in canvas_dict.items():
    canvas.cd()
    # maximum_data = hist_data_dict[var].GetMaximum()
    # maximum_mc = hist_data_dict[var].GetMaximum()
    # maximum = max(maximum_data, maximum_mc)
    left_edge = hist_data_dict[var].GetXaxis().GetBinLowEdge(1)
    right_edge = hist_data_dict[var].GetXaxis().GetBinUpEdge(
        hist_data_dict[var].GetNbinsX())
    x_title = hist_data_dict[var].GetXaxis().GetTitle()
    canvas.DrawFrame(left_edge, 0., right_edge, y_limit_dict[var], f';{x_title};')
    hist_data_dict[var].Draw('HISTO SAME')
    hist_mc_dict[var].Draw('HISTO SAME')

    legend = ROOT.TLegend(0.4, 0.7, 0.6, 0.8, '', 'brNDC')
    legend.SetLineWidth(0)
    legend.AddEntry(hist_data_dict[var], 'Data', 'L')
    legend.AddEntry(hist_mc_dict[var], 'MC', 'L')
    legend.Draw()
    output_file.cd()
    canvas.Write()
