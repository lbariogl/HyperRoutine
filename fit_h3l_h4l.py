from hipe4ml.tree_handler import TreeHandler
import yaml
import argparse
import uproot
import numpy as np
import os
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
from signal_extraction import SignalExtraction

kOrangeC  = ROOT.TColor.GetColor('#ff7f00')


import sys
sys.path.append('utils')
import utils as utils


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
input_file_name_data = config['input_files_data']
input_file_name_mc_h3l = config['input_files_mc_h3l']
input_file_name_mc_h4l = config['input_files_mc_h4l']
output_dir = config['output_dir']
output_file_name = config['output_file']
colliding_system = config['colliding_system']

selections = config['selection']
is_matter = config['is_matter']
calibrate_he_momentum = config['calibrate_he_momentum']

selections_string = utils.convert_sel_to_string(selections)

if is_matter == 'matter':
    selections_string += ' and fIsMatter == True'
elif is_matter == 'antimatter':
    selections_string += ' and fIsMatter == False'


tree_names = ['O2datahypcands','O2hypcands', 'O2hypcandsflow']
tree_keys = uproot.open(input_file_name_data[0]).keys()
for tree in tree_names:
    for key in tree_keys:
        if tree in key:
            tree_name = tree
            break

data_hdl = TreeHandler(input_file_name_data, tree_name, folder_name='DF*')
mc_hdl_h3l_full = TreeHandler(input_file_name_mc_h3l, 'O2mchypcands', folder_name='DF*')
mc_hdl_h4l_full = TreeHandler(input_file_name_mc_h4l, 'O2mchypcands', folder_name='DF*')

# Add columns to the handlers
utils.correct_and_convert_df(data_hdl, calibrate_he3_pt=calibrate_he_momentum, isMC=False)
utils.correct_and_convert_df(mc_hdl_h3l_full, calibrate_he3_pt=calibrate_he_momentum, isMC=True)
utils.correct_and_convert_df(mc_hdl_h4l_full, calibrate_he3_pt=calibrate_he_momentum, isMC=True)

## reweight the pt spectrum of the MCs  
if colliding_system == 'pp':
    spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
    he3_spectrum = ROOT.TF1('mtexpo', '[2]*x*exp(-TMath::Sqrt([0]*[0]+x*x)/[1])', 0.1, 6)
    he3_spectrum.FixParameter(0, 2.99131)
    he3_spectrum.FixParameter(1, 0.5199)
    he3_spectrum.FixParameter(2, 1.0)
    he4_spectrum = he3_spectrum.Clone('he4_spectrum')
    he4_spectrum.FixParameter(0, 3.72738)
    # print('he3_spectrum parameters: ', he3_spectrum.GetParameter(0), he3_spectrum.GetParameter(1), he3_spectrum.GetParameter(2), he3_spectrum.GetParameter(3))

else:
    spectra_file = ROOT.TFile.Open('utils/bw.root')
    he3_spectrum = spectra_file.Get('helium3')
    he4_spectrum = spectra_file.Get('helium4')
    spectra_file.Close()

th1_pt_h3l = ROOT.TH1D('pt_h3l_mc', 'pt_h3l_mc', 100, 0, 10)
th2_n_sigma_he3_pt_mc = ROOT.TH2D('n_sigma_he3_pt_mc', 'n_sigma_he3_pt_mc', 100, 0, 10, 100, -5, 5)
th1_pt_h4l = ROOT.TH1D('pt_h4l_mc', 'pt_h4l_mc', 100, 0, 10)
th2_n_sigma_he4_pt_mc = ROOT.TH2D('n_sigma_he4_pt_mc', 'n_sigma_he4_pt_mc', 100, 0, 10, 100, -5, 5)

th2_wrong_mass_h3l_pt_mc = ROOT.TH2D('wrong_mass_h3l_pt_mc', 'wrong_mass_h3l_pt_mc', 100, 0, 10, 100, 2.96, 3.04)
th2_wrong_mass_h4l_pt_mc = ROOT.TH2D('wrong_mass_h4l_pt_mc', 'wrong_mass_h4l_pt_mc', 100, 0, 10, 100, 3.89, 3.97)

th2_n_sigma_he3_pt_data = ROOT.TH2D('n_sigma_he3_pt_data', 'n_sigma_he3_pt_data', 100, 0, 10, 100, -5, 5)
th2_n_sigma_he4_pt_data = ROOT.TH2D('n_sigma_he4_pt_data', 'n_sigma_he4_pt_data', 100, 0, 10, 100, -5, 5)

utils.reweight_pt_spectrum(mc_hdl_h3l_full, 'fAbsGenPt', he3_spectrum)
mc_hdl_h3l_full.apply_preselections('rej==True')
mc_hdl_h3l = mc_hdl_h3l_full.apply_preselections('fIsReco==True', inplace=False)

utils.reweight_pt_spectrum(mc_hdl_h4l_full, 'fAbsGenPt', he4_spectrum)
mc_hdl_h4l_full.apply_preselections('rej==True')
mc_hdl_h4l = mc_hdl_h4l_full.apply_preselections('fIsReco==True', inplace=False)



## apply the selections
data_hdl.apply_preselections(selections_string)
mc_hdl_h3l.apply_preselections(selections_string)
mc_hdl_h4l.apply_preselections(selections_string)

utils.fill_th1_hist(th1_pt_h3l, mc_hdl_h3l, 'fAbsGenPt')
utils.fill_th1_hist(th1_pt_h4l, mc_hdl_h4l, 'fAbsGenPt')
utils.fill_th2_hist(th2_n_sigma_he3_pt_mc, mc_hdl_h3l, 'fAbsGenPt', 'fNSigmaHe3')
utils.fill_th2_hist(th2_n_sigma_he4_pt_mc, mc_hdl_h4l, 'fAbsGenPt', 'fNSigmaHe4')
utils.fill_th2_hist(th2_n_sigma_he3_pt_data, data_hdl, 'fPt', 'fNSigmaHe3')
utils.fill_th2_hist(th2_n_sigma_he4_pt_data, data_hdl, 'fPt', 'fNSigmaHe4')
utils.fill_th2_hist(th2_wrong_mass_h3l_pt_mc, mc_hdl_h3l, 'fAbsGenPt', 'fMassH4L')
utils.fill_th2_hist(th2_wrong_mass_h4l_pt_mc, mc_hdl_h4l, 'fAbsGenPt', 'fMassH3L')

data_hdl.apply_preselections('fMassH4L >  3.89 and fMassH4L <  3.97 and fMassH3L > 2.96 and fMassH3L < 3.04')
mc_hdl_h3l.apply_preselections('fMassH4L > 3.89 and fMassH4L <  3.97 and fMassH3L > 2.96 and fMassH3L < 3.04')
mc_hdl_h4l.apply_preselections('fMassH4L > 3.89 and fMassH4L <  3.97 and fMassH3L > 2.96 and fMassH3L < 3.04')


## get the invariant mass template from MC, use a rookeyspdf
inv_mass_string_H4L = '#it{M}_{^{4}He+#pi^{-}}' if is_matter else '#it{M}_{^{4}#bar{He}+#pi^{+}}'
inv_mass_string_H3L = '#it{M}_{^{3}He+#pi^{-}}' if is_matter else '#it{M}_{^{3}#bar{He}+#pi^{+}}'


mass3HL = ROOT.RooRealVar('mass3HL', inv_mass_string_H3L, data_hdl['fMassH3L'].min(), data_hdl['fMassH3L'].max(), 'GeV/c^{2}')
mass4HL = ROOT.RooRealVar('mass4HL', inv_mass_string_H4L, data_hdl['fMassH4L'].min(), data_hdl['fMassH4L'].max(), 'GeV/c^{2}')


## first extract h3l and h4l templates from MC with a double sided crystal ball
mass_roo_mc_h3l = utils.ndarray2roo(np.array(mc_hdl_h3l['fMassH3L'].values, dtype=np.float64), mass3HL, 'histo_mc_h3l')
mass_roo_mc_h4l = utils.ndarray2roo(np.array(mc_hdl_h4l['fMassH4L'].values, dtype=np.float64), mass4HL, 'histo_mc_h4l')

mu3HL = ROOT.RooRealVar('mu_h3l', 'hypernucl mass', data_hdl['fMassH3L'].min(), data_hdl['fMassH3L'].max(), 'GeV/c^{2}')
sigma_h3l = ROOT.RooRealVar('sigma_h3l', 'hypernucl width', 0.001, 0.0024, 'GeV/c^{2}')
a1_h3l = ROOT.RooRealVar('a1_h3l', 'a1_h3l', 0., 5.)
a2_h3l = ROOT.RooRealVar('a2_h3l', 'a2_h3l', 0., 5.)
n1_h3l = ROOT.RooRealVar('n1_h3l', 'n1_h3l', 0., 5.)
n2_h3l = ROOT.RooRealVar('n2_h3l', 'n2_h3l', 0., 5.)
pars_h3l = [mu3HL, sigma_h3l, a1_h3l, a2_h3l, n1_h3l, n2_h3l]
signal_h3l = ROOT.RooCrystalBall('cb_h3l', 'cb_h3l_cl', mass3HL, mu3HL, sigma_h3l, a1_h3l, n1_h3l, a2_h3l, n2_h3l)
signal_h3l.fitTo(mass_roo_mc_h3l, ROOT.RooFit.Extended(True))
frame_h3l = mass3HL.frame()
frame_h3l.SetName('frame_h3l_mc')
mass_roo_mc_h3l.plotOn(frame_h3l)
signal_h3l.plotOn(frame_h3l)

mu4HL = ROOT.RooRealVar('mu_h4l', 'hypernucl mass', 3.91, 3.95, 'GeV/c^{2}')
sigma_h4l = ROOT.RooRealVar('sigma_h4l', 'hypernucl width', 0.001, 0.004, 'GeV/c^{2}')
a1_h4l = ROOT.RooRealVar('a1_h4l', 'a1_h4l', 0., 5.)
a2_h4l = ROOT.RooRealVar('a2_h4l', 'a2_h4l', 0., 5.)
n1_h4l = ROOT.RooRealVar('n1_h4l', 'n1_h4l', 0., 5.)
n2_h4l = ROOT.RooRealVar('n2_h4l', 'n2_h4l', 0., 5.)
pars_h4l = [mu4HL, sigma_h4l, a1_h4l, a2_h4l, n1_h4l, n2_h4l]
signal_h4l = ROOT.RooCrystalBall('cb_h4l', 'cb_h4l_cl', mass4HL, mu4HL, sigma_h4l, a1_h4l, n1_h4l, a2_h4l, n2_h4l)
signal_h4l.fitTo(mass_roo_mc_h4l, ROOT.RooFit.Extended(True))
frame_h4l = mass4HL.frame()
frame_h4l.SetName('frame_h4l_mc')
mass_roo_mc_h4l.plotOn(frame_h4l)
signal_h4l.plotOn(frame_h4l)

## fix all the params except for the mu
for par in pars_h3l:
    if par.GetName() == 'mu_h3l':
        continue
    if par.GetName() == 'sigma_h3l':
        print('Range H3L: ', par.getVal(), par.getVal()*1.5)
        par.setRange(par.getVal(), par.getVal()*1.5)
        continue
    par.setConstant(True)

for par in pars_h4l:
    if par.GetName() == 'mu_h4l':
        continue
    if par.GetName() == 'sigma_h4l':
        print('Range H4L: ', par.getVal(), par.getVal()*1.5)
        par.setRange(par.getVal(), par.getVal()*1.5)
        continue
    par.setConstant(True)

### now get the mc template of h4l for real h3l candidates and vice versa
mass3HL.setBins(60)
mass4HL.setBins(60)
mass_roo_mc_h3l_wrong_mass = utils.ndarray2roo(np.array(mc_hdl_h4l['fMassH3L'].values, dtype=np.float64), mass3HL, 'histo_mc_h3l_wrong_mass')
mass_roo_mc_h4l_wrong_mass = utils.ndarray2roo(np.array(mc_hdl_h3l['fMassH4L'].values, dtype=np.float64), mass4HL, 'histo_mc_h4l_wrong_mass')
roo_hist_h3l_wrong_mass = ROOT.RooDataHist('mc_datahist_h3l_wrong_mass', 'mc_datahist_h3l_wrong_mass', ROOT.RooArgList(mass3HL), mass_roo_mc_h3l_wrong_mass)
roo_hist_h4l_wrong_mass = ROOT.RooDataHist('mc_datahist_h4l_wrong_mass', 'mc_datahist_h4l_wrong_mass', ROOT.RooArgList(mass4HL), mass_roo_mc_h4l_wrong_mass)
pdf_h3l_wrong_mass = ROOT.RooHistPdf('mc_pdf_h3l_wrong_mass', 'mc_pdf_h3l_wrong_mass', ROOT.RooArgSet(mass3HL), roo_hist_h3l_wrong_mass)
pdf_h4l_wrong_mass = ROOT.RooHistPdf('mc_pdf_h4l_wrong_mass', 'mc_pdf_h4l_wrong_mass', ROOT.RooArgSet(mass4HL), roo_hist_h4l_wrong_mass)

### plot the histograms
frame_h3l_wrong_mass = mass3HL.frame()
frame_h3l_wrong_mass.SetName('frame_h3l_wrong_mass')
mass_roo_mc_h3l_wrong_mass.plotOn(frame_h3l_wrong_mass)

frame_h4l_wrong_mass = mass4HL.frame()
frame_h4l_wrong_mass.SetName('frame_h4l_wrong_mass')
mass_roo_mc_h4l_wrong_mass.plotOn(frame_h4l_wrong_mass)


## building models for h3l and h4l
c0_bkg_h3l = ROOT.RooRealVar('c0_bkg_h3l', 'c0_bkg_h3l', -1, 1.)
c1_bkg_h3l = ROOT.RooRealVar('c1_bkg_h3l', 'c1_bkg_h3l', -1, 1.)
bkg_h3l = ROOT.RooChebychev('bkg_h3l', 'bkg_h3l', mass3HL, ROOT.RooArgList(c0_bkg_h3l, c1_bkg_h3l))

c0_bkg_h4l = ROOT.RooRealVar('c0_bkg_h4l', 'c0_bkg_h4l', -1, 1.)
c1_bkg_h4l = ROOT.RooRealVar('c1_bkg_h4l', 'c1_bkg_h4l', -1, 1.)
bkg_h4l = ROOT.RooChebychev('bkg_h4l', 'bkg_h4l', mass4HL, ROOT.RooArgList(c0_bkg_h4l, c1_bkg_h4l))

nsig_h3l = ROOT.RooRealVar('nsig_h3l', 'signal events', 100, 0, 1e5)
n_bkg = ROOT.RooRealVar('nbkg', 'background events', 100, 0, 1e5)
n_sig_h4l = ROOT.RooRealVar('nsig_h4l', 'signal events', 100, 0, 1e5)

model_h3l = ROOT.RooAddPdf('model_h3l', 'model_h3l', ROOT.RooArgList(signal_h3l, pdf_h3l_wrong_mass, bkg_h3l), ROOT.RooArgList(nsig_h3l, n_sig_h4l, n_bkg))
model_h4l = ROOT.RooAddPdf('model_h4l', 'model_h4l', ROOT.RooArgList(signal_h4l, pdf_h4l_wrong_mass, bkg_h4l), ROOT.RooArgList(n_sig_h4l, nsig_h3l, n_bkg))
mass_roo_data_h3l = utils.ndarray2roo(np.array(data_hdl['fMassH3L'].values, dtype=np.float64), mass3HL, 'histo_data_h3l')
mass_roo_data_h4l = utils.ndarray2roo(np.array(data_hdl['fMassH4L'].values, dtype=np.float64), mass4HL, 'histo_data_h4l')


# simultaneous fit
categories = ROOT.RooCategory('categories', 'categories')
categories.defineType('h3l')
categories.defineType('h4l')
data = ROOT.RooDataSet('data', 'data', ROOT.RooArgSet(mass3HL, mass4HL, categories), ROOT.RooFit.Index(categories), ROOT.RooFit.Import('h3l', mass_roo_data_h3l), ROOT.RooFit.Import('h4l', mass_roo_data_h4l))
roosim = ROOT.RooSimultaneous('roosim', 'roosim', categories)
roosim.addPdf(model_h3l, 'h3l')
roosim.addPdf(model_h4l, 'h4l')
fit_results = roosim.fitTo(data, ROOT.RooFit.Extended(True), ROOT.RooFit.Save(True))

# model_h3l.fitTo(mass_roo_data_h3l, ROOT.RooFit.Extended(True))
# model_h4l.fitTo(mass_roo_data_h4l, ROOT.RooFit.Extended(True))


mass3HL.setRange('signal', mu3HL.getVal()-3*sigma_h3l.getVal(), mu3HL.getVal()+3*sigma_h3l.getVal())
mass4HL.setRange('signal', mu4HL.getVal()-3*sigma_h4l.getVal(), mu4HL.getVal()+3*sigma_h4l.getVal())
signal_h3l_int = signal_h3l.createIntegral(ROOT.RooArgSet(mass3HL), ROOT.RooArgSet(mass3HL), 'signal')
signal_h3l_int_val_3s = signal_h3l_int.getVal()*nsig_h3l.getVal()
signal_h3l_int_val_3s_error = signal_h3l_int_val_3s*nsig_h3l.getError()/nsig_h3l.getVal()
signal_h4l_int = signal_h4l.createIntegral(ROOT.RooArgSet(mass4HL), ROOT.RooArgSet(mass4HL), 'signal')
signal_h4l_int_val_3s = signal_h4l_int.getVal()*n_sig_h4l.getVal()
signal_h4l_int_val_3s_error = signal_h4l_int_val_3s*n_sig_h4l.getError()/n_sig_h4l.getVal()

##TPave infos
bkg_int_h3l = bkg_h3l.createIntegral(ROOT.RooArgSet(mass3HL), ROOT.RooArgSet(mass3HL))
bkg_int_h3l_val_3s = bkg_int_h3l.getVal()*n_bkg.getVal()
bkg_int_h3l_val_3s_error = bkg_int_h3l_val_3s*n_bkg.getError()/n_bkg.getVal()
wrong_mass_h3l_int = pdf_h3l_wrong_mass.createIntegral(ROOT.RooArgSet(mass3HL), ROOT.RooArgSet(mass3HL))
wrong_mass_h3l_int_val_3s = wrong_mass_h3l_int.getVal()*n_sig_h4l.getVal()
wrong_mass_h3l_int_val_3s_error = wrong_mass_h3l_int_val_3s*n_sig_h4l.getError()/n_sig_h4l.getVal()
s_b_ratio_h3l = signal_h3l_int_val_3s / (bkg_int_h3l_val_3s + wrong_mass_h3l_int_val_3s)
s_b_ratio_h3l_error = s_b_ratio_h3l*np.sqrt((signal_h3l_int_val_3s_error/signal_h3l_int_val_3s)**2 + (bkg_int_h3l_val_3s_error/bkg_int_h3l_val_3s)**2 + (wrong_mass_h3l_int_val_3s_error/wrong_mass_h3l_int_val_3s)**2)

bkg_int_h4l = bkg_h4l.createIntegral(ROOT.RooArgSet(mass4HL), ROOT.RooArgSet(mass4HL))
bkg_int_h4l_val_3s = bkg_int_h4l.getVal()*n_bkg.getVal()
bkg_int_h4l_val_3s_error = bkg_int_h4l_val_3s*n_bkg.getError()/n_bkg.getVal()
wrong_mass_h4l_int = pdf_h4l_wrong_mass.createIntegral(ROOT.RooArgSet(mass4HL), ROOT.RooArgSet(mass4HL))
wrong_mass_h4l_int_val_3s = wrong_mass_h4l_int.getVal()*nsig_h3l.getVal()
wrong_mass_h4l_int_val_3s_error = wrong_mass_h4l_int_val_3s*nsig_h3l.getError()/nsig_h3l.getVal()
s_b_ratio_h4l = signal_h4l_int_val_3s / (bkg_int_h4l_val_3s + wrong_mass_h4l_int_val_3s)
s_b_ratio_h4l_error = s_b_ratio_h4l*np.sqrt((signal_h4l_int_val_3s_error/signal_h4l_int_val_3s)**2 + (bkg_int_h4l_val_3s_error/bkg_int_h4l_val_3s)**2 + (wrong_mass_h4l_int_val_3s_error/wrong_mass_h4l_int_val_3s)**2)

pinfo_h3l = ROOT.TPaveText(0.632, 0.5, 0.932, 0.85, 'NDC')
pinfo_h3l.SetBorderSize(0)
pinfo_h3l.SetFillStyle(0)
pinfo_h3l.SetTextAlign(11)
pinfo_h3l.SetTextFont(42)
pinfo_h3l.AddText(f'Signal (S): {signal_h3l_int_val_3s:.0f} #pm {signal_h3l_int_val_3s_error:.0f}')
pinfo_h3l.AddText(f'S/B (3 #sigma): {s_b_ratio_h3l:.1f} #pm {s_b_ratio_h3l_error:.1f}')
pinfo_h3l.AddText('#mu = ' + f'{mu3HL.getVal()*1e3:.2f} #pm {mu3HL.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
pinfo_h3l.AddText('#sigma = ' + f'{sigma_h3l.getVal()*1e3:.2f} #pm {sigma_h3l.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')

pinfo_h4l = ROOT.TPaveText(0.632, 0.5, 0.932, 0.85, 'NDC')
pinfo_h4l.SetBorderSize(0)
pinfo_h4l.SetFillStyle(0)
pinfo_h4l.SetTextAlign(11)
pinfo_h4l.SetTextFont(42)
pinfo_h4l.AddText(f'Signal (S): {signal_h4l_int_val_3s:.0f} #pm {signal_h4l_int_val_3s_error:.0f}')
pinfo_h4l.AddText(f'S/B (3 #sigma): {s_b_ratio_h4l:.1f} #pm {s_b_ratio_h4l_error:.1f}')
pinfo_h4l.AddText('#mu = ' + f'{mu4HL.getVal()*1e3:.2f} #pm {mu4HL.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
pinfo_h4l.AddText('#sigma = ' + f'{sigma_h4l.getVal()*1e3:.2f} #pm {sigma_h4l.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')


## print the values of all the parameters
br_h3l = 0.25
br_h4l = 0.55
eff_3hl_mc = len(mc_hdl_h3l) / len(mc_hdl_h3l_full)
eff_4hl_mc = len(mc_hdl_h4l) / len(mc_hdl_h4l_full)

mass3HL.setBins(50)
mass4HL.setBins(30)
## plot the results
frame_data_h3l = mass3HL.frame()
frame_data_h3l.SetName('frame_data_h3l')
mass_roo_data_h3l.plotOn(frame_data_h3l)
model_h3l.plotOn(frame_data_h3l)
model_h3l.plotOn(frame_data_h3l, ROOT.RooFit.Components('cb_h3l'), ROOT.RooFit.LineColor(kOrangeC), ROOT.RooFit.LineStyle(2))
model_h3l.plotOn(frame_data_h3l, ROOT.RooFit.Components('mc_pdf_h3l_wrong_mass'), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(2))
model_h3l.plotOn(frame_data_h3l, ROOT.RooFit.Components('bkg_h3l'), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(2))
frame_data_h3l.addObject(pinfo_h3l)

frame_data_h4l = mass4HL.frame()
frame_data_h4l.SetName('frame_data_h4l')
mass_roo_data_h4l.plotOn(frame_data_h4l)
model_h4l.plotOn(frame_data_h4l)
model_h4l.plotOn(frame_data_h4l, ROOT.RooFit.Components('cb_h4l'), ROOT.RooFit.LineColor(kOrangeC), ROOT.RooFit.LineStyle(2))
model_h4l.plotOn(frame_data_h4l, ROOT.RooFit.Components('mc_pdf_h4l_wrong_mass'), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(2))
model_h4l.plotOn(frame_data_h4l, ROOT.RooFit.Components('bkg_h4l'), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(2))
frame_data_h4l.addObject(pinfo_h4l)



print ('----------------- Signal extraction -----------------')
print('efficiency h3l: ', len(mc_hdl_h3l) / len(mc_hdl_h3l_full))
print('efficiency h4l: ', len(mc_hdl_h4l) / len(mc_hdl_h4l_full))
print('nsig_h3l: ', nsig_h3l.getVal(), nsig_h3l.getError())
print('nsig_h4l: ', n_sig_h4l.getVal(), n_sig_h4l.getError())
print('h3l / h4l: ', nsig_h3l.getVal() / n_sig_h4l.getVal() / (eff_3hl_mc / eff_4hl_mc) / (br_h3l / br_h4l))

print ('----------------- Fit parameters -----------------')
print('nbkg: ', n_bkg.getVal(), n_bkg.getError())
print('c0_bkg_h3l: ', c0_bkg_h3l.getVal(), c0_bkg_h3l.getError())
print('c1_bkg_h3l: ', c1_bkg_h3l.getVal(), c1_bkg_h3l.getError())
print('c0_bkg_h4l: ', c0_bkg_h4l.getVal(), c0_bkg_h4l.getError())
print('c1_bkg_h4l: ', c1_bkg_h4l.getVal(), c1_bkg_h4l.getError())
print('mu_h3l: ', mu3HL.getVal(), mu3HL.getError())
print('sigma_h3l: ', sigma_h3l.getVal(), sigma_h3l.getError())
print('a1_h3l: ', a1_h3l.getVal(), a1_h3l.getError())
print('a2_h3l: ', a2_h3l.getVal(), a2_h3l.getError())
print('n1_h3l: ', n1_h3l.getVal(), n1_h3l.getError())
print('n2_h3l: ', n2_h3l.getVal(), n2_h3l.getError())
print('mu_h4l: ', mu4HL.getVal(), mu4HL.getError())
print('sigma_h4l: ', sigma_h4l.getVal(), sigma_h4l.getError())
print('a1_h4l: ', a1_h4l.getVal(), a1_h4l.getError())
print('a2_h4l: ', a2_h4l.getVal(), a2_h4l.getError())
print('n1_h4l: ', n1_h4l.getVal(), n1_h4l.getError())
print('n2_h4l: ', n2_h4l.getVal(), n2_h4l.getError())






tfile = ROOT.TFile.Open(f'{output_dir}/{output_file_name}', 'recreate')
tfile.cd()

he3_spectrum.Write()
he4_spectrum.Write()

th1_pt_h3l.Write()
th1_pt_h4l.Write()
th2_n_sigma_he3_pt_mc.Write()
th2_n_sigma_he4_pt_mc.Write()
th2_n_sigma_he3_pt_data.Write()
th2_n_sigma_he4_pt_data.Write()
th2_wrong_mass_h3l_pt_mc.Write()
th2_wrong_mass_h4l_pt_mc.Write()

frame_h3l.Write()
frame_h4l.Write()
frame_h3l_wrong_mass.Write()
frame_h4l_wrong_mass.Write()

frame_data_h3l.Write()
frame_data_h4l.Write()



tfile.Close()

exit()
