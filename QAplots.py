from hipe4ml.tree_handler import TreeHandler
import pandas as pd
import ROOT
import uproot
import argparse
import os
import yaml

import sys
sys.path.append('utils')
import utils as utils

ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser(
        description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file',
                    help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
input_file_name_data = config['input_files_data']
input_file_name_mc = config['input_files_mc']
input_analysis_results_file = config['input_analysis_results_file']
input_analysis_results_file_mc = config['input_analysis_results_file_mc']
output_dir_name = config['output_dir']
qa_dir_name = output_dir_name + 'QAplots'
output_file_name_qa = config['output_file_qa']
pt_bins = config['pt_bins']
selections_std = config['selection']
is_matter = config['is_matter']

matter_options = ['matter', 'antimatter', 'both']
if is_matter not in matter_options:
    raise ValueError(
        f'Invalid is-matter option. Expected one of: {matter_options}')

args = parser.parse_args()

if not os.path.exists(qa_dir_name):
   os.makedirs(qa_dir_name)

print('**********************************')
print('    Running QAplots.py')
print('**********************************\n')
print("----------------------------------")
print("** Loading data and apply preselections **")


tree_names = ['O2datahypcands','O2hypcands', 'O2hypcandsflow']
tree_keys = uproot.open(input_file_name_data[0]).keys()
for tree in tree_names:
    for key in tree_keys:
        if tree in key:
            tree_name = tree
            break
print(f"Data tree found: {tree_name}")
data_hdl = TreeHandler(input_file_name_data, tree_name, folder_name='DF*')
mc_hdl = TreeHandler(input_file_name_mc, 'O2mchypcands', folder_name='DF*')

# declare output file
output_file = ROOT.TFile.Open(
    f'{output_dir_name}/{output_file_name_qa}.root', 'recreate')

# Add columns to the handlers
utils.correct_and_convert_df(data_hdl, calibrate_he3_pt=False)
utils.correct_and_convert_df(mc_hdl, calibrate_he3_pt=False, isMC=True)

# apply preselections
# selections = 'fCosPA > 0.99 and fNSigmaHe > -2 and fTPCsignalPi < 1000'
# data_hdl.apply_preselections(selections)
# data_hdl.apply_preselections(selections)

matter_sel = ''
mc_matter_sel = ''
if is_matter == 'matter':
    matter_sel = 'fIsMatter == True'
    mc_matter_sel = 'fGenPt > 0'

elif is_matter == 'antimatter':
    matter_sel = 'fIsMatter == False'
    mc_matter_sel = 'fGenPt < 0'

if matter_sel != '':
    data_hdl.apply_preselections(matter_sel)
    mc_hdl.apply_preselections(mc_matter_sel)

# reweight MC pT spectrum
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()
utils.reweight_pt_spectrum(mc_hdl, 'fAbsGenPt', he3_spectrum)

mc_hdl.apply_preselections('rej==True')
# Needed to remove the peak at 28.5 cm in the anchored MC
mc_hdl.apply_preselections('fGenCt < 28.5 or fGenCt > 28.6 and fIsReco == 1')
mc_reco_hdl = mc_hdl.apply_preselections('fIsReco == 1', inplace=False)

print("----------------------------------")
print("** Creating and filling histograms **")

############# Create histograms #############

# Data
hCosPA = ROOT.TH1F("hCosPA", r";cos(#theta_{PA})", 50, 0.95, 1)
utils.setHistStyle(hCosPA, ROOT.kRed+1)
hNTPCclus = ROOT.TH1F("hNTPCclus", r";n TPC clusters", 50, 60, 200)
utils.setHistStyle(hNTPCclus, ROOT.kRed+1)
hMass3LH = ROOT.TH1F("hMass3LH", r"; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 40, 2.96, 3.04)
utils.setHistStyle(hMass3LH, ROOT.kRed+1)
hMass4LH = ROOT.TH1F("hMass4LH", r";  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})", 30, 3.89, 3.97)
utils.setHistStyle(hMass4LH, ROOT.kRed+1)
hPtRec = ROOT.TH1F("hPtRec", r";#it{p}_{T} (GeV/#it{c})", 50, 0, 5)
utils.setHistStyle(hPtRec, ROOT.kRed+1)
hCtRec = ROOT.TH1F("hCtRec", r";#it{c#tau} (cm)", 50, 0, 40)
utils.setHistStyle(hCtRec, ROOT.kRed+1)
hRadius = ROOT.TH1F("hRadius", r";Radius (cm)", 100, 0, 40)
utils.setHistStyle(hRadius, ROOT.kRed+1)
hDecLen = ROOT.TH1F("hDecLen", r";Decay length (cm)", 100, 0, 40)
utils.setHistStyle(hDecLen, ROOT.kRed+1)
hNSigHe = ROOT.TH1F("hNSigmaHe", r";n_{#sigma}^{TPC}({}^{3}He)", 50, -3, 3)
utils.setHistStyle(hNSigHe, ROOT.kRed+1)
h2MassCosPA = ROOT.TH2F("h2MassCosPA", r";cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0.99, 1, 50, 2.96, 3.04)
h2MassDecLen = ROOT.TH2F("h2MassDecLen", r";Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 40, 50, 2.96, 3.04)
h2MassDCADaughters = ROOT.TH2F("h2MassDCADaughters", r";DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 200, 0, 0.3, 50, 2.96, 3.04)
h2MassDCAHePv = ROOT.TH2F("h2MassDCAHe", r";DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 2, 50, 2.96, 3.04)
h2MassPt = ROOT.TH2F("h2MassPt", r";#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 0, 7, 50, 2.96, 3.04)
h2Mass4LHnSigmaHe = ROOT.TH2F("h2Mass4LHnSigmaHe", r";n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})", 50, -4, 4, 30, 3.89, 3.97)
# MC
hCosPA_MC = ROOT.TH1F("hCosPA_MC", r";cos(#theta_{PA})", 50, 0.95, 1)
utils.setHistStyle(hCosPA_MC, ROOT.kAzure+2)
hNTPCclus_MC = ROOT.TH1F("hNTPCclus_MC", r";n TPC clusters", 50, 60, 200)
utils.setHistStyle(hNTPCclus_MC, ROOT.kAzure+2)
hMass3LH_MC = ROOT.TH1F("hMass3LH_MC", r"; m({}^{3}_{#Lambda}H) (GeV/#it{c})", 40, 2.96, 3.04)
utils.setHistStyle(hMass3LH_MC, ROOT.kAzure+2)
hMass4LH_MC = ROOT.TH1F("hMass4LH_MC", r";  m({}^{4}_{#Lambda}H) (GeV/#it{c^{2}})", 30, 3.89, 3.97)
utils.setHistStyle(hMass4LH_MC, ROOT.kAzure+2)
hPtRec_MC = ROOT.TH1F("hPtRec_MC", r";#it{p}_{T} (GeV/#it{c})", 50, 0, 5)
utils.setHistStyle(hPtRec_MC, ROOT.kAzure+2)
hCtRec_MC = ROOT.TH1F("hCtRec_MC", r";#it{c#tau} (cm)", 50, 0, 40)
utils.setHistStyle(hCtRec_MC, ROOT.kAzure+2)
hRadius_MC = ROOT.TH1F("hRadiu_MCs", r";Radius (cm)", 100, 0, 40)
utils.setHistStyle(hRadius_MC, ROOT.kAzure+2)
hDecLen_MC = ROOT.TH1F("hDecLen_MC", r";Decay length (cm)", 100, 0, 40)
utils.setHistStyle(hDecLen_MC, ROOT.kAzure+2)
hNSigHe_MC = ROOT.TH1F("hNSigmaHe_MC", r";n_{#sigma}^{TPC}({}^{3}He)", 50, -3, 3)
utils.setHistStyle(hNSigHe_MC, ROOT.kAzure+2)
h2MassCosPA_MC = ROOT.TH2F("h2MassCosPA_MC", r";cos(#theta_{PA}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0.99, 1, 50, 2.96, 3.04)
h2MassDecLen_MC = ROOT.TH2F("h2MassDecLen_MC", r";Decay length (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 40, 50, 2.96, 3.04)
h2MassDCADaughters_MC = ROOT.TH2F("h2MassDCADaughters_MC", r";DCA daughters (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 200, 0, 0.3, 50, 2.96, 3.04)
h2MassDCAHePv_MC = ROOT.TH2F("h2MassDCAHe_MC", r";DCA He3 PVs (cm); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 100, 0, 2, 50, 2.96, 3.04)
h2MassPt_MC = ROOT.TH2F("h2MassPt_MC", r";#it{p}_{T} (GeV/#it{c}); m({}^{3}_{#Lambda}H) (GeV/#it{c})", 50, 0, 7, 50, 2.96, 3.04)
h2Mass4LHnSigmaHe_MC = ROOT.TH2F("h2Mass4LHnSigmaHe_MC", r";n_{#sigma}^{TPC}({}^{3}He); m({}^{4}_{#Lambda}H) (GeV/#it{c})", 50, -4, 4, 30, 3.89, 3.97)

# for MC only
hPtGen = ROOT.TH1F("hPtGen", r";#it{p}_{T}^{gen} (GeV/#it{c})", 50, 0, 5)
utils.setHistStyle(hPtGen, ROOT.kAzure+2)
hCtGen = ROOT.TH1F("hCtGen", r";#it{c#tau}^{gen} (cm)", 50, 0, 40)
utils.setHistStyle(hCtGen, ROOT.kAzure+2)
hResolutionPt = ROOT.TH1F("hResolutionPt", r";(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, -0.2, 0.2)
utils.setHistStyle(hResolutionPt, ROOT.kAzure+2)
hResolutionPtvsPt = ROOT.TH2F("hResolutionPtvsPt", r";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionHe3PtvsPt = ROOT.TH2F("hResolutionHe3PtvsPt", r";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 5, 50, -0.2, 0.2)
hResolutionPiPtvsPt = ROOT.TH2F("hResolutionPiPtvsPt", r";#it{p}_{T}^{gen} (GeV/#it{c});(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}", 50, 0, 1, 50, -0.2, 0.2)
hResolutionDecVtxX = ROOT.TH1F("hResolutionDecVtxX", r"; Resolution Dec X", 50, -0.2, 0.2)
utils.setHistStyle(hResolutionDecVtxX, ROOT.kAzure+2)
hResolutionDecVtxY = ROOT.TH1F("hResolutionDecVtxY", r"; Resolution Dec Y", 50, -0.2, 0.2)
utils.setHistStyle(hResolutionDecVtxY, ROOT.kAzure+2)
hResolutionDecVtxZ = ROOT.TH1F("hResolutionDecVtxZ", r"; Resolution Dec Z", 50, -0.2, 0.2)
utils.setHistStyle(hResolutionDecVtxZ, ROOT.kAzure+2)

# filling histograms

# data
utils.fill_th1_hist(hPtRec, data_hdl, 'fPt')
utils.fill_th1_hist(hCtRec, data_hdl, 'fCt')
utils.fill_th1_hist(hCosPA, data_hdl, 'fCosPA')
utils.fill_th1_hist(hRadius, data_hdl, 'fDecRad')
utils.fill_th1_hist(hDecLen, data_hdl, 'fDecLen')
utils.fill_th1_hist(hNTPCclus, data_hdl, 'fNTPCclusHe')
utils.fill_th1_hist(hNSigHe, data_hdl, 'fNSigmaHe')
utils.fill_th1_hist(hMass3LH, data_hdl, 'fMassH3L')
utils.fill_th1_hist(hMass4LH, data_hdl, 'fMassH4L')
utils.fill_th2_hist(h2MassCosPA, data_hdl, 'fCosPA', 'fMassH3L')
utils.fill_th2_hist(h2MassDecLen, data_hdl, 'fDecLen', 'fMassH3L')
utils.fill_th2_hist(h2MassDCADaughters, data_hdl, 'fDcaV0Daug', 'fMassH3L')
utils.fill_th2_hist(h2MassDCAHePv, data_hdl, 'fDcaHe', 'fMassH3L')
utils.fill_th2_hist(h2MassPt, data_hdl, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHnSigmaHe, data_hdl, 'fNSigmaHe', 'fMassH4L')

# MC
utils.fill_th1_hist(hPtRec_MC, mc_hdl, 'fPt')
utils.fill_th1_hist(hCtRec_MC, mc_hdl, 'fCt')
utils.fill_th1_hist(hCosPA_MC, mc_hdl, 'fCosPA')
utils.fill_th1_hist(hRadius_MC, mc_hdl, 'fDecRad')
utils.fill_th1_hist(hDecLen_MC, mc_hdl, 'fDecLen')
utils.fill_th1_hist(hNTPCclus_MC, mc_hdl, 'fNTPCclusHe')
utils.fill_th1_hist(hNSigHe_MC, mc_hdl, 'fNSigmaHe')
utils.fill_th1_hist(hMass3LH_MC, mc_hdl, 'fMassH3L')
utils.fill_th1_hist(hMass4LH_MC, mc_hdl, 'fMassH4L')
utils.fill_th2_hist(h2MassCosPA_MC, mc_hdl, 'fCosPA', 'fMassH3L')
utils.fill_th2_hist(h2MassDecLen_MC,mc_hdl, 'fDecLen', 'fMassH3L')
utils.fill_th2_hist(h2MassDCADaughters_MC, mc_hdl, 'fDcaV0Daug', 'fMassH3L')
utils.fill_th2_hist(h2MassDCAHePv_MC, mc_hdl, 'fDcaHe', 'fMassH3L')
utils.fill_th2_hist(h2MassPt_MC, mc_hdl, 'fPt', 'fMassH3L')
utils.fill_th2_hist(h2Mass4LHnSigmaHe_MC, mc_hdl, 'fNSigmaHe', 'fMassH4L')

# MC truth
mc_hdl._full_data_frame.eval('resPt = (fPt - fAbsGenPt)/fAbsGenPt', inplace=True)
mc_hdl._full_data_frame.eval('ResDecX = (fXDecVtx - fGenXDecVtx)/fGenXDecVtx', inplace=True)
mc_hdl._full_data_frame.eval('ResDecY = (fYDecVtx - fGenYDecVtx)/fGenYDecVtx', inplace=True)
mc_hdl._full_data_frame.eval('ResDecZ = (fZDecVtx - fGenZDecVtx)/fGenZDecVtx', inplace=True)

utils.fill_th1_hist(hResolutionPt, mc_hdl, 'resPt')
utils.fill_th1_hist(hResolutionDecVtxX, mc_hdl, 'ResDecX')
utils.fill_th1_hist(hResolutionDecVtxY, mc_hdl, 'ResDecY')
utils.fill_th1_hist(hResolutionDecVtxZ, mc_hdl, 'ResDecZ')
utils.fill_th2_res_hist(hResolutionPtvsPt, mc_hdl, 'fAbsGenPt', 'fPt')


print("----------------------------------")
print("** Saving histograms **")

output_file.cd()


# Data

cCosPA = ROOT.TCanvas('cCosPA', 'cCosPA', 800, 600)
hCosPA.Draw('HISTO')
hCosPA.Write()
cCosPA.SaveAs(f'{qa_dir_name}/cCosPA.pdf')

cNTPCclus = ROOT.TCanvas('cNTPCclus', 'cNTPCclus', 800, 600)
hNTPCclus.Draw('HISTO')
hNTPCclus.Write()
cNTPCclus.SaveAs(f'{qa_dir_name}/cNTPCclus.pdf')

cMass3LH = ROOT.TCanvas('cMass3LH', 'cMass3LH', 800, 600)
hMass3LH.Draw('HISTO')
hMass3LH.Write()
cMass3LH.SaveAs(f'{qa_dir_name}/cMass3LH.pdf')

cMass4LH = ROOT.TCanvas('cMass4LH', 'cMass4LH', 800, 600)
hMass4LH.Draw('HISTO')
hMass4LH.Write()
cMass4LH.SaveAs(f'{qa_dir_name}/cMass4LH.pdf')

cPtRec = ROOT.TCanvas('cPtRec', 'cPtRec', 800, 600)
hPtRec.Draw('HISTO')
hPtRec.Write()
cPtRec.SaveAs(f'{qa_dir_name}/cPtRec.pdf')

cRadius = ROOT.TCanvas('cRadius', 'cRadius', 800, 600)
hRadius.Draw('HISTO')
hRadius.Write()
cRadius.SaveAs(f'{qa_dir_name}/cRadius.pdf')

cDecLen = ROOT.TCanvas('cDecLen', 'cDecLen', 800, 600)
hDecLen.Draw('HISTO')
hDecLen.Write()
cDecLen.SaveAs(f'{qa_dir_name}/cDecLen.pdf')

cNSigHe = ROOT.TCanvas('cNSigHe', 'cNSigHe', 800, 600)
hNSigHe.Draw('HISTO')
hNSigHe.Write()
cNSigHe.SaveAs(f'{qa_dir_name}/cNSigHe.pdf')

c2MassCosPA = ROOT.TCanvas('c2MassCosPA', 'c2MassCosPA', 800, 600)
c2MassCosPA.SetRightMargin(1.2)
h2MassCosPA.Draw('COLZ')
h2MassCosPA.Write()
c2MassCosPA.SaveAs(f'{qa_dir_name}/c2MassCosPA.pdf')

c2MassDecLen = ROOT.TCanvas('c2MassDecLen', 'c2MassDecLen', 800, 600)
c2MassDecLen.SetRightMargin(1.2)
h2MassDecLen.Draw('COLZ')
h2MassDecLen.Write()
c2MassDecLen.SaveAs(f'{qa_dir_name}/c2MassDecLen.pdf')

c2MassDCADaughters = ROOT.TCanvas('c2MassDCADaughters', 'c2MassDCADaughters', 800, 600)
c2MassDCADaughters.SetRightMargin(1.2)
h2MassDCADaughters.Draw('COLZ')
h2MassDCADaughters.Write()
c2MassDCADaughters.SaveAs(f'{qa_dir_name}/c2MassDCADaughters.pdf')

c2MassDCAHePv = ROOT.TCanvas('c2MassDCAHePv', 'c2MassDCAHePv', 800, 600)
c2MassDCAHePv.SetRightMargin(1.2)
h2MassDCAHePv.Draw('COLZ')
c2MassDCAHePv.Write()
c2MassDCAHePv.SaveAs(f'{qa_dir_name}/c2MassDCAHePv.pdf')

c2MassPt = ROOT.TCanvas('c2MassPt', 'c2MassPt', 800, 600)
c2MassPt.SetRightMargin(1.2)
h2MassPt.Draw('COLZ')
h2MassPt.Write()
c2MassPt.SaveAs(f'{qa_dir_name}/c2MassPt.pdf')

# MC

cCosPA_MC = ROOT.TCanvas('cCosPA_MC', 'cCosPA_MC', 800, 600)
hCosPA_MC.Draw('HISTO')
hCosPA_MC.Write()
cCosPA_MC.SaveAs(f'{qa_dir_name}/cCosPA_MC.pdf')

cNTPCclus_MC = ROOT.TCanvas('cNTPCclus_MC', 'cNTPCclus_MC', 800, 600)
hNTPCclus_MC.Draw('HISTO')
hNTPCclus_MC.Write()
cNTPCclus_MC.SaveAs(f'{qa_dir_name}/cNTPCclus_MC.pdf')

cMass3LH_MC = ROOT.TCanvas('cMass3LH_MC', 'cMass3LH_MC', 800, 600)
hMass3LH_MC.Draw('HISTO')
hMass3LH_MC.Write()
cMass3LH_MC.SaveAs(f'{qa_dir_name}/cMass3LH_MC.pdf')

cMass4LH_MC = ROOT.TCanvas('cMass4LH_MC', 'cMass4LH_MC', 800, 600)
hMass4LH_MC.Draw('HISTO')
hMass4LH_MC.Write()
cMass4LH_MC.SaveAs(f'{qa_dir_name}/cMass4LH_MC.pdf')

cPtRec_MC = ROOT.TCanvas('cPtRec_MC', 'cPtRec_MC', 800, 600)
hPtRec_MC.Draw('HISTO')
hPtRec_MC.Write()
cPtRec_MC.SaveAs(f'{qa_dir_name}/cPtRec_MC.pdf')

cRadius_MC = ROOT.TCanvas('cRadius_MC', 'cRadius_MC', 800, 600)
hRadius_MC.Draw('HISTO')
hRadius_MC.Write()
cRadius_MC.SaveAs(f'{qa_dir_name}/cRadius_MC.pdf')

cDecLen_MC = ROOT.TCanvas('cDecLen_MC', 'cDecLen_MC', 800, 600)
hDecLen_MC.Draw('HISTO')
hDecLen_MC.Write()
cDecLen_MC.SaveAs(f'{qa_dir_name}/cDecLen_MC.pdf')

cNSigHe_MC = ROOT.TCanvas('cNSigHe_MC', 'cNSigHe_MC', 800, 600)
hNSigHe_MC.Draw('HISTO')
hNSigHe_MC.Write()
cNSigHe_MC.SaveAs(f'{qa_dir_name}/cNSigHe_MC.pdf')

c2MassCosPA_MC = ROOT.TCanvas('c2MassCosPA_MC', 'c2MassCosPA_MC', 800, 600)
c2MassCosPA_MC.SetRightMargin(0.12)
c2MassCosPA_MC.SetLeftMargin(0.12)
h2MassCosPA_MC.Draw('COLZ')
h2MassCosPA_MC.Write()
c2MassCosPA_MC.SaveAs(f'{qa_dir_name}/c2MassCosPA_MC.pdf')

c2MassDecLen_MC = ROOT.TCanvas('c2MassDecLen_MC', 'c2MassDecLen_MC', 800, 600)
c2MassDecLen_MC.SetRightMargin(0.12)
c2MassDecLen_MC.SetLeftMargin(0.12)
h2MassDecLen_MC.Draw('COLZ')
h2MassDecLen_MC.Write()
c2MassDecLen_MC.SaveAs(f'{qa_dir_name}/c2MassDecLen_MC.pdf')

c2MassDCADaughters_MC = ROOT.TCanvas('c2MassDCADaughters_MC', 'c2MassDCADaughters_MC', 800, 600)
c2MassDCADaughters_MC.SetRightMargin(0.12)
c2MassDCADaughters_MC.SetLeftMargin(0.12)
h2MassDCADaughters_MC.Draw('COLZ')
h2MassDCADaughters_MC.Write()
c2MassDCADaughters_MC.SaveAs(f'{qa_dir_name}/c2MassDCADaughters_MC.pdf')

c2MassDCAHePv_MC = ROOT.TCanvas('c2MassDCAHePv_MC', 'c2MassDCAHePv_MC', 800, 600)
c2MassDCAHePv_MC.SetRightMargin(0.12)
c2MassDCAHePv_MC.SetLeftMargin(0.12)
h2MassDCAHePv_MC.Draw('COLZ')
c2MassDCAHePv_MC.Write()
c2MassDCAHePv_MC.SaveAs(f'{qa_dir_name}/c2MassDCAHePv_MC.pdf')

c2MassPt_MC = ROOT.TCanvas('c2MassPt_MC', 'c2MassPt_MC', 800, 600)
c2MassPt.SetRightMargin(0.12)
c2MassPt.SetLeftMargin(0.12)
h2MassPt_MC.Draw('COLZ')
h2MassPt_MC.Write()
c2MassPt_MC.SaveAs(f'{qa_dir_name}/c2MassPt_MC.pdf')

cPtGen = ROOT.TCanvas('cPtGen', 'cPtGen', 800, 600)
hPtGen.Draw('HISTO')
hPtGen.Write()
cPtGen.SaveAs(f'{qa_dir_name}/cPtGen.pdf')

cResolutionPt = ROOT.TCanvas('cResolutionPt', 'cResolutionPt', 800, 600)
hResolutionPt.Draw('HISTO')
hResolutionPt.Write()
cResolutionPt.SaveAs(f'{qa_dir_name}/cResolutionPt.pdf')

cResolutionPtvsPt = ROOT.TCanvas('cResolutionPtvsPt', 'cResolutionPtvsPt', 800, 600)
cResolutionPtvsPt.SetRightMargin(0.12)
cResolutionPtvsPt.SetLeftMargin(0.12)
hResolutionPtvsPt.Draw('COLZ')
hResolutionPtvsPt.Write()
cResolutionPtvsPt.SaveAs(f'{qa_dir_name}/cResolutionPtvsPt.pdf')

cResolutionDecVtxX = ROOT.TCanvas('cResolutionDecVtxX', 'cResolutionDecVtxX', 800, 600)
hResolutionDecVtxX.Draw('HISTO')
hResolutionDecVtxX.Write()
cResolutionDecVtxX.SaveAs(f'{qa_dir_name}/cResolutionDecVtxX.pdf')

cResolutionDecVtxY = ROOT.TCanvas('cResolutionDecVtxY', 'cResolutionDecVtxY', 800, 600)
hResolutionDecVtxY.Draw('HISTO')
hResolutionDecVtxY.Write()
cResolutionDecVtxY.SaveAs(f'{qa_dir_name}/cResolutionDecVtxY.pdf')

cResolutionDecVtxZ = ROOT.TCanvas('cResolutionDecVtxZ', 'cResolutionDecVtxZ', 800, 600)
hResolutionDecVtxZ.Draw('HISTO')
hResolutionDecVtxZ.Write()
cResolutionDecVtxZ.SaveAs(f'{qa_dir_name}/cResolutionDecVtxZ.pdf')

# general plots

input_analysis_results = ROOT.TFile(input_analysis_results_file)

hZvtx = input_analysis_results.Get('hyper-reco-task/hZvtx')
hZvtx.SetDirectory(0)
utils.setHistStyle(hZvtx, ROOT.kRed+1)
cZvtx = ROOT.TCanvas('cZvtx', 'cZvtx', 800, 600)
hZvtx.Draw('HISTO')
cZvtx.SaveAs(f'{qa_dir_name}/cZvtx.pdf')
output_file.cd()
hZvtx.Write()

hEvents = input_analysis_results.Get('hyper-reco-task/hEvents')
hEvents.SetDirectory(0)
utils.setHistStyle(hEvents, ROOT.kRed+1)
cEvents = ROOT.TCanvas('cEvents', 'cEvents', 800, 600)
hEvents.Draw('HISTO')
cEvents.SaveAs(f'{qa_dir_name}/cEvents.pdf')
output_file.cd()
hEvents.Write()

input_analysis_results_mc = ROOT.TFile(input_analysis_results_file_mc)

hZvtx_MC = input_analysis_results_mc.Get('hyper-reco-task/hZvtx')
hZvtx_MC.SetDirectory(0)
hZvtx_MC.SetName('hZvtx_MC')
utils.setHistStyle(hZvtx_MC, ROOT.kAzure+2)
cZvtx_MC = ROOT.TCanvas('cZvtx_MC', 'cZvtx_MC', 800, 600)
hZvtx_MC.Draw('HISTO')
cZvtx_MC.SaveAs(f'{qa_dir_name}/cZvtx_MC.pdf')
output_file.cd()
hZvtx_MC.Write()

hEvents_MC = input_analysis_results_mc.Get('hyper-reco-task/hEvents')
hEvents_MC.SetDirectory(0)
hEvents_MC.SetName('hEvents_MC')
utils.setHistStyle(hEvents_MC, ROOT.kAzure+2)
cEvents_MC = ROOT.TCanvas('cEvents_MC', 'cEvents_MC', 800, 600)
hEvents_MC.Draw('HISTO')
cEvents_MC.SaveAs(f'{qa_dir_name}/cEvents_MC.pdf')
output_file.cd()
hEvents_MC.Write()

