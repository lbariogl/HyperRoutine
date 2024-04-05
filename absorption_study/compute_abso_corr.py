#!/usr/bin/env python3
import argparse
import yaml
import uproot
import ROOT
import numpy as np

import sys
sys.path.append('../utils')
import utils as utils


parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
pt_bins = config['pt_bins']
pt_bins = np.array(pt_bins, dtype=np.float32)

HE_3_MASS = 2.809230089
reweight_pt = False


h_he_rad = {"mat": ROOT.TH1F('h_he_rad_mat', 'He3 radius (cm)', 400, 0, 100),
                "antimat": ROOT.TH1F('h_he_rad_antimat', 'Absorption radius (cm)', 400, 0, 100)}
h_he_ct = {"mat": ROOT.TH1F('h_he_ct_mat', 'He3 ct (cm)', 20, 0, 40),
                "antimat": ROOT.TH1F('h_he_ct_antimat', 'Absorption ct (cm)', 20, 0, 40)}
h_he_pt = {"mat": ROOT.TH1F('h_he_pt_mat', 'He3 p_{T} (GeV/c)', len(pt_bins)-1, pt_bins),
                "antimat": ROOT.TH1F('h_he_pt_antimat', 'Absorption p_{T} (GeV/c)', len(pt_bins)-1, pt_bins)}


## histos for only the absorbed particles before the decay with h3l lifetime
h_abso_ct = {"mat": ROOT.TH1F('h_abso_ct_mat', '; ct (cm)', 20, 0, 40),
                "antimat": ROOT.TH1F('h_abso_ct_antimat', '; ct (cm)', 20, 0, 40)}

h_abso_pt = {"mat": ROOT.TH1F('h_abso_pt_mat', '; p_{T} (GeV/c)', len(pt_bins)-1, pt_bins),
                "antimat": ROOT.TH1F('h_abso_pt_antimat', '; p_{T} (GeV/c)', len(pt_bins)-1, pt_bins)}    


# mc input file
mc_file = '../../results/absorption/absorption_tree_x5.root'
tree = uproot.open(mc_file)['he3candidates'].arrays(library='pd')
# reweight MC pT spectrum
spectra_file = ROOT.TFile.Open('../utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()

utils.reweight_pt_spectrum(tree, 'pt', he3_spectrum)
# tree = tree.query('rej>0')

for (pt, eta, phi, absoX, absoY, absoZ, process, pdg) in zip(tree['pt'], tree['eta'], tree['phi'], tree['absoX'], tree['absoY'], tree['absoZ'], tree['process'], tree['pdg']):
    
    mat_type = 'mat' if pdg > 0 else 'antimat'
    he3TLV = ROOT.TLorentzVector()
    he3TLV.SetPtEtaPhiM(pt, eta, phi, HE_3_MASS)
    he3p = he3TLV.P()
    absoR = np.sqrt(absoX**2 + absoY**2)
    absoL = np.sqrt(absoX**2 + absoY**2 + absoZ**2)
    absoCt = absoL * HE_3_MASS / he3p
    decCt = ROOT.gRandom.Exp(7.6)

    pt_weight = 1 if not reweight_pt else he3_spectrum.Eval(pt)

    h_he_rad[mat_type].Fill(absoR)
    h_he_ct[mat_type].Fill(decCt)
    h_he_pt[mat_type].Fill(pt, pt_weight)

    if absoCt > decCt:
        h_abso_ct[mat_type].Fill(decCt)
        h_abso_pt[mat_type].Fill(pt, pt_weight)


outfile = ROOT.TFile('../../results/absorption/absorption_histos_x5.root', 'RECREATE')

for mat_type in ['mat', 'antimat']:
    h_he_rad[mat_type].Write()
    h_he_ct[mat_type].Write()
    h_he_pt[mat_type].Write()
    h_abso_ct[mat_type].Write()
    h_abso_pt[mat_type].Write()

    ## compute the correction factor
    h_abso_frac_ct = h_abso_ct[mat_type].Clone(f'h_abso_frac_ct_{mat_type}')
    h_abso_frac_ct.Divide(h_he_ct[mat_type])



    h_abso_frac_ct.GetYaxis().SetTitle('1 - f_{abs}')
    h_abso_frac_ct.Write()

    h_abso_frac_pt = h_abso_pt[mat_type].Clone(f'h_abso_frac_pt_{mat_type}')
    h_abso_frac_pt.Divide(h_he_pt[mat_type])
    h_abso_frac_pt.GetYaxis().SetTitle('1 - f_{abs}')
    h_abso_frac_pt.GetYaxis().SetTitleSize(0.06)


    ## binomial errors bin by bin
    for ibin in range(1, h_abso_frac_pt.GetNbinsX()+1):
        n = h_he_pt[mat_type].GetBinContent(ibin)
        k = h_abso_pt[mat_type].GetBinContent(ibin)
        p = h_abso_frac_pt.GetBinContent(ibin)
        err = np.sqrt(p*(1-p)/n)
        h_abso_frac_pt.SetBinError(ibin, err)

    h_abso_frac_pt.Write()

outfile.Close()

    





