import ROOT
import numpy as np

h4l_mass = 3.921
h3l_mass = 1.991
pi_mass = 0.139570
he3_mass = 2.8083916
he4_mass = 3.7273794
daug_masses = np.array([pi_mass, he4_mass], dtype=np.float64)

debug = False

def print_vector(vec):
  print(f'Pt: {vec.Pt()} Eta: {vec.Eta()} Phi: {vec.Phi()} M: {vec.M()}')

h2MassMass = ROOT.TH2F('h2MassMass',';m({}^{3}_{#Lambda}H) (GeV/#it{c}^{2});m({}^{4}_{#Lambda}H) (GeV/#it{c}^{2})', 80, 2.92, 3.08, 60, 3.85, 4.01)

# get file for pT rejection
spectra_file = ROOT.TFile.Open('utils/heliumSpectraMB.root')
he3_spectrum = spectra_file.Get('fCombineHeliumSpecLevyFit_0-100')
spectra_file.Close()

for i in range(0, 100000):
  h4l_pt = 0.
  while(h4l_pt < 1.2) :
    h4l_pt = he3_spectrum.GetRandom()
  h4l_eta = ROOT.gRandom.Uniform(-0.8, 0.8)
  h4l_phi = ROOT.gRandom.Uniform(0., 2 * np.pi)

  h4l_vec = ROOT.TLorentzVector()
  h4l_vec.SetPtEtaPhiM(h4l_pt, h4l_eta, h4l_phi, h4l_mass)

  event = ROOT.TGenPhaseSpace()
  event.SetDecay(h4l_vec, 2, daug_masses)

  event.Generate()
  pion_vec = event.GetDecay(0)
  he4_vec = event.GetDecay(1)

  if debug:
    print('************************')
    print('************************')
    print(f'Event {i}')
    print('************************')
    print('h4l')
    print_vector(h4l_vec)
    print('***********************')
    print('Before smearing')
    print('pion')
    print_vector(pion_vec)
    print('he4')
    print_vector(he4_vec)

  # smearing
  pion_vec.SetPtEtaPhiM(pion_vec.Pt()*ROOT.gRandom.Gaus(1, 0.05), pion_vec.Eta(), pion_vec.Phi(), pi_mass)

  he4_vec.SetPtEtaPhiM(he4_vec.Pt()*ROOT.gRandom.Gaus(1, 0.05), he4_vec.Eta(), he4_vec.Phi(), he4_mass)

  if debug:
    print('***********************')
    print('After smearing')
    print('pion')
    print_vector(pion_vec)
    print('he4')
    print_vector(he4_vec)

  he3_vec = ROOT.TLorentzVector()
  he3_vec.SetPtEtaPhiM(he4_vec.Pt(), he4_vec.Eta(), he4_vec.Phi(), he3_mass)

  h4l_vec_rec = pion_vec + he4_vec
  h3l_vec_rec = pion_vec + he3_vec

  if debug:
    print('************************')
    print('rec 4lh')
    print_vector(h4l_vec_rec)
    print('rec 3lh')
    print_vector(h3l_vec_rec)
    print('\n\n')

  h2MassMass.Fill(h3l_vec_rec.M(), h4l_vec_rec.M())

output_file = ROOT.TFile('../results/check_he_mass.root', 'recreate')
h2MassMass.Write()
