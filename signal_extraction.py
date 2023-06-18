import ROOT
import uproot
import numpy as np
import pandas as pd

import argparse
import yaml

import sys
sys.path.append('utils')
import utils as utils



ROOT.gROOT.LoadMacro('utils/RooCustomPdfs/RooDSCBShape.cxx++')
from ROOT import RooDSCBShape

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.")
parser.set_defaults(config_file='')
args = parser.parse_args()


config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
input_parquet_data = config['input_parquet_data']
input_parquet_mc = config['input_parquet_mc']
output_dir = config['output_dir']
output_file = config['output_file']
matter_type = config['matter_type']
histo_maximum = config['histo_maximum']

if matter_type=="matter":
    inv_mass_string = "#it{M}_{^{3}He+#pi^{-}}"

elif matter_type=="antimatter":
    inv_mass_string = "#it{M}_{^{3}#bar{He}+#pi^{+}}"

else:
    inv_mass_string = "#it{M}_{^{3}He+#pi^{-}} + c.c."





mass = ROOT.RooRealVar('m', inv_mass_string, 2.96, 3.04, 'GeV/c^{2}')
mu = ROOT.RooRealVar('mu', 'hypernucl mass', 2.98, 3.0, 'GeV/c^{2}')
sigma = ROOT.RooRealVar('sigma', 'hypernucl width', 0.001, 0.004, 'GeV/c^{2}')
a1 = ROOT.RooRealVar('a1', 'a1', 0, 5.)
a2 = ROOT.RooRealVar('a2', 'a2', 0, 10.)
n1 = ROOT.RooRealVar('n1', 'n1', 1, 10.)
n2 = ROOT.RooRealVar('n2', 'n2', 1, 10.)
# signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)
signal = ROOT.RooDSCBShape('cb', 'cb', mass, mu, sigma, a1, n1, a2, n2)
# signal = ROOT.RooGaussian('gaus', 'gaus', mass, mu, sigma)
c0 = ROOT.RooRealVar('c0', 'constant c0', -1., 1)
c1 = ROOT.RooRealVar('c1', 'constant c1', -1., 1)
c2 = ROOT.RooRealVar('c2', 'constant c2', -1., 1)
c3 = ROOT.RooRealVar('c3', 'constant c3', -1., 1)



# ### fix DSCB parameters to MC
df_mc = pd.read_parquet(input_parquet_mc)
mass_roo_mc = utils.ndarray2roo(np.array(df_mc['fMassH3L'].values, dtype=np.float64), mass, "histo_mc")
signal.fitTo(mass_roo_mc, ROOT.RooFit.Range(2.96, 3.01))
##set a1, a2, n1, n2 to constant
a1.setConstant()
a2.setConstant()
n1.setConstant()
n2.setConstant()

output_file = ROOT.TFile(output_dir + "/" + output_file, 'recreate')
## save histo with mc invariant mass
frame = mass.frame(80)
mass_roo_mc.plotOn(frame)
signal.plotOn(frame)
## add fit parameters to the frame (in a TPaveText)
fit_param = ROOT.TPaveText(0.6, 0.6, 0.9, 0.9, 'NDC')
fit_param.SetBorderSize(0)
fit_param.SetFillStyle(0)
fit_param.SetTextAlign(12)

fit_param.AddText('#mu = ' + f'{mu.getVal()*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
fit_param.AddText('#sigma = ' + f'{sigma.getVal()*1e3:.2f} #pm {sigma.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')

frame.addObject(fit_param)
frame.Write("histo_mc")






background = ROOT.RooChebychev('bkg', 'pol1 bkg', mass, ROOT.RooArgList(c0))
n = ROOT.RooRealVar('n', 'n const', 0.01, 0.4)
# define the fit funciton and perform the actual fit
fit_function = ROOT.RooAddPdf('total_pdf', 'signal + background', ROOT.RooArgList(signal, background), ROOT.RooArgList(n))

### if input_parquet_data is a list of files, loop over them
if type(input_parquet_data) == list:
    mass_array = np.array([])
    for file in input_parquet_data:
        df = pd.read_parquet(file)
        mass_array = np.append(mass_array, df['fMassH3L'].values)
        mass_array = np.array(mass_array, dtype=np.float64)
else:
    df = pd.read_parquet(input_parquet_data)
    mass_array = np.array(df['fMassH3L'].values, dtype=np.float64)


mass_roo_data = utils.ndarray2roo(mass_array, mass)


utils.fit_and_plot(mass_roo_data, mass, fit_function, signal, background, sigma, mu, n, n_ev=330, matter_type=matter_type, histo_maximum=histo_maximum)