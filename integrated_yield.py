## script to get the integrated yield of the hypertriton spectrum
import ROOT
import numpy as np
import yaml
import argparse

## levy-tsallis is defined in the file AdditionalFunctions.h
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine('.L utils/AdditionalFunctions.h++')
from ROOT import LevyTsallis

parser = argparse.ArgumentParser(description='Configure the parameters of the script.')
parser.add_argument('--config-file', dest='config_file', help="path to the YAML file with configuration.", default='')
args = parser.parse_args()
if args.config_file == "":
    print('** No config file provided. Exiting. **')
    exit()

config_file = open(args.config_file, 'r')
config = yaml.full_load(config_file)
input_file = ROOT.TFile(config['input_file_spectrum'])
h3l_spectrum = input_file.Get(config['histo_stat_name'])
h3l_spectrum_syst = input_file.Get(config['histo_syst_name'])
syst_sig_extr = input_file.Get('h_yields')
h3l_spectrum.SetDirectory(0)
syst_sig_extr.SetDirectory(0)
syst_sig_extr.Fit('gaus', 'R')

h3l_mass = 2.99131
mt_expo = ROOT.TF1('mtexpo', '[2]*x*exp(-TMath::Sqrt(([0]*[0]+x*x))/[1])', 0., 10)
mt_expo.FixParameter(0, h3l_mass)
mt_expo.SetParLimits(1, 0.1, 1)
mt_expo.SetParLimits(2, 1.e-08, 1)

pt_expo = ROOT.TF1('ptexpo', '[1]*x*exp(-x/[0])', 0., 10)
pt_expo.SetParLimits(0, 0.1, 1)
pt_expo.SetParLimits(1, 1.e-08, 1)
levy = LevyTsallis('levy', h3l_mass)
levy.SetParLimits(1, 5, 10)
levy.SetParLimits(3, 1e-08, 4e-08)
## fit the spectrum with all the functions and get the integral of the fit functions

integral_histo = 0
integral_histo_error = 0
for i in range(1, h3l_spectrum.GetNbinsX()+1):
    bin_width = h3l_spectrum.GetXaxis().GetBinWidth(i)
    bin_content = h3l_spectrum.GetBinContent(i)
    integral_histo += bin_content * bin_width
    bin_error = h3l_spectrum.GetBinError(i)
    integral_histo_error += (bin_error * bin_width)**2
integral_histo_error = np.sqrt(integral_histo_error)
print(f'Integral of the histogram: {integral_histo} +/- {integral_histo_error}')

lowest_pt_edge = h3l_spectrum.GetXaxis().GetBinLowEdge(1)

h3l_spectrum.Fit(mt_expo, 'R')
mt_expo_integral = [mt_expo.Integral(0., lowest_pt_edge), mt_expo.IntegralError(0., lowest_pt_edge)]
h3l_spectrum.Fit(pt_expo, 'R')
pt_expo_integral = [pt_expo.Integral(0., lowest_pt_edge), pt_expo.IntegralError(0., lowest_pt_edge)]
h3l_spectrum.Fit(levy, 'R')
levy_integral = [levy.Integral(0, lowest_pt_edge), levy.IntegralError(0, lowest_pt_edge)]

rms_extr= np.std([mt_expo_integral[0], pt_expo_integral[0], levy_integral[0]])
rel_unc_extr = rms_extr / np.mean([mt_expo_integral[0], pt_expo_integral[0], levy_integral[0]])
syst_sig_extr_rms = syst_sig_extr.GetFunction('gaus').GetParameter(2)
syst_sig_extr_mean = syst_sig_extr.GetFunction('gaus').GetParameter(1)


## do a gaussian sampling of the statistical distribution, for each toy compute the yield with the three functions and get the RMS of the distribution

n_toys = 10000
yield_toy_mt = []
yield_toy_pt = []
yield_toy_levy = []
for i in range(n_toys):
    # sample the histogram
    h3l_spectrum_toy = h3l_spectrum.Clone('h3l_spectrum_toy')
    for j in range(1, h3l_spectrum_toy.GetNbinsX()+1):
        content = np.random.normal(h3l_spectrum.GetBinContent(j), h3l_spectrum.GetBinError(j))
        if content < 0:
            content = 0.
        h3l_spectrum_toy.SetBinContent(j, content)
    # fit with the three functions
    h3l_spectrum_toy.Fit(mt_expo, 'SRQ')
    mt_yield = mt_expo.Integral(0., lowest_pt_edge)
    h3l_spectrum_toy.Fit(pt_expo, 'RSQ')
    pt_yield = pt_expo.Integral(0., lowest_pt_edge)
    h3l_spectrum_toy.Fit(levy, 'RSQ')
    levy_yield = levy.Integral(0., lowest_pt_edge)

    yield_toy_mt.append(mt_yield)
    yield_toy_pt.append(pt_yield)
    yield_toy_levy.append(levy_yield)

yield_toy_mt = np.array(yield_toy_mt)
yield_toy_pt = np.array(yield_toy_pt)
yield_toy_levy = np.array(yield_toy_levy)

histo_toys = []
for yields, name in zip([yield_toy_mt, yield_toy_pt, yield_toy_levy], ['mt_expo', 'pt_expo', 'levy']):
    histo = ROOT.TH1F(f'histo_toys_{name}', f'histo_toys_{name}', 100, np.min(yields), np.max(yields))
    for y in yields:
        histo.Fill(y)
    histo_toys.append(histo)

yield_final = levy_integral[0] + integral_histo
stat_unc = integral_histo_error

fit_function_syst = np.std([np.mean(yield_toy_mt), np.mean(yield_toy_pt), np.mean(yield_toy_levy)])
histo_fit_func_syst = ROOT.TH1F('histo_fit_func_syst', 'histo_fit_func_syst', 1, 0, 1)
histo_fit_func_syst.SetBinContent(1, fit_function_syst)
absorption_relative_syst = 0.03
br_relative_syst = 0.08
extrapolation_syst = yield_toy_levy.std()

normalistion_relative_unc = 0.1

total_syst = np.sqrt((br_relative_syst * yield_final)**2 + (absorption_relative_syst * yield_final)**2 + (syst_sig_extr_rms)**2 + (fit_function_syst)**2 + (extrapolation_syst)**2 + (normalistion_relative_unc * yield_final)**2)
relative_syst = total_syst / yield_final


print('--------------------------------------------')
print('Final result for integrated yield:')
print(f'  - dN/dy = {yield_final:.6e} +/- {stat_unc:.6e} (stat) +/- {total_syst:.6e} (syst)')
print('--------------------------------------------')

print('Breakdown of systematic uncertainties:')
print(f'  - signal selection and extraction: {syst_sig_extr_rms:.6e} ({syst_sig_extr_rms / yield_final * 100:.2f} %)')
print(f'  - fit function choice: {fit_function_syst:.6e} ({fit_function_syst / yield_final * 100:.2f} %)')
print(f'  - absorption correction: {(absorption_relative_syst * yield_final):.6e} ({absorption_relative_syst * 100:.2f} %)')
print(f'  - branching ratio: {(br_relative_syst * yield_final):.6e} ({br_relative_syst * 100:.2f} %)')
print(f'  - extrapolation to zero pT: {extrapolation_syst:.6e} ({extrapolation_syst / yield_final * 100:.2f} %)')
print('--------------------------------------------')

## plot all the fit functions and the datapoint into a single canvas
canvas = ROOT.TCanvas('canvas', 'canvas', 800, 600)
canvas.SetLogy()
## draw a new frame
canvas.DrawFrame(0, h3l_spectrum.GetMinimum(), 10, h3l_spectrum.GetMaximum()*10, f';{h3l_spectrum.GetXaxis().GetTitle()};{h3l_spectrum.GetYaxis().GetTitle()}')
## set x-axis range from 0 to 6
# h3l_spectrum.GetXaxis().SetRangeUser(0., 6.)
h3l_spectrum.SetMarkerStyle(20)
h3l_spectrum.SetMarkerSize(0.5)
h3l_spectrum.SetMarkerColor(ROOT.kBlack)
h3l_spectrum.SetLineColor(ROOT.kBlack)
## remove fit function attached to the histogram
h3l_spectrum.SetStats(0)
h3l_spectrum.GetListOfFunctions().Delete()
h3l_spectrum.Draw("same")
mt_expo.SetLineColor(ROOT.kRed)
mt_expo.SetLineWidth(2)
mt_expo.Draw('same')
pt_expo.SetLineColor(ROOT.kBlue)
pt_expo.SetLineWidth(2)
pt_expo.Draw('same')
levy.SetLineColor(ROOT.kGreen)
levy.SetLineWidth(2)
levy.Draw('same')
leg_canvas = ROOT.TLegend(0.15, 0.6, 0.4, 0.85)
leg_canvas.SetFillStyle(0)
leg_canvas.SetBorderSize(0)
leg_canvas.SetTextFont(42)
leg_canvas.SetMargin(0.1)
leg_canvas.SetTextSize(0.037)
leg_canvas.AddEntry(levy, 'Levy-Tsallis fit', 'L')
leg_canvas.AddEntry(mt_expo, '#it{m}_{T} exponential fit', 'L')
leg_canvas.AddEntry(pt_expo, '#it{p}_{T} exponential fit', 'L')
leg_canvas.Draw()



outfile = ROOT.TFile(config['output_file'], 'RECREATE')
h3l_spectrum.Write('hStat')
h3l_spectrum_syst.Write('hSyst')
syst_sig_extr.Write('hSystDistr')

mt_expo.Write('mt_expo')
pt_expo.Write('pt_expo')
levy.Write('levy')
canvas.Write('canvas')

histo_fit_func_syst.Write('histo_fit_func_syst')
for histo in histo_toys:
    histo.Write()