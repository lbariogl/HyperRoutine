import ROOT


ffile_input = ROOT.TFile("pt_analysis_antimat_sigloss_levy.root")
func_levy = ffile_input.Get('std/fitLT3')

ffile2_input = ROOT.TFile("pt_analysis_antimat_sigloss_mt.root")
yield_stat = ffile2_input.Get('std/hStat')
yield_stat.SetDirectory(0)
yield_syst = ffile2_input.Get('std/hSyst')
yield_syst.SetDirectory(0)

func_mt = ffile2_input.Get('std/mtexpo')
func_levy.SetLineColor(ROOT.kRed)
## create canvas



pinfo_alice = ROOT.TPaveText(0.5, 0.7, 0.9, 0.85, 'NDC')
pinfo_alice.SetBorderSize(0)
pinfo_alice.SetFillStyle(0)
pinfo_alice.SetTextAlign(11)
pinfo_alice.SetTextFont(42)
pinfo_alice.AddText('ALICE Preliminary')
pinfo_alice.AddText('Run 3, pp #sqrt{#it{s}} = 13.6 TeV')
pinfo_alice.Draw()

cFinalSpectrum = ROOT.TCanvas('cFinalSpectrum', 'cFinalSpectrum', 800, 600)
cFinalSpectrum.SetLogy()
cFinalSpectrum.DrawFrame(0.5,0.1e-09, 6, 1.5 * yield_stat.GetMaximum(), r';#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}} #frac{#it{d}N}{#it{d}y#it{d}#it{p}_{T}} (GeV/#it{c})^{-1}')
yield_stat.Draw('PEX0 same')
yield_syst.Draw('PE2 same')
func_mt.Draw('same')
pinfo_alice.Draw()
leg = ROOT.TLegend(0.5, 0.5, 0.8, 0.7)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.AddEntry(func_mt, '#it{m}_{T} - exponential', 'L')
leg.AddEntry(yield_stat, '{}^{3}_{#bar{#Lambda}}#bar{H}', 'PE')
leg.Draw()

cFuncComparison = ROOT.TCanvas('cFuncComparison', 'cFuncComparison', 800, 600)
cFuncComparison.DrawFrame(0.1e-09, 0, 6, 1.8 * yield_stat.GetMaximum(), r';#it{p}_{T} (GeV/#it{c});#frac{1}{N_{ev}}#frac{#it{d}N}{#it{d}y#it{d}#it{p}_{T}} (GeV/#it{c})^{-1}')
func_mt.Draw('same')
func_levy.Draw('same')
yield_stat.Draw('PEX0 same')
yield_syst.Draw('PE2 same')



leg2 = ROOT.TLegend(0.6, 0.6, 0.9, 0.85)
leg2.SetFillStyle(0)
leg2.SetBorderSize(0)
leg2.SetTextFont(42)
leg2.AddEntry(func_mt, 'MT exponential', 'L')
leg2.AddEntry(func_levy, 'Levy-Tsallis', 'L')
leg2.Draw()




outfile = ROOT.TFile("spectra_prel.root", "RECREATE")
cFinalSpectrum.Write()
cFuncComparison.Write()
outfile.Close()