import ROOT
import argparse
import os

ROOT.SetOptStat(0)

parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--mc', dest='mc', action='store_true',
                    help='if True MC information is stored.', default=False)
parser.add_argument('--input-file', dest='input_file',
                    help='path to the input file.', default='../results/HypertritonResults_tot.root')
parser.add_argument('--output-dir', dest='output_dir',
                    help='path to the output directory.', default='../results/plots')
args = parser.parse_args()

def setHistStyle(histo, color = ROOT.kBlack, line_width = 1, marker_style = 20):
    histo.SetLineColor(color)
    histo.SetLineWidth(line_width)
    histo.SetMarkerColor(color)
    histo.SetMarkerStyle(marker_style)
    histo.GetXaxis().SetTitleOffset(1.2)

mc = args.mc
input_file_name = args.input_file
outout_dir = args.output_dir

if not os.path.exists(outout_dir):
   os.makedirs(outout_dir)

input_file = ROOT.TFile(input_file_name)

hCosPA = input_file.Get('hCosPA')
setHistStyle(hCosPA)
cCosPA = ROOT.TCanvas('cCosPA', 'cCosPA', 800, 600)
hCosPA.Draw('HISTO')
cCosPA.SaveAs(f'{outout_dir}/cCosPA.pdf')

hNTPCclus = input_file.Get('hNTPCclus')
setHistStyle(hNTPCclus)
cNTPCclus = ROOT.TCanvas('cNTPCclus', 'cNTPCclus', 800, 600)
hNTPCclus.Draw('HISTO')
cNTPCclus.SaveAs(f'{outout_dir}/cNTPCclus.pdf')

hMass3LH = input_file.Get('h_3lh_mass')
setHistStyle(hMass3LH)
cMass3LH = ROOT.TCanvas('cMass3LH', 'cMass3LH', 800, 600)
hMass3LH.Draw('HISTO')
cMass3LH.SaveAs(f'{outout_dir}/cMass3LH.pdf')

hMass4LH = input_file.Get('h_4lh_mass')
setHistStyle(hMass4LH)
cMass4LH = ROOT.TCanvas('cMass4LH', 'cMass4LH', 800, 600)
hMass4LH.Draw('HISTO')
cMass4LH.SaveAs(f'{outout_dir}/cMass4LH.pdf')

hPtRec = input_file.Get('hPtRec')
setHistStyle(hPtRec)
cPtRec = ROOT.TCanvas('cPtRec', 'cPtRec', 800, 600)
hPtRec.Draw('HISTO')
cPtRec.SaveAs(f'{outout_dir}/cPtRec.pdf')

hRadius = input_file.Get('hRadius')
setHistStyle(hRadius)
cRadius = ROOT.TCanvas('cRadius', 'cRadius', 800, 600)
hRadius.Draw('HISTO')
cRadius.SaveAs(f'{outout_dir}/cRadius.pdf')

hDecLen = input_file.Get('hDecLen')
setHistStyle(hDecLen)
cDecLen = ROOT.TCanvas('cDecLen', 'cDecLen', 800, 600)
hDecLen.Draw('HISTO')
cDecLen.SaveAs(f'{outout_dir}/cDecLen.pdf')

hNSigHe = input_file.Get('hNSigmaHe')
setHistStyle(hNSigHe)
cNSigHe = ROOT.TCanvas('cNSigHe', 'cNSigHe', 800, 600)
hNSigHe.Draw('HISTO')
cNSigHe.SaveAs(f'{outout_dir}/cNSigHe.pdf')

h2MassCosPA = input_file.Get('h2MassCosPA')
c2MassCosPA = ROOT.TCanvas('c2MassCosPA', 'c2MassCosPA', 800, 600)
c2MassCosPA.SetRightMargin(1.2)
h2MassCosPA.Draw('COLZ')
c2MassCosPA.SaveAs(f'{outout_dir}/c2MassCosPA.pdf')

h2MassDecLen = input_file.Get('h2MassDecLen')
c2MassDecLen = ROOT.TCanvas('c2MassDecLen', 'c2MassDecLen', 800, 600)
c2MassDecLen.SetRightMargin(1.2)
h2MassDecLen.Draw('COLZ')
c2MassDecLen.SaveAs(f'{outout_dir}/c2MassDecLen.pdf')

h2MassDCADaughters = input_file.Get('h2MassDCADaughters')
c2MassDCADaughters = ROOT.TCanvas('c2MassDCADaughters', 'c2MassDCADaughters', 800, 600)
c2MassDCADaughters.SetRightMargin(1.2)
h2MassDCADaughters.Draw('COLZ')
c2MassDCADaughters.SaveAs(f'{outout_dir}/c2MassDCADaughters.pdf')

h2MassDCAHePv = input_file.Get('h2MassDCAHe')
c2MassDCAHePv = ROOT.TCanvas('c2MassDCAHePv', 'c2MassDCAHePv', 800, 600)
c2MassDCAHePv.SetRightMargin(1.2)
h2MassDCAHePv.Draw('COLZ')
c2MassDCAHePv.SaveAs(f'{outout_dir}/c2MassDCAHePv.pdf')

h2MassPt = input_file.Get('h2MassPt')
c2MassPt = ROOT.TCanvas('c2MassPt', 'c2MassPt', 800, 600)
c2MassPt.SetRightMargin(1.2)
h2MassPt.Draw('COLZ')
c2MassPt.SaveAs(f'{outout_dir}/c2MassPt.pdf')

if mc:
    hPtGen = input_file.Get('MC/hPtGen')
    setHistStyle(hPtGen)
    cPtGen = ROOT.TCanvas('cPtGen', 'cPtGen', 800, 600)
    hPtGen.Draw('HISTO')
    cPtGen.SaveAs(f'{outout_dir}/cPtGen.pdf')

    hResolutionPt = input_file.Get('MC/hResolutionPt')
    setHistStyle(hResolutionPt)
    cResolutionPt = ROOT.TCanvas('cResolutionPt', 'cResolutionPt', 800, 600)
    hResolutionPt.Draw('HISTO')
    cResolutionPt.SaveAs(f'{outout_dir}/cResolutionPt.pdf')

    hResolutionPtvsPt = input_file.Get('MC/hResolutionPtvsPt')
    cResolutionPtvsPt = ROOT.TCanvas('cResolutionPtvsPt', 'cResolutionPtvsPt', 800, 600)
    cResolutionPtvsPt.SetRightMargin(1.2)
    hResolutionPtvsPt.Draw('COLZ')
    cResolutionPtvsPt.SaveAs(f'{outout_dir}/cResolutionPtvsPt.pdf')

    hResolutionP = input_file.Get('MC/hResolutionP')
    setHistStyle(hResolutionP)
    cResolutionP = ROOT.TCanvas('cResolutionP', 'cResolutionP', 800, 600)
    hResolutionP.Draw('HISTO')
    cResolutionP.SaveAs(f'{outout_dir}/cResolutionP.pdf')

    hResolutionPvsP = input_file.Get('MC/hResolutionPvsP')
    cResolutionPvsP = ROOT.TCanvas('cResolutionPvsP', 'cResolutionPvsP', 800, 600)
    cResolutionPvsP.SetRightMargin(1.2)
    hResolutionPvsP.Draw('COLZ')
    cResolutionPvsP.SaveAs(f'{outout_dir}/cResolutionPvsP.pdf')

    hResolutionDecVtxX = input_file.Get('MC/hResolutionDecVtxX')
    setHistStyle(hResolutionDecVtxX)
    cResolutionDecVtxX = ROOT.TCanvas('cResolutionDecVtxX', 'cResolutionDecVtxX', 800, 600)
    hResolutionDecVtxX.Draw('HISTO')
    cResolutionDecVtxX.SaveAs(f'{outout_dir}/cResolutionDecVtxX.pdf')

    hResolutionDecVtxY = input_file.Get('MC/hResolutionDecVtxY')
    setHistStyle(hResolutionDecVtxY)
    cResolutionDecVtxY = ROOT.TCanvas('cResolutionDecVtxY', 'cResolutionDecVtxY', 800, 600)
    hResolutionDecVtxY.Draw('HISTO')
    cResolutionDecVtxY.SaveAs(f'{outout_dir}/cResolutionDecVtxY.pdf')

    hResolutionDecVtxZ = input_file.Get('MC/hResolutionDecVtxZ')
    setHistStyle(hResolutionDecVtxZ)
    cResolutionDecVtxZ = ROOT.TCanvas('cResolutionDecVtxZ', 'cResolutionDecVtxZ', 800, 600)
    hResolutionDecVtxZ.Draw('HISTO')
    cResolutionDecVtxZ.SaveAs(f'{outout_dir}/cResolutionDecVtxZ.pdf')
