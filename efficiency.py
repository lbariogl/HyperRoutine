import ROOT
import utils
import argparse
import yaml


parser = argparse.ArgumentParser(
    description='Configure the parameters of the script.')
parser.add_argument('--input-file-matter', dest='input_file_matter',
                    help='path to the input file for matter.', default='../results/HypertritonResults_inj_gap_matter.root')
parser.add_argument('--input-file-antimatter', dest='input_file_antimatter',
                    help='path to the input file for antimatter.', default='../results/HypertritonResults_inj_gap_antimatter.root')
parser.add_argument('--output-dir', dest='output_dir',
                    help='path to the directory in which the output is stored.', default='../results')
parser.add_argument('--output-file', dest='output_file',
                    help='path to the output file.', default='../results/Efficiency_antimatter.root')
parser.add_argument('--config-file', dest='config_file',
                    help='path to the YAML file with configuration.', default='')
args = parser.parse_args()

input_file_matter_name = args.input_file_matter
input_file_antimatter_name = args.input_file_antimatter
output_dir_name = args.output_dir
output_file_name = args.output_file

if args.config_file != '':
    config_file = open(args.config_file, 'r')
    config = yaml.full_load(config_file)
    input_file_matter_name = config['input_file_matter']
    input_file_antimatter_name = config['input_file_antimatter']
    output_dir_name = config['output_dir']
    output_file_name = config['output_file']

print('**********************************')
print('    Running efficiency.py')
print('**********************************')

# matter
input_file_matter = ROOT.TFile(input_file_matter_name)
hPtRecMatter = input_file_matter.Get('hPtRec')
hPtRecMatter.SetDirectory(0)
hPtGenMatter = input_file_matter.Get('MC/hPtGen')
hPtGenMatter.SetDirectory(0)

hEffMatter = utils.computeEfficiency(
    hPtGenMatter, hPtRecMatter, 'hEffMatter', 2)
utils.setHistStyle(hEffMatter, ROOT.kBlack)

# antimatter
input_file_antimatter = ROOT.TFile(input_file_antimatter_name)
hPtRecAntimatter = input_file_antimatter.Get('hPtRec')
hPtRecAntimatter.SetDirectory(0)
hPtGenAntimatter = input_file_antimatter.Get('MC/hPtGen')
hPtGenAntimatter.SetDirectory(0)

hEffAntimatter = utils.computeEfficiency(
    hPtGenAntimatter, hPtRecAntimatter, 'hEffAntimatter', 2)
utils.setHistStyle(hEffAntimatter, ROOT.kRed)

# joined results
hPtRecTot = hPtRecMatter.Clone("hPtRecTot")
hPtRecTot.Add(hPtRecAntimatter)
hPtGenTot = hPtGenMatter.Clone("hPtGenTot")
hPtGenTot.Add(hPtGenAntimatter)

hEffTot = utils.computeEfficiency(hPtGenTot, hPtRecTot, 'hEffTot')
utils.setHistStyle(hEffTot, ROOT.kBlue)

cEff = ROOT.TCanvas('cEff', 'cEff', 800, 600)
hFrame = cEff.DrawFrame(0.001, 0.001, 10, 0.30,
                        r';#it{p}_{T} (GeV/#it{c});#epsilon #times Acc.')
cEff.SetLeftMargin(1.8)
hFrame.GetXaxis().SetTitleOffset(1.2)
hEffMatter.Draw('PE SAME')
hEffAntimatter.Draw('PE SAME')
lEff = ROOT.TLegend(0.21, 0.74, 0.48, 0.86, '', 'brNDC')
lEff.SetBorderSize(0)
lEff.AddEntry(hEffMatter, 'Matter', 'PE')
lEff.AddEntry(hEffAntimatter, 'Antimatter', 'PE')
lEff.Draw()

rec_file = ROOT.TFile('../results/fmazzasc_sim/hyp_study.root')
hEffRec = rec_file.Get('histoEffvsPt')
hEffRec.SetDirectory(0)
utils.setHistStyle(hEffRec, ROOT.kBlack)
hEffRecNoIts = rec_file.Get('histoEffvsPtNoITS')
hEffRecNoIts.SetDirectory(0)
utils.setHistStyle(hEffRecNoIts, ROOT.kRed)

cEffComp = ROOT.TCanvas('cEffComp', 'cEffComp', 800, 600)
hFrameComp = cEffComp.DrawFrame(0.001, 0.001, 10, 0.40,
                        r';#it{p}_{T} (GeV/#it{c});#epsilon #times Acc.')
cEffComp.SetLeftMargin(1.8)
hFrameComp.GetXaxis().SetTitleOffset(1.2)
hEffRec.Draw('HISTO SAME')
hEffRecNoIts.Draw('HISTO SAME')
hEffTot.Draw('PE SAME')
lEffComp = ROOT.TLegend(0.21, 0.74, 0.48, 0.86, '', 'brNDC')
lEffComp.SetBorderSize(0)
lEffComp.AddEntry(hEffRec, 'Reconstruction', 'L')
lEffComp.AddEntry(hEffRecNoIts, 'Reconstruction - no ITS-only tracks', 'L')
lEffComp.AddEntry(hEffTot, 'Task', 'PE')
lEffComp.Draw()

# save outputs
output_file = ROOT.TFile(f'{output_dir_name}/{output_file_name}', 'RECREATE')
hEffMatter.Write()
hEffAntimatter.Write()
hEffTot.Write()
cEff.Write()
hEffRec.Write()
hEffRecNoIts.Write()
cEffComp.Write()
