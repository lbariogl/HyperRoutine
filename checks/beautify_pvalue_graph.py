import ROOT
ROOT.gROOT.SetBatch(True)
def compute_pvalue_from_sign(significance):
    return ROOT.Math.chisquared_cdf_c(significance**2, 1) / 2

# get TGraphErrors from file

ff = ROOT.TFile("../../results/fit_results_antimat_4lh.root")
gr = ff.Get("p0_values")
gr.SetMaximum(1)
gr.SetMinimum(1e-6)



## draw horizontal lines for significance levels
pval_1sigma = compute_pvalue_from_sign(1)
pval_2sigma = compute_pvalue_from_sign(2)
pval_3sigma = compute_pvalue_from_sign(3)
pval_4sigma = compute_pvalue_from_sign(4)
pval_5sigma = compute_pvalue_from_sign(5)

line_1sigma = ROOT.TLine(gr.GetXaxis().GetXmin(), pval_1sigma, gr.GetXaxis().GetXmax(), pval_1sigma)
line_2sigma = ROOT.TLine(gr.GetXaxis().GetXmin(), pval_2sigma, gr.GetXaxis().GetXmax(), pval_2sigma)
line_3sigma = ROOT.TLine(gr.GetXaxis().GetXmin(), pval_3sigma, gr.GetXaxis().GetXmax(), pval_3sigma)
line_4sigma = ROOT.TLine(gr.GetXaxis().GetXmin(), pval_4sigma, gr.GetXaxis().GetXmax(), pval_4sigma)
line_5sigma = ROOT.TLine(gr.GetXaxis().GetXmin(), pval_5sigma, gr.GetXaxis().GetXmax(), pval_5sigma)

## dashed lines
line_1sigma.SetLineStyle(2)
line_2sigma.SetLineStyle(2)
line_3sigma.SetLineStyle(2)
line_4sigma.SetLineStyle(2)
line_5sigma.SetLineStyle(2)
#save to file
c = ROOT.TCanvas("c", "c", 800, 600)
gr.Draw()
line_1sigma.Draw("same")
line_2sigma.Draw("same")
line_3sigma.Draw("same")
line_4sigma.Draw("same")

## write nsigma on the left side of the canvas
text_1sigma = ROOT.TLatex(gr.GetXaxis().GetXmin() + 0.003, pval_1sigma + 0.08, "1#sigma")
text_2sigma = ROOT.TLatex(gr.GetXaxis().GetXmin() + 0.003, pval_2sigma + 0.15*1e-1, "2#sigma")
text_3sigma = ROOT.TLatex(gr.GetXaxis().GetXmin() + 0.003, pval_3sigma + 0.066*1e-2, "3#sigma")
text_4sigma = ROOT.TLatex(gr.GetXaxis().GetXmin() + 0.003, pval_4sigma + 0.02*1e-3, "4#sigma")
text_5sigma = ROOT.TLatex(gr.GetXaxis().GetXmin() + 0.003, pval_5sigma + 0.05, "5#sigma")

text_1sigma.SetTextAlign(12)
text_2sigma.SetTextAlign(12)
text_3sigma.SetTextAlign(12)
text_4sigma.SetTextAlign(12)
text_5sigma.SetTextAlign(12)

text_1sigma.Draw("same")
text_2sigma.Draw("same")
text_3sigma.Draw("same")
text_4sigma.Draw("same")


## set log scale for y axis
c.SetLogy()

outf = ROOT.TFile("pvalue_graph.root", "recreate")
c.Write()
gr.Write()
outf.Close()



