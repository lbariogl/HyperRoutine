import ROOT
import numpy as np

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC  = ROOT.TColor.GetColor("#ff7f00")

def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)


def fill_th1_hist_abs(h, df, var):
    for var_val in df[var]:
        h.Fill(abs(var_val))


def fill_th2_hist(h, df, var1, var2):
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(var1_val, var2_val)


def fill_th2_hist_abs(h, df, var1, var2):
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(abs(var1_val), var2_val)


def computeEfficiency(gen_hist, rec_hist, name, rebin=0):
    if rebin > 1:
        gen_hist.Rebin(rebin)
        rec_hist.Rebin(rebin)
    eff_hist = gen_hist.Clone(name)
    eff_hist.Reset()
    eff_hist.GetYaxis().SetTitle(r'#epsilon #times Acc')
    eff_hist.GetYaxis().SetRangeUser(0., 1.1)
    for iPt in range(1, rec_hist.GetNbinsX() + 1):
        gen_val = gen_hist.GetBinContent(iPt)
        if gen_val < 1e-24:
            continue
        rec_val = rec_hist.GetBinContent(iPt)
        eff_val = rec_val / gen_val
        eff_err = np.sqrt(eff_val * (1 - eff_val) / gen_val)
        # print('iPt: ', iPt, ' eff: ', eff_val, ' +- ', eff_err)
        eff_hist.SetBinContent(iPt, eff_val)
        eff_hist.SetBinError(iPt, eff_err)
    return eff_hist


def setHistStyle(hist, colour, marker=20, fillstyle=0, linewidth=1):
    hist.SetMarkerColor(colour)
    hist.SetLineColor(colour)
    hist.SetMarkerStyle(marker)
    hist.SetFillStyle(fillstyle)
    hist.SetLineWidth(linewidth)


def fill_th1_hist(h, df, var):
    for var_val in df[var]:
        h.Fill(var_val)

def fill_th2_hist(h, df, var1, var2):
    for i in range(df.shape[0]):
        h.Fill(df[var1].iloc[i], df[var2].iloc[i])


def fill_res_hist(h, df, var1, var2):
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill((var_val1 - var_val2)/var_val1)

def fill_res_hist_th2(h, df, var1, var2):
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill(var_val1,(var_val1 - var_val2)/var_val1)


def fill_mass_weighted_hist(h, df, var, weight = [1,1]):

    for var_val, w in zip(df[var], df['isSignal']):
        if w == 1:
            h.Fill(var_val, weight[0])
        else:
            h.Fill(var_val, weight[1])

def significance_error(signal, background, signal_error, background_error):

    sb = signal + background + 1e-10
    sb_sqrt = np.sqrt(sb)

    s_propag = (sb_sqrt + signal / (2 * sb_sqrt))/sb * signal_error
    b_propag = signal / (2 * sb_sqrt)/sb * background_error

    if signal+background == 0:
        return 0

    return np.sqrt(s_propag * s_propag + b_propag * b_propag)

def scale_hist_content(h, scale):
    ## generate poissonian counts
    for i in range(1, h.GetNbinsX()+1):
        pois = ROOT.gRandom.Poisson(scale)
        pois_sqrt = np.sqrt(pois)
        h.SetBinContent(i, h.GetBinContent(i)+pois)
        h.SetBinError(i, np.sqrt(pois_sqrt*pois_sqrt + h.GetBinError(i)*h.GetBinError(i)))


def set_style():
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptDate(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetLabelSize(0.04,"xyz")
    ROOT.gStyle.SetTitleSize(0.05,"xyz")
    ROOT.gStyle.SetTitleFont(42,"xyz")
    ROOT.gStyle.SetLabelFont(42,"xyz")
    ROOT.gStyle.SetTitleOffset(1.05,"x")
    ROOT.gStyle.SetTitleOffset(1.1,"y")
    ROOT.gStyle.SetCanvasDefW(800)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadGridX(0)
    ROOT.gStyle.SetPadGridY(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetPaperSize(20,24)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetEndErrorSize(0.)
    ROOT.gStyle.SetMarkerSize(1)




def fit_and_plot(dataset, var, fit_function, signal, background, sigma, mu, f, n_ev=300, matter_type="both", bdt_eff=None):

    fit_function.fitTo(dataset, ROOT.RooFit.Extended(False), ROOT.RooFit.Save(True))
    frame = var.frame(30)
    frame.SetName(f'data_extr_{bdt_eff}')
    frame.SetTitle('')
    set_style()
    frame.GetYaxis().SetTitleSize(0.06)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetYaxis().SetMaxDigits(2)
    frame.GetXaxis().SetTitleOffset(1.1)

    dataset.plotOn(frame, ROOT.RooFit.Name('data'))
    fit_function.plotOn(frame, ROOT.RooFit.LineColor(kBlueC), ROOT.RooFit.Name('fit_func'))
    fit_function.plotOn(frame, ROOT.RooFit.Components('bkg'), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(kOrangeC))
    # fit_function.plotOn(frame, ROOT.RooFit.Components('cb'), ROOT.RooFit.LineStyle(ROOT.kDashed))

    sigma_val = sigma.getVal()
    mu_val = mu.getVal()

    signal_counts = f.getVal()*dataset.sumEntries()
    signal_counts_error = (f.getError()/f.getVal())*f.getVal()*dataset.sumEntries()


    background_counts = (1-f.getVal())*dataset.sumEntries()
    background_counts_error = (1-f.getVal())*dataset.sumEntries()*f.getError()/f.getVal()

    #signal within 3 sigma
    var.setRange('signal', mu_val-3*sigma_val, mu_val+3*sigma_val)
    signal_int = signal.createIntegral(ROOT.RooArgSet(var), ROOT.RooArgSet(var), 'signal')
    signal_int_val_3s = signal_int.getVal()*signal_counts
    signal_int_val_3s_error = signal_int_val_3s*signal_counts_error/signal_counts
    #background within 3 sigma
    var.setRange('bkg', mu_val-3*sigma_val, mu_val+3*sigma_val)
    bkg_int = background.createIntegral(ROOT.RooArgSet(var), ROOT.RooArgSet(var), 'bkg')
    bkg_int_val_3s = bkg_int.getVal()*background_counts
    bkg_int_val_3s_error = bkg_int_val_3s*background_counts_error/background_counts
    significance = signal_int_val_3s/np.sqrt(signal_int_val_3s + bkg_int_val_3s)
    significance_err = significance_error(signal_int_val_3s, bkg_int_val_3s, signal_int_val_3s_error, bkg_int_val_3s_error)
    s_b_ratio_err = np.sqrt((signal_int_val_3s_error/signal_int_val_3s)**2 + (bkg_int_val_3s_error/bkg_int_val_3s)**2)*signal_int_val_3s/bkg_int_val_3s


    chi2 = frame.chiSquare("fit_func", "data", 6)
    fit_probability = ROOT.TMath.Prob(chi2*(frame.GetNbinsX() -6), frame.GetNbinsX() -6)

    print('chi2: ', chi2)
    print('fit probability: ', fit_probability)
    print('muv: ', mu_val)


    pinfo = ROOT.TPaveText(0.632, 0.5, 0.932, 0.85, 'NDC')
    pinfo.SetBorderSize(0)
    pinfo.SetFillStyle(0)
    pinfo.SetTextAlign(11)
    pinfo.SetTextFont(42)
    string_list = []

    # string_list.append('Fit Probability: ' + f'{fit_probability:.2f}')
    string_list.append(f'Signal (S): {signal_counts:.0f} #pm {signal_counts_error:.0f}')
    string_list.append(f'S/B (3 #sigma): {signal_int_val_3s/bkg_int_val_3s:.1f} #pm {s_b_ratio_err:.1f}')
    string_list.append('S/#sqrt{S+B} (3 #sigma): ' + f'{significance:.1f} #pm {significance_err:.1f}')
    string_list.append('#mu = ' + f'{mu_val*1e3:.2f} #pm {mu.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
    string_list.append('#sigma = ' + f'{sigma_val*1e3:.2f} #pm {sigma.getError()*1e3:.2f}' + ' MeV/#it{c}^{2}')
    for s in string_list:
        pinfo.AddText(s)

    string_list = []
    pinfo2 = ROOT.TPaveText(0.17, 0.6, 0.45, 0.85, "NDC")
    pinfo2.SetBorderSize(0)
    pinfo2.SetFillStyle(0)
    pinfo2.SetTextAlign(11)
    pinfo2.SetTextFont(42)

    if matter_type == "matter":
        matter_string = "{}^{3}_{#Lambda}H #rightarrow ^{3}He+#pi^{-}"

    elif matter_type == "antimatter":
        matter_string = "{}^{3}_{#bar{#Lambda}}#bar{H} #rightarrow ^{3}#bar{He}+#pi^{+}"

    else:
        matter_string = "{}^{3}_{#Lambda}H #rightarrow ^{3}He+#pi^{-} + c.c."

    string_list.append("ALICE Performance")
    string_list.append("Run 3, pp #sqrt{#it{s}} = 13.6 TeV")
    string_list.append("N_{ev} = " f"{n_ev:.0f} "  "#times 10^{9}")
    string_list.append(matter_string)

    if bdt_eff != None:
        string_list.append(f"BDT Efficiency: {bdt_eff:.2f}")

    for s in string_list:
        pinfo2.AddText(s)


    frame.addObject(pinfo)
    frame.addObject(pinfo2)

    fit_stats = {"signal": [signal_counts, signal_counts_error],
    "significance": [significance, significance_err], "s_b_ratio": [signal_int_val_3s/bkg_int_val_3s, s_b_ratio_err]}

    return frame




def ndarray2roo(ndarray, var, name='data'):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'

    name = var.GetName()
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{name}', x ,f'{name}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()

    array_roo = ROOT.RooDataSet(name, 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo
