import ROOT
import numpy as np
import pandas as pd

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC = ROOT.TColor.GetColor('#ff7f00')

## set numpy seed
np.random.seed(42)


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
        h.Fill(var_val1, (var_val1 - var_val2)/var_val1)


def fill_mass_weighted_hist(h, df, var, weight=[1, 1]):

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
    # generate poissonian counts
    for i in range(1, h.GetNbinsX()+1):
        pois = ROOT.gRandom.Poisson(scale)
        pois_sqrt = np.sqrt(pois)
        h.SetBinContent(i, h.GetBinContent(i)+pois)
        h.SetBinError(i, np.sqrt(pois_sqrt*pois_sqrt +
                      h.GetBinError(i)*h.GetBinError(i)))


def set_style():
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptDate(0)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetLabelSize(0.04, 'xyz')
    ROOT.gStyle.SetTitleSize(0.05, 'xyz')
    ROOT.gStyle.SetTitleFont(42, 'xyz')
    ROOT.gStyle.SetLabelFont(42, 'xyz')
    ROOT.gStyle.SetTitleOffset(1.05, 'x')
    ROOT.gStyle.SetTitleOffset(1.1, 'y')
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
    ROOT.gStyle.SetPaperSize(20, 24)
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetEndErrorSize(0.)
    ROOT.gStyle.SetMarkerSize(1)


def ndarray2roo(ndarray, var, name='data'):
    if isinstance(ndarray, ROOT.RooDataSet):
        print('Already a RooDataSet')
        return ndarray

    assert isinstance(ndarray, np.ndarray), 'Did not receive NumPy array'
    assert len(ndarray.shape) == 1, 'Can only handle 1d array'
    x = np.zeros(1, dtype=np.float64)

    tree = ROOT.TTree('tree', 'tree')
    tree.Branch(f'{var.GetName()}', x, f'{var.GetName()}/D')

    for i in ndarray:
        x[0] = i
        tree.Fill()
    array_roo = ROOT.RooDataSet(name, 'dataset from tree', tree, ROOT.RooArgSet(var))
    return array_roo


### reweight a distribution with rejection sampling
def reweight_pt_spectrum(df, var, distribution):
    rej_flag = np.ones(len(df))
    random_arr = np.random.rand(len(df))
    max_bw = distribution.GetMaximum()

    for ind, (val, rand) in enumerate(zip(df[var],random_arr)):
        frac = distribution.Eval(val)/max_bw
        if rand > frac:
            rej_flag[ind] = -1
    ## check if it is a pandas dataframe
    if isinstance(df, pd.DataFrame):
        df['rej'] = rej_flag
        return
    df._full_data_frame['rej'] = rej_flag

# create histogram for momentum correction

def create_pt_shift_histo(df):
    h2MomResoVsPtHe3 = ROOT.TH2F('h2MomResoVsPtHe3', ';{}^{3}He #it{p}_{T} (GeV/#it{c});{}^{3}He #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c})', 50, 1., 5, 50, -0.4, 0.4)
    df.eval('PtResHe3 = (fPtHe3 - fGenPtHe3)', inplace=True)
    fill_th2_hist(h2MomResoVsPtHe3, df, 'fPtHe3', 'PtResHe3')
    h2MomResoVsPtHe3.FitSlicesY()
    hShiftVsPtHe3 = ROOT.gDirectory.Get('h2MomResoVsPtHe3_1')
    hShiftVsPtHe3.SetName('hShiftVsPtHe3')
    return h2MomResoVsPtHe3, hShiftVsPtHe3

# put dataframe in the correct format

def correct_and_convert_df(df, calibrate_he3_pt = False, isMC=False, isH4L=False):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    # correct 3He momentum
    if calibrate_he3_pt:
        df["fPtHe3"] += 2.98019e-02 + 7.66100e-01 * np.exp(-1.31641e+00 * df["fPtHe3"]) ### functional form given by mpuccio

    # 3He momentum
    df.eval('fPxHe3 = fPtHe3 * cos(fPhiHe3)', inplace=True)
    df.eval('fPyHe3 = fPtHe3 * sin(fPhiHe3)', inplace=True)
    df.eval('fPzHe3 = fPtHe3 * sinh(fEtaHe3)', inplace=True)
    df.eval('fPHe3 = fPtHe3 * cosh(fEtaHe3)', inplace=True)
    df.eval('fEnHe3 = sqrt(fPHe3**2 + 2.8083916**2)', inplace=True)
    df.eval('fEnHe4 = sqrt(fPHe3**2 + 3.7273794**2)', inplace=True)
    # pi momentum
    df.eval('fPxPi = fPtPi * cos(fPhiPi)', inplace=True)
    df.eval('fPyPi = fPtPi * sin(fPhiPi)', inplace=True)
    df.eval('fPzPi = fPtPi * sinh(fEtaPi)', inplace=True)
    df.eval('fPPi = fPtPi * cosh(fEtaPi)', inplace=True)
    df.eval('fEnPi = sqrt(fPPi**2 + 0.139570**2)', inplace=True)
    # hypertriton momentum
    df.eval('fPx = fPxHe3 + fPxPi', inplace=True)
    df.eval('fPy = fPyHe3 + fPyPi', inplace=True)
    df.eval('fPz = fPzHe3 + fPzPi', inplace=True)
    df.eval('fP = sqrt(fPx**2 + fPy**2 + fPz**2)', inplace=True)
    df.eval('fEn = fEnHe3 + fEnPi', inplace=True)
    df.eval('fEn4 = fEnHe4 + fEnPi', inplace=True)
    # Momentum variables to be stored
    df.eval('fPt = sqrt(fPx**2 + fPy**2)', inplace=True)
    df.eval('fEta = arccosh(fP/fPt)', inplace=True)
    df.eval('fPhi = arctan(fPy/fPx)', inplace=True)
    # Variables of interest
    df.eval('fDecLen = sqrt(fXDecVtx**2 + fYDecVtx**2 + fZDecVtx**2)', inplace=True)
    if not isH4L:
        df.eval('fCt = fDecLen * 2.99131 / fP', inplace=True)
    else:
        print('Using H4L decay length')
        df.eval('fCt = fDecLen * 3.922 / fP', inplace=True)

    df.eval('fDecRad = sqrt(fXDecVtx**2 + fYDecVtx**2)', inplace=True)
    df.eval('fCosPA = (fPx * fXDecVtx + fPy * fYDecVtx + fPz * fZDecVtx) / (fP * fDecLen)', inplace=True)
    df.eval('fMassH3L = sqrt(fEn**2 - fP**2)', inplace=True)
    df.eval('fMassH4L = sqrt(fEn4**2 - fP**2)', inplace=True)

    if isMC:
        df.eval('fGenDecLen = sqrt(fGenXDecVtx**2 + fGenYDecVtx**2 + fGenZDecVtx**2)', inplace=True)
        df.eval('fGenPz = fGenPt * sinh(fGenEta)', inplace=True)
        df.eval('fGenP = sqrt(fGenPt**2 + fGenPz**2)', inplace=True)
        df.eval("fAbsGenPt = abs(fGenPt)", inplace=True)

        if not isH4L:
            df.eval('fGenCt = fGenDecLen * 2.99131 / fGenP', inplace=True)
        else:
            df.eval('fGenCt = fGenDecLen * 3.922 / fGenP', inplace=True)
            

    # remove useless columns
    df.drop(columns=['fPxHe3', 'fPyHe3', 'fPzHe3', 'fPHe3', 'fEnHe3', 'fPxPi', 'fPyPi', 'fPzPi', 'fPPi', 'fEnPi', 'fPx', 'fPy', 'fPz', 'fP', 'fEn'])


def compute_pvalue_from_sign(significance):
    return ROOT.Math.chisquared_cdf_c(significance**2, 1) / 2

def convert_sel_to_string(selection):
    sel_string = ''
    conj = ' and '
    for _, val in selection.items():
        sel_string = sel_string + val + conj
    return sel_string[:-len(conj)]
