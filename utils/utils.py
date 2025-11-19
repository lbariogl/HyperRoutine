import ROOT
import numpy as np
import pandas as pd

kBlueC = ROOT.TColor.GetColor('#1f78b4')
kOrangeC = ROOT.TColor.GetColor('#ff7f00')

## set numpy seed
np.random.seed(1995)

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
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val in df[var]:
        h.Fill(var_val)


def fill_th1_hist_abs(h, df, var):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val in df[var]:
        h.Fill(abs(var_val))


def fill_th2_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(var1_val, var2_val)


def fill_th2_hist_abs(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var1_val, var2_val in zip(df[var1], df[var2]):
        h.Fill(abs(var1_val), var2_val)


def fill_res_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill((var_val1 - var_val2)/var_val1)


def fill_th2_res_hist(h, df, var1, var2):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
    for var_val1, var_val2 in zip(df[var1], df[var2]):
        h.Fill(var_val1, (var_val2 - var_val1)/var_val1)

def fill_mass_weighted_hist(h, df, var, weight=[1, 1]):
    if not type(df) == pd.DataFrame:
        df = df._full_data_frame
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
    array_roo = ROOT.RooDataSet(ROOT.RooStringView(name), ROOT.RooStringView('dataset from tree'), ROOT.RooArgSet(var), ROOT.RooFit.Import(tree))
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
    h2MomResoVsPtHe3 = ROOT.TH2F('h2MomResoVsPtHe3', ';{}^{3}He #it{p}_{T} (GeV/#it{c});{}^{3}He #it{p}_{T}^{gen} - #it{p}_{T}^{reco} (GeV/#it{c})', 30, 1.3, 5, 50, -0.4, 0.4)
    df.eval('PtResHe3 = (fGenPtHe3 - fPtHe3)', inplace=True)
    fill_th2_hist(h2MomResoVsPtHe3, df, 'fPtHe3', 'PtResHe3')
    h2MomResoVsPtHe3.FitSlicesY()
    hShiftVsPtHe3 = ROOT.gDirectory.Get('h2MomResoVsPtHe3_1')
    hShiftVsPtHe3.SetName('hShiftVsPtHe3')
    return h2MomResoVsPtHe3, hShiftVsPtHe3

def heBB(rigidity, mass, isMC=False):

    if isMC:
        p1 = -611.592275784047
        p2 = 1.7160058505736953
        p3 = 3.2070341629338754
        p4 = -0.7453309037015425
        p5 = 4.107037208529488

    else:
        p1 = -321.34
        p2 = 0.6539
        p3 = 1.591
        p4 = 0.8225
        p5 = 2.363

    # p1 = -241.4902
    # p2 = 0.374245
    # p3 = 1.397847
    # p4 = 1.0782504
    # p5 = 2.048336

    betagamma = rigidity * 2 / mass
    beta = betagamma / np.sqrt(1 + betagamma**2)
    aa = beta**p4
    bb = np.log(p3 + (1 / betagamma)**p5)
    return (p2 - aa - bb) * p1 / aa

def computeNSigmaHe4(df, isMC=False):
    expBB = heBB(df['fTPCmomHe'], 3.727, isMC)
    nSigma = (df['fTPCsignalHe'] - expBB) / (0.09*df['fTPCsignalHe'])
    return nSigma

def computeNSigmaHe3(df, isMC=False):
    expBB = heBB(df['fTPCmomHe'], 2.80839160743, isMC)
    nSigma = (df['fTPCsignalHe'] - expBB) / (0.09*df['fTPCsignalHe'])
    return nSigma


def get_bin(h, label):
    for ibin in range(1, h.GetNbinsX() + 1):
        if h.GetXaxis().GetBinLabel(ibin) == label:
            return ibin
    return -1


def getNEvents(an_files, is_trigger=False):
    n_ev = 0
    if isinstance(an_files, str):
        an_files = [an_files]

    for an_file in an_files:
        ff = ROOT.TFile(an_file)
        if not ff or ff.IsZombie():
            print(f"Could not open file: {an_file}")
            continue
        if is_trigger:
            # === Replacement logic starts here ===
            ff.cd('hyper-reco-task/Zorro')
            total_sum = 0.0
            for key in ROOT.gDirectory.GetListOfKeys():
                key_name = key.GetName()
                # Only process integer-named directories
                if not key_name.isdigit():
                    continue
                # Get histograms
                h_selections = ff.Get(f"hyper-reco-task/Zorro/{key_name}/Selections")
                h_insp_tvx  = ff.Get(f"hyper-reco-task/Zorro/{key_name}/InspectedTVX")
                h_anal_trig = ff.Get(f"hyper-reco-task/Zorro/{key_name}/AnalysedTriggersOfInterest")

                if not (h_selections and h_insp_tvx and h_anal_trig):
                    print(f"⚠️ Missing histograms in {key_name}, skipping.")
                    continue
                ibin_fHe = get_bin(h_selections, "fHe")
                if ibin_fHe == -1:
                    print(f"'fHe' bin not found in {key_name}, skipping.")
                    continue
                n_sel = h_selections.GetBinContent(ibin_fHe)
                n_trig = h_anal_trig.GetBinContent(get_bin(h_anal_trig, "fHe"))
                n_insp = h_insp_tvx.Integral()
                if n_sel == 0:
                    print(f"Zero selections in {key_name}, skipping to avoid division by zero.")
                    continue

                ntvx = n_insp * n_trig / n_sel
                total_sum += ntvx
            n_ev += total_sum
            # === Replacement logic ends here ===
        else:
            # Default (non-trigger) behavior
            h_zvtx = ff.Get('hyper-reco-task/hZvtx')
            if h_zvtx:
                n_ev += h_zvtx.Integral()
            else:
                print(f"⚠️ Missing hZvtx in {an_file}")

        ff.Close()

    return n_ev




def correct_and_convert_df(df, calibrate_he3_pt = False, isMC=False, isH4L=False):

    kDefaultPID = 15
    kPionPID = 2
    kTritonPID = 6

    if not type(df) == pd.DataFrame:
        df = df._full_data_frame

    if 'fFlags' in df.columns:
        df['fHePIDHypo'] = np.right_shift(df['fFlags'], 4)
        df['fPiPIDHypo'] = np.bitwise_and(df['fFlags'], 0b1111)
    
    if not 'fTPCChi2He' in df.columns or isMC: ## chi2 for MC is not reliable
        ## set dummy column to one
        df['fTPCChi2He'] = 1
    
    ## set Nan chi2 values to one
    df['fTPCChi2He'] = df['fTPCChi2He'].replace(np.nan, 1)

    # correct 3He momentum    

    if calibrate_he3_pt:
        # print(df.query('fIsReco==True')['fHePIDHypo'])
        no_pid_mask = np.logical_and(df['fHePIDHypo'] != kDefaultPID, df['fHePIDHypo'] != kPionPID)

        if (no_pid_mask.sum() == 0):
            print("PID in tracking not detected, using old momentum re-calibration")
            df.loc[:, "fPtHe3"] = df["fPtHe3"] + 2.98019e-02 + 7.66100e-01 * np.exp(-1.31641e+00 * df["fPtHe3"]) ### functional form given by mpuccio
        else:
            print("PID in tracking detected, using new momentum re-calibration")
            df_Trit_PID = df.query('fHePIDHypo==6')
            df_else = df.query('fHePIDHypo!=6')
            ##pt_new = pt + kp0 + kp1 * pt + kp2 * pt^2 curveParams = {'kp0': -0.200281,'kp1': 0.103039,'kp2': -0.012325}, functional form given by G.A. Lucia
            df_Trit_PID.loc[:, "fPtHe3"] = df_Trit_PID["fPtHe3"] + -0.1286 - 0.1269 * df_Trit_PID["fPtHe3"] + 0.06 * df_Trit_PID["fPtHe3"]**2
            df_new = pd.concat([df_Trit_PID, df_else])
            
            ## assign the new dataframe to the original one
            df[:] = df_new
    
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
    df.eval('fCosLambda = fPt/fP', inplace=True)
    df.eval('fCosLambdaHe = fPtHe3/fPHe3', inplace=True)
    df['fNSigmaHe4'] = computeNSigmaHe4(df, isMC)
    df['fNSigmaHe3'] = computeNSigmaHe3(df, isMC)

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

    ## signed TPC mom
    df.eval('fTPCSignMomHe3 = fTPCmomHe * (-1 + 2*fIsMatter)', inplace=True)
    df.eval('fGloSignMomHe3 = fPHe3 / 2 * (-1 + 2*fIsMatter)', inplace=True)

    if "fITSclusterSizesHe" in df.columns:
    ## loop over the candidates and compute the average cluster size
        clSizesHe = df['fITSclusterSizesHe'].to_numpy()
        clSizesPi = df['fITSclusterSizesPi'].to_numpy()
        clSizeHeAvg = np.zeros(len(clSizesHe))
        clSizePiAvg = np.zeros(len(clSizesPi))
        nHitsHe = np.zeros(len(clSizesHe))
        nHitsPi = np.zeros(len(clSizesPi))
        for iLayer in range(7):
            clSizeHeAvg += np.right_shift(clSizesHe, 4*iLayer) & 0b1111
            clSizePiAvg += np.right_shift(clSizesPi, 4*iLayer) & 0b1111
            nHitsHe += np.right_shift(clSizesHe, 4*iLayer) & 0b1111 > 0
            nHitsPi += np.right_shift(clSizesPi, 4*iLayer) & 0b1111 > 0

        nHitsHe += 1e-10
        nHitsPi += 1e-10
        clSizeHeAvg /= nHitsHe
        clSizePiAvg /= nHitsPi
        df['fAvgClusterSizeHe'] = clSizeHeAvg
        df['fAvgClusterSizePi'] = clSizePiAvg
        df['nITSHitsHe'] = nHitsHe
        df['nITSHitsPi'] = nHitsPi
        df.eval('fAvgClSizeCosLambda = fAvgClusterSizeHe * fCosLambdaHe', inplace=True)

    if "fPsiFT0C" in df.columns:
        df.eval('fPhi = arctan2(fPy, fPx)', inplace=True)
        df.eval('fV2 = cos(2*(fPhi - fPsiFT0C))', inplace=True)


    if isMC:
        df.eval('fGenDecLen = sqrt(fGenXDecVtx**2 + fGenYDecVtx**2 + fGenZDecVtx**2)', inplace=True)
        df.eval('fGenDecRad = sqrt(fGenXDecVtx**2 + fGenYDecVtx**2)', inplace=True)
        df.eval('fGenPz = fGenPt * sinh(fGenEta)', inplace=True)
        df.eval('fGenP = sqrt(fGenPt**2 + fGenPz**2)', inplace=True)
        df.eval("fAbsGenPt = abs(fGenPt)", inplace=True)

        if not isH4L:
            df.eval('fGenCt = fGenDecLen * 2.99131 / fGenP', inplace=True)
        else:
            df.eval('fGenCt = fGenDecLen * 3.922 / fGenP', inplace=True)
    
        if 'fIsTwoBodyDecay' in df.columns:
            df.query('fIsTwoBodyDecay==True', inplace=True)


    # remove useless columns
    df.drop(columns=['fPxHe3', 'fPyHe3', 'fPzHe3', 'fEnHe3', 'fPxPi', 'fPyPi', 'fPzPi', 'fPPi', 'fEnPi', 'fPx', 'fPy', 'fPz', 'fP', 'fEn'])


def compute_pvalue_from_sign(significance):
    return ROOT.Math.chisquared_cdf_c(significance**2, 1) / 2

def convert_sel_to_string(selection):
    sel_string = ''
    conj = ' and '
    for _, val in selection.items():
        sel_string = sel_string + val + conj
    return sel_string[:-len(conj)]

def saveCanvasAsPDF(histo, plots_dir, is2D=False):
    histo_name = histo.GetName()
    canvas_name = histo_name.replace('h', 'c', 1)
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 800, 600)
    canvas.SetBottomMargin(0.13)
    canvas.SetLeftMargin(0.13)
    if not is2D:
        histo.Draw('histo')
    else:
        histo.Draw('colz')
    canvas.SaveAs(f'{plots_dir}/{canvas_name}.pdf')
