import ROOT
import numpy as np

def computeEfficiency(gen_hist, rec_hist, name):
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
