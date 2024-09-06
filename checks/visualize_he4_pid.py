import ROOT
import numpy as np
from hipe4ml.tree_handler import TreeHandler
import sys
sys.path.append('../utils')
import utils as utils


# template <typename T>
# GPUdi() T BetheBlochAleph(T bg, T kp1, T kp2, T kp3, T kp4, T kp5)
# {
#   T beta = bg / o2::gpu::GPUCommonMath::Sqrt(static_cast<T>(1.) + bg * bg);

#   T aa = o2::gpu::GPUCommonMath::Pow(beta, kp4);
#   T bb = o2::gpu::GPUCommonMath::Pow(static_cast<T>(1.) / bg, kp5);
#   bb = o2::gpu::GPUCommonMath::Log(kp3 + bb);

#   return (kp2 - aa - bb) * kp1 / aa;
# }


def heBB(rigidity, mass):
    p1 = -321.34
    p2 = 0.6539
    p3 = 1.591
    p4 = 0.8225
    p5 = 2.363

    betagamma = rigidity * 2 / mass
    beta = betagamma / np.sqrt(1 + betagamma**2)
    aa = beta**p4
    bb = np.log(p3 + (1 / betagamma)**p5)
    return (p2 - aa - bb) * p1 / aa

tree_list = [
    "/data3/fmazzasc/hyp_run_3/pp2024/ag/AO2D_custom.root",
    "/data3/fmazzasc/hyp_run_3/pp2024/aj/AO2D_custom.root",
    "/data3/fmazzasc/hyp_run_3/pp2024/af/AO2D_custom.root",
    "/data3/fmazzasc/hyp_run_3/pp2024/al/AO2D_custom.root"
  ]    

hdl = TreeHandler(tree_list, 'O2hypcands', folder_name='DF*')
utils.correct_and_convert_df(hdl, True, False, True)

hdl.apply_preselections('fAvgClusterSizeHe > 6 and fNSigmaHe > -3 and fNTPCclusHe > 90 and fIsMatter==False')

h2TPCSigClusSize = ROOT.TH2F('h2TPCSigMomHeTPC', r'; p_{TPC}; TPC signal', 50, 0.5, 5, 200, 0.5, 2000)
utils.fill_th2_hist(h2TPCSigClusSize, hdl, 'fTPCmomHe', 'fTPCsignalHe')

rigidity = np.linspace(0.5, 5, 100)
bbHe3 = heBB(rigidity, 2.809)
bbHe4 = heBB(rigidity, 3.727)

## create 1sigma bands, resolution: 8%
bbHe4_low = bbHe4 - 2*0.08 * bbHe4
bbHe4_high = bbHe4 + 2*0.08 * bbHe4
bbHe3_low = bbHe3 - 2*0.08 * bbHe3
bbHe3_high = bbHe3 + 2*0.08 * bbHe3


print(bbHe3)

## create 2 tgraphs for He3 and He4
grHe3 = ROOT.TGraph(len(rigidity), rigidity, bbHe3)
grHe4 = ROOT.TGraph(len(rigidity), rigidity, bbHe4)

grHe3_low = ROOT.TGraph(len(rigidity), rigidity, bbHe3_low)
grHe4_low = ROOT.TGraph(len(rigidity), rigidity, bbHe4_low)

grHe3_high = ROOT.TGraph(len(rigidity), rigidity, bbHe3_high)
grHe4_high = ROOT.TGraph(len(rigidity), rigidity, bbHe4_high)


##Plot th2 and tgraphs in the same canvas
c = ROOT.TCanvas('c', 'c', 800, 600)
h2TPCSigClusSize.Draw('colz')
grHe3.SetLineColor(ROOT.kRed)
grHe3_low.SetLineColor(ROOT.kRed)
grHe3_high.SetLineColor(ROOT.kRed)

grHe4.SetLineColor(ROOT.kBlue)
grHe4_low.SetLineColor(ROOT.kBlue)
grHe4_high.SetLineColor(ROOT.kBlue)

grHe3.Draw('same')
grHe3_low.Draw('same')
grHe3_high.Draw('same')
grHe4.Draw('same')
grHe4_low.Draw('same')
grHe4_high.Draw('same')


outfile = ROOT.TFile('he4_pid_calibration.root', 'recreate')
h2TPCSigClusSize.Write()
c.Write()
outfile.Close()


