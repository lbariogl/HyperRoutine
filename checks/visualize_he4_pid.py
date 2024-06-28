import ROOT
import numpy as np
from hipe4ml.tree_handler import TreeHandler
import sys
sys.path.append('utils')
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
    

hdl = TreeHandler("/data3/fmazzasc/hyp_run_3/pbpb/pass3/AO2D.root", 'O2hypcandsflow', folder_name='DF*')
utils.correct_and_convert_df(hdl, True, False, True)

hdl.apply_preselections('fAvgClusterSizeHe > 4.5 and fNSigmaHe > -2 and fNTPCclusHe > 90')

h2TPCSigClusSize = ROOT.TH2F('h2TPCSigMomHeTPC', r'; p_{TPC}; TPC signal', 50, 0.5, 5, 200, 0.5, 2000)
utils.fill_th2_hist(h2TPCSigClusSize, hdl, 'fTPCmomHe', 'fTPCsignalHe')

rigidity = np.linspace(0.5, 5, 100)
bbHe3 = heBB(rigidity, 2.809)
bbHe4 = heBB(rigidity, 3.727)

print(bbHe3)

## create 2 tgraphs for He3 and He4
grHe3 = ROOT.TGraph(len(rigidity), rigidity, bbHe3)
grHe4 = ROOT.TGraph(len(rigidity), rigidity, bbHe4)


##Plot th2 and tgraphs in the same canvas
c = ROOT.TCanvas('c', 'c', 800, 600)
h2TPCSigClusSize.Draw('colz')
grHe3.SetLineColor(ROOT.kRed)
grHe4.SetLineColor(ROOT.kBlue)
grHe3.Draw('same')
grHe4.Draw('same')


outfile = ROOT.TFile('he4_pid_calibration.root', 'recreate')
h2TPCSigClusSize.Write()
c.Write()
outfile.Close()


