import ROOT


## static const std::vector<int> pdgCodes{211, 321, 2212, 1000010020, 1000020030, 310, 3122, 1010010030};

input_file = ROOT.TFile("AnalysisResults.root", "READ")

mc_gen_all = input_file.Get("mcsignalloss/mcGenAll")
mc_gen_reco_coll = input_file.Get("mcsignalloss/mcGenAllRecoColl")
mc_gen_v0s = input_file.Get("mcsignalloss/mcGenV0s")
mc_gen_reco_v0s = input_file.Get("mcsignalloss/mcGenV0sRecoColl")

mc_gen_coll = input_file.Get("mcsignalloss/mcGenCollisionVtx")
mc_reco_coll = input_file.Get("mcsignalloss/recCollisionVtx")
## divide for the number of entries of the two histograms
event_loss = mc_reco_coll.GetEntries() / mc_gen_coll.GetEntries()
print("Event loss: ", event_loss)



## divide the two histograms
signal_loss_all = mc_gen_reco_coll.Clone("signal_loss_all")
signal_loss_all.Divide(mc_gen_all)
## put the pdg codes in the y axis, and the pT in the x axis
signal_loss_all.GetYaxis().SetBinLabel(1, "#pi^{#pm}")
signal_loss_all.GetYaxis().SetBinLabel(2, "K^{#pm}")
signal_loss_all.GetYaxis().SetBinLabel(3, "p")
signal_loss_all.GetYaxis().SetBinLabel(4, "d")
signal_loss_all.GetYaxis().SetBinLabel(5, "^{3}He")
signal_loss_all.GetYaxis().SetBinLabel(6, "K^{0}_{S}")
signal_loss_all.GetYaxis().SetBinLabel(7, "#Lambda")
signal_loss_all.GetYaxis().SetBinLabel(8, "^{3}_{#Lambda}H")
signal_loss_all.GetXaxis().SetTitle("p_{T} [GeV/c]")
signal_loss_all.GetZaxis().SetTitle("Signal loss")


## get lambdas (second bin of the z axis), project to xy plane
mc_gen_v0s.GetZaxis().SetRange(7, 7)
mc_gen_v0s_xy = mc_gen_v0s.Project3D("xy")
mc_gen_reco_v0s.GetZaxis().SetRange(7, 7)
mc_gen_reco_v0s_xy = mc_gen_reco_v0s.Project3D("xy")
## divide the two histograms
signal_loss_lambda = mc_gen_reco_v0s_xy.Clone("signal_loss_lambda")
signal_loss_lambda.Divide(mc_gen_v0s_xy)

### get only the x axis for lambdas
mc_gen_v0s_x = mc_gen_v0s.Project3D("x")
mc_gen_reco_v0s_x = mc_gen_reco_v0s.Project3D("x")
signal_loss_lambda_pt = mc_gen_reco_v0s_x.Clone("signal_loss_lambda_pt")
signal_loss_lambda_pt.Divide(mc_gen_v0s_x)

mc_gen_v0s_y = mc_gen_v0s.Project3D("y")
mc_gen_reco_v0s_y = mc_gen_reco_v0s.Project3D("y")
signal_loss_lambda_dL = mc_gen_reco_v0s_y.Clone("signal_loss_lambda_dL")
signal_loss_lambda_dL.Divide(mc_gen_v0s_y)


## get signal loss for protons
signal_loss_pion_pt = signal_loss_all.ProjectionX("signal_loss_pion", 1, 1)
signal_loss_kaon_pt = signal_loss_all.ProjectionX("signal_loss_kaon", 2, 2)
signal_loss_proton_pt = signal_loss_all.ProjectionX("signal_loss_proton", 3, 3)
signal_loss_deuteron_pt = signal_loss_all.ProjectionX("signal_loss_deuteron", 4, 4)
signal_loss_helium3_pt = signal_loss_all.ProjectionX("signal_loss_helium3", 5, 5)


### get ratio between lambda and proton signal losses
signal_loss_lambda_proton_pt = signal_loss_lambda_pt.Clone("signal_loss_lambda_proton_pt")
signal_loss_lambda_proton_pt.Divide(signal_loss_proton_pt)
signal_loss_lambda_proton_pt.SetTitle("; p_{T} [GeV/c]; Signal loss ratio #Lambda/p")

## get lambda signal loss divided for event loss
sigloss_evloss = signal_loss_lambda_pt.Clone("sigloss_evloss")
sigloss_evloss.Scale(1/event_loss)
sigloss_evloss.SetDrawOption("hist")
sigloss_evloss.GetYaxis().SetTitle("Signal loss / event loss")
sigloss_evloss.GetXaxis().SetTitle("p_{T} [GeV/c]")


outfile = ROOT.TFile("signal_loss_23k2f_rofcut.root", "RECREATE")
signal_loss_all.Write()
signal_loss_lambda.Write()
signal_loss_lambda_pt.Write()
signal_loss_lambda_dL.Write()

signal_loss_pion_pt.Write()
signal_loss_kaon_pt.Write()
signal_loss_proton_pt.Write()
signal_loss_deuteron_pt.Write()
signal_loss_helium3_pt.Write()
signal_loss_lambda_proton_pt.Write()

sigloss_evloss.Write()

outfile.Close()
