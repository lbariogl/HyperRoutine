#include <iostream>

#include <TDatabasePDG.h>
#include <TFile.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include "Pythia8/Pythia.h"

//__________________________________________________________________________________________________
void simulateK0sDecay(int nEvents = 30000000, int seed = 42, bool linearSmearing = true);
Pythia8::Particle GetParticle(int pdgCode = 310);

//__________________________________________________________________________________________________
void simulateK0sDecay(int nEvents, int seed, bool linearSmearing)
{
    //__________________________________________________________
    // create and configure pythia generator

    Pythia8::Pythia pythia;
    pythia.readString("SoftQCD:inelastic = on");

    // keep only interesting decay to charged pions
    pythia.readString("310:onMode = off");
    pythia.readString("310:onIfMatch = 211 211");

    // init
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed %d", seed));
    pythia.init();

    gRandom->SetSeed(seed);

    //__________________________________________________________
    // create output TH1D

    TH2D *hDauPtRes = new TH2D("hDauPtRes", ";#it{p}^{gen}_{T} (+) (GeV/#it{c}); (#it{p}^{gen}_{T} - #it{p}_{T}) / #it{p}^{gen}_{T} (+)", 100, 0, 5, 100, -0.2, 0.2);
    TH2D *hDauPtInvMass = new TH2D("hDauPtInvMass", ";#it{p}^{+}_{T} (GeV/#it{c}); (#it{m} (#pi^{+}#pi^{-}) - #it{m} (k^{0}_{s})  / #it{m} (k^{0}_{s})", 100, 0, 5, 100, -0.2, 0.2);

    TH2D *hDauPtDiffVsPt = new TH2D("hDauPtDiffVsPt", ";#it{p}^{gen}_{T} (k^{0}_{s}) (GeV/#it{c}); #it{p}_{T} (#pi^{+}) - #it{p}_{T} (#pi^{-}) (GeV/#it{c})", 100, 0, 5, 400, -3, 3);
    TH2D *hDauPtVsPt = new TH2D("hDauPtVsPt", ";#it{p}^{gen}_{T} (k^{0}_{s}) (GeV/#it{c}); #it{p}_{T} (#pi^{+}) (GeV/#it{c})", 100, 0, 5, 100, 0, 5);
    TH2D *hPtInvMass = new TH2D("hPtInvMass", ";#it{p}^{gen}_{T} (k^{0}_{s}) (GeV/#it{c}); #it{m} (#pi^{+}#pi^{-}) - #it{m} (K^{0}_{s}) (GeV/#it{c}^{2})", 100, 0, 10, 100, -0.1, 0.1);
    TH1D *hSigmaMassVsPt = new TH1D("hSigmaMassVsPt", ";#it{p}^{gen}_{T} (k^{0}_{s}) (GeV/#it{c}); #sigma_{#it{m}} (#pi^{+}#pi^{-}) (GeV/#it{c}^{2})", 100, 0, 10);

    TH2D *hDeltaPtCM = new TH2D("hDeltaPtCM", ";#it{p}^{+}_{T} (GeV/#it{c}); (#it{p}^{+}_{T} - #it{p}^{-}_{T}) / (#sqrt{2} #it{p}^{+}_{T}) [C.M.]", 100, 0, 5, 100, -0.2, 0.2);
    TH2D *hDeltaPtCMTrueBoost = new TH2D("hDeltaPtCMTrueBoost", ";#it{p}^{+}_{T} (GeV/#it{c});  (#it{p}^{+}_{T} - #it{p}^{-}_{T}) / (#sqrt{2} #it{p}^{+}_{T}) [C.M.]", 100, 0, 5, 100, -0.2, 0.2);

    //__________________________________________________________
    // perform the simulation

    for (auto iEvent{1}; iEvent < nEvents; ++iEvent)
    {
        // reset pythia event and put resonance only
        pythia.event.reset();
        pythia.event.append(GetParticle(310));
        int idPart = pythia.event[1].id();
        pythia.particleData.mayDecay(idPart, true);
        pythia.moreDecays();

        for (auto iPart{1}; iPart < pythia.event.size(); ++iPart)
        {

            int pdg = pythia.event[iPart].id();
            if (pdg != 310)
            {
                continue;
            }

            auto dauIds = pythia.event[iPart].daughterList();
            if (dauIds.size() != 2)
            { // something went wrong with forced decays
                continue;
            }

            auto k0s = pythia.event[iPart];
            auto piPos = pythia.event[dauIds[0]];
            auto piNeg = pythia.event[dauIds[1]];

            TLorentzVector k0sTrueVec, k0sRecoVec, piPosVec, piNegVec;
            k0sTrueVec.SetPtEtaPhiM(k0s.pT(), k0s.eta(), k0s.phi(), k0s.m());

            // smear daughter pTs and set TLorentzVectors.
            float ptResNeg = 0., ptResPos = 0.;
            if (linearSmearing)
            {
                // Linear pT smearing from 0.1 and 0.01 at 2.5 GeV/c and vice versa at 5 GeV/c
                ptResNeg = piNeg.pT() < 2.5 ? -0.09 / 2.5 * piNeg.pT() + 0.1 : 0.09 / 2.5 * piNeg.pT() - 0.08;
                ptResPos = piPos.pT() < 2.5 ? -0.09 / 2.5 * piPos.pT() + 0.1 : 0.09 / 2.5 * piPos.pT() - 0.08;
            }
            else
            {
                // Flat pT smearing with 3% resolution
                ptResNeg = 0.03;
                ptResPos = 0.03;
            }

            float ptPiNeg = piNeg.pT() * (1. + gRandom->Gaus(0., ptResNeg));
            float ptPiPos = piPos.pT() * (1. + gRandom->Gaus(0., ptResPos));
            piNegVec.SetPtEtaPhiM(ptPiNeg, piNeg.eta(), piNeg.phi(), piNeg.m());
            piPosVec.SetPtEtaPhiM(ptPiPos, piPos.eta(), piPos.phi(), piPos.m());

            // momenta
            TVector3 piNegMomVec = piNegVec.Vect();
            TVector3 piPosMomVec = piPosVec.Vect();

            TVector3 k0sRecoMomVec = piNegMomVec + piPosMomVec;

            // skip events with pT imbalance, otherwise unfolding should be used
            if (abs(piNegVec.Pt() - piPosVec.Pt()) > 0.1)
                continue;

            // k0sRecoVec.SetPtEtaPhiM(ptPiNeg + ptPiPos, k0s.eta(), k0s.phi(), k0s.m());
            // k0sRecoVec.SetPtEtaPhiM(TMath::Hypot(piNegVec.Px() + piPosVec.Px(), piNegVec.Py() + piPosVec.Py()), k0s.eta(), k0s.phi(), k0s.m());
            k0sRecoVec.SetVectM(k0sRecoMomVec, k0s.m());

            hDauPtRes->Fill(piPos.pT(), (ptPiPos - piPos.pT()) / piPos.pT());
            hDauPtInvMass->Fill(ptPiPos, ((piNegVec + piPosVec).M() - k0s.m()) / k0s.m());

            hDauPtDiffVsPt->Fill(k0s.pT(), ptPiNeg - ptPiPos);
            hDauPtVsPt->Fill(k0s.pT(), ptPiPos);
            hPtInvMass->Fill(k0s.pT(), (piNegVec + piPosVec).M() - k0s.m());

            // boost the pions in the reco mother CM
            TLorentzVector piPosVecCM = piPosVec;
            TLorentzVector piNegVecCM = piNegVec;
            piPosVecCM.Boost(-k0sRecoVec.BoostVector());
            piNegVecCM.Boost(-k0sRecoVec.BoostVector());
            float deltaPtCM = (piPosVecCM.Pt() - piNegVecCM.Pt()) / piPosVecCM.Pt() / TMath::Sqrt(2.);
            hDeltaPtCM->Fill(piPosVec.Pt(), deltaPtCM);

            // boost the pions in the true mother CM
            TLorentzVector piPosVecCMTrueBoost = piPosVec;
            TLorentzVector piNegVecCMTrueBoost = piNegVec;
            piPosVecCMTrueBoost.Boost(-k0sTrueVec.BoostVector());
            piNegVecCMTrueBoost.Boost(-k0sTrueVec.BoostVector());

            float deltaPtCMTrue = (piPosVecCMTrueBoost.Pt() - piNegVecCMTrueBoost.Pt()) / piPosVecCMTrueBoost.Pt() / TMath::Sqrt(2.);
            hDeltaPtCMTrueBoost->Fill(piPosVec.Pt(), deltaPtCMTrue);
        }
    }

    // fit slices of momentum and mass resolutions
    TF1 *f = new TF1("f", "gaus", -0.1, 0.1);
    hDauPtInvMass->FitSlicesY(f, 0, -1, 0, "QNR");
    TH1D *hMean = (TH1D *)gDirectory->Get("hDauPtInvMass_1");
    TH1D *hSigmaMass = (TH1D *)gDirectory->Get("hDauPtInvMass_2");
    hSigmaMass->SetName("hSigmaMass");

    TF1 *fRes = new TF1("fRes", "gaus", -0.1, 0.1);
    hDauPtRes->FitSlicesY(fRes, 0, -1, 0, "QNR");
    TH1D *hSigmaPt = (TH1D *)gDirectory->Get("hDauPtRes_2");
    hSigmaPt->SetName("hDauSigmaPt");

    // put sigmas in the same canvas
    TCanvas *cSigma = new TCanvas("cSigma", "cSigma", 800, 600);
    hSigmaMass->SetLineColor(kRed);
    hSigmaPt->Draw();
    hSigmaMass->Draw("same");
    hSigmaPt->GetYaxis()->SetTitle("Relative Resolutions");
    hSigmaPt->GetXaxis()->SetTitle("p^{gen}_{T} (#pi^{+}) [GeV/c]");
    // build legend
    TLegend *leg = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg->AddEntry(hSigmaPt, "#sigma_{p_{T}} / p_{T} (#pi^{+})", "l");
    leg->AddEntry(hSigmaMass, "#sigma_{M} / M (#pi^{+} + #pi^{-})", "l");
    leg->Draw();

    // Compare the ratios
    TH1D *hRatio = (TH1D *)hSigmaMass->Clone("hRatio");
    hRatio->Divide(hSigmaPt);
    hRatio->SetLineColor(kRed);
    hRatio->SetMarkerColor(kRed);
    hRatio->SetMarkerStyle(20);
    hRatio->GetYaxis()->SetTitle("[#sigma_{M} / M (#pi^{+} + #pi^{-})] / [#sigma_{p_{T}} / p_{T} (#pi^{+})]");

    //__________________________________________________________

    // Fit slices of deltaPtCM and deltaPtCMTrue
    TF1 *fDeltaCM = new TF1("fDeltaCM", "gaus", -0.05, 0.05);
    hDeltaPtCM->FitSlicesY(fDeltaCM, 0, -1, 0, "QNR");
    TH1D *hDeltaPtCMMean = (TH1D *)gDirectory->Get("hDeltaPtCM_1");
    TH1D *hDeltaPtCMSigma = (TH1D *)gDirectory->Get("hDeltaPtCM_2");
    hDeltaPtCMSigma->SetName("hDeltaPtCMSigma");

    TF1 *fDeltaCMTrue = new TF1("fDeltaCMTrue", "gaus", -0.2, 0.2);
    hDeltaPtCMTrueBoost->FitSlicesY(fDeltaCMTrue, 0, -1, 0, "QNR");
    TH1D *hDeltaPtCMTrueMean = (TH1D *)gDirectory->Get("hDeltaPtCMTrueBoost_1");
    TH1D *hDeltaPtCMTrueSigma = (TH1D *)gDirectory->Get("hDeltaPtCMTrueBoost_2");
    hDeltaPtCMTrueSigma->SetName("hDeltaPtCMTrueSigma");

    // put sigmas in the same canvas and compare with momentum resolution
    TCanvas *cDeltaPtCM = new TCanvas("cDeltaPtCM", "cDeltaPtCM", 800, 600);
    cDeltaPtCM->DrawFrame(0., 0., 5., 0.12, ";#it{p}_{T}^{gen}(#pi^{+}) (GeV/#it{c}); #sigma_{p_{T}} / p_{T} (#pi^{+})");
    hDeltaPtCMSigma->SetLineColor(kRed);
    hDeltaPtCMTrueSigma->SetLineColor(kBlue);
    hSigmaPt->SetLineColor(kGreen);
    hSigmaPt->Draw("same");
    hDeltaPtCMSigma->Draw("same");
    hDeltaPtCMTrueSigma->Draw("same");

    // build legend
    TLegend *leg2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    leg2->AddEntry(hSigmaPt, "Real value", "l");
    leg2->AddEntry(hDeltaPtCMSigma, "CM method", "l");
    leg2->AddEntry(hDeltaPtCMTrueSigma, "CM method (True boost)", "l");
    leg2->Draw();

    hSigmaPt->GetYaxis()->SetTitle("#sigma_{p_{T}} / p_{T} (#pi^{+})");
    hSigmaPt->GetXaxis()->SetTitle("p^{gen}_{T} (#pi^{+}) [GeV/c]");

    // ratio to momentum resolution
    TH1D *hRatioDeltaCM = (TH1D *)hDeltaPtCMSigma->Clone("hRatioDeltaCM");
    hRatioDeltaCM->Divide(hSigmaPt);
    hRatioDeltaCM->SetLineColor(kRed);

    TH1D *hRatioDeltaCMTrue = (TH1D *)hDeltaPtCMTrueSigma->Clone("hRatioDeltaCMTrue");
    hRatioDeltaCMTrue->Divide(hSigmaPt);
    hRatioDeltaCMTrue->SetLineColor(kBlue);

    // put ratios in the same canvas
    TCanvas *cRatioDeltaCM = new TCanvas("cRatioDeltaCM", "cRatioDeltaCM", 800, 600);
    hRatioDeltaCM->GetYaxis()->SetTitle("Ratio to #sigma_{p_{T}} / p_{T} (#pi^{+})");

    hRatioDeltaCM->Draw();
    hRatioDeltaCMTrue->Draw("same");

    // build legend
    TLegend *legRatio2 = new TLegend(0.6, 0.6, 0.9, 0.9);
    legRatio2->AddEntry(hRatioDeltaCM, "CM method", "l");
    legRatio2->AddEntry(hRatioDeltaCMTrue, "CM method (True boost)", "l");
    legRatio2->Draw();

    //__________________________________________________________
    // save output
    std::string linearSm = linearSmearing ? "linear" : "flat";
    TFile *outFile = new TFile(Form("results_%s_vect.root", linearSm.c_str()), "RECREATE");
    hDauPtRes->Write();
    hDauPtInvMass->Write();

    hDauPtDiffVsPt->Write();
    hDauPtVsPt->Write();
    hDeltaPtCM->Write();
    hDeltaPtCMTrueBoost->Write();
    hPtInvMass->Write();

    hSigmaMass->Write();
    hSigmaPt->Write();
    cSigma->Write();
    hRatio->Write();

    cDeltaPtCM->Write();
    cRatioDeltaCM->Write();

    outFile->Close();
}

//__________________________________________________________
Pythia8::Particle GetParticle(int pdg)
{
    double mass = gRandom->BreitWigner(TDatabasePDG::Instance()->GetParticle(pdg)->Mass(), TDatabasePDG::Instance()->GetParticle(pdg)->Width()); // not needed for K0S, but ok
    double phi = gRandom->Uniform(2 * TMath::Pi());                                                                                              // flat distribution between 0. and 2Pi
    double pt = gRandom->Uniform(15.);                                                                                                           // flat distribution between 0. and 5. GeV/c
    double mt = TMath::Sqrt(mass * mass + pt * pt);
    double y = -1. + gRandom->Uniform(2.); // flat distribution between -1 and 1

    auto fourMom = ROOT::Math::PxPyPzMVector(pt * TMath::Cos(phi), pt * TMath::Sin(phi), mt * TMath::SinH(y), mass);

    Pythia8::Particle part;
    part.id(pdg);
    part.status(11);
    part.xProd(0.);
    part.yProd(0.);
    part.zProd(0.);
    part.tProd(0.);
    part.m(fourMom.M());
    part.e(fourMom.E());
    part.px(fourMom.Px());
    part.py(fourMom.Py());
    part.pz(fourMom.Pz());

    return part;
}
