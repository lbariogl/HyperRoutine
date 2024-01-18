#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ITSMFTSimulation/Hit.h"

#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "ITSBase/GeometryTGeo.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITStracking/IOUtils.h"
#include <gsl/gsl>
#include <TLorentzVector.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TSystemDirectory.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TLegend.h"
#include "CommonDataFormat/RangeReference.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

const int motherPDG = 1010010040;
const int firstDaughterPDG = 1000020040;
const int secondDaughterPDG = 211;

// const int motherPDG = 3122;
// const int firstDaughterPDG = 2212;
// const int secondDaughterPDG = 211;

// const int motherPDG = 310;
// const int firstDaughterPDG = 211;
// const int secondDaughterPDG = -211;

o2::its::TrackITS *getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack);
double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2]);
double calcPt(const V0 &v0, double dauMass[2], int dauCharges[2]);
double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);
double calcV0alpha(const V0 &v0);
bool checkV0Decay(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int firstDauPDG, int secDauPDG);
double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2]);

void v0Study(std::string path = "/data/fmazzasc/its_data/sim/hyp_gap_trig/h4l/")
{

    double bins[2] = {2.96, 3.04};
    double motherMass = 2.99131;
    double dauMass[2] = {2.80839160743, 0.13957};
    int dauCharges[2] = {2, 1};

    if (std::abs(motherPDG) == 3122)
    {
        motherMass = 1.115683;
        bins[0] = 1.0;
        bins[1] = 1.2;
        dauMass[0] = 0.938272;
        dauMass[1] = 0.13957;
        dauCharges[0] = 1;
        dauCharges[1] = 1;
    }
    if (motherPDG == 310)
    {
        motherMass = 0.493677;
        bins[0] = 0.4;
        bins[1] = 0.6;
        dauMass[0] = 0.13957;
        dauMass[1] = 0.13957;
        dauCharges[0] = 1;
        dauCharges[1] = 1;
    }

    if (std::abs(motherPDG) == 1010010040)
    {
        motherMass = 3.929;
        bins[0] = 3.8;
        bins[1] = 4.2;
        dauMass[0] = 3.727;
        dauMass[1] = 0.13957;
        dauCharges[0] = 2;
        dauCharges[1] = 1;
    }

    int injectedParticles = 0;

    std::vector<TH1D *> hists(5);
    hists[0] = new TH1D("recoPDGits", "Reconstructed ITS PDG;;Efficiency", 3, 0, 3);
    hists[1] = new TH1D("recoPDGtpc", "Reconstructed TPC PDG;;Efficiency", 3, 0, 3);
    hists[2] = new TH1D("recoPDGitsTPC", "Reconstructed ITS-TPC PDG;;Efficiency", 3, 0, 3);
    hists[3] = new TH1D("recoPDGtpcTOF", "Reconstructed TPC-TOF PDG;;Efficiency", 3, 0, 3);
    hists[4] = new TH1D("recoPDGitsTPCTOF", "Reconstructed ITS-TPC-TOF PDG;;Efficiency", 3, 0, 3);

    TH1D *histInvMass = new TH1D("V0 invariant mass", "; V0 Mass (GeV/c^{2}); Counts", 30, bins[0], bins[1]);

    TH1D *histGenRadius = new TH1D("Gen Radius", "; Gen Radius (cm); Counts", 300, 0, 90);
    TH1D *histRecRadius = new TH1D("Rec Radius", "; Rec Radius (cm); Counts", 300, 0, 90);
    TH1D *histoGenPt = new TH1D("Gen Pt", "; Gen Pt (GeV/c); Counts", 50, 0, 10);
    TH1D *histoRecPt = new TH1D("Rec Pt", "; Rec Pt (GeV/c); Counts", 50, 0, 10);
    TH1D *histoRecPtITS = new TH1D("Rec Pt ITS", "; Rec Pt (GeV/c); Counts", 50, 0, 10);
    TH1D *histoRecPtNoITS = new TH1D("Rec Pt NoITS", "; Rec Pt (GeV/c); Counts", 50, 0, 10);

    TH1D *histGenDecLength = new TH1D("Gen Dec Length", "; Gen Dec Length (cm); Counts", 300, 0, 90);
    TH1D *histGenLifetime = new TH1D("Gen ct", "; Gen ct (cm); Counts", 300, 0, 90);
    TH1D *histRecDecLength = new TH1D("Rec Dec Length", "; Gen Dec Length (cm); Counts", 300, 0, 90);

    TH1D *histGeneratedV0s = new TH1D("# of generated V0s", ";; Counts", 1, 0, 1);

    TH2D *histV0radiusRes = new TH2D("V0 radius resolution", "; Gen Radius (cm); Gen - Rec / Gen; Counts", 400, 0, 90, 20, -1, 1);

    TH1D *histITSHits = new TH1D("V0 candidate ITS hits", " Number of ITS hits; #Hits; Counts", 7, 0.5, 7.5);
    TH1D *histITScounter = new TH1D("V0 candidate ITS hits and tracks counter", ";; Counts/(# of V0s)", 3, 0, 3);

    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            int dirnum = stoi(file.substr(2, file.size()));
            // if (dirnum != 4)
            //     continue;
            LOG(info) << "Processing " << dirnum;
            dirs.push_back(path + file);
            auto innerdir = (TSystemDirectory *)fileObj;
            auto innerfiles = innerdir->GetListOfFiles();
            for (auto innerfileObj : *innerfiles)
            {
                TString innerfile = ((TSystemFile *)innerfileObj)->GetName();
                if (innerfile.EndsWith("Kine.root") && innerfile.Contains("sgn"))
                {
                    kine_files.push_back(innerfile);
                }
            }
        }
    }
    int counter = 0;
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        counter++;
        LOG(info) << "Processing " << dirs[i] << "  #: " << counter << "/" << dirs.size() << " files";
        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
        auto fSecondaries = TFile::Open((dir + "/o2_secondary_vertex.root").data());
        auto fITSTPC = TFile::Open((dir + "/o2match_itstpc.root").data());
        auto fTPCTOF = TFile::Open((dir + "/o2match_tof_tpc.root").data());
        auto fTPCTRD = TFile::Open((dir + "/trdmatches_tpc.root").data());
        auto fITSTPCTOF = TFile::Open((dir + "/o2match_tof_itstpc.root").data());
        auto fITS = TFile::Open((dir + "/o2trac_its.root").data());
        auto fClusITS = TFile::Open((dir + "/o2clus_its.root").data());
        auto fTPC = TFile::Open((dir + "/tpctracks.root").data());
        auto fITSTPCTRD = TFile::Open((dir + "/trdmatches_itstpc.root").data());
        auto fTPCTRDTOF = TFile::Open((dir + "/o2match_tof_tpctrd.root").data());
        auto fITSTPCTRDTOF = TFile::Open((dir + "/o2match_tof_itstpctrd.root").data());

        // Trees
        auto treeMCTracks = (TTree *)fMCTracks->Get("o2sim");
        auto treeSecondaries = (TTree *)fSecondaries->Get("o2sim");
        auto treeITSTPC = (TTree *)fITSTPC->Get("matchTPCITS");
        auto treeITSTPCTOF = (TTree *)fITSTPCTOF->Get("matchTOF");
        auto treeTPCTOF = (TTree *)fTPCTOF->Get("matchTOF");
        auto treeTPCTRD = (TTree *)fTPCTRD->Get("tracksTRD");
        auto treeITSTPCTRD = (TTree *)fITSTPCTRD->Get("tracksTRD");
        auto treeTPCTRDTOF = (TTree *)fTPCTRDTOF->Get("matchTOF");
        auto treeITSTPCTRDTOF = (TTree *)fITSTPCTRDTOF->Get("matchTOF");

        auto treeITS = (TTree *)fITS->Get("o2sim");
        auto treeITSclus = (TTree *)fClusITS->Get("o2sim");
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;
        std::vector<o2::itsmft::Hit> *ITSHits = nullptr;

        // Secondary Vertices
        std::vector<Cascade> *cascVec = nullptr;
        std::vector<V0> *v0Vec = nullptr;

        // ITS tracks
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;

        // Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDTOFvec = nullptr;

        // Clusters
        std::vector<CompClusterExt> *ITSclus = nullptr;
        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;
        std::vector<unsigned char> *ITSpatt = nullptr;

        treeSecondaries->SetBranchAddress("Cascades", &cascVec);
        treeSecondaries->SetBranchAddress("V0s", &v0Vec);

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);
        treeITSTPCTRD->SetBranchAddress("labelsTRD", &labITSTPCTRDvec);
        treeTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTRDTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITSTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTRDTOFvec);

        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);
        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);

        // define detector map
        std::map<std::string, std::vector<o2::MCCompLabel> *> map{{"ITS", labITSvec}, {"TPC", labTPCvec}, {"ITS-TPC", labITSTPCvec}, {"TPC-TOF", labTPCTOFvec}, {"TPC-TRD", labTPCTRDvec}, {"ITS-TPC-TOF", labITSTPCTOFvec}, {"ITS-TPC-TRD", labITSTPCTRDvec}, {"TPC-TRD-TOF", labTPCTRDTOFvec}, {"ITS-TPC-TRD-TOF", labITSTPCTRDTOFvec}};

        // fill MC matrix
        std::vector<std::vector<o2::MCTrack>> mcTracksMatrix;
        auto nev = treeMCTracks->GetEntriesFast();
        mcTracksMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);

            mcTracksMatrix[n].resize(MCtracks->size());
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                mcTracksMatrix[n][mcI] = MCtracks->at(mcI);
                if (std::abs(MCtracks->at(mcI).GetPdgCode()) == motherPDG && MCtracks->at(mcI).isPrimary())
                {
                    auto &mcTrack = mcTracksMatrix[n][mcI];
                    if (!checkV0Decay(MCtracks, mcTrack, firstDaughterPDG, secondDaughterPDG))
                    {
                        continue;
                    }
                    injectedParticles++;
                    histGeneratedV0s->Fill(0.5);
                    double L = calcDecLength(MCtracks, mcTrack, firstDaughterPDG);
                    histGenDecLength->Fill(L);
                    histGenLifetime->Fill(L * motherMass / mcTrack.GetP());
                    histGenRadius->Fill(calcRadius(MCtracks, mcTrack, firstDaughterPDG));
                    histoGenPt->Fill(mcTrack.GetPt());
                }
            }
        }

        for (int frame = 0; frame < treeSecondaries->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame))
                continue;
            for (unsigned int iV0 = 0; iV0 < v0Vec->size(); iV0++)
            {
                auto &v0 = v0Vec->at(iV0);
                std::vector<int> motherIDvec;
                std::vector<int> daughterIDvec;
                std::vector<int> evIDvec;
                std::vector<bool> isITSvec;

                for (int iDaugh = 0; iDaugh < 2; iDaugh++)
                {

                    if (map[v0.getProngID(iDaugh).getSourceName()])
                    {
                        auto source = v0.getProngID(iDaugh).getSourceName();
                        auto labTrackType = map[source];
                        auto lab = labTrackType->at(v0.getProngID(iDaugh).getIndex());

                        int trackID, evID, srcID;
                        bool fake;
                        lab.get(trackID, evID, srcID, fake);
                        if (!lab.isNoise() && lab.isValid() && lab.isCorrect())
                        {
                            auto motherID = mcTracksMatrix[evID][trackID].getMotherTrackId();
                            motherIDvec.push_back(motherID);
                            daughterIDvec.push_back(trackID);
                            evIDvec.push_back(evID);
                            isITSvec.push_back(source == "ITS");
                        }
                    }
                }

                if (motherIDvec.size() < 2)
                    continue;
                if (motherIDvec[0] != motherIDvec[1] || evIDvec[0] != evIDvec[1])
                    continue;
                if (motherIDvec[0] < 0 || motherIDvec[0] > 10000)
                    continue;

                int pdg0 = mcTracksMatrix[evIDvec[0]][daughterIDvec[0]].GetPdgCode();
                int pdg1 = mcTracksMatrix[evIDvec[0]][daughterIDvec[1]].GetPdgCode();

                if (!(std::abs(pdg0) == firstDaughterPDG && std::abs(pdg1) == secondDaughterPDG) && !(std::abs(pdg0) == secondDaughterPDG && std::abs(pdg1) == firstDaughterPDG))
                    continue;

                auto motherTrack = mcTracksMatrix[evIDvec[0]][motherIDvec[0]];

                if (!motherTrack.isPrimary() || std::abs(motherTrack.GetPdgCode()) != motherPDG)
                    continue;

                auto genRad = calcRadius(&mcTracksMatrix[evIDvec[0]], motherTrack, firstDaughterPDG);

                double dauMassTemp[2] = {dauMass[0], dauMass[1]};
                int dauChargesTemp[2] = {dauCharges[0], dauCharges[1]};

                if (calcV0alpha(v0) < 0)
                {
                    std::swap(dauMassTemp[0], dauMassTemp[1]);
                    std::swap(dauChargesTemp[0], dauChargesTemp[1]);
                }

                histInvMass->Fill(calcMass(v0, dauMassTemp, dauChargesTemp));
                histRecDecLength->Fill(TMath::Sqrt(v0.calcR2() + v0.getZ() * v0.getZ()));

                auto recRad = TMath::Sqrt(v0.calcR2());
                histRecRadius->Fill(recRad);
                histoRecPt->Fill(calcPt(v0, dauMassTemp, dauChargesTemp));
                if ((std::abs(pdg0) == firstDaughterPDG && isITSvec[0]) || (std::abs(pdg1) == firstDaughterPDG && isITSvec[1]))
                {
                    histoRecPtITS->Fill(calcPt(v0, dauMassTemp, dauChargesTemp));
                }
                else
                {
                    histoRecPtNoITS->Fill(calcPt(v0, dauMassTemp, dauChargesTemp));
                }

                histV0radiusRes->Fill(genRad, (recRad - genRad) / genRad);

                bool he3index = pdg1 == 1000020030;
                auto he3track = he3index ? v0.getProng(1) : v0.getProng(0);
                auto he3trackID = he3index ? v0.getProngID(1) : v0.getProngID(0);

                if (motherTrack.leftTrace())
                {
                    histITScounter->Fill(0);
                    std::cout << "ITS sees mother hits! " << std::endl;
                    o2::its::TrackITS *motherITStrack = getITSTrack(evIDvec[0], motherIDvec[0], treeITS, labITSvec, ITStracks);

                    if (motherITStrack != nullptr)
                    {
                        histITSHits->Fill(motherITStrack->getNClusters());
                        motherITStrack->getNFakeClusters() > 0 ? histITScounter->Fill(2) : histITScounter->Fill(1);
                    }
                }
            }
        }
    }
    auto outFile = TFile("hyp_study.root", "recreate");

    TH1D *histoEffvsRadius = (TH1D *)histRecRadius->Clone("histoEffvsRadius");
    histoEffvsRadius->Divide(histGenRadius);
    histoEffvsRadius->GetYaxis()->SetTitle("Efficiency");
    histoEffvsRadius->GetXaxis()->SetTitle("V0 Radius (cm)");
    histoEffvsRadius->Write();

    TH1D *histoEffvsPt = (TH1D *)histoRecPt->Clone("histoEffvsPt");
    histoEffvsPt->Divide(histoGenPt);
    histoEffvsPt->GetYaxis()->SetTitle("#epsilon #times Acc.");
    histoEffvsPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    histoEffvsPt->Write();

    TH1D *histoEffvsPtITS = (TH1D *)histoRecPtITS->Clone("histoEffvsPtITS");
    histoEffvsPtITS->Divide(histoGenPt);
    histoEffvsPtITS->GetYaxis()->SetTitle("#epsilon #times Acc.");
    histoEffvsPtITS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    histoEffvsPtITS->Write();

    TH1D *histoEffvsPtNoITS = (TH1D *)histoRecPtNoITS->Clone("histoEffvsPtNoITS");
    histoEffvsPtNoITS->Divide(histoGenPt);
    histoEffvsPtNoITS->GetYaxis()->SetTitle("#epsilon #times Acc.");
    histoEffvsPtNoITS->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    histoEffvsPtNoITS->Write();

    histInvMass->Write();
    histGenLifetime->Write();
    histRecRadius->Write();
    histoGenPt->Write();
    histoRecPt->Write();
    histoRecPtITS->Write();
    histoRecPtNoITS->Write();
    histRecDecLength->Write();
    histGenDecLength->Write();
    histGenRadius->Write();
    histV0radiusRes->Write();
    histITSHits->Write();
    histITScounter->Write();
    histGeneratedV0s->Write();
    outFile.Close();
}

o2::its::TrackITS *getITSTrack(int motherEvID, int motherTrackID, TTree *ITStree, std::vector<o2::MCCompLabel> *ITSlabel, std::vector<o2::its::TrackITS> *ITStrack)
{
    o2::its::TrackITS *motherTrack{nullptr};

    for (int frame = 0; frame < ITStree->GetEntriesFast(); frame++)
    {
        if (!ITStree->GetEvent(frame) || !ITStree->GetEvent(frame))
            continue;
        if (!ITStree->GetEvent(frame))
        {
            continue;
        }
        for (unsigned int iTrack{0}; iTrack < ITSlabel->size(); ++iTrack)
        {
            auto lab = ITSlabel->at(iTrack);
            int trackID, evID, srcID;
            bool fake;
            lab.get(trackID, evID, srcID, fake);
            if (!lab.isNoise() && lab.isValid())
            {
                if (evID == motherEvID and trackID == motherTrackID)
                {
                    std::cout << "Matching indexes: " << evID << "   " << motherTrackID << std::endl;
                    motherTrack = &ITStrack->at(iTrack);
                    std::cout << "ITS sees mother track! " << std::endl;
                    return motherTrack;
                }
            }
        }
    }
    return motherTrack;
};

double calcV0alpha(const V0 &v0)
{
    std::array<float, 3> fV0mom, fPmom, fNmom = {0, 0, 0};
    v0.getProng(0).getPxPyPzGlo(fPmom);
    v0.getProng(1).getPxPyPzGlo(fNmom);
    v0.getPxPyPzGlo(fV0mom);

    TVector3 momNeg(fNmom[0], fNmom[1], fNmom[2]);
    TVector3 momPos(fPmom[0], fPmom[1], fPmom[2]);
    TVector3 momTot(fV0mom[0], fV0mom[1], fV0mom[2]);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (std::abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto decLength = (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) *
                                 (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) +
                             (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) *
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY());
            return sqrt(decLength);
        }
    }
    return -1;
}

bool checkV0Decay(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int firstDauPDG, int secDauPDG)
{
    bool dau1 = false;
    bool dau2 = false;
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    if (idStart == -1 || idStop == -1)
    {
        return false;
    }
    for (auto iD{idStart}; iD <= idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        // LOG(info) << "Dau PDG: " << dauTrack.GetPdgCode();
        if (std::abs(dauTrack.GetPdgCode()) == firstDauPDG)
        {
            dau1 = true;
        }
        if (std::abs(dauTrack.GetPdgCode()) == secDauPDG)
        {
            dau2 = true;
        }
    }
    return dau1 && dau2;
}

double calcDecLength(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD < idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (std::abs(dauTrack.GetPdgCode()) == dauPDG)
        {
            auto decLength = (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) *
                                 (dauTrack.GetStartVertexCoordinatesX() - motherTrack.GetStartVertexCoordinatesX()) +
                             (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) *
                                 (dauTrack.GetStartVertexCoordinatesY() - motherTrack.GetStartVertexCoordinatesY()) +
                             (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ()) *
                                 (dauTrack.GetStartVertexCoordinatesZ() - motherTrack.GetStartVertexCoordinatesZ());
            return sqrt(decLength);
        }
    }
    return -1;
}

double calcMass(const V0 &v0, double dauMass[2], int dauCharges[2])
{
    std::vector<o2::dataformats::V0::Track> dauTracks = {v0.getProng(0), v0.getProng(1)};
    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMass[i];
        track.getPxPyPzGlo(p);
        int charge = dauCharges[i];
        prong.SetVectM({charge * p[0], charge * p[1], charge * p[2]}, mass);
        moth += prong;
    }
    return moth.M();
}

double calcPt(const V0 &v0, double dauMass[2], int dauCharges[2])
{
    std::vector<o2::dataformats::V0::Track> dauTracks = {v0.getProng(0), v0.getProng(1)};
    TLorentzVector moth, prong;
    std::array<float, 3> p;
    for (int i = 0; i < 2; i++)
    {
        auto &track = dauTracks[i];
        auto &mass = dauMass[i];
        track.getPxPyPzGlo(p);
        int charge = dauCharges[i];
        prong.SetVectM({charge * p[0], charge * p[1], charge * p[2]}, mass);
        moth += prong;
    }
    return moth.Pt();
}
