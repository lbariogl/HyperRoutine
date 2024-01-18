// simple counter for checking the number of reconstructed V0s

#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/DecayNBodyIndex.h"

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

const int motherPDG = 1010010030;
const int firstDaughterPDG = 1000020030;
const int secondDaughterPDG = 211;

bool checkV0Decay(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int firstDauPDG, int secDauPDG);

void findableV0s(std::string path = "/data/fmazzasc/its_data/sim/hyp_gap_trig/test_relval")
{

    TH1D* hGen2B = new TH1D("hGen2B", "hGen2B", 1, 0., 1.);
    TH1D* hRec2B = new TH1D("hRec2B", "hRec2B", 1, 0., 1.);

    int injectedParticles = 0;

    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<int> dirnums;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            int dirnum = stoi(file.substr(2, file.size()));
            if (dirnum > 3)
                continue;
            LOG(info) << "Processing " << dirnum;
            dirnums.push_back(dirnum);
        }
    }
    int counter = 0;
    for (unsigned int i = 0; i < dirnums.size(); i++)
    {
        counter++;
        LOG(info) << "Processing " << dirnums[i] << "  #: " << counter << "/" << dirnums.size() << " files";
        auto dir = path + "/tf" + std::to_string(dirnums[i]);
        auto kine_file = dir + "/" + "sgn_" + std::to_string(dirnums[i]) + "_Kine.root";
        // Files
        auto fMCTracks = TFile::Open(kine_file.data());
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
        auto treeTPC = (TTree *)fTPC->Get("tpcrec");

        // MC Tracks
        std::vector<o2::MCTrack> *MCtracks = nullptr;

        std::vector<o2::dataformats::V0Index> *v0Vec = nullptr;
        // MC Labels
        std::vector<o2::MCCompLabel> *labITSvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDvec = nullptr;
        std::vector<o2::MCCompLabel> *labTPCTRDTOFvec = nullptr;
        std::vector<o2::MCCompLabel> *labITSTPCTRDTOFvec = nullptr;

        treeSecondaries->SetBranchAddress("V0sID", &v0Vec);

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);
        treeITSTPCTRD->SetBranchAddress("labelsTRD", &labITSTPCTRDvec);
        treeTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTRDTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITSTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTRDTOFvec);
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
                    hGen2B->Fill(0.5);
                }
            }
        }

        LOG(info) << "Injected particles with PDG " << motherPDG << " : " << injectedParticles;

        for (int frame = 0; frame < treeSecondaries->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeSecondaries->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame))
                continue;

            for (unsigned int iV0 = 0; iV0 < v0Vec->size(); iV0++)
            {
                auto &v0 = v0Vec->at(iV0);
                std::vector<int> motherIDvec;
                std::vector<int> daughterIDvec;
                std::vector<int> evIDvec;

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
                LOG(info) << "--------------------------------";
                LOG(info) << "Found injected V0 with PDG " << motherPDG << " and daughters " << firstDaughterPDG << " and " << secondDaughterPDG;
                LOG(info) << "Daughter 0 PDG: " << pdg0 << "  Daughter 1 PDG: " << pdg1;
                LOG(info) << "Track sources, Dau 0: " << v0.getProngID(0).getSourceName() << "  Dau 1: " << v0.getProngID(1).getSourceName();
                hRec2B->Fill(0.5);
            }
        }
    }

    TFile *outFile = new TFile("findableV0s.root", "RECREATE");
    outFile->cd();
    hGen2B->Write();
    hRec2B->Write();
    outFile->Close();
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