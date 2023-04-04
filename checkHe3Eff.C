#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

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
#include "DetectorsVertexing/DCAFitterN.h"
#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;

const int hypPDG = 1010010030;
const int he3PDG = 1000020030;

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

void checkHe3Eff()
{

    struct GPart
    {
        float genRad = -1, genPt = -1;
        int detectorBMap = 0;
        bool isSecHe3 = false;
        bool isReco = false;
        bool isAB = false;
    };

    // write the output tree
    TFile outFile = TFile("He3TreeMC.root", "recreate");
    TTree *He3Tree = new TTree("He3TreeMC", "He3TreeMC");
    float genRad, genPt;
    // source bitmap
    // 0: ITS, 1: TPC, 2: ITS-TPC, 3: TPC-TOF, 4: TPC-TRD, 5: ITS-TPC-TOF, 6: ITS-TPC-TRD, 7: TPC-TRD-TOF, 8: ITS-TPC-TRD-TOF
    int detBMap = 0;
    bool isReco = false;
    bool isAB = false;

    He3Tree->Branch("genRad", &genRad);
    He3Tree->Branch("genPt", &genPt);
    He3Tree->Branch("isReco", &isReco);
    He3Tree->Branch("detBMap", &detBMap);
    He3Tree->Branch("isAB", &isAB);

    // Path to the directory with the files
    std::string path = "/data/fmazzasc/its_data/sim/hyp_2_body/";
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
            // if (dirnum > 3)
            //     continue;
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

        auto &dir = dirs[i];
        auto &kine_file = kine_files[i];
        LOG(info) << "Processing " << dir;
        LOG(info) << "File # " << counter << " of " << dirs.size();

        // Files
        auto fMCTracks = TFile::Open((TString(dir + "/") + kine_file));
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
        std::vector<o2::dataformats::TrackTPCITS> *TPCITStracks = nullptr;

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

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITSTPC->SetBranchAddress("TPCITS", &TPCITStracks);

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
        std::map<std::string, unsigned int> bmap{{"ITS", 0}, {"TPC", 1}, {"ITS-TPC", 2}, {"TPC-TOF", 3}, {"TPC-TRD", 4}, {"ITS-TPC-TOF", 5}, {"ITS-TPC-TRD", 6}, {"TPC-TRD-TOF", 7}, {"ITS-TPC-TRD-TOF", 8}};

        // fill MC matrix
        std::vector<std::vector<GPart>> GPartMatrix;

        auto nev = treeMCTracks->GetEntriesFast();
        GPartMatrix.resize(nev);
        for (int n = 0; n < nev; n++)
        { // loop over MC events
            treeMCTracks->GetEvent(n);
            GPartMatrix[n].resize(MCtracks->size());
            for (unsigned int mcI{0}; mcI < MCtracks->size(); ++mcI)
            {
                GPart he3Part;
                auto &mcTrack = MCtracks->at(mcI);
                auto pdg = mcTrack.GetPdgCode();
                if (abs(pdg) == he3PDG)
                {
                    auto motherID = mcTrack.getMotherTrackId();
                    if (motherID < 0)
                        continue;
                    auto motherPDG = MCtracks->at(motherID).GetPdgCode();
                    if (abs(motherPDG) != hypPDG)
                        continue;
                    auto &hypTrack = MCtracks->at(motherID);
                    he3Part.isSecHe3 = true;
                    he3Part.genRad = calcRadius(MCtracks, hypTrack, he3PDG);
                    he3Part.genPt = mcTrack.GetPt();
                }

                GPartMatrix[n][mcI] = he3Part;
            }
        }

        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame))
                continue;

            // loop over all the map entries
            for (auto const &[detector, labelVector] : map)
            {
                // loop over all the labels
                for (unsigned int i = 0; i < labelVector->size(); i++)
                {
                    auto &label = labelVector->at(i);
                    // check if the label is valid
                    if (!label.isValid())
                        continue;
                    // check if the label is a primary particle
                    auto &gPart = GPartMatrix[label.getEventID()][label.getTrackID()];
                    if (!gPart.isSecHe3)
                        continue;
                    gPart.isReco = true;
                    // update the bit map
                    gPart.detectorBMap |= (1 << bmap[detector]);
                    if (detector == "ITS-TPC")
                    {
                        auto &track = TPCITStracks->at(i);
                        // LOG(info) << "ITS Ref: " << track.getRefITS().getIndex() << " TPC Ref: " << track.getRefTPC().getIndex();
                        if (track.getRefITS().getSource() == 24)
                            gPart.isAB = true;
                    }
                }
            }
        }

        // loop over gPartMatrix, save He3 particles
        for (auto &gPartVec : GPartMatrix)
        {
            for (auto &gPart : gPartVec)
            {
                if (!gPart.isSecHe3)
                    continue;

                genRad = gPart.genRad;
                genPt = gPart.genPt;
                detBMap = gPart.detectorBMap;
                isReco = gPart.isReco;
                isAB = gPart.isAB;
                He3Tree->Fill();
            }
        }
    }
    outFile.cd();
    He3Tree->Write();
    outFile.Close();
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