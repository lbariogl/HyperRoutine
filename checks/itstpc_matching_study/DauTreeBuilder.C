#if !defined(CLING) || defined(ROOTCLING)
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/V0.h"
#include "ReconstructionDataFormats/Cascade.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTrack.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsTPC/TrackTPC.h"

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

#include "GlobalTracking/MatchTPCITS.h"

#endif

using GIndex = o2::dataformats::VtxTrackIndex;
using V0 = o2::dataformats::V0;
using MCTrack = o2::MCTrack;
using Cascade = o2::dataformats::Cascade;
using RRef = o2::dataformats::RangeReference<int, int>;
using VBracket = o2::math_utils::Bracket<int>;
using namespace o2::itsmft;
using Vec3 = ROOT::Math::SVector<double, 3>;
using TrackLocTPC = o2::globaltracking::TrackLocTPC;
using TrackLocITS = o2::globaltracking::TrackLocITS;
using rejMap = o2::globaltracking::TrackRejFlag;

// source detBMap
// 0: ITS, 1: TPC, 2: ITS-TPC, 3: TPC-TOF, 4: TPC-TRD, 5: ITS-TPC-TOF, 6: ITS-TPC-TRD, 7: TPC-TRD-TOF, 8: ITS-TPC-TRD-TOF

struct GPart
{
    float genRad = -1, genPt = -1, genRap = -5, genEta = -5;
    int detectorBMap = 0;
    int clusLayerMap = 0;
    int itsRef = -1;
    int tpcRef = -1;
    float itsPt = -1;
    float tpcPt = -1;
    float itsTPCPt = -1;
    float chi2Match = -1;
    bool isSecDau = false;
    bool isReco = false;
    bool isAB = false;
    bool isITSfake = false;
    bool isTPCfake = false;
    bool isITSTPCfake = false;
    int rejFlag = -1;
    int tfNum = -1;
    int pdg = -1;
    int nRefs = 0;
    int nClus = 0;
    int clRefL5 = -1;
    int clRefL6 = -1;
    bool clL5tracked = false;
    bool clL6tracked = false;
};

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG);

void DauTreeBuilder(int dau0PDG = 211, int dau1PDG = 1000020030, int mothPDG = 1010010030,
                    bool debug = false, std::string path = "/data3/fmazzasc/sim/hyp_sim/", std::string outsuffix = "no_mat")
{

    if(outsuffix != "")
        outsuffix = "_" + outsuffix;

    std::string outFileName = "../../DauTreeMC" + outsuffix + ".root";
    std::string treeName = "DauTreeMC";
    TFile outFile = TFile(outFileName.data(), "recreate");
    TTree *DauTree = new TTree(treeName.data(), treeName.data());
    GPart outPart;

    // write the struct to the tree
    DauTree->Branch("DauTree", &outPart);

    // Path to the directory with the files
    TSystemDirectory dir("MyDir", path.data());
    auto files = dir.GetListOfFiles();
    std::vector<std::string> dirs;
    std::vector<int> dirnums;
    std::vector<TString> kine_files;

    for (auto fileObj : *files)
    {
        std::string file = ((TSystemFile *)fileObj)->GetName();
        if (file.substr(0, 2) == "tf")
        {
            int dirnum = stoi(file.substr(2, file.size()));
            // if (dirnum > 1)
            //     continue;
            dirs.push_back(path + file);
            dirnums.push_back(dirnum);
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
    // first we load the geometry
    LOG(info) << "Loading geometry...";
    o2::base::GeometryManager::loadGeometry(dirs[0] + "/o2sim_geometry.root");
    auto gman = o2::its::GeometryTGeo::Instance();
    gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));
    LOG(info) << "Geometry loaded.";

    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        counter++;

        auto &dir = dirs[i];
        auto &dirnum = dirnums[i];
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
        std::vector<o2::its::TrackITS> *ITStracks = nullptr;
        std::vector<int> *ITSTrackClusIdx = nullptr;

        std::vector<o2::tpc::TrackTPC> *TPCtracks = nullptr;

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

        o2::dataformats::MCTruthContainer<o2::MCCompLabel> *clusLabArr = nullptr;
        std::vector<CompClusterExt> *ITSclus = nullptr;

        treeMCTracks->SetBranchAddress("MCTrack", &MCtracks);
        treeITSTPC->SetBranchAddress("TPCITS", &TPCITStracks);
        treeITS->SetBranchAddress("ITSTrack", &ITStracks);
        treeITS->SetBranchAddress("ITSTrackClusIdx", &ITSTrackClusIdx);
        treeTPC->SetBranchAddress("TPCTracks", &TPCtracks);

        treeITS->SetBranchAddress("ITSTrackMCTruth", &labITSvec);
        treeTPC->SetBranchAddress("TPCTracksMCTruth", &labTPCvec);
        treeITSTPC->SetBranchAddress("MatchMCTruth", &labITSTPCvec);
        treeTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTOFvec);
        treeTPCTRD->SetBranchAddress("labels", &labTPCTRDvec);
        treeITSTPCTRD->SetBranchAddress("labelsTRD", &labITSTPCTRDvec);
        treeTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labTPCTRDTOFvec);
        treeITSTPCTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTOFvec);
        treeITSTPCTRDTOF->SetBranchAddress("MatchTOFMCTruth", &labITSTPCTRDTOFvec);

        treeITSclus->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
        treeITSclus->SetBranchAddress("ITSClusterComp", &ITSclus);

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
                GPart dauPart;
                auto &mcTrack = MCtracks->at(mcI);
                auto mcPdg = mcTrack.GetPdgCode();

                if (abs(mcPdg) == dau0PDG || abs(mcPdg) == dau1PDG)
                {
                    auto motherID = mcTrack.getMotherTrackId();
                    if (motherID < 0)
                        continue;
                    auto motherPDG = MCtracks->at(motherID).GetPdgCode();
                    if (abs(motherPDG) != mothPDG)
                        continue;
                    
                    LOG(info) << "Found daughter with PDG " << mcPdg << " from mother with PDG " << motherPDG << " in event " << n << " (global track ID " << mcI << ", mother ID " << motherID << ")";

                    auto &mothTrack = MCtracks->at(motherID);
                    dauPart.isSecDau = true;
                    dauPart.genRad = calcRadius(MCtracks, mothTrack, mcPdg);
                    dauPart.genPt = mcTrack.GetPt();
                    dauPart.genRap = mcTrack.GetRapidity();
                    dauPart.genEta = mcTrack.GetEta();
                    dauPart.tfNum = dirnum;
                    dauPart.pdg = mcPdg;
                }

                GPartMatrix[n][mcI] = dauPart;
            }
        }

        for (int frame = 0; frame < treeITS->GetEntriesFast(); frame++)
        {
            if (!treeITS->GetEvent(frame) || !treeITS->GetEvent(frame) || !treeITSTPC->GetEvent(frame) || !treeTPC->GetEvent(frame) ||
                !treeITSTPCTOF->GetEvent(frame) || !treeTPCTOF->GetEvent(frame) || !treeITSclus->GetEvent(frame) || !treeTPCTRD->GetEvent(frame) ||
                !treeITSTPCTRD->GetEvent(frame) || !treeTPCTRDTOF->GetEvent(frame) || !treeITSTPCTRDTOF->GetEvent(frame))
                continue;

            std::vector<bool> isClusTracked(ITSclus->size(), false);

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
                    if (!gPart.isSecDau)
                        continue;
                    gPart.isReco = true;
                    // update the bit map
                    gPart.detectorBMap |= (1 << bmap[detector]);
                    if (detector == "ITS")
                    {
                        gPart.itsRef = i;
                        gPart.isITSfake = label.isFake();
                        
                        auto &ITStrack = ITStracks->at(i);
                        gPart.itsPt = ITStrack.getPt();
                        
                        // flag tracked clusters
                        auto firstClus = ITStrack.getFirstClusterEntry();
                        auto ncl = ITStrack.getNumberOfClusters();
                        for (unsigned int icl = 0; icl < ncl; icl++)
                        {
                            int clInd = ITSTrackClusIdx->at(firstClus + icl);
                            isClusTracked[clInd] = true;
                        }
                    }

                    if (detector == "TPC")
                    {
                        gPart.tpcRef = i;
                        gPart.isTPCfake = label.isFake();
                        gPart.tpcPt = TPCtracks->at(i).getPt();
                    }

                    if (detector == "ITS-TPC")
                    {
                        gPart.nRefs += 1;
                        auto &track = TPCITStracks->at(i);

                        if (track.getRefITS().getSource() == 24)
                            gPart.isAB = true;
                        gPart.isITSTPCfake = label.isFake();
                        gPart.itsTPCPt = track.getPt();
                        gPart.chi2Match = track.getChi2Match();

                    }
                }
            }

            for (unsigned int iClus{0}; iClus < ITSclus->size(); iClus++)
            {
                auto &clus = ITSclus->at(iClus);
                auto &clusLab = (clusLabArr->getLabels(iClus))[0];
                int layer = gman->getLayer(clus.getSensorID());

                if (!clusLab.isValid())
                    continue;
                auto &gPart = GPartMatrix[clusLab.getEventID()][clusLab.getTrackID()];
                if (!gPart.isSecDau)
                    continue;
                gPart.clusLayerMap |= (1 << layer);

                gPart.nClus += 1;
                if (layer == 5)
                {
                    gPart.clRefL5 = iClus;
                    gPart.clL5tracked = isClusTracked[iClus];
                }
                if (layer == 6)
                {
                    gPart.clRefL6 = iClus;
                    gPart.clL6tracked = isClusTracked[iClus];
                }

                // LOG(info) << "TPC ref: " << gPart.tpcRef;
                // LOG(info) << "Found cluster " << iClus << " in layer " << layer;
                // LOG(info) << "Cluster layer map: " << gPart.clusLayerMap;
            }
        }

        if (debug)
        {
            LOG(info) << "Start to access debug info";
            // loop over debug file
            auto fDebug = TFile::Open((dir + "/dbg_TPCITSmatch.root").data());
            auto fTreeDebug = (TTree *)fDebug->Get("match");
            LOG(info) << "Found " << fTreeDebug->GetEntriesFast() << " entries in the debug tree";
            TrackLocITS *trackLocITS = nullptr;
            TrackLocTPC *trackLocTPC = nullptr;
            o2::MCCompLabel *itsLab = nullptr;
            o2::MCCompLabel *tpcLab = nullptr;
            float chi2match = -1;
            int rejflag = -1;

            fTreeDebug->SetBranchAddress("its", &trackLocITS);
            fTreeDebug->SetBranchAddress("tpc", &trackLocTPC);
            fTreeDebug->SetBranchAddress("itsLbl", &itsLab);
            fTreeDebug->SetBranchAddress("tpcLbl", &tpcLab);
            fTreeDebug->SetBranchAddress("chi2Match", &chi2match);
            fTreeDebug->SetBranchAddress("rejFlag", &rejflag);

            for (int frame = 0; frame < fTreeDebug->GetEntriesFast(); frame++)
            {

                fTreeDebug->GetEntry(frame);
                if (!itsLab->isValid() || !tpcLab->isValid())
                    continue;
                if (itsLab->getEventID() != tpcLab->getEventID() || itsLab->getTrackID() != tpcLab->getTrackID())
                    continue;

                auto &gPart = GPartMatrix[itsLab->getEventID()][itsLab->getTrackID()];
                if (!gPart.isSecDau)
                    continue;
                gPart.rejFlag = rejflag;
                gPart.chi2Match = chi2match;
            }
        }

        // loop over gPartMatrix, save He3 particles
        for (auto &gPartVec : GPartMatrix)
        {
            for (auto &gPart : gPartVec)
            {
                if (!gPart.isSecDau)
                    continue;

                outPart = gPart;
                DauTree->Fill();
            }
        }
    }
    outFile.cd();
    DauTree->Write();
    outFile.Close();
}

double calcRadius(std::vector<MCTrack> *MCTracks, const MCTrack &motherTrack, int dauPDG)
{
    auto idStart = motherTrack.getFirstDaughterTrackId();
    auto idStop = motherTrack.getLastDaughterTrackId();
    for (auto iD{idStart}; iD <= idStop; ++iD)
    {
        auto dauTrack = MCTracks->at(iD);
        if (std::abs(dauTrack.GetPdgCode()) == abs(dauPDG))
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