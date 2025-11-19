#include <TChain.h>
#include <TSystemDirectory.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "SimulationDataFormat/MCTrack.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>

void buildAbsoTree(std::string path = "/data3/fmazzasc/sim/xs_studies/he3_x1.5_new/")
{
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
      // if (dirnum > 10)
      //   continue;
      dirs.push_back(path + file + "/" + "sgn_Kine.root");
      std::cout << dirs.back() << std::endl;
    }
  }

  // create a new tree
  float pt, eta, phi, absoX, absoY, absoZ;
  int process, pdg;

  TFile *fout = new TFile("absorption_tree_x1.5.root", "recreate");
  TTree *t = new TTree("he3candidates", "he3candidates");
  t->Branch("pt", &pt, "pt/F");
  t->Branch("eta", &eta, "eta/F");
  t->Branch("phi", &phi, "phi/F");
  t->Branch("absoX", &absoX, "absoX/F");
  t->Branch("absoY", &absoY, "absoY/F");
  t->Branch("absoZ", &absoZ, "absoZ/F");
  t->Branch("process", &process, "process/I");
  t->Branch("pdg", &pdg, "pdg/I");

  for (auto filename : dirs)
  {
    TTreeReader fReader;
    TFile file(filename.data());
    TTree *tree = (TTree *)file.Get("o2sim");
    TTreeReaderArray<o2::MCTrackT<float>> MCTracks = {fReader, "MCTrack"};
    fReader.SetTree(tree);

    while (fReader.Next())
    {
      for (auto &part : MCTracks)
      {
        if (std::abs(part.GetPdgCode()) != 1000020030 || !part.isPrimary())
        {
          continue;
        }

        bool isTreeFilled = false;
        pt = part.GetPt();
        eta = part.GetEta();
        phi = part.GetPhi();
        pdg = part.GetPdgCode();

        for (int j{part.getFirstDaughterTrackId()}; j <= part.getLastDaughterTrackId(); ++j)
        {
          auto &daughter = MCTracks.At(j);
          if (abs(daughter.GetPdgCode()) == 11 || daughter.GetPdgCode() == 22)
            continue;

          if (isTreeFilled)
            continue;

          absoX = daughter.GetStartVertexCoordinatesX();
          absoY = daughter.GetStartVertexCoordinatesY();
          absoZ = daughter.GetStartVertexCoordinatesZ();
          process = daughter.getProcess();
          t->Fill();
          isTreeFilled = true;
        }

        if(!isTreeFilled){
          process = -1;
          absoX = -999;
          absoY = -999;
          absoZ = -999;
          t->Fill();
        }
      }
    }
  }

  fout->cd();
  t->Write();
  std::cout << "End of the loop" << std::endl;
}