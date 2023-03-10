#include <TFile.h>
#include <TTree.h>

// loop over AO2D directories and merge the trees into a single one

void MergeTrees(const std::string &inputFileName, const std::string &treeName, const std::string &outputFileName)
{
    TFile inputFile(inputFileName.c_str());
    TDirectory *dir = gDirectory;
    TList *keys = dir->GetListOfKeys();
    TList *outlist = new TList;
    for (int i = 0; i < keys->GetEntries(); i++)
    {
        TKey *key = (TKey *)keys->At(i);
        // if(i>3) continue;
        TString keyName = key->GetName();
        key->Print();
        if (key->IsFolder() && keyName.BeginsWith("DF"))
        {
            TDirectory *dir = (TDirectory *)inputFile.Get(keyName);
            TTree *tree = (TTree *)dir->Get(treeName.c_str());

            if (!tree)
            {
                std::cerr << "Error: could not get TTree from directory " << keyName << std::endl;
                continue;
            }
            else
            {
                // tree->SetDirectory(0);
                outlist->Add(tree);
                std::cout << "Tree Added!" << std::endl;
            }
        }
    }

    //   std::cout << "Loop ended" << std::endl;

    //   // Merge the trees into a single output file
    TFile outputFile(outputFileName.c_str(), "RECREATE");
    TTree *newtree = TTree::MergeTrees(outlist);
    newtree->SetName(treeName.c_str());
    newtree->Write();

    //   outputTree->Write();
    //   outputFile.Close();
}