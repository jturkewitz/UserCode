#include "RooDataSet.h"
#include "TFile.h"
#include <vector>
#include <string>

using namespace std;

int main()
{
  string directory =
    "/local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_dataWithTightRPCTrig_absEta_ptErrorPresel_feb6/outputs/";
  vector<string> fileNames;
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_0.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_1.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_2.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_3.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_4.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_5.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_6.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_7.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_8.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_9.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_10.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_11.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_12.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_13.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_14.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_15.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_16.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_17.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_18.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_19.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_feb6_Data2011_20.root");

  RooDataSet* rooDataSetOneCandidatePerEvent;
  RooDataSet* rooDataSetCandidates;
  bool initEvent = true;
  bool initCand = true;

  for(vector<string>::const_iterator fileNameItr = fileNames.begin();
      fileNameItr != fileNames.end(); ++fileNameItr)
  {
    std::cout << "Merging in file: " << *fileNameItr << std::endl;
    TFile* thisFile = TFile::Open(fileNameItr->c_str());
    if(initEvent)
    {
      rooDataSetOneCandidatePerEvent = (RooDataSet*) thisFile->Get("rooDataSetOneCandidatePerEvent");
      initEvent = false;
    }
    else
      rooDataSetOneCandidatePerEvent->append(*((RooDataSet*)thisFile->Get("rooDataSetOneCandidatePerEvent")));

    if(initCand)
    {
      rooDataSetCandidates = (RooDataSet*) thisFile->Get("rooDataSetCandidates");
      initCand = false;
    }
    else
      rooDataSetCandidates->append(*((RooDataSet*)thisFile->Get("rooDataSetCandidates")));

    thisFile->Close();
  }

  TFile* outputFile = new TFile((directory+"makeHSCParticlePlots_feb1_Data2011_all.root").c_str(),"recreate");
  outputFile->cd();
  rooDataSetOneCandidatePerEvent->Write();
  rooDataSetCandidates->Write();
  outputFile->Close();

  return 0;
}
