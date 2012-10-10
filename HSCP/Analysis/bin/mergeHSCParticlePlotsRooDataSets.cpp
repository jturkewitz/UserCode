#include "RooDataSet.h"
#include "TFile.h"
#include <vector>
#include <string>

using namespace std;

int main()
{
  string directory =
    "/local/cms/user/turkewitz/HSCP/CMSSW_5_2_3_patch4_HSCP/src/HSCP/Analysis/bin/FARM_MakeHSCParticlePlots_data_2012_Aug15/outputs/";
  vector<string> fileNames;
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_0.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_1.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_2.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_3.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_4.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_5.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_6.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_7.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_8.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_9.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_10.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_11.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_12.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_13.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_14.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_15.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_16.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_17.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_18.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_19.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_20.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_21.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Aug15_Data2012_22.root");

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

  TFile* outputFile = new TFile((directory+"makeHSCParticlePlots_Aug15_Data2012_all.root").c_str(),"recreate");
  outputFile->cd();
  rooDataSetOneCandidatePerEvent->Write();
  rooDataSetCandidates->Write();
  outputFile->Close();

  return 0;
}
