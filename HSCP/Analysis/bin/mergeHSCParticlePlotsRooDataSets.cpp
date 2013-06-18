#include "RooDataSet.h"
#include "TFile.h"
#include <vector>
#include <string>

using namespace std;

int main()
{
  string directory =
    "/local/cms/user/turkewitz/HSCP/CMSSW_5_3_8_HSCP/src/HSCP/Analysis/bin/FARM_MakeHSCParticlePlots_data_2012_Mar10/outputs/";
  vector<string> fileNames;
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_0.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_1.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_2.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_3.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_4.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_5.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_6.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_7.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_8.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_9.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_10.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_11.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_12.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_13.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_14.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_15.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_16.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_17.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_18.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_19.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_20.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_21.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_22.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_23.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_24.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_25.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_26.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_27.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_28.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_29.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_30.root");
  fileNames.push_back(directory+"makeHSCParticlePlots_Mar10_Data2012_31.root");

  RooDataSet* rooDataSetOneCandidatePerEvent;
//  RooDataSet* rooDataSetOneCandidatePerEventReduced;
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
    rooDataSetOneCandidatePerEvent->Print();
    thisFile->Close();
  }
  rooDataSetOneCandidatePerEvent->Print();
  //RooDataSet* rooDataSetOneCandidatePerEventReduced = (RooDataSet*) rooDataSetOneCandidatePerEvent->reduce("rooVarPt>50.0");

//  RooDataSet* rooDataSetOneCandidatePerEventReduced = (RooDataSet*) rooDataSetOneCandidatePerEvent->reduce("rooVarIas>0.1&&rooVarPt>70.0");
//  rooDataSetOneCandidatePerEventReduced->Print();
  TFile* outputFile = new TFile((directory+"makeHSCParticlePlots_Mar10_Data2012_some.root").c_str(),"recreate");
  outputFile->cd();
//  rooDataSetOneCandidatePerEventReduced->Write();
  rooDataSetOneCandidatePerEvent->Write();
  rooDataSetCandidates->Write();
  outputFile->Close();

  return 0;
}
