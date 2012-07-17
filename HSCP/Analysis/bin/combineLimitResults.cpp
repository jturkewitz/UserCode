//#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "RooStats/HypoTestInverterResult.h"

using namespace std;
using namespace RooStats;

vector<string> getFilesInDir(string dir, string signal)
{
  vector<string> files;
  DIR *dp;
  struct dirent *dirp;
  if((dp  = opendir(dir.c_str())) == NULL)
  {
    cout << "Error(" << errno << ") opening " << dir << endl;
    return vector<string>();
  }

  while ((dirp = readdir(dp)) != NULL) {
    string fileName = string(dirp->d_name);
    if(fileName.find("expected") != string::npos)
      continue;
    if(fileName.find("doLimits") == string::npos)
      continue;
    if(fileName.find(signal) == string::npos)
      continue;
    if(fileName.find("combined") != string::npos)
      continue;
    if(signal.find("N") == string::npos && fileName.find("N") != string::npos)
      continue;
    if(signal.find("f5") == string::npos && fileName.find("f5") != string::npos)
      continue;

    int dotRootPos = fileName.find(".root");
    int indexPos = fileName.find("index");
    string indexStr = fileName.substr(indexPos+5,dotRootPos-1);
    int index = atoi(indexStr.c_str());
    //if(fileName.find("GMStau100") != string::npos)
    //{
    //  if(index > 12 || index < 5)
    //    continue;
    //}
    //else if(fileName.find("Gluino1000") != string::npos)
    //{
    //  if(index > 12 || index < 5)
    //    continue;
    //}
    //if(fileName.find("Gluino1200") != string::npos)
    //{
    //  if(index > 14 || index < 7)
    //    continue;
    //}
    //if(fileName.find("GMStau") != string::npos)
    //{
    //  if(index > 7 || index < 1)
    //    continue;
    //}
    //else if(fileName.find("Gluino") != string::npos)
    //{
    //  if(index > 20 || index < 17)
    //    continue;
    //}
    //else
    //{
    //  if(index > 20 || index < 1)
    //    continue;
    //}
    
    string toAdd = dir;
    toAdd+="/";
    toAdd+=fileName;
    files.push_back(toAdd);
  }
  closedir(dp);

  return files;
}

int main(int argc, char* argv[])
{
  if(argc < 2)
  {
    cout << "ERROR: Argument count less than 2." << endl;
    cout << "Usage: combineLimitResults DirectoryOfFiles SignalName" << endl;
    return -1;
  }

  string outputDir(argv[1]);
  string signalName(argv[2]);
  if(outputDir.length() < 1)
  {
    cout << "Found outputDir = " << outputDir << " which has length < 1; exiting" << endl;
    return -1;
  }
  if(signalName.length() < 1)
  {
    cout << "Found signalName = " << signalName << " which has length < 1; exiting" << endl;
    return -1;
  }

  //'doLimits_'+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)+'.root'


  // kill RooFit/RooStats output
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  vector<double> trialCrossSections;
  int fileCounter = 0;
  cout << "Doing " << signalName << "..." << endl;
  vector<string> thisSignalFiles = getFilesInDir(outputDir, signalName);
  HypoTestInverterResult* r = 0;
  HypoTestInverterResult* thisResult = 0;
  TFile* thisTFile = 0;
  for(vector<string>::const_iterator fileItr = thisSignalFiles.begin();
      fileItr != thisSignalFiles.end(); ++fileItr)
  {
    //cout << "Examining file: " << *fileItr << endl;
    thisTFile = new TFile(fileItr->c_str(),"read");
    if(!r)
    {
      r =  dynamic_cast<HypoTestInverterResult*>(thisTFile->Get("result_SigXsec"));
      ++fileCounter;
    }
    else
    {
      delete thisResult;
      thisResult = dynamic_cast<HypoTestInverterResult*>(thisTFile->Get("result_SigXsec"));
      //cout << "Number of x-values in result: " << thisResult->ArraySize() << endl;
      //for(int i=0; i < thisResult->ArraySize(); ++i)
      //  cout << "\t" << thisResult->GetXValue(i) << endl;
      //if(find(trialCrossSections.begin(),trialCrossSections.end(),thisResult->GetXValue(0)) == trialCrossSections.end())
      bool add = true;
      //if(signalName=="GMStau200" && thisResult->GetXValue(0) < 0.0014)
      //  add = false;
      //if(signalName=="Gluino1000" && (thisResult->GetXValue(0) < 0.0025 || thisResult->GetXValue(0) > 0.004))
      //  add = false;

      if(add)
      {
        trialCrossSections.push_back(thisResult->GetXValue(0));
        r->Add(*thisResult);
        //cout << "Add cross section: " << thisResult->GetXValue(0) << endl;
        ++fileCounter;
      }
      //else
      //  cout << "Skipping already-seen cross section: " << thisResult->GetXValue(0) << endl;
    }
    //cout << endl;
  }

  string outputPath = outputDir;
  outputPath+=signalName;
  outputPath+="_doLimits_combined.root";
  cout << "Done merging " << fileCounter << " files, writing file " << outputPath << endl;

  TFile* thisSignalCombinedTFile = new TFile(outputPath.c_str(),"recreate");
  thisSignalCombinedTFile->cd();
  r->SetName("result_SigXsec");
  r->Write();
  thisSignalCombinedTFile->Close();

  double upperLimit = r->UpperLimit();
  double ulError = r->UpperLimitEstimatedError();
  std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;
  // compute expected limit
  std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
  std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
  std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
  std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
  std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;

  return 0;
}
