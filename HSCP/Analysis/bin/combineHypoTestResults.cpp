//#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "RooStats/HypoTestResult.h"

using namespace std;
using namespace RooStats;

std::string floatToString(float num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}

vector<string> getFilesInDir(string dir, string signal, float poi)
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
    if(fileName.find("asymptotic") != string::npos)
      continue;
    if(fileName.find("root") == string::npos)
      continue;
    if(fileName.find(signal) == string::npos)
      continue;
    if(fileName.find(floatToString(poi)) == string::npos)
      continue;
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
  if(argc < 3)
  {
    cout << "ERROR: Argument count less than 3." << endl;
    cout << "Usage: combineHypoTestResults DirectoryOfFiles SignalName poiVal" << endl;
    return -1;
  }

  // directory with outputs of toys
  string outputDir = string(argv[1]);
  if(outputDir.length() < 1)
  {
    cout << "Found outputDir = " << outputDir << " which has length < 1; exiting" << endl;
    return -1;
  }
  string signalName = string(argv[2]);
  if(signalName.length() < 1)
  {
    cout << "Found signalName = " << signalName << " which has length < 1; exiting" << endl;
    return -1;
  }
  float poi = atof(argv[3]);
  if(poi < 0 || poi == 0)
  {
    cout << "Found poi < 0 || poi == 0: " << poi << " exiting" << endl;
    return -1;
  }

  int fileCounter = 0;
  vector<string> thisSignalFiles = getFilesInDir(outputDir, signalName, poi);
  HypoTestResult* r = 0;
  HypoTestResult* thisResult = 0;
  TFile* thisTFile = 0;
  for(vector<string>::const_iterator fileItr = thisSignalFiles.begin();
      fileItr != thisSignalFiles.end(); ++fileItr)
  {
    thisTFile = new TFile(fileItr->c_str(),"read");
    if(!r)
      r =  dynamic_cast<HypoTestResult*>(thisTFile->Get("HypoTestCalculator_result"));
    else
    {
      std::cout << "merging in result from file number: " << fileCounter << " named: " << fileItr->c_str() << std::endl;
      thisResult = dynamic_cast<HypoTestResult*>(thisTFile->Get("HypoTestCalculator_result"));
      r->Append(thisResult);
      delete thisResult;
    }
    ++fileCounter;
  }

  cout << "Done merging " << fileCounter << " files, writing file " << signalName+"_poi"+floatToString(poi)+"_doLimits_combined.root"
    << endl;
  r->Print();
  string outputPath = outputDir;
  outputPath+="_";
  outputPath+=signalName;
  outputPath+="_poi";
  outputPath+=floatToString(poi);
  outputPath+="_combined.root";
  TFile* thisSignalCombinedTFile = new TFile(outputPath.c_str(),"recreate");
  thisSignalCombinedTFile->cd();
  r->SetName("result");
  r->Write();
  thisSignalCombinedTFile->Close();

  return 0;
}
