#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TF1.h"
#include "TProfile.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "HFitInterface.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>


std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}


// ********* MAIN
int main(int argc, char** argv)
{
  using namespace std;

  char* inputFile = argv[1];
  if(!inputFile)
  {
    cout << "Error: Missing input file." << endl;
    return -1;
  }

  gROOT->Macro("../rootstyle.C");

  TFile* inputTFile = new TFile(inputFile);
  inputTFile->cd();

  TFile* outputTFile = new TFile("plotDeDxSingleHitEtaSlices.root","recreate");

  //int nPBins = 200;
  const int nEtaBins = 25;

  TH1F* dedxNomCombined320Hists[nEtaBins];
  TH1F* dedxNomCombined500Hists[nEtaBins];
  TH1F* dedxNomCombinedAllHists[nEtaBins];

  string dir = "dedxHistsInNomBins";
  inputTFile->cd();
  string name320 = "dedxSingleHit320";
  string fullpath320 = dir;
  fullpath320+="/";
  fullpath320+=name320;
  TH3F* rawHist320 = (TH3F*) inputTFile->Get(fullpath320.c_str());
  string name500 = "dedxSingleHit500";
  string fullpath500 = dir;
  fullpath500+="/";
  fullpath500+=name500;
  TH3F* rawHist500 = (TH3F*) inputTFile->Get(fullpath500.c_str());
  string nameAll = "dedxSingleHitAll";
  string fullpathAll = dir;
  fullpathAll+="/";
  fullpathAll+=nameAll;
  TH3F* rawHistAll = (TH3F*) inputTFile->Get(fullpathAll.c_str());

  // loop over eta bins
  for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
  {
    rawHist320->GetZaxis()->SetRange(etaBin,etaBin);
    //XXX Set P range to low P sideband (0-->100 GeV)
    //rawHist->GetYaxis()->SetRange(1,20);
    string projName = "x";
    projName+="_EtaBin";
    projName+=intToString(etaBin);
    dedxNomCombined320Hists[etaBin-1] = (TH1F*)rawHist320->Project3D(projName.c_str());

    rawHist500->GetZaxis()->SetRange(etaBin,etaBin);
    //XXX Set P range to low P sideband (0-->100 GeV)
    //rawHist->GetYaxis()->SetRange(1,20);
    string projName2 = "x";
    projName2+="_EtaBin";
    projName2+=intToString(etaBin);
    dedxNomCombined500Hists[etaBin-1] = (TH1F*)rawHist500->Project3D(projName2.c_str());

    rawHistAll->GetZaxis()->SetRange(etaBin,etaBin);
    //XXX Set P range to low P sideband (0-->100 GeV)
    //rawHist->GetYaxis()->SetRange(1,20);
    string projName3 = "x";
    projName3+="_EtaBin";
    projName3+=intToString(etaBin);
    dedxNomCombinedAllHists[etaBin-1] = (TH1F*)rawHistAll->Project3D(projName3.c_str());
  }


  outputTFile->cd();
  for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
  {
    dedxNomCombined500Hists[etaBin-1]->Write();
    dedxNomCombined320Hists[etaBin-1]->Write();
    dedxNomCombinedAllHists[etaBin-1]->Write();
  }

}

