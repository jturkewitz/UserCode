#include "TFile.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooPlot.h"


#include <set>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace RooFit;
using namespace std;

//
std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
std::string floatToString(float num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
std::string doubleToString(double num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
void setHistAttributes(TH2F* hist, int colorCounter, int markerStyleCounter)
{
  hist->SetLineColor(colorCounter);
  hist->SetMarkerColor(colorCounter);
  hist->SetMarkerStyle(markerStyleCounter);
  hist->SetMarkerSize(0.5);
}


int rooDatasetSignalPlots(bool doIhVsP = false)
{

  cout << "Run using the + option in root, and root 5.32: rooDatasetSignalPlots.C+" << endl;

  vector<string> signalNames;
  vector<string> fileNames;
  const int numModels = 4;
  RooDataSet* signalCandidatesRooDataSets[numModels];
  RooDataSet* signalCandidatesEtaCutRooDataSets[numModels];
  signalNames.push_back("Gluino800");
  signalNames.push_back("Gluino600");
  signalNames.push_back("Gluino400");
  signalNames.push_back("Gluino300");
  fileNames.push_back("/data/whybee1a/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_Signals_absEta_ptErrorPresel_Feb21/outputs/makeHSCParticlePlots_Feb21_Gluino800.root");
  fileNames.push_back("/data/whybee1a/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_Signals_absEta_ptErrorPresel_Feb21/outputs/makeHSCParticlePlots_Feb21_Gluino600.root");
  fileNames.push_back("/data/whybee1a/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_Signals_absEta_ptErrorPresel_Feb21/outputs/makeHSCParticlePlots_Feb21_Gluino400.root");
  fileNames.push_back("/data/whybee1a/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_Signals_absEta_ptErrorPresel_Feb21/outputs/makeHSCParticlePlots_Feb21_Gluino300.root");

  //RooRealVar rooVarNumGenHSCPEvents("rooVarNumGenHSCPEvents","numGenHSCPEvents",0,5e6);
  //RooRealVar rooVarNumGenHSCPTracks("rooVarNumGenHSCPTracks","numGenHSCPTracks",0,5e6);
  //RooRealVar rooVarNumGenChargedHSCPTracks("rooVarNumGenChargedHSCPTracks","numGenChargedHSCPTracks",0,5e6);
  //RooRealVar rooVarSignalEventCrossSection("rooVarSignalEventCrossSection","signalEventCrossSection",0,100); // pb
  //// get roodataset with gen information from signal file
  //RooDataSet* rooDataSetGenSignal = (RooDataSet*) _file0->Get("rooDataSetGenHSCPTracks");
  //if(!rooDataSetGenSignal)
  //{
  //  cout << "Problem with RooDataSet named rooDataSetGenSignal in signal file " << endl;
  //  return -3;
  //}

  //const RooArgSet* genArgSet = rooDataSetGenSignal->get();
  //RooRealVar* numGenHSCPTracksRooVar = (RooRealVar*)genArgSet->find(rooVarNumGenHSCPTracks.GetName());
  //RooRealVar* numGenHSCPEventsRooVar = (RooRealVar*)genArgSet->find(rooVarNumGenHSCPEvents.GetName());
  //rooDataSetGenSignal->get(0);
  //cout << " total signal tracks gen: " <<  numGenHSCPTracksRooVar->getVal()
  //  << " total signal events gen: " << numGenHSCPEventsRooVar->getVal()
  // << endl; 

  // define observables
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",-2.5,2.5);
  RooRealVar rooVarRun("rooVarRun","run",0,4294967295);
  RooRealVar rooVarLumiSection("rooVarLumiSection","lumiSection",0,4294967295);
  RooRealVar rooVarEvent("rooVarEvent","event",0,4294967295);
  
  int fileNameCount = 0;
  TFile* tFiles[numModels];
  for(vector<string>::const_iterator fileNameItr = fileNames.begin();
      fileNameItr != fileNames.end(); fileNameItr++)
  {
    tFiles[fileNameCount] = new TFile(fileNameItr->c_str(),"r");
    // get roodataset from signal file
    signalCandidatesRooDataSets[fileNameCount] = (RooDataSet*)tFiles[fileNameCount] ->Get("rooDataSetCandidates");
    //RooDataSet* rooDataSetAllSignal = (RooDataSet*) _file0->Get("rooDataSetOneCandidatePerEvent");
    if(!signalCandidatesRooDataSets[fileNameCount])
    {
      cout << "Problem with RooDataSet named rooDataSetCandidates in file: " << endl;
      return -3;
    }
    string etaCutString = "rooVarEta<1.5&&rooVarEta>-1.5";
    signalCandidatesEtaCutRooDataSets[fileNameCount] =
      (RooDataSet*) signalCandidatesRooDataSets[fileNameCount]->reduce(Cut(etaCutString.c_str()));
    fileNameCount++;
  }

  // dE/dx calib
  double dEdx_k = 2.529;
  double dEdx_c = 2.772;

  //string pSearchCutString = "rooVarP>100";
  ////string ptSearchCutString = "rooVarPt>50&&rooVarPt<65";
  //string ptSearchCutString = "rooVarPt>50";
  //string ptSidebandCutString = "rooVarPt<50";
  //string iasSearchCutString = "rooVarIas>0.1";
  //string iasSidebandCutString = "rooVarIas<0.1";
  //string ihSearchCutString = "rooVarIh>3";
  //string etaCutString = "rooVarEta<1.5&&rooVarEta>-1.5";
  //string ptMinCutString = "rooVarPt>45";
  //cout << "Thresholds: ";
  //if(usePtForSideband)
  //{
  //  cout << ptSearchCutString;
  //}
  //else
  //  cout << pSearchCutString;
  //cout << " ias > " << iasCut <<
  //  " mass > " << massCut <<
  //  //" Ih > 3.5 MeV/cm" << 
  //cout << iasSearchCutString << 
  cout <<  " |eta| < 1.5" << endl;
  
  //RooDataSet* ihDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(ihSearchCutString.c_str()));
  //RooDataSet* ptCutDataSet = (RooDataSet*) ihDataSet->reduce(Cut(ptMinCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) ptCutDataSet->reduce(Cut(etaCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(etaCutString.c_str()));

  //RooDataSet* iasSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(iasSidebandCutString.c_str()));
  //RooDataSet* aRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSidebandCutString.c_str()));
  //RooDataSet* cRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSearchCutString.c_str()));
  //RooDataSet* ptSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(ptSidebandCutString.c_str()));
  //RooDataSet* bRegionDataSet = (RooDataSet*) ptSidebandDataSet->reduce(Cut(iasSearchCutString.c_str()));
  //RooDataSet* cRegionPCutDataSet = (RooDataSet*) cRegionDataSet->reduce(Cut("rooVarP<1200"));
  //RooDataSet* ptSearchDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(ptSearchCutString.c_str()));
  //RooDataSet* dRegionDataSet = (RooDataSet*) ptSearchDataSet->reduce(Cut(iasSearchCutString.c_str()));
  //RooDataSet* dRegionPCutDataSet = (RooDataSet*) dRegionDataSet->reduce(Cut("rooVarP<1200"));

  //std::cout << "Using file: " << _file0->GetName() << " and dataset: " << rooDataSetAllSignal->GetName() << std::endl;
  //std::cout << "Entries in dataset: " << rooDataSetAllSignal->numEntries() << std::endl;
  ////std::cout << "Entries in dataset with ih > 3: " << ihDataSet->numEntries() << std::endl;
  ////std::cout << "Entries in dataset with ih > 3 and pt > 45: " << ptCutDataSet->numEntries() << std::endl;
  ////std::cout << "Entries in dataset with ih > 3 and pt > 45 and eta < 1.5: " << etaCutDataSet->numEntries() << std::endl;
  //std::cout << "Entries in dataset with |eta| < 1.5: " << etaCutDataSet->numEntries() << std::endl;
  //std::cout << "Region\tCuts\t\t\t\t\tEntries" <<
  //  std::endl << "A: " << ptSidebandCutString << "(+)" << iasSidebandCutString << "\t\t\t\t" << aRegionDataSet->numEntries() << 
  //  std::endl << "B: " << ptSidebandCutString << "(+)" << iasSearchCutString << "\t\t\t\t" << bRegionDataSet->numEntries() << 
  //  std::endl << "C: " << ptSearchCutString << "(+)" << iasSidebandCutString << "\t\t" << cRegionDataSet->numEntries() << 
  //  std::endl << "C with P < 1200 GeV\t\t\t\t\t\t" << cRegionPCutDataSet->numEntries() << 
  //  std::endl;

  //double bRegionErr = sqrt(bRegionDataSet->numEntries());
  //double cRegionErr = sqrt(cRegionDataSet->numEntries());
  //double aRegionErr = sqrt(aRegionDataSet->numEntries());
  //double bcOverA = bRegionDataSet->numEntries()*cRegionDataSet->numEntries()/(double)aRegionDataSet->numEntries();
  //double bcOverAErr = bcOverA*sqrt(pow(bRegionErr/bRegionDataSet->numEntries(),2)+pow(cRegionErr/cRegionDataSet->numEntries(),2)+pow(aRegionErr/aRegionDataSet->numEntries(),2));
  //std::cout << "BC/A = " << bcOverA << " +/- " << bcOverAErr << std::endl;

  // plots
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  gStyle->SetTitleFontSize(0.025);
  gStyle->SetTitleX(0.6);
  gStyle->SetTitleY(0.975);
  rooVarIh.setBins(140);
  rooVarIas.setBins(100);
  rooVarP.setBins(200);

  //TH1F* ihVsPHist = new TH1F("ihVsP",";P [GeV];Ih [MeV/cm]",80,0,15);
  //ihVsPHist->Draw();
  TLegend ihLegend(0.7,0.7,0.9,0.9);
  TH2F* ihVsPHists[numModels];
  int colorCounter = 1;
  int markerStyleCounter = 20;
  for(int i = 0; i < numModels; i++)
  {
    cout << "NumEntries this dataset (" << signalNames[i] << ") : " << signalCandidatesEtaCutRooDataSets[i]->numEntries()
         << endl;

    ihVsPHists[i] = signalCandidatesEtaCutRooDataSets[i]->createHistogram(rooVarP,rooVarIh,
        "",("ihVsP"+signalNames[i]).c_str());
    if(colorCounter==5)
      colorCounter+=2;
    setHistAttributes(ihVsPHists[i],colorCounter,markerStyleCounter);


    if(i==0)
      ihVsPHists[i]->SetTitle(";P [GeV];Ih [MeV/cm]");
    //string histName = "ihVsP_gluinos";
    //histName+=intToString(i);
    //histName+=".C";
    //t.Print(histName.c_str());
    ++colorCounter;
    //--markerStyleCounter;
    //cout << "Tlegend rows: " << ihLegend->GetNRows() << endl;
  }

  if(doIhVsP)
  {
    TCanvas t;
    t.cd();
    for(int i = 0; i < numModels; i++)
    {
      if(i==0)
      {
        ihVsPHists[i]->GetXaxis()->SetRangeUser(0,1500);
        ihVsPHists[i]->Draw();
      }
      else
        ihVsPHists[i]->Draw("same");

      ihLegend.AddEntry(ihVsPHists[i],signalNames[i].c_str(),"p");
    }

    ihLegend.SetBorderSize(0);
    ihLegend.Draw();
    t.Print("ihVsP_gluinos.png");
    t.Print("ihVsP_gluinos.eps");
  }

  // mass
  TH1F* massHists[numModels];
  TLegend massLegend(0.7,0.7,0.9,0.9);
  const RooArgSet* argSets[numModels];
  RooRealVar* ihData[numModels];
  RooRealVar* pData[numModels];
  for(int i = 0; i < numModels; i++)
  {
    argSets[i] = signalCandidatesEtaCutRooDataSets[i]->get();
    ihData[i] = (RooRealVar*)argSets[i]->find(rooVarIh.GetName());
    pData[i] = (RooRealVar*)argSets[i]->find(rooVarP.GetName());
    massHists[i] = new TH1F(("massHist"+signalNames[i]).c_str(),";Mass [GeV];arbitrary units",200,0,2000);
    // loop over entries
    for(int index=0; index < signalCandidatesEtaCutRooDataSets[i]->numEntries(); ++index)
    {
      signalCandidatesEtaCutRooDataSets[i]->get(index);
      float ih = ihData[i]->getVal();
      float trackP = pData[i]->getVal();
      float massSqr = (ih-dEdx_c)*pow(trackP,2)/dEdx_k;
      if(massSqr >= 0)
        massHists[i]->Fill(sqrt(massSqr));
    }
  }

  TCanvas t;
  t.cd();
  colorCounter = 1;
  markerStyleCounter = 23;
  for(int i = 0; i < numModels; i++)
  {
    if(colorCounter==5)
      colorCounter+=3;
    if(colorCounter==3)
      colorCounter++;

    massHists[i]->SetLineColor(colorCounter);
    //massHists[i]->SetMarkerColor(colorCounter);
    //massHists[i]->SetMarkerStyle(markerStyleCounter);
    //massHists[i]->SetMarkerSize(0.5);
    if(i==0)
    {
      massHists[i]->GetYaxis()->SetRangeUser(5e-1,5e3);
      massHists[i]->GetXaxis()->SetRangeUser(0,1500);
      massHists[i]->Draw();
    }
    else
      massHists[i]->Draw("same");

    massLegend.AddEntry(massHists[i],signalNames[i].c_str(),"l");
    ++colorCounter;
    markerStyleCounter--;
  }

  massLegend.SetBorderSize(0);
  massLegend.Draw();
  t.SetLogy();
  t.Print("mass_gluinos.png");
  t.Print("ihVsP_gluinos.eps");


  return 0;

}


