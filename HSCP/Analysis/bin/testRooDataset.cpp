#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

#include <set>
#include <iostream>

using namespace RooFit;
using namespace std;

struct EventInfo
{
  double runNumber;
  double eventNumber;
  double lumiSection;

  EventInfo()
  {
    runNumber = 0;
    eventNumber = 0;
    lumiSection = 0;
  }

  EventInfo(double runNum, double lumiNum, double eventNum) :
    runNumber(runNum),  eventNumber(eventNum), lumiSection(lumiNum)
  {
  }

  bool operator==(const EventInfo& comp) const
  {
    return (comp.runNumber==runNumber) && (comp.eventNumber==eventNumber)
      && (comp.lumiSection==lumiSection);
  }

  bool operator<(const EventInfo& comp) const
  {
    // run number < ; run numbers equal but lumi < ; run, lumis equal but event <
    return (runNumber<comp.runNumber) ||
      (runNumber==comp.runNumber && lumiSection < comp.lumiSection) ||
      (runNumber==comp.runNumber && lumiSection==comp.lumiSection && eventNumber < comp.eventNumber);
  }
};

double getRatioError(int entries1, int entries2)
{
  double err1 = sqrt(entries1);
  double err2 = sqrt(entries2);
  return (entries1/(double)entries2)*sqrt(pow(err1/entries1,2)+pow(err2/entries2,2));
}

int main()
{
  bool usePtForSideband = true;
  float iasCut = 0.4;
  float massCut = 190;
  //string fileName = "/local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_feb1/outputs/makeHSCParticlePlots_feb1_Data2011_all.root";
  //string fileName = "/local/cms/user/cooper/data/hscp/428/makeHSCParticlePlotsOutput/MakeHSCParticlePlots_2011data_jan13/outputs/makeHSCParticlePlots_jan13_Data2011.root";
  //string fileName = "/local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_dataWithTightRPCTrig_feb2/outputs/makeHSCParticlePlots_feb1_Data2011_all.root";
  string fileName = "/local/cms/user/cooper/cmssw/428/HSCPAnalysis/src/HSCP2011/ToyMC/bin/FARM_MakeHSCParticlePlots_dataWithTightRPCTrig_absEta_ptErrorPresel_feb6/outputs/makeHSCParticlePlots_feb1_Data2011_all.root";
  TFile* _file0 = TFile::Open(fileName.c_str());

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
  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",0,2.5);
  RooRealVar rooVarRun("rooVarRun","run",0,4294967295);
  RooRealVar rooVarLumiSection("rooVarLumiSection","lumiSection",0,4294967295);
  RooRealVar rooVarEvent("rooVarEvent","event",0,4294967295);
  // get roodataset from signal file
  RooDataSet* rooDataSetAllSignal = (RooDataSet*) _file0->Get("rooDataSetCandidates");
  //RooDataSet* rooDataSetAllSignal = (RooDataSet*) _file0->Get("rooDataSetOneCandidatePerEvent");
  if(!rooDataSetAllSignal)
  {
    cout << "Problem with RooDataSet named rooDataSetCandidates in file: " << endl;
    return -3;
  }

  // dE/dx calib
  double dEdx_k = 2.529;
  double dEdx_c = 2.772;

  //// construct D region for signal
  //cout << "Thresholds: ";
  //if(usePtForSideband)
  //{
  //  cout << "pt > 70";
  //}
  //else
  //  cout << "p > 100";
  //cout << " ias > " << iasCut <<
  //  " mass > " << massCut <<
  //  //" Ih > 3.5 MeV/cm" << 
  //  " |eta| < 1.5" << endl;

  //string pSearchCutString = "rooVarP>100";
  //string ptSearchCutString = "rooVarPt>50&&rooVarPt<65";
  string ptSearchCutString = "rooVarPt>50";
  string ptSidebandCutString = "rooVarPt<50";
  string iasSearchCutString = "rooVarIas>0.1";
  string iasSidebandCutString = "rooVarIas<0.1";
  string ihSearchCutString = "rooVarIh>3";
  string etaOverallCutString = "rooVarEta<1.5&&rooVarEta>-1.5";
  string ptMinCutString = "rooVarPt>45";
  
  //RooDataSet* ihDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(ihSearchCutString.c_str()));
  //RooDataSet* ptCutDataSet = (RooDataSet*) ihDataSet->reduce(Cut(ptMinCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) ptCutDataSet->reduce(Cut(etaCutString.c_str()));
  std::string nomCutString = "rooVarNoMias==9";
  nomCutString+="||rooVarNoMias==10";
  std::string etaCutString = "(rooVarEta>1.0";
  etaCutString+="&&rooVarEta<1.2";
  etaCutString+=")||(rooVarEta>-1.2";
  etaCutString+="&&rooVarEta<-1.0)";

  RooDataSet* etaCutDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(etaCutString.c_str()));
  RooDataSet* etaCutNoMCutDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(nomCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) etaCutNoMCutDataSet->reduce(Cut(etaOverallCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(etaOverallCutString.c_str()));

  RooDataSet* iasSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(iasSidebandCutString.c_str()));
  RooDataSet* aRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSidebandCutString.c_str()));
  RooDataSet* cRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSearchCutString.c_str()));
  RooDataSet* ptSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(ptSidebandCutString.c_str()));
  RooDataSet* bRegionDataSet = (RooDataSet*) ptSidebandDataSet->reduce(Cut(iasSearchCutString.c_str()));
  RooDataSet* cRegionPCutDataSet = (RooDataSet*) cRegionDataSet->reduce(Cut("rooVarP<1200"));
  RooDataSet* ptSearchDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(ptSearchCutString.c_str()));
  RooDataSet* dRegionDataSet = (RooDataSet*) ptSearchDataSet->reduce(Cut(iasSearchCutString.c_str()));
  RooDataSet* dRegionPCutDataSet = (RooDataSet*) dRegionDataSet->reduce(Cut("rooVarP<1200"));

  std::cout << "Using file: " << _file0->GetName() << " and dataset: " << rooDataSetAllSignal->GetName() << std::endl;
  std::cout << "Entries in dataset: " << rooDataSetAllSignal->numEntries() << std::endl;
  //std::cout << "Entries in dataset with ih > 3: " << ihDataSet->numEntries() << std::endl;
  //std::cout << "Entries in dataset with ih > 3 and pt > 45: " << ptCutDataSet->numEntries() << std::endl;
  //std::cout << "Entries in dataset with ih > 3 and pt > 45 and eta < 1.5: " << etaCutDataSet->numEntries() << std::endl;
  std::cout << "Entries in dataset with |eta| < 1.5: " << etaCutDataSet->numEntries() << std::endl;
  std::cout << "Region\tCuts\t\t\t\t\tEntries" <<
    std::endl << "A: " << ptSidebandCutString << "(+)" << iasSidebandCutString << "\t\t\t\t" << aRegionDataSet->numEntries() << 
    std::endl << "B: " << ptSidebandCutString << "(+)" << iasSearchCutString << "\t\t\t\t" << bRegionDataSet->numEntries() << 
    std::endl << "C: " << ptSearchCutString << "(+)" << iasSidebandCutString << "\t\t" << cRegionDataSet->numEntries() << 
    std::endl << "C with P < 1200 GeV\t\t\t\t\t\t" << cRegionPCutDataSet->numEntries() << 
    std::endl << "C/A\t\t\t\t\t\t" << cRegionDataSet->numEntries()/(double)aRegionDataSet->numEntries() << 
    std::endl;

  double bRegionErr = sqrt(bRegionDataSet->numEntries());
  double cRegionErr = sqrt(cRegionDataSet->numEntries());
  double aRegionErr = sqrt(aRegionDataSet->numEntries());
  double bcOverA = bRegionDataSet->numEntries()*cRegionDataSet->numEntries()/(double)aRegionDataSet->numEntries();
  double bcOverAErr = bcOverA*sqrt(pow(bRegionErr/bRegionDataSet->numEntries(),2)+pow(cRegionErr/cRegionDataSet->numEntries(),2)+pow(aRegionErr/aRegionDataSet->numEntries(),2));
  std::cout << "BC/A = " << bcOverA << " +/- " << bcOverAErr << std::endl;

  // make hists
  TH1F* iasInARegion = (TH1F*) aRegionDataSet->createHistogram("iasInARegion",rooVarIas);
  TCanvas t;
  t.cd();
  t.SetLogy();
  iasInARegion->Draw();
  t.Print("iasInARegion.png");
  TH1F* iasInCRegion = (TH1F*) cRegionDataSet->createHistogram("iasInCRegion",rooVarIas);
  iasInCRegion->Draw();
  t.Print("iasInCRegion.png");
  TH1F* ptInARegion = (TH1F*) aRegionDataSet->createHistogram("ptInARegion",rooVarPt);
  ptInARegion->Draw();
  t.Print("ptInARegion.png");
  TH1F* ptInCRegion = (TH1F*) cRegionDataSet->createHistogram("ptInCRegion",rooVarPt);
  ptInCRegion->Draw();
  t.Print("ptInCRegion.png");
  
  //std::cout << "D: " << ptSearchCutString << "(+)" << iasSearchCutString << "\t\t\t" << dRegionDataSet->numEntries() << std::endl;
  //std::cout << "D with P < 1200 GeV: " << ptSearchCutString << "(+)" << iasSearchCutString << "\t\t" << dRegionDataSet->numEntries() << std::endl;
  //// do ratios and errors
  //std::cout << std::endl;
  //std::cout << "Ratios and Errors" << std::endl;
  //double dOverB = dRegionDataSet->numEntries()/(double)bRegionDataSet->numEntries();
  //double dOverBError = getRatioError(dRegionDataSet->numEntries(),bRegionDataSet->numEntries());
  //std::cout << "D/B = " << dOverB << " +/- " << dOverBError << std::endl;
  //double cOverA = cRegionDataSet->numEntries()/(double)aRegionDataSet->numEntries();
  //double cOverAError = getRatioError(cRegionDataSet->numEntries(),aRegionDataSet->numEntries());
  //std::cout << "C/A = " << cOverA << " +/- " << cOverAError << std::endl;
  //std::cout << "-->with P < 1200 GeV cut" << std::endl;
  //double dOverBPCut = dRegionPCutDataSet->numEntries()/(double)bRegionDataSet->numEntries();
  //double dOverBPCutError = getRatioError(dRegionPCutDataSet->numEntries(),bRegionDataSet->numEntries());
  //std::cout << "D/B = " << dOverBPCut << " +/- " << dOverBPCutError << std::endl;
  //double cOverAPCut = cRegionPCutDataSet->numEntries()/(double)aRegionDataSet->numEntries();
  //double cOverAPCutError = getRatioError(cRegionPCutDataSet->numEntries(),aRegionDataSet->numEntries());
  //std::cout << "C/A = " << cOverAPCut << " +/- " << cOverAPCutError << std::endl;


  //RooDataSet* regionD1DataSetSignal;
  //if(usePtForSideband)
  //  regionD1DataSetSignal = (RooDataSet*)rooDataSetAllSignal->reduce(Cut(ptSearchCutString.c_str()));
  //else
  //  regionD1DataSetSignal = (RooDataSet*)rooDataSetAllSignal->reduce(Cut(pSearchCutString.c_str()));
  //RooDataSet* ihCutRegionDDataSetSignal = (RooDataSet*)regionD1DataSetSignal->reduce(Cut(ihSearchCutString.c_str()));
  ////RooDataSet* regionDDataSetSignal = (RooDataSet*)ihCutRegionDDataSetSignal->reduce(Cut(etaCutString.c_str()));
  //// below for no Ih cut
  //RooDataSet* regionDDataSetSignal = (RooDataSet*)regionD1DataSetSignal->reduce(Cut(etaCutString.c_str()));

  //// construct signal prediction hist (with mass cut)
  //set<EventInfo> selectedSignalEventsSet;
  //int numSignalTracksAboveMassCut = 0;
  //int numSignalTracksAboveMassCutAndIasCut = 0;
  //int numSignalEventsOverIasCut = 0;
  //const RooArgSet* argSet = regionDDataSetSignal->get();
  //RooRealVar* iasData = (RooRealVar*)argSet->find(rooVarIas.GetName());
  //RooRealVar* ihData = (RooRealVar*)argSet->find(rooVarIh.GetName());
  //RooRealVar* pData = (RooRealVar*)argSet->find(rooVarP.GetName());
  //RooRealVar* eventNumData = (RooRealVar*)argSet->find(rooVarEvent.GetName());
  //RooRealVar* lumiSecData = (RooRealVar*)argSet->find(rooVarLumiSection.GetName());
  //RooRealVar* runNumData = (RooRealVar*)argSet->find(rooVarRun.GetName());
  //for(int evt=0; evt < regionDDataSetSignal->numEntries(); ++evt)
  //{
  //  regionDDataSetSignal->get(evt);
  //  // apply mass cut
  //  float massSqr = (ihData->getVal()-dEdx_c)*pow(pData->getVal(),2)/dEdx_k;
  //  if(massSqr < 0)
  //    continue;
  //  else if(sqrt(massSqr) >= massCut)
  //  {
  //    numSignalTracksAboveMassCut++;
  //    if(iasData->getVal() > iasCut)
  //    {
  //      numSignalTracksAboveMassCutAndIasCut++;
  //      // if track over ias cut, put the event in the set (one track per event kept)
  //      EventInfo evtInfo(runNumData->getVal(),lumiSecData->getVal(),eventNumData->getVal());
  //      pair<set<EventInfo>::iterator,bool> ret;
  //      ret = selectedSignalEventsSet.insert(evtInfo);
  //      if(ret.second)
  //        numSignalEventsOverIasCut++;
  //    }
  //  }
  //}

  //cout << "signal events with at least one track > mass,ias cuts: " << numSignalEventsOverIasCut
  //  << " = " << numSignalEventsOverIasCut/numGenHSCPEventsRooVar->getVal() << endl;

  ////cout << " total signal tracks gen: " <<  
  ////  << " total signal events gen: " << numGenHSCPEventsRooVar->getVal()

  return 0;
}
