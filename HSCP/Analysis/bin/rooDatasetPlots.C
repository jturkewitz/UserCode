#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooPlot.h"

#include "TFile.h"

#include <set>
#include <iostream>

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
void setHistAttributes(TH1* hist, int colorCounter, int markerStyleCounter)
{
  hist->SetLineColor(colorCounter);
  hist->SetMarkerColor(colorCounter);
  hist->SetMarkerStyle(markerStyleCounter);
  if(markerStyleCounter > 21)
    hist->SetMarkerSize(1.0);
}


rooDatasetPlots(bool doNoMplots = false, bool doEtaPlots = false, bool doPtPlots = false, bool doPtVsIasPlot = true)
{
  bool usePtForSideband = true;
  float iasCut = 0.4;

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
  RooRealVar rooVarEta("rooVarEta","eta",-2.5,2.5);
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


  string pSearchCutString = "rooVarP>100";
  //string ptSearchCutString = "rooVarPt>50&&rooVarPt<65";
  string ptSearchCutString = "rooVarPt>50";
  string ptSidebandCutString = "rooVarPt<50";
  string iasSearchCutString = "rooVarIas>0.4";
  string iasSidebandCutString = "rooVarIas<0.4";
  string ihSearchCutString = "rooVarIh>3";
  string etaCutString = "rooVarEta<1.5&&rooVarEta>-1.5";
  string ptMinCutString = "rooVarPt>45";
  cout << "Thresholds: ";
  if(usePtForSideband)
  {
    cout << ptSearchCutString;
  }
  else
    cout << pSearchCutString;
  //cout << " ias > " << iasCut <<
  //  " mass > " << massCut <<
  //  //" Ih > 3.5 MeV/cm" << 
  cout << iasSearchCutString << 
    " |eta| < 1.5" << endl;
  
  //RooDataSet* ihDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(ihSearchCutString.c_str()));
  //RooDataSet* ptCutDataSet = (RooDataSet*) ihDataSet->reduce(Cut(ptMinCutString.c_str()));
  //RooDataSet* etaCutDataSet = (RooDataSet*) ptCutDataSet->reduce(Cut(etaCutString.c_str()));
  RooDataSet* etaCutDataSet = (RooDataSet*) rooDataSetAllSignal->reduce(Cut(etaCutString.c_str()));

  RooDataSet* iasSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(iasSidebandCutString.c_str()));
  //RooDataSet* aRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSidebandCutString.c_str()));
  //RooDataSet* cRegionDataSet = (RooDataSet*) iasSidebandDataSet->reduce(Cut(ptSearchCutString.c_str()));
  RooDataSet* ptSidebandDataSet = (RooDataSet*) etaCutDataSet->reduce(Cut(ptSidebandCutString.c_str()));
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
  rooVarPt.setBins(1000);
  // plot of Ih in the AB (pt sideband) region
  if(doNoMplots)
  {
    int nomMin_=5;
    int nomMax_=20;
    TH1F* ihInNoMSlices = new TH1F("ihInNoMSlices",
        "CMS Preliminary    #sqrt{s} = 7 TeV;Ih [MeV/cm];arbitrary units",80,0,15);
    TLegend* ihLegend = new TLegend(0.13,0.72,0.51,0.92);
    TH1F* iasInNoMSlices = new TH1F("iasInNoMSlices",
        "CMS Preliminary    #sqrt{s} = 7 TeV;Ias;arbitrary units",100,0,1);
    TLegend* ihLegend = new TLegend(0.13,0.72,0.51,0.92);
    TLegend* iasLegend = new TLegend(0.65,0.7,0.92,0.92);
    TCanvas c1;
    c1.cd();
    ihInNoMSlices->Draw();
    TCanvas c2;
    c2.cd();
    iasInNoMSlices->Draw();
    TH1* ihInNoMSlicesHists[9];
    TH1* iasInNoMSlicesHists[9];
    int nomSlice = -1;
    int colorCounter = 1;
    int markerStyleCounter = 23;
    // loop over NoM 
    for(int nom = nomMin_; nom < nomMax_; nom+=2)
    {
      ++nomSlice;
      if(nomSlice % 2 != 0) continue;
      //std::cout << "Doing slice: " << nomSlice+etaSlice << std::endl;
      std::string nomCutString = "(rooVarNoMias==";
      nomCutString+=intToString(nom);
      nomCutString+="||rooVarNoMias==";
      nomCutString+=intToString(nom+1);
      nomCutString+=")";
      std::string thisCut = nomCutString;
      RooDataSet* nomPtSBCutDataSet = (RooDataSet*) ptSidebandDataSet->reduce(thisCut.c_str());
      //std::cout << "cut: " << thisCut << " entries: " << nomEtaPtSBCutDataSet->numEntries() << std::endl;
      ihInNoMSlicesHists[nomSlice] = nomPtSBCutDataSet->createHistogram(("ihHistNoM"+intToString(nom)).c_str(),rooVarIh);
      iasInNoMSlicesHists[nomSlice] = nomPtSBCutDataSet->createHistogram(("iasHistNoM"+intToString(nom)).c_str(),rooVarIas);
      if(colorCounter==5)
        colorCounter+=2;
      setHistAttributes(ihInNoMSlicesHists[nomSlice],colorCounter,markerStyleCounter);
      setHistAttributes(iasInNoMSlicesHists[nomSlice],colorCounter,markerStyleCounter);
      c1.cd();
      ihInNoMSlicesHists[nomSlice]->DrawNormalized("same");
      c2.cd();
      iasInNoMSlicesHists[nomSlice]->DrawNormalized("same");
      std::string legString = "NoM ";
      legString+=intToString(nom);
      legString+="-";
      legString+=intToString(nom+1);
      ihLegend->AddEntry(ihInNoMSlicesHists[nomSlice],legString.c_str(),"p");
      iasLegend->AddEntry(iasInNoMSlicesHists[nomSlice],legString.c_str(),"p");
      delete nomPtSBCutDataSet;
      ++colorCounter;
      markerStyleCounter--;
    }

    c1.cd();
    ihLegend->SetBorderSize(0);
    ihLegend->Draw("same");
    c1.SetLogy();
    ihInNoMSlices->GetXaxis()->SetRangeUser(0.6,5.2);
    ihInNoMSlices->GetYaxis()->SetRangeUser(1e-6,5e-1);
    c1.Print("ihInNoMSlices.png");
    c1.Print("ihInNoMSlices.eps");
    c2.cd();
    iasLegend->SetBorderSize(0);
    iasLegend->Draw("same");
    c2.SetLogy();
    iasInNoMSlices->GetYaxis()->SetRangeUser(1e-6,5e-1);
    c2.Print("iasInNoMSlices.png");
    c2.Print("iasInNoMSlices.eps");

  } // if doNoMplots



  if(doEtaPlots)
  {
    // do plots of ih/ias in eta slices
    float etaMin_=0.0;
    float etaMax_=1.3;
    TH1F* ihInEtaSlices = new TH1F("ihInEtaSlices",
        "CMS Preliminary    #sqrt{s} = 7 TeV;Ih [MeV/cm];arbitrary units",80,0,15);
    TLegend* ihLegend = new TLegend(0.13,0.72,0.51,0.92);
    TH1F* iasInEtaSlices = new TH1F("iasInEtaSlices",
        "CMS Preliminary    #sqrt{s} = 7 TeV;Ias;arbitrary units",100,0,1);
    TLegend* ihLegend = new TLegend(0.13,0.72,0.51,0.92);
    TLegend* iasLegend = new TLegend(0.55,0.7,0.9,0.92);
    TCanvas c1;
    c1.cd();
    ihInEtaSlices->Draw();
    TCanvas c2;
    c2.cd();
    iasInEtaSlices->Draw();
    TH1* ihInEtaSlicesHists[9];
    TH1* iasInEtaSlicesHists[9];
    int etaSlice = -1;
    int colorCounter = 1;
    int markerStyleCounter = 23;
    for(float lowerEta = etaMin_; lowerEta < etaMax_; lowerEta+=0.2)
    {
      ++etaSlice;
      if(etaSlice % 2 != 0) continue;
      std::string etaCutString = "(rooVarEta>";
      etaCutString+=floatToString(lowerEta);
      etaCutString+="&&rooVarEta<";
      std::string upperEtaLimit;
      upperEtaLimit=floatToString(lowerEta+0.2);
      etaCutString+=upperEtaLimit;
      etaCutString+=")||(rooVarEta>-";
      etaCutString+=upperEtaLimit;
      etaCutString+="&&rooVarEta<-";
      etaCutString+=floatToString(lowerEta);
      etaCutString+=")";

      std::string thisCut = etaCutString;
      RooDataSet* etaPtSBCutDataSet = (RooDataSet*) ptSidebandDataSet->reduce(thisCut.c_str());
      //std::cout << "cut: " << thisCut << " entries: " << nomEtaPtSBCutDataSet->numEntries() << std::endl;
      ihInEtaSlicesHists[etaSlice] = etaPtSBCutDataSet->createHistogram(
          ("ihHistEta"+floatToString(lowerEta)).c_str(),rooVarIh);
      iasInEtaSlicesHists[etaSlice] = etaPtSBCutDataSet->createHistogram(
          ("iasHistEta"+floatToString(lowerEta)).c_str(),rooVarIas);
      if(colorCounter==5)
        colorCounter+=2;
      setHistAttributes(ihInEtaSlicesHists[etaSlice],colorCounter,markerStyleCounter);
      setHistAttributes(iasInEtaSlicesHists[etaSlice],colorCounter,markerStyleCounter);
      c1.cd();
      ihInEtaSlicesHists[etaSlice]->DrawNormalized("same");
      c2.cd();
      iasInEtaSlicesHists[etaSlice]->DrawNormalized("same");
      std::string legString = floatToString(lowerEta);
      legString+=" < |#eta| < ";
      legString+=floatToString(lowerEta+0.2);
      ihLegend->AddEntry(ihInEtaSlicesHists[etaSlice],legString.c_str(),"p");
      iasLegend->AddEntry(iasInEtaSlicesHists[etaSlice],legString.c_str(),"p");
      delete etaPtSBCutDataSet;
      ++colorCounter;
      markerStyleCounter--;
    }

    c1.cd();
    ihLegend->SetBorderSize(0);
    ihLegend->Draw("same");
    c1.SetLogy();
    ihInEtaSlices->GetXaxis()->SetRangeUser(0.6,5.2);
    ihInEtaSlices->GetYaxis()->SetRangeUser(1e-6,5e-1);
    c1.Print("ihInEtaSlices.png");
    c1.Print("ihInEtaSlices.eps");
    c2.cd();
    iasLegend->SetBorderSize(0);
    iasLegend->Draw("same");
    c2.SetLogy();
    iasInEtaSlices->GetXaxis()->SetRangeUser(0,1);
    iasInEtaSlices->GetYaxis()->SetRangeUser(1e-6,5e-1);
    c2.Print("iasInEtaSlices.png");
    c2.Print("iasInEtaSlices.eps");

  } // ifdoEtaPlots

  if(doPtPlots)
  {
    // do plots of pt in ias slices
    float iasMin=0.0;
    float iasMax=0.4;
    TH1F* ptInIasSlices = new TH1F("ptInIasSlices",
        "CMS Preliminary    #sqrt{s} = 7 TeV;Pt [GeV];arbitrary units",1000,0,5000);
    TLegend* iasLegend = new TLegend(0.55,0.7,0.9,0.92);
    TCanvas c1;
    c1.cd();
    ptInIasSlices->Draw();
    TH1* ptInIasSlicesHists[9];
    int iasSlice = -1;
    int colorCounter = 1;
    int markerStyleCounter = 23;
    for(float lowerIas = iasMin; lowerIas < iasMax; lowerIas+=0.1)
    {
      ++iasSlice;
      std::string iasCutString = "rooVarIas>";
      iasCutString+=floatToString(lowerIas);
      iasCutString+="&&rooVarIas<";
      std::string upperIasLimit;
      upperIasLimit=floatToString(lowerIas+0.1);
      iasCutString+=upperIasLimit;

      std::string thisCut = iasCutString;
      RooDataSet* cutDataSet = (RooDataSet*) iasSidebandDataSet->reduce(thisCut.c_str());
      //std::cout << "cut: " << thisCut << " entries: " << nomEtaPtSBCutDataSet->numEntries() << std::endl;
      ptInIasSlicesHists[iasSlice] = cutDataSet->createHistogram(
          ("ptHistIas"+floatToString(lowerIas)).c_str(),rooVarPt);
      if(colorCounter==5)
        colorCounter+=2;
      setHistAttributes(ptInIasSlicesHists[iasSlice],colorCounter,markerStyleCounter);
      c1.cd();
      ptInIasSlicesHists[iasSlice]->DrawNormalized("same");

      std::string legString = floatToString(lowerIas);
      legString+=" < Ias < ";
      legString+=floatToString(lowerIas+0.1);
      iasLegend->AddEntry(ptInIasSlicesHists[iasSlice],legString.c_str(),"p");
      delete cutDataSet;
      ++colorCounter;
      markerStyleCounter--;
    }

    c1.cd();
    iasLegend->SetBorderSize(0);
    iasLegend->Draw("same");
    c1.SetLogy();
    ptInIasSlices->GetXaxis()->SetRangeUser(0,1000);
    ptInIasSlices->GetYaxis()->SetRangeUser(1e-8,1);
    c1.Print("ptInIasSlices.png");
    c1.Print("ptInIasSlices.eps");

  } // ifdoPtPlots

  if(doPtVsIasPlot)
  {
    TCanvas c1;
    c1.cd();
    TH2F* ptVsIasHist = etaCutDataSet->createHistogram(rooVarIas,rooVarPt);
    ptVsIasHist->GetYaxis()->SetRangeUser(0,500);
    ptVsIasHist->GetXaxis()->SetRangeUser(0,0.5);
    ptVsIasHist->SetTitle("CMS Preliminary    #sqrt{s} = 7 TeV;Ias;Pt [GeV]");
    ptVsIasHist->Draw("colz");
    c1.SetLogz();
    c1.Print("ptVsIas.png");
    c1.Print("ptVsIas.eps");

  } // if doPtVsIasPlot

  // overall dists?



}
