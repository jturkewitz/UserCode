// -*- C++ -*-
//
// Package:    HSCP
// Class:      GenerateHSCPToyMC
// 
/**\class GenerateHSCPToyMC GenerateHSCPToyMC.cc HSCPValidation/GenerateHSCPToyMC/src/GenerateHSCPToyMC.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seth Cooper,27 1-024,+41227672342,
//         Created:  Wed Apr 14 14:27:52 CEST 2010
// $Id: GenerateHSCPToyMC.cc,v 1.3 2011/09/13 11:20:08 scooper Exp $
//
//


// system include files
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TCut.h"
#include "TGraphErrors.h"
#include "TPaveLabel.h"

#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooLandau.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooFFTConvPdf.h"
#include "RooBinning.h"
#include "RooMCStudy.h"
#include "RooHist.h"
#include "RooTFnBinding.h"
#include "RooTFnPdfBinding.h"
#include "RooRandom.h"
#include "RooMsgService.h"
#include "RooBinning.h"

#include "HSCP2011/ToyMC/interface/GenerateHSCPToyMC.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//
//edm::Service<TFileService> fileService;

//
// constructors and destructor
//
GenerateHSCPToyMC::GenerateHSCPToyMC(const edm::ParameterSet& iConfig) :
  genSignalAndBackground_ (iConfig.getParameter<bool>("GenerateSignalAndBackground")),
  signalFraction_ (iConfig.getParameter<double>("SignalFraction")),
  numberOfTrials_ (iConfig.getParameter<int>("NumberOfTrials")),
  numberOfEventsPerTrial_ (iConfig.getParameter<int>("EventsPerTrial"))
{
  // Generation params
  numSignalEvents_ = 0;
  numBackgroundEvents_ = numberOfEventsPerTrial_;
  if(genSignalAndBackground_)
  {
    numSignalEvents_ = round(numberOfEventsPerTrial_*signalFraction_);
    numBackgroundEvents_ = numberOfEventsPerTrial_-numSignalEvents_;
  }
  edm::LogInfo("GenerateHSCPToyMC") << "numberOfEventsPerTrial*signalFraction = " << numberOfEventsPerTrial_ << 
    " * " << signalFraction_ << " = " << numberOfEventsPerTrial_*signalFraction_ << std::endl;

}


GenerateHSCPToyMC::~GenerateHSCPToyMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GenerateHSCPToyMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace RooFit;
  using namespace std;
  SetEStyle();


  // ************************ Define observables *****************************************
  RooRealVar TrackdEdxE1("TrackdEdxE1","dEdxE1",3,25,"MeV/cm"); // consider events only in the search region
  RooRealVar Trackpt("Trackpt","pt",100,2000,"GeV");
  RooRealVar Trackp("Trackp","p",50,5000,"GeV");
  RooRealVar Tracketa("Tracketa","eta",-2.5,2.5);
  RooRealVar Tracknom("Tracknom","nom",1,50);
  RooRealVar sigFrac("sigFrac","sigFrac",0.01,0,1);

  TrackdEdxE1.setBins(88);
  double nomBinBounds[6] = {1,6,11,16,21,50};
  RooBinning nomBinning(5,nomBinBounds);
  Tracknom.setBinning(nomBinning);
  double pBinBounds[6] = {15,50,150,300,450,5000};
  RooBinning pBinning(5,pBinBounds);
  Trackp.setBinning(pBinning);

  // **************** Create GEN PDFs (use muon ToF sideband) ****************************
  muonGenInputFile_ = new TFile("1.may17.readHSCPTreeTKplots.muTimingSbPnom.mu2010abTrigSkim.root");
  hscpGenInputFile_ = new TFile("3.may11.readHSCPTreeTKplots.PnomSlices.gluino500TrigSkim.root");
  muonEtaVsPtHist_ = (TH2F*) muonGenInputFile_->Get("trackEtaVsPtMuTimingSB");
  hscpEtaVsPtHist_ = (TH2F*) hscpGenInputFile_->Get("trackEtaVsPtMatched");
  hscpDeDxVsPHist_ = (TH2F*) hscpGenInputFile_->Get("trackDeDxE1VsPMatched");
  // hists in slices of P and NoM -- P sideband
  muonGenDeDxPNomSliceHists_[0] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsIndex2");
  muonGenDeDxPNomSliceHists_[1] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsIndex3");
  muonGenDeDxPNomSliceHists_[2] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsIndex4");
  muonGenDeDxPNomSliceHists_[3] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsIndex5");
  // hists in slices of P and NoM -- mu SB sideband
  //muonGenDeDxPNomSliceHists_[0] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlicesMuTimingSB/histForPDFNoMPBinsMuSBIndex2");
  //muonGenDeDxPNomSliceHists_[1] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlicesMuTimingSB/histForPDFNoMPBinsMuSBIndex3");
  //muonGenDeDxPNomSliceHists_[2] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlicesMuTimingSB/histForPDFNoMPBinsMuSBIndex4");
  //muonGenDeDxPNomSliceHists_[3] = (TH1F*) muonGenInputFile_->Get("trackDeDxNoMAndPSlicesMuTimingSB/histForPDFNoMPBinsMuSBIndex5");
  // hscp
  hscpGenDeDxPNomSliceHists_[0] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex1");
  hscpGenDeDxPNomSliceHists_[1] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex2");
  hscpGenDeDxPNomSliceHists_[2] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex3");
  //hscpGenDeDxPNomSliceHists_[3] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex4");
  //hscpGenDeDxPNomSliceHists_[4] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex5");
  hscpGenDeDxPNomSliceHists_[5] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex6");
  hscpGenDeDxPNomSliceHists_[6] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex7");
  hscpGenDeDxPNomSliceHists_[7] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex8");
  hscpGenDeDxPNomSliceHists_[8] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex9");
  hscpGenDeDxPNomSliceHists_[9] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex10");
  hscpGenDeDxPNomSliceHists_[10] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex11");
  hscpGenDeDxPNomSliceHists_[11] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex12");
  hscpGenDeDxPNomSliceHists_[12] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex13");
  hscpGenDeDxPNomSliceHists_[13] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex14");
  hscpGenDeDxPNomSliceHists_[14] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex15");
  hscpGenDeDxPNomSliceHists_[15] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex16");
  hscpGenDeDxPNomSliceHists_[16] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex17");
  hscpGenDeDxPNomSliceHists_[17] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex18");
  hscpGenDeDxPNomSliceHists_[18] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex19");
  hscpGenDeDxPNomSliceHists_[19] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex20");
  hscpGenDeDxPNomSliceHists_[20] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex21");
  hscpGenDeDxPNomSliceHists_[21] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex22");
  hscpGenDeDxPNomSliceHists_[22] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex23");
  hscpGenDeDxPNomSliceHists_[23] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex24");
  hscpGenDeDxPNomSliceHists_[24] = (TH1F*) hscpGenInputFile_->Get("trackDeDxNoMAndPSlices/histForPDFNoMPBinsMatchedIndex25");
  // below now only does nom vs. eta hist
  initGenHists();
  //XXX debug
  //std::cout << "Entries in muonEtaVsPtHist: " << muonEtaVsPtHist_->GetEntries() << std::endl;
  //std::cout << "Entries in hscpEtaVsPtHist: " << hscpEtaVsPtHist_->GetEntries() << std::endl;
  //std::cout << "Entries in avgTrackNoMvsEtaProf:" << avgTrackNoMvsEtaProf_->GetEntries() << std::endl;
  //std::cout << "Entries in muonGenDeDxPNomSliceHists_[0]: " << muonGenDeDxPNomSliceHists_[0]->GetEntries() << std::endl;
  //std::cout << "Entries in muonGenDeDxPNomSliceHists_[1]: " << muonGenDeDxPNomSliceHists_[1]->GetEntries() << std::endl;
  //std::cout << "Entries in muonGenDeDxPNomSliceHists_[2]: " << muonGenDeDxPNomSliceHists_[2]->GetEntries() << std::endl;
  //std::cout << "Entries in muonGenDeDxPNomSliceHists_[3]: " << muonGenDeDxPNomSliceHists_[3]->GetEntries() << std::endl;
  //std::cout << "Entries in hscpGenDeDxPNomSliceHists_[0]: " << hscpGenDeDxPNomSliceHists_[0]->GetEntries() << std::endl;
  //std::cout << "Entries in hscpGenDeDxPNomSliceHists_[1]: " << hscpGenDeDxPNomSliceHists_[1]->GetEntries() << std::endl;
  //std::cout << "Entries in hscpGenDeDxPNomSliceHists_[2]: " << hscpGenDeDxPNomSliceHists_[2]->GetEntries() << std::endl;
  //std::cout << "Entries in hscpGenDeDxPNomSliceHists_[3]: " << hscpGenDeDxPNomSliceHists_[3]->GetEntries() << std::endl;
  // \debug
  RooDataHist muonEtaPtGenDataHist("muonEtaPtGenDataHist","muonEtaPtGenDataHist",RooArgList(Trackpt,Tracketa),muonEtaVsPtHist_);
  RooHistPdf muonEtaPtGenHistPdf("muonEtaPtGenHistPdf","muonEtaPtGenHistPdf",RooArgSet(Trackpt,Tracketa),muonEtaPtGenDataHist);
  RooDataHist* muonDeDxGenPNomSliceDataHists[4];
  RooHistPdf* muonDeDxGenPNomSliceHistPdfs[4];
  muonDeDxGenPNomSliceDataHists[0] = new RooDataHist("muonDeDxGenPNomSliceDataHist0","muonDeDxGenPNomSliceDataHist0",TrackdEdxE1,muonGenDeDxPNomSliceHists_[0]);
  muonDeDxGenPNomSliceDataHists[1] = new RooDataHist("muonDeDxGenPNomSliceDataHist1","muonDeDxGenPNomSliceDataHist1",TrackdEdxE1,muonGenDeDxPNomSliceHists_[1]);
  muonDeDxGenPNomSliceDataHists[2] = new RooDataHist("muonDeDxGenPNomSliceDataHist2","muonDeDxGenPNomSliceDataHist2",TrackdEdxE1,muonGenDeDxPNomSliceHists_[2]);
  muonDeDxGenPNomSliceDataHists[3] = new RooDataHist("muonDeDxGenPNomSliceDataHist3","muonDeDxGenPNomSliceDataHist3",TrackdEdxE1,muonGenDeDxPNomSliceHists_[3]);
  muonDeDxGenPNomSliceHistPdfs[0] = new RooHistPdf("muonDeDxGenPNomSliceHistPdf0","muonDeDxGenPNomSliceHistPdf0",TrackdEdxE1,*muonDeDxGenPNomSliceDataHists[0]);
  muonDeDxGenPNomSliceHistPdfs[1] = new RooHistPdf("muonDeDxGenPNomSliceHistPdf1","muonDeDxGenPNomSliceHistPdf1",TrackdEdxE1,*muonDeDxGenPNomSliceDataHists[1]);
  muonDeDxGenPNomSliceHistPdfs[2] = new RooHistPdf("muonDeDxGenPNomSliceHistPdf2","muonDeDxGenPNomSliceHistPdf2",TrackdEdxE1,*muonDeDxGenPNomSliceDataHists[2]);
  muonDeDxGenPNomSliceHistPdfs[3] = new RooHistPdf("muonDeDxGenPNomSliceHistPdf3","muonDeDxGenPNomSliceHistPdf3",TrackdEdxE1,*muonDeDxGenPNomSliceDataHists[3]);

  RooDataHist hscpEtaPtGenDataHist("hscpEtaPtGenDataHist","hscpEtaPtGenDataHist",RooArgList(Trackpt,Tracketa),hscpEtaVsPtHist_);
  RooHistPdf hscpEtaPtGenHistPdf("hscpEtaPtGenHistPdf","hscpEtaPtGenHistPdf",RooArgSet(Trackpt,Tracketa),hscpEtaPtGenDataHist);
  RooDataHist* hscpDeDxGenPNomSliceDataHists[25];
  RooHistPdf* hscpDeDxGenPNomSliceHistPdfs[25];
  hscpDeDxGenPNomSliceDataHists[0] = new RooDataHist("hscpDeDxGenPNomSliceDataHist0","hscpDeDxGenPNomSliceDataHist0",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[0]);
  hscpDeDxGenPNomSliceDataHists[1] = new RooDataHist("hscpDeDxGenPNomSliceDataHist1","hscpDeDxGenPNomSliceDataHist1",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[1]);
  hscpDeDxGenPNomSliceDataHists[2] = new RooDataHist("hscpDeDxGenPNomSliceDataHist2","hscpDeDxGenPNomSliceDataHist2",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[2]);
  //hscpDeDxGenPNomSliceDataHists[3] = new RooDataHist("hscpDeDxGenPNomSliceDataHist3","hscpDeDxGenPNomSliceDataHist3",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[3]);
  //hscpDeDxGenPNomSliceDataHists[4] = new RooDataHist("hscpDeDxGenPNomSliceDataHist4","hscpDeDxGenPNomSliceDataHist4",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[4]);
  hscpDeDxGenPNomSliceDataHists[5] = new RooDataHist("hscpDeDxGenPNomSliceDataHist5","hscpDeDxGenPNomSliceDataHist5",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[5]);
  hscpDeDxGenPNomSliceDataHists[6] = new RooDataHist("hscpDeDxGenPNomSliceDataHist6","hscpDeDxGenPNomSliceDataHist6",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[6]);
  hscpDeDxGenPNomSliceDataHists[7] = new RooDataHist("hscpDeDxGenPNomSliceDataHist7","hscpDeDxGenPNomSliceDataHist7",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[7]);
  hscpDeDxGenPNomSliceDataHists[8] = new RooDataHist("hscpDeDxGenPNomSliceDataHist8","hscpDeDxGenPNomSliceDataHist8",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[8]);
  hscpDeDxGenPNomSliceDataHists[9] = new RooDataHist("hscpDeDxGenPNomSliceDataHist9","hscpDeDxGenPNomSliceDataHist9",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[9]);
  hscpDeDxGenPNomSliceDataHists[10] = new RooDataHist("hscpDeDxGenPNomSliceDataHist10","hscpDeDxGenPNomSliceDataHist10",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[10]);
  hscpDeDxGenPNomSliceDataHists[11] = new RooDataHist("hscpDeDxGenPNomSliceDataHist11","hscpDeDxGenPNomSliceDataHist11",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[11]);
  hscpDeDxGenPNomSliceDataHists[12] = new RooDataHist("hscpDeDxGenPNomSliceDataHist12","hscpDeDxGenPNomSliceDataHist12",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[12]);
  hscpDeDxGenPNomSliceDataHists[13] = new RooDataHist("hscpDeDxGenPNomSliceDataHist13","hscpDeDxGenPNomSliceDataHist13",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[13]);
  hscpDeDxGenPNomSliceDataHists[14] = new RooDataHist("hscpDeDxGenPNomSliceDataHist14","hscpDeDxGenPNomSliceDataHist14",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[14]);
  hscpDeDxGenPNomSliceDataHists[15] = new RooDataHist("hscpDeDxGenPNomSliceDataHist15","hscpDeDxGenPNomSliceDataHist15",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[15]);
  hscpDeDxGenPNomSliceDataHists[16] = new RooDataHist("hscpDeDxGenPNomSliceDataHist16","hscpDeDxGenPNomSliceDataHist16",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[16]);
  hscpDeDxGenPNomSliceDataHists[17] = new RooDataHist("hscpDeDxGenPNomSliceDataHist17","hscpDeDxGenPNomSliceDataHist17",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[17]);
  hscpDeDxGenPNomSliceDataHists[18] = new RooDataHist("hscpDeDxGenPNomSliceDataHist18","hscpDeDxGenPNomSliceDataHist18",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[18]);
  hscpDeDxGenPNomSliceDataHists[19] = new RooDataHist("hscpDeDxGenPNomSliceDataHist19","hscpDeDxGenPNomSliceDataHist19",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[19]);
  hscpDeDxGenPNomSliceDataHists[20] = new RooDataHist("hscpDeDxGenPNomSliceDataHist20","hscpDeDxGenPNomSliceDataHist20",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[20]);
  hscpDeDxGenPNomSliceDataHists[21] = new RooDataHist("hscpDeDxGenPNomSliceDataHist21","hscpDeDxGenPNomSliceDataHist21",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[21]);
  hscpDeDxGenPNomSliceDataHists[22] = new RooDataHist("hscpDeDxGenPNomSliceDataHist22","hscpDeDxGenPNomSliceDataHist22",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[22]);
  hscpDeDxGenPNomSliceDataHists[23] = new RooDataHist("hscpDeDxGenPNomSliceDataHist23","hscpDeDxGenPNomSliceDataHist23",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[23]);
  hscpDeDxGenPNomSliceDataHists[24] = new RooDataHist("hscpDeDxGenPNomSliceDataHist24","hscpDeDxGenPNomSliceDataHist24",TrackdEdxE1,hscpGenDeDxPNomSliceHists_[24]);
  hscpDeDxGenPNomSliceHistPdfs[0] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf0","hscpDeDxGenPNomSliceHistPdf0",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[0]);
  hscpDeDxGenPNomSliceHistPdfs[1] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf1","hscpDeDxGenPNomSliceHistPdf1",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[1]);
  hscpDeDxGenPNomSliceHistPdfs[2] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf2","hscpDeDxGenPNomSliceHistPdf2",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[2]);
  //hscpDeDxGenPNomSliceHistPdfs[3] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf3","hscpDeDxGenPNomSliceHistPdf3",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[3]);
  //hscpDeDxGenPNomSliceHistPdfs[4] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf4","hscpDeDxGenPNomSliceHistPdf4",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[4]);
  hscpDeDxGenPNomSliceHistPdfs[5] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf5","hscpDeDxGenPNomSliceHistPdf5",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[5]);
  hscpDeDxGenPNomSliceHistPdfs[6] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf6","hscpDeDxGenPNomSliceHistPdf6",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[6]);
  hscpDeDxGenPNomSliceHistPdfs[7] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf7","hscpDeDxGenPNomSliceHistPdf7",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[7]);
  hscpDeDxGenPNomSliceHistPdfs[8] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf8","hscpDeDxGenPNomSliceHistPdf8",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[8]);
  hscpDeDxGenPNomSliceHistPdfs[9] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf9","hscpDeDxGenPNomSliceHistPdf9",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[9]);
  hscpDeDxGenPNomSliceHistPdfs[10] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf10","hscpDeDxGenPNomSliceHistPdf10",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[10]);
  hscpDeDxGenPNomSliceHistPdfs[11] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf11","hscpDeDxGenPNomSliceHistPdf11",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[11]);
  hscpDeDxGenPNomSliceHistPdfs[12] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf12","hscpDeDxGenPNomSliceHistPdf12",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[12]);
  hscpDeDxGenPNomSliceHistPdfs[13] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf13","hscpDeDxGenPNomSliceHistPdf13",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[13]);
  hscpDeDxGenPNomSliceHistPdfs[14] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf14","hscpDeDxGenPNomSliceHistPdf14",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[14]);
  hscpDeDxGenPNomSliceHistPdfs[15] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf15","hscpDeDxGenPNomSliceHistPdf15",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[15]);
  hscpDeDxGenPNomSliceHistPdfs[16] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf16","hscpDeDxGenPNomSliceHistPdf16",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[16]);
  hscpDeDxGenPNomSliceHistPdfs[17] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf17","hscpDeDxGenPNomSliceHistPdf17",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[17]);
  hscpDeDxGenPNomSliceHistPdfs[18] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf18","hscpDeDxGenPNomSliceHistPdf18",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[18]);
  hscpDeDxGenPNomSliceHistPdfs[19] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf19","hscpDeDxGenPNomSliceHistPdf19",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[19]);
  hscpDeDxGenPNomSliceHistPdfs[20] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf20","hscpDeDxGenPNomSliceHistPdf20",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[20]);
  hscpDeDxGenPNomSliceHistPdfs[21] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf21","hscpDeDxGenPNomSliceHistPdf21",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[21]);
  hscpDeDxGenPNomSliceHistPdfs[22] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf22","hscpDeDxGenPNomSliceHistPdf22",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[22]);
  hscpDeDxGenPNomSliceHistPdfs[23] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf23","hscpDeDxGenPNomSliceHistPdf23",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[23]);
  hscpDeDxGenPNomSliceHistPdfs[24] = new RooHistPdf("hscpDeDxGenPNomSliceHistPdf24","hscpDeDxGenPNomSliceHistPdf24",TrackdEdxE1,*hscpDeDxGenPNomSliceDataHists[24]);

  // ******** Create CALC PDFs (use muon dE/dx E1 from P sideband; Nom slices) ***************
  dEdxE1PAndNoMSliceFuncs_[0] = new TF1("dEdxE1PAndNoMSliceFunc0",fitTailSlice2,3,25,7);
  dEdxE1PAndNoMSliceFuncs_[1] = new TF1("dEdxE1PAndNoMSliceFunc1",fitTailSlice3,3,25,7);
  dEdxE1PAndNoMSliceFuncs_[2] = new TF1("dEdxE1PAndNoMSliceFunc2",fitTailSlice4,3,25,7);
  dEdxE1PAndNoMSliceFuncs_[3] = new TF1("dEdxE1PAndNoMSliceFunc3",fitTailSlice5,3,25,7);
  initCalcPDFTF1s();

  RooAbsPdf* dEdxE1PAndNoMSlicePdfs[4];
  dEdxE1PAndNoMSlicePdfs[0] = bindPdf(dEdxE1PAndNoMSliceFuncs_[0],TrackdEdxE1);
  dEdxE1PAndNoMSlicePdfs[1] = bindPdf(dEdxE1PAndNoMSliceFuncs_[1],TrackdEdxE1);
  dEdxE1PAndNoMSlicePdfs[2] = bindPdf(dEdxE1PAndNoMSliceFuncs_[2],TrackdEdxE1);
  dEdxE1PAndNoMSlicePdfs[3] = bindPdf(dEdxE1PAndNoMSliceFuncs_[3],TrackdEdxE1);


  // set random seed to clock time
  RooRandom::randomGenerator()->SetSeed(time(0));

  // turn off progress messages for below
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  // Background
  RooMCStudy mcStudyEtaPtMuon(muonEtaPtGenHistPdf,RooArgSet(Trackpt,Tracketa));
  RooMCStudy* mcStudyDeDxSlicesMuon[4];
  //mcStudyDeDxSlicesMuon[0] = new RooMCStudy(*muonDeDxGenPNomSliceHistPdfs[0],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //mcStudyDeDxSlicesMuon[1] = new RooMCStudy(*muonDeDxGenPNomSliceHistPdfs[1],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //mcStudyDeDxSlicesMuon[2] = new RooMCStudy(*muonDeDxGenPNomSliceHistPdfs[2],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //mcStudyDeDxSlicesMuon[3] = new RooMCStudy(*muonDeDxGenPNomSliceHistPdfs[3],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //XXX -- use the analytical function PDF to gen (AND calc)
  mcStudyDeDxSlicesMuon[0] = new RooMCStudy(*dEdxE1PAndNoMSlicePdfs[0],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesMuon[1] = new RooMCStudy(*dEdxE1PAndNoMSlicePdfs[1],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesMuon[2] = new RooMCStudy(*dEdxE1PAndNoMSlicePdfs[2],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesMuon[3] = new RooMCStudy(*dEdxE1PAndNoMSlicePdfs[3],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));

  // Signal
  RooMCStudy mcStudyEtaPtHscp(hscpEtaPtGenHistPdf,RooArgSet(Trackpt,Tracketa),Extended(true));
  RooMCStudy* mcStudyDeDxSlicesHscp[25];
  mcStudyDeDxSlicesHscp[0] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[0],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[1] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[1],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[2] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[2],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //TODO: eventually fill these slices?
  //mcStudyDeDxSlicesHscp[3] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[3],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  //mcStudyDeDxSlicesHscp[4] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[4],RooArgSet(TrackdEdxE1),
  //    FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[5] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[5],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[6] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[6],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[7] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[7],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[8] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[8],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[9] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[9],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[10] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[10],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[11] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[11],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[12] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[12],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[13] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[13],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[14] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[14],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[15] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[15],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[16] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[16],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[17] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[17],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[18] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[18],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[19] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[19],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[20] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[20],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[21] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[21],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[22] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[22],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[23] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[23],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));
  mcStudyDeDxSlicesHscp[24] = new RooMCStudy(*hscpDeDxGenPNomSliceHistPdfs[24],RooArgSet(TrackdEdxE1),
      FitOptions(Save(true),PrintEvalErrors(0),Hesse(true)));


  TFile* outFile;
  if(genSignalAndBackground_)
    outFile = new TFile("generateToyMcRooFit.signalAndBackground.root","recreate");
  else
    outFile = new TFile("generateToyMcRooFit.backgroundOnly.root","recreate");

  outFile->cd();

  TH1F* nllSignalHist = new TH1F("nllSignalHist","NLL signal model",2000,-12000,-10000);
  TH1F* nllBackHist = new TH1F("nllBackHist","NLL background model",20000,0,20000);
  TH1F* nllSBHist = new TH1F("nllSBHist","NLL signal+background model",2000,-12000,-10000);
  TH1F* dedxOfHighestNLLHist = new TH1F("dedxOfHighestNLL","dE/dx of highest 5 NLL points in each trial",500,0,100);
  TH1F* ptGenHistMuon = new TH1F("ptGenMuon","Track pt of generated muons;GeV",1000,0,1000);
  TH1F* etaGenHistMuon = new TH1F("etaGenMuon","Track #eta of generated muons;GeV",600,-3,3);
  TH1F* dEdxGenHistMuon = new TH1F("dEdxGenMuon","dE/dx of generated muons;MeV/cm",250,0,50);
  TH2F* etaVsPtGenHistMuon = new TH2F("etaVsPtGenMuon","Track #eta vs. pt of generated muons;Pt [GeV];#eta",1000,0,1000,600,-3,3);
  TH2F* nomVsEtaGenHistMuon = new TH2F("nomVsEtaGenMuon","Track nom vs. #eta of generated muons;#eta",600,-3,3,35,0,35);
  TH1F* ptGenHistHscp = new TH1F("ptGenHscp","Track pt of generated hscps;GeV",1000,0,1000);
  TH1F* etaGenHistHscp = new TH1F("etaGenHscp","Track #eta of generated hscps;GeV",600,-3,3);
  TH1F* dEdxGenHistHscp = new TH1F("dEdxGenHscp","dE/dx of generated hscps;MeV/cm",250,0,50);
  TH2F* etaVsPtGenHistHscp = new TH2F("etaVsPtGenHscp","Track #eta vs. pt of generated hscps;Pt [GeV];#eta",1000,0,1000,600,-3,3);
  TH2F* nomVsEtaGenHistHscp = new TH2F("nomVsEtaGenHscp","Track nom vs. #eta of generated hscps;#eta",600,-3,3,35,0,35);
  TH1F* nllPerHscpEventHist = new TH1F("nllPerHscpEventHist","NLL per HSCP event",100,0,100);
  TH1F* nllPerMuonEventHist = new TH1F("nllPerMuonEventHist","NLL per muon event",100,0,100);
  TH2F* dEdxVsPtGenHscpHist = new TH2F("dEdxVsPtGenHscp","dE/dx vs. Pt gen HSCP;GeV;MeV/cm",500,0,1000,250,0,50);
  TH2F* dEdxVsPtGenMuonHist = new TH2F("dEdxVsPtGenMuon","dE/dx vs. Pt gen muon;GeV;MeV/cm",500,0,1000,250,0,50);
  TH2F* nllPerEventVsDeDxHist = new TH2F("nllPerEventVsDeDx","NLL/event vs. dE/dx",50,0,25,100,-50,50);
  TH2F* nllPerEventVsDeDxSlice0Hist = new TH2F("nllPerEventVsDeDxSlice0","NLL/event vs. dE/dx slice 0",50,0,25,100,-50,50);
  TH2F* nllPerEventVsDeDxSlice1Hist = new TH2F("nllPerEventVsDeDxSlice1","NLL/event vs. dE/dx slice 1",50,0,25,100,-50,50);
  TH2F* nllPerEventVsDeDxSlice2Hist = new TH2F("nllPerEventVsDeDxSlice2","NLL/event vs. dE/dx slice 2",50,0,25,100,-50,50);
  TH2F* nllPerEventVsDeDxSlice3Hist = new TH2F("nllPerEventVsDeDxSlice3","NLL/event vs. dE/dx slice 3",50,0,25,100,-50,50);
  TH1F* slicesUsedHscpHist = new TH1F("slicesUsedHscpHist","PDF sSlice2lices used HSCP",4,0,4);
  TH1F* slicesUsedMuonHist = new TH1F("slicesUsedMuonHist","PDF slices used muon",4,0,4);
  TH1F* chi2Hist = new TH1F("chi2Hist","Chi2",5000,0,5000);
  TH1F* nllSamplesLowDeDxHist = new TH1F("nllSamplesLowDeDx","NLL for samples with no track with dE/dx > 20 MeV/cm",20000,0,20000);
  TH2F* chi2VsNLLBinnedHist = new TH2F("chi2VsNLLBinnedHist","chi2 vs. nll binned",20000,0,20000,5000,0,5000);

  //// Generation params -- TODO make these configurable (or when we move to EDAnalyzer)
  //int numberOfTrials = 1000;
  //int numberOfEventsPerTrial = 27893+51;
  //float signalFraction = 51./27893; // # HSCP tracks / # muon tracks in ~ 35/pb, Pt > 100 GeV && dE/dx > 3 MeV/cm
  ////XXX for testing ONLY -- compare to gluino500 xc limit
  ////signalFraction = 14./27893;
  ////numberOfEventsPerTrial = 27893+14;
  ////XXX end testing
  //int numSignalEvents = 0;
  //int numBackgroundEvents = numberOfEventsPerTrial;
  //if(genSignalAndBackground_)
  //{
  //  numSignalEvents = round(numberOfEventsPerTrial*signalFraction);
  //  numBackgroundEvents = numberOfEventsPerTrial-numSignalEvents;
  //}

  edm::LogInfo("GenerateHSCPToyMC") << "Generating " << numBackgroundEvents_
    << " background events and " << numSignalEvents_
    << " signal events with Poisson fluctuations.";

  if(!genSignalAndBackground_)
  {
    mcStudyEtaPtMuon.generate(numberOfTrials_,numBackgroundEvents_,true);
  }
  else
  {
    mcStudyEtaPtMuon.generate(numberOfTrials_,numBackgroundEvents_,true);
    mcStudyEtaPtHscp.generate(numberOfTrials_,numSignalEvents_,true);
  }

  


  //Loop over toys, generate the dE/dx, and then calculate the likelihood
  for(int i=0; i<numberOfTrials_; ++i)
  {
    bool hasHighDeDxTrack = false;
    double NLLthisSample = 0;
    map<float,float> nllsAndDeDxs;

    const RooDataSet* muonData = ( RooDataSet*) (mcStudyEtaPtMuon.genData(i));
    RooDataSet muonDeDxData("muonDeDxData","muonDeDxData",RooArgSet(TrackdEdxE1,Tracknom,Trackp));
    RooDataSet hscpDeDxData("hscpDeDxData","hscpDeDxData",RooArgSet(TrackdEdxE1,Tracknom,Trackp));
    double chi2 = 0;
    RooArgSet* tempMuon = new RooArgSet();
    RooArgSet* tempHscp = new RooArgSet();
    RooArgSet* temp = new RooArgSet();
    for(int entry = 0; entry < muonData->sumEntries(); ++entry)
    {
      tempMuon = const_cast<RooArgSet*>(muonData->get(entry));
      Tracketa = tempMuon->getRealValue("Tracketa");
      Trackpt = tempMuon->getRealValue("Trackpt");
      int nomMuon = getNomFromEta(Tracketa.getVal());
      Tracknom = nomMuon;
      double theta = 2*TMath::ATan(TMath::Exp(-1*Tracketa.getVal()));
      if(nomMuon < 6) //XXX cutting on tracks with < 6 measurements
        continue;
      int slice = getSliceFromNom(nomMuon);
      int index = getSliceFromNomAndPHscp(nomMuon,Trackpt.getVal()/TMath::Sin(theta));
      if(index < 1)
        continue;
      mcStudyDeDxSlicesMuon[slice-2]->generate(1,1,true);
      RooDataSet* deDxDataSetMuon(const_cast<RooDataSet*>((RooDataSet*)mcStudyDeDxSlicesMuon[slice-2]->genData(0)));
      const RooArgSet* temp2muon = deDxDataSetMuon->get(0);
      TrackdEdxE1 = temp2muon->getRealValue("TrackdEdxE1");
      Trackp = Trackpt.getVal()/TMath::Sin(theta);
      muonDeDxData.add(RooArgSet(TrackdEdxE1,Tracknom,Trackp));
      // Fill hists
      etaVsPtGenHistMuon->Fill(Trackpt.getVal(),Tracketa.getVal());
      nomVsEtaGenHistMuon->Fill(Tracketa.getVal(),Tracknom.getVal());
      ptGenHistMuon->Fill(Trackpt.getVal());
      dEdxGenHistMuon->Fill(TrackdEdxE1.getVal());
      etaGenHistMuon->Fill(Tracketa.getVal());
      dEdxVsPtGenMuonHist->Fill(Trackpt.getVal()/TMath::Sin(theta),TrackdEdxE1.getVal());
      if(TrackdEdxE1.getVal() > 20)
        hasHighDeDxTrack = true;
    }
    if(genSignalAndBackground_)
    {
      const RooDataSet* hscpData = (RooDataSet*) mcStudyEtaPtHscp.genData(i);
      for(int entry = 0; entry < hscpData->sumEntries(); ++entry)
      {
        tempHscp = const_cast<RooArgSet*>(hscpData->get(entry));
        Tracketa = tempHscp->getRealValue("Tracketa");
        Trackpt = tempHscp->getRealValue("Trackpt");
        int nomHscp = getNomFromEta(Tracketa.getVal());
        Tracknom = nomHscp;
        double theta = 2*TMath::ATan(TMath::Exp(-1*Tracketa.getVal()));
        if(nomHscp < 6) //XXX cutting on tracks with < 6 measurements
          continue;
        int index = getSliceFromNomAndPHscp(nomHscp,Trackpt.getVal()/TMath::Sin(theta));
        if(index < 1)
          continue;
        mcStudyDeDxSlicesHscp[index-1]->generate(1,1,true);
        RooDataSet* deDxDataSetHscp(const_cast<RooDataSet*>((RooDataSet*) mcStudyDeDxSlicesHscp[index-1]->genData(0)));
        const RooArgSet* temp2hscp = deDxDataSetHscp->get(0);
        TrackdEdxE1 = temp2hscp->getRealValue("TrackdEdxE1");
        Trackp = Trackpt.getVal()/TMath::Sin(theta);
        hscpDeDxData.add(RooArgSet(TrackdEdxE1,Tracknom,Trackp));
        // Fill hists
        etaVsPtGenHistHscp->Fill(Trackpt.getVal(),Tracketa.getVal());
        nomVsEtaGenHistHscp->Fill(Tracketa.getVal(),Tracknom.getVal());
        ptGenHistHscp->Fill(Trackpt.getVal());
        dEdxGenHistHscp->Fill(TrackdEdxE1.getVal());
        etaGenHistHscp->Fill(Tracketa.getVal());
        dEdxVsPtGenHscpHist->Fill(Trackpt.getVal(),TrackdEdxE1.getVal());
        if(TrackdEdxE1.getVal() > 20)
          hasHighDeDxTrack = true;
      }
    }

    // Make combined datasets and datahists
    RooDataSet combinedDataSet = muonDeDxData;
    combinedDataSet.append(hscpDeDxData);
    RooDataHist* deDxDataHistCombined = combinedDataSet.binnedClone();
    RooDataHist* deDxDataHistCombinedNomOnly = (RooDataHist*)deDxDataHistCombined->reduce(Tracknom);
    RooDataHist* deDxDataHistCombinedDeDxNom = (RooDataHist*)deDxDataHistCombined->reduce(RooArgSet(TrackdEdxE1,Tracknom));


    //double nllFromRooFit = 0;
    //double llrFromRooFit = 0;
    double llrByHand = 0;
    RooDataHist* deDxDataHistCombinedTrackNomCut = (RooDataHist*) deDxDataHistCombinedDeDxNom->reduce("Tracknom>5");
    for(int nomBin=1; nomBin < Tracknom.numBins(); ++nomBin)
    {
      RooDataHist* deDxDataHistCombinedThisBin = new RooDataHist();
      if(nomBin==1)
      {
        deDxDataHistCombinedThisBin = (RooDataHist*) deDxDataHistCombinedTrackNomCut->reduce(RooArgSet(TrackdEdxE1),"Tracknom<11");
        //TCanvas* dedxDataHistCanvas = new TCanvas("dedxDataHistCanvas","dedxDataHistCanvas",1,1,2500,1000);
        //dedxDataHistCanvas->cd();
        //RooPlot* dedxFrame = TrackdEdxE1.frame();
        //deDxDataHistCombinedThisBin->plotOn(dedxFrame);
        //dEdxE1PAndNoMSlicePdfs[0]->plotOn(dedxFrame);
        //dedxFrame->Draw("e0");
        //dedxDataHistCanvas->Write();
        //dedxFrame->Write();
      }
      else if(nomBin==2)          
        deDxDataHistCombinedThisBin = (RooDataHist*) deDxDataHistCombinedTrackNomCut->reduce(RooArgSet(TrackdEdxE1),"Tracknom<16&&Tracknom>=11");
      else if(nomBin==3)          
        deDxDataHistCombinedThisBin = (RooDataHist*) deDxDataHistCombinedTrackNomCut->reduce(RooArgSet(TrackdEdxE1),"Tracknom<21&&Tracknom>=16");
      else if(nomBin==4)
        deDxDataHistCombinedThisBin = (RooDataHist*) deDxDataHistCombinedTrackNomCut->reduce(RooArgSet(TrackdEdxE1),"Tracknom<50&&Tracknom>=21");

      if(deDxDataHistCombinedThisBin->sumEntries() < 1)
        continue;

      //RooRealVar norm("norm","norm",deDxDataHistCombinedThisBin->sumEntries());
      ////std::cout << "nomBin=" << nomBin << " #entries in dataset=" << deDxDataHistCombinedThisBin->sumEntries() << std::endl;
      ////double nllValModel = dEdxE1PAndNoMSlicePdfs[nomBin-1]->createNLL(*deDxDataHistCombinedThisBin,Extended(true))->getVal();
      //RooExtendPdf extendPdf("extendPdf","extendPdf",*dEdxE1PAndNoMSlicePdfs[nomBin-1],norm);
      //double nllValModel = extendPdf.createNLL(*deDxDataHistCombinedThisBin,Extended(true))->getVal();
      ////std::cout << "NLLValModel: " << nllValModel << std::endl;
      ////std::cout << "ExtendedTermModel: (expt-obsv*log(expt)= " << extendPdf.extendedTerm((int)norm.getVal(),&RooArgSet(TrackdEdxE1)) << std::endl;
      //RooHistPdf thisBinDataHistPdf("thisBinDataHistPdf","thisBinDataHistPdf",TrackdEdxE1,*deDxDataHistCombinedThisBin);
      ////double nllValMax = thisBinDataHistPdf.createNLL(*deDxDataHistCombinedThisBin,Extended(true))->getVal();
      //RooExtendPdf extend2Pdf("extend2Pdf","extend2Pdf",thisBinDataHistPdf,norm);
      //double nllValMax = extend2Pdf.createNLL(*deDxDataHistCombinedThisBin,Extended(true))->getVal();
      ////std::cout << "NLLValMax: " << nllValMax << std::endl;
      ////std::cout << "ExtendedTermModel: (expt-obsv*log(expt)= " << extend2Pdf.extendedTerm((int)norm.getVal(),&RooArgSet(TrackdEdxE1)) << std::endl;
      //double llr = 2*(nllValModel-nllValMax);
      ////std::cout << "This LLR roofit ( =2*(NLLValModel-NLLValMax) ): " << llr << std::endl;
      //llrFromRooFit+=llr;
      //nllFromRooFit+=nllValModel;

      //calc by hand
      double llrByHandThisBin = 0;
      double chi2ThisBin  = 0;
      for(int dedxBin=0; dedxBin < TrackdEdxE1.numBins(); ++dedxBin)
      {
        temp = const_cast<RooArgSet*>(deDxDataHistCombinedThisBin->get(dedxBin));
        if(!deDxDataHistCombinedThisBin->valid())
          continue;
        double norm = deDxDataHistCombinedThisBin->sumEntries()*(TrackdEdxE1.getMax()-TrackdEdxE1.getMin())/TrackdEdxE1.getBins();
        TrackdEdxE1 = temp->getRealValue("TrackdEdxE1");
        double pdfVal = dEdxE1PAndNoMSlicePdfs[nomBin-1]->getVal(new RooArgSet(TrackdEdxE1));
        pdfVal*=norm;
        double binc = deDxDataHistCombinedThisBin->weight();
        //std::cout << "BY HAND2 DEBUG: trackdedx: " << TrackdEdxE1.getVal() << " nom: " << nomBin-1 << " binc= " << binc << " pdfVal=" << pdfVal << " norm=" << norm << std::endl << std::endl;
        double thisLLR = 0;
        if(binc==0)
        {
          thisLLR = pdfVal;
        }
        else
        {
          thisLLR = pdfVal-binc+binc*log(binc/pdfVal);
        }
        thisLLR*=2;
        llrByHand+=thisLLR;
        llrByHandThisBin+=thisLLR;
        //chi2
        if(pdfVal > 0)
        {
          chi2ThisBin+=pow(pdfVal-binc,2)/pdfVal;
          chi2+=pow(pdfVal-binc,2)/pdfVal;
        }
      }
      //std::cout << "This bin 2*LLR model by hand: " << llrByHandThisBin << std::endl << std::endl;
      //std::cout << "This bin chi2 model: " << chi2ThisBin << std::endl;

      delete deDxDataHistCombinedThisBin;
    }
    //std::cout << "-2*log(lamba) roofit = " << llrFromRooFit << std::endl;
    edm::LogInfo("GenerateHSCPToyMC") << "2*Nll by hand this sample: " << llrByHand
      << "chi2 this sample: " << chi2 << std::endl;

    delete deDxDataHistCombined;
    delete deDxDataHistCombinedNomOnly;
    delete deDxDataHistCombinedDeDxNom;

    //// plot an example -- slice 0
    //TCanvas* eventsCanvas = new TCanvas("eventsCanvas","eventsCanvas",1,1,2500,1000);
    //eventsCanvas->cd();
    //RooPlot* dedxFrame = TrackdEdxE1.frame();
    //deDxDataHistMuon->plotOn(dedxFrame);
    //dEdxE1PAndNoMSlicePdfs[0]->plotOn(dedxFrame);
    //dedxFrame->Draw("e0");
    //eventsCanvas->Write();
    //dedxFrame->Write();

    chi2VsNLLBinnedHist->Fill(llrByHand,chi2);
    nllBackHist->Fill(llrByHand);
    //nllSignalHist->Fill(nllBinnedSignalOnly);
    //nllSBHist->Fill(nllBinnedPlusSignal);
    chi2Hist->Fill(chi2);
    //std::cout << "-2*LogLikelihoodBack = " << nllBinned << std::endl;
    //std::cout << "-2*LogLikelihoodSig = " << nllBinnedSignalOnly << std::endl;
    //std::cout << "chi2 = " << chi2 << std::endl;
    //nllSignalHist->Fill(logLikelihoodSignal);
    //std::cout << "-LogLikelihood = " << logLikelihood << std::endl;
    //std::cout << "-LogLikelihoodBack = " << logLikelihoodBack << std::endl;
    //std::cout << "-LogLikelihoodSig = " << logLikelihoodSignal << std::endl;
    if(!hasHighDeDxTrack)
      nllSamplesLowDeDxHist->Fill(NLLthisSample);

    //if((i+1) % 100 == 0)
    edm::LogInfo("GenerateHSCPToyMC") << "Done with trial: " << i+1 << std::endl;
  }

  //// One example to plot
  //mcStudy.generate(1,numberOfEventsPerTrial,true);
  //RooDataSet* data(const_cast<RooDataSet*>(mcStudy.genData(0)));
  //TCanvas* eventsCanvas = new TCanvas("eventsCanvas","eventsCanvas",1,1,2500,1000);
  //eventsCanvas->cd();
  //RooPlot* dedxFrame = TrackdEdxE1.frame();
  //data->plotOn(dedxFrame);
  ////langausDeDxMuonPt100to200.plotOn(dedxFrame);
  //sigAndBackPt100to200ModelGen.plotOn(dedxFrame);
  //sigAndBackPt100to200ModelGen.plotOn(dedxFrame,Components("langausDeDxHSCP"),LineStyle(kDashed),LineColor(2));
  //dedxFrame->Draw("e0");
  //eventsCanvas->Write();
  //dedxFrame->Write();
  //RooHist* hist = dedxFrame->getHist();
  //hist->Write();

  //gen
  muonEtaVsPtHist_->Write();
  hscpEtaVsPtHist_->Write();

  // close input files
  muonGenInputFile_->Close();
  hscpGenInputFile_->Close();

  nllPerEventVsDeDxSlice0Hist->Write();
  nllPerEventVsDeDxSlice1Hist->Write();
  nllPerEventVsDeDxSlice2Hist->Write();
  nllPerEventVsDeDxSlice3Hist->Write();
  nllPerEventVsDeDxHist->Write();
  nllPerMuonEventHist->Write();
  nllPerHscpEventHist->Write();
  nllSignalHist->Write();
  nllBackHist->Write();
  nllSBHist->Write();
  dedxOfHighestNLLHist->Write();
  nomVsEtaGenHistMuon->Write();
  etaVsPtGenHistMuon->Write();
  ptGenHistMuon->Write();
  dEdxGenHistMuon->Write();
  etaGenHistMuon->Write();
  dEdxVsPtGenMuonHist->Write();
  slicesUsedMuonHist->Write();
  chi2Hist->Write();

  //hscp
  nomVsEtaGenHistHscp->Write();
  etaVsPtGenHistHscp->Write();
  ptGenHistHscp->Write();
  dEdxGenHistHscp->Write();
  etaGenHistHscp->Write();
  dEdxVsPtGenHscpHist->Write();
  slicesUsedHscpHist->Write();

  nllSamplesLowDeDxHist->Write();

  TDirectory* muonCalcFuncDir = outFile->mkdir("muonCalcPDFsPNoMslices");
  muonCalcFuncDir->cd();
  for(int i=0; i <4; ++i)
    dEdxE1PAndNoMSliceFuncs_[i]->Write();

  outFile->Close();
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenerateHSCPToyMC::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenerateHSCPToyMC::endJob()
{
}

// ------------- Convert int to string for printing -------------------------------------
std::string GenerateHSCPToyMC::intToString(int num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}


void GenerateHSCPToyMC::initGenHists()
{
  //Profile
  avgTrackNoMvsEtaProf_ = new TProfile("avgTrackNoMvsEtaProf","Track number of dE/dx measurements vs. #eta",600,-3,3,"");
  initTrackEtaVsNomProfile(avgTrackNoMvsEtaProf_);
}

// **************************************************************************************
// Lookup for eta--><NoM>
int GenerateHSCPToyMC::getNomFromEta(float eta)
{
  if(fabs(eta) > 2.5)
    return -1; // definitely outside tracker fiducial acceptance

  // Use average profile histogram
  int binx = avgTrackNoMvsEtaProf_->GetXaxis()->FindBin(eta);
  return avgTrackNoMvsEtaProf_->GetBinContent(binx);
}

// **************************************************************************************
// Lookup for NoM->Slice
int GenerateHSCPToyMC::getSliceFromNom(int nom)
{
  if(nom==-1)
    return -1; // generated outside TK fiducial acceptance

  if(nom >=1 && nom <= 5)
    return 1;
  else if(nom >= 6 && nom <= 10)
    return 2;
  else if(nom >= 11 && nom <= 15)
    return 3;
  else if(nom >= 16 && nom <= 20)
    return 4;
  else if(nom > 20)
    return 5;
  else
    return -1;  // something strange happened (nom==0, etc.)
}

// **************************************************************************************
// Lookup for NoM,P->Slice Hscp
int GenerateHSCPToyMC::getSliceFromNomAndPHscp(int nom, float p)
{
  if(nom==-1)
    return -1; // generated outside TK fiducial acceptance

  int column = -1;
  if(nom >=1 && nom <= 5)
    column = 1;
  else if(nom >= 6 && nom <= 10)
    column = 2;
  else if(nom >= 11 && nom <= 15)
    column = 3;
  else if(nom >= 16 && nom <= 20)
    column = 4;
  else if(nom > 20)
    column = 5;
  else
    return -1;  // something strange happened (nom==0, etc.)

  int row = -1;
  if(p > 15 && p < 50)
    row = 1;
  else if(p < 150)
    row = 2;
  else if(p < 300)
    row = 3;
  else if(p < 450)
    row = 4;
  else if(p >= 450)
    row = 5;
  else
    return -1; // something strange (p < 15 somehow)

  return (row-1)*5+column;
}

// **************************************************************************************
// Define Calc PDF TF1s
void GenerateHSCPToyMC::initCalcPDFTF1s()
{
  // index 0 --> bin 2 in root file aka slice as above (NoM 6-10)
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(0,0.1096);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(1,-3.017);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(2,17.26);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(3,4.545e4);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(4,3.14);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(5,0.3092);
  dEdxE1PAndNoMSliceFuncs_[0]->FixParameter(6,19.5);

  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(0,0.1391);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(1,-3.52);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(2,18.76);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(3,6.371e5);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(4,3.069);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(5,0.2745);
  dEdxE1PAndNoMSliceFuncs_[1]->FixParameter(6,25.31);

  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(3,0);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(0,2.076);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(1,-23.54);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(2,67.66);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(4,3);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(5,0.25);
  dEdxE1PAndNoMSliceFuncs_[2]->FixParameter(6,24.04);

  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(3,0);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(0,2.833);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(1,-30.74);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(2,80.94);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(4,3);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(5,0.25);
  dEdxE1PAndNoMSliceFuncs_[3]->FixParameter(6,23.67);
}

// **************************************************************************************
// *******************  Function to set up my nice TStyle  ******************************
void GenerateHSCPToyMC::SetEStyle()
{
  TStyle* EStyle = new TStyle("EStyle", "E's not Style");

  //set the background color to white
  EStyle->SetFillColor(10);
  EStyle->SetFrameFillColor(10);
  EStyle->SetFrameFillStyle(0);
  EStyle->SetFillStyle(0);
  EStyle->SetCanvasColor(10);
  EStyle->SetPadColor(10);
  EStyle->SetTitleFillColor(0);
  EStyle->SetStatColor(10);

  //dont put a colored frame around the plots
  EStyle->SetFrameBorderMode(0);
  EStyle->SetCanvasBorderMode(0);
  EStyle->SetPadBorderMode(0);

  //use the primary color palette
  EStyle->SetPalette(1,0);
  EStyle->SetNumberContours(255);

  //set the default line color for a histogram to be black
  EStyle->SetHistLineColor(kBlack);

  //set the default line color for a fit function to be red
  EStyle->SetFuncColor(kRed);

  //make the axis labels black
  EStyle->SetLabelColor(kBlack,"xyz");

  //set the default title color to be black
  EStyle->SetTitleColor(kBlack);

  // Sizes

  //For Small Plot needs
  //set the margins
  //  EStyle->SetPadBottomMargin(.2);
  //   EStyle->SetPadTopMargin(0.08);
  //   EStyle->SetPadLeftMargin(0.12);
  //   EStyle->SetPadRightMargin(0.12);

  //   //set axis label and title text sizes
  //   EStyle->SetLabelSize(0.06,"xyz");
  //   EStyle->SetTitleSize(0.06,"xyz");
  //   EStyle->SetTitleOffset(1.,"x");
  //   EStyle->SetTitleOffset(.9,"yz");
  //   EStyle->SetStatFontSize(0.04);
  //   EStyle->SetTextSize(0.06);
  //   EStyle->SetTitleBorderSize(0.5);


  //set the margins
  EStyle->SetPadBottomMargin(.15);
  EStyle->SetPadTopMargin(0.08);
  EStyle->SetPadLeftMargin(0.14);
  EStyle->SetPadRightMargin(0.12);

  //set axis label and title text sizes
  EStyle->SetLabelSize(0.04,"xyz");
  EStyle->SetTitleSize(0.06,"xyz");
  EStyle->SetTitleOffset(1.,"x");
  EStyle->SetTitleOffset(1.1,"yz");
  EStyle->SetStatFontSize(0.04);
  EStyle->SetTextSize(0.04);
  EStyle->SetTitleBorderSize(Width_t(0.5));
  //EStyle->SetTitleY(0.5);

  //set line widths
  EStyle->SetHistLineWidth(1);
  EStyle->SetFrameLineWidth(2);
  EStyle->SetFuncWidth(2);

  //Paper Size
  EStyle->SetPaperSize(TStyle::kUSLetter);

  // Misc

  //align the titles to be centered
  //Style->SetTextAlign(22);

  //set the number of divisions to show
  EStyle->SetNdivisions(506, "xy");

  //turn off xy grids
  EStyle->SetPadGridX(0);
  EStyle->SetPadGridY(0);

  //set the tick mark style
  EStyle->SetPadTickX(1);
  EStyle->SetPadTickY(1);

  //show the fit parameters in a box
  EStyle->SetOptFit(111111);

  //turn on all other stats
  //EStyle->SetOptStat(0000000);
  //EStyle->SetOptStat(1111111);
  EStyle->SetOptStat(2222);

  //Move stats box
  EStyle->SetStatX(0.75);

  //marker settings
  EStyle->SetMarkerStyle(8);
  EStyle->SetMarkerSize(0.8);

  // Fonts
  EStyle->SetStatFont(42);
  EStyle->SetLabelFont(42,"xyz");
  EStyle->SetTitleFont(42,"xyz");
  EStyle->SetTextFont(42);
  //  EStyle->SetStatFont(82);
  //   EStyle->SetLabelFont(82,"xyz");
  //   EStyle->SetTitleFont(82,"xyz");
  //   EStyle->SetTextFont(82);


  //done
  EStyle->cd();
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenerateHSCPToyMC);
