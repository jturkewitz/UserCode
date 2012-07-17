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
//         Created:  Mon May 23 14:02:32 CEST 2010
// $Id: GenerateHSCPToyMC.h,v 1.2 2011/06/30 13:51:12 scooper Exp $
//
//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFile.h"

#include "HSCP2011/ToyMC/interface/GenerateToyMcRooFit.calcPdfTF1Defs.h"
#include "HSCP2011/ToyMC/interface/GenerateToyMcRooFit.trackNomVsEta.mu2010abTrigSkimPlot.h"

//
// class declaration
//

class GenerateHSCPToyMC : public edm::EDAnalyzer {
   public:
      explicit GenerateHSCPToyMC(const edm::ParameterSet&);
      ~GenerateHSCPToyMC();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      std::string intToString(int num);

      // ----------member data ---------------------------
      void initGenHists();
      int getNomFromEta(float eta);
      int getSliceFromNom(int nom);
      int getSliceFromNomAndPHscp(int nom, float p);
      void initCalcPDFTF1s();
      void SetEStyle();

      TFile* muonGenInputFile_;
      TFile* hscpGenInputFile_;
      TH2F* muonEtaVsPtHist_;
      TH2F* hscpEtaVsPtHist_;
      TH2F* hscpDeDxVsPHist_;
      TProfile* avgTrackNoMvsEtaProf_;
      TH1F* muonGenDeDxPNomSliceHists_[4];
      TH1F* hscpGenDeDxPNomSliceHists_[25];
      TF1* dEdxE1PAndNoMSliceFuncs_[4];
      bool genSignalAndBackground_;
      double signalFraction_;
      int numberOfTrials_;
      int numberOfEventsPerTrial_;
      int numSignalEvents_;
      int numBackgroundEvents_;

};
