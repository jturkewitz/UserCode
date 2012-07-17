#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooWorkspace.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooProduct.h"
#include "RooRealSumPdf.h"
#include "RooProdPdf.h"
#include "RooMsgService.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodTestStat.h"

using namespace RooFit;
using namespace RooStats;

analyticalWithoutHistFactory()
{
  RooWorkspace* wspace = new RooWorkspace("combined");

  // lumi
  RooRealVar Lumi("Lumi","Lumi",1.0,0,2);
  RooRealVar nominalLumi("nominalLumi","nominalLumi",1.0);
  nominalLumi.setConstant();
  RooGaussian lumiConstraint("lumiConstraint","lumiConstraint",Lumi,nominalLumi,RooRealConstant::value(0.022));

  const int numIasBins = 50;
  // observable
  RooRealVar ias("ias","ias",0,1);
  ias.setBins(numIasBins);
  // poi
  RooRealVar sigXSec("sigXSec","sigXSec",0,1);
  RooArgSet poi("poi");
  poi.add(sigXSec);

  // background pdf -- gluino1000 exp
  RooRealVar expoSlope("expoSlope","expoSlope",-2.38);
  expoSlope.setConstant();
  RooExponential backgroundExp("backgroundExp","backgroundExp",ias,expoSlope);
  RooRealVar backNorm("backNorm","backNorm",2.151); // normalize to exp events
  backNorm.setConstant();
  RooExtendPdf backgroundOnly("backgroundOnly","backgroundOnly",backgroundExp,backNorm);

  // signal pdf -- taken from Gluino1000, from apr16 directory on laptop (nom11-12,eta0-0.2)
  RooRealVar gausMean("gausMean","gausMean",0.7481);
  RooRealVar gausWidth("gausWidth","gausWidth",0.04995);
  gausMean.setConstant();
  gausWidth.setConstant();
  RooGaussian signalGaus("signalGaus","signalGaus",ias,gausMean,gausWidth);
  RooRealVar sigNorm("sigNorm","sigNorm",868.440); // normalize to exp events with 1 pb cross section
  sigNorm.setConstant();
  RooProduct lumiXCrossSection("lumiXCrossSection","lumiXCrossSection",RooArgSet(Lumi,sigXSec));
  RooProduct sigNormXLumiXCrossSection("sigNormXLumiXCrossSection","sigNormXLumiXCrossSection",RooArgSet(sigNorm,lumiXCrossSection));
  //RooExtendPdf signal("signal","signal",signalGaus,lumiXCrossSection);

  RooAddPdf combinedModel("combinedModel","combinedModel",RooArgList(backgroundExp,signalGaus),
      RooArgList(backNorm,sigNormXLumiXCrossSection));

  RooProdPdf model("model","model",RooArgSet(lumiConstraint,combinedModel));
  // bg only
  RooProdPdf bgOnlyModel("bgOnlyModel","bgOnlyModel",RooArgSet(lumiConstraint,backgroundOnly));

  // generate background-only data
  RooDataSet* asimovDataSet = bgOnlyModel.generate(ias,RooFit::ExpectedData());
  wspace->import(*asimovDataSet,Rename("AsimovData"));

  // create set of global observables (need to be defined as constants)
  RooArgSet globalObs("global_obs");
  globalObs.add(nominalLumi);

  // create set of nuisance parameters
  RooArgSet nuis("nuis");
  nuis.add(Lumi);

  // create set of observables
  RooArgSet obs("observables");
  obs.add(ias);

  ModelConfig* modelConfig = new ModelConfig("modelConfig");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf(model);
  modelConfig->SetParametersOfInterest(poi);
  modelConfig->SetGlobalObservables(globalObs);
  modelConfig->SetNuisanceParameters(nuis);
  modelConfig->SetObservables(obs);

  // set parameter snapshot that corresponds to the best fit to data
  RooAbsReal* pNll = modelConfig->GetPdf()->createNLL( *asimovDataSet );
  RooAbsReal* pProfile = pNll->createProfile( globalObs ); // do not profile global observables
  pProfile->getVal(); // this will do fit and set POI and nuisance parameters to fitted values
  RooArgSet* pPoiAndNuisance = new RooArgSet("poiAndNuisance");
  pPoiAndNuisance->add(*modelConfig->GetParametersOfInterest());
  modelConfig->SetSnapshot(*pPoiAndNuisance);
  delete pProfile;
  delete pNll;
  delete pPoiAndNuisance;

  wspace->import(*modelConfig);
  wspace->Print();
  modelConfig->Print();
  
  TFile* myFile = new TFile("myFile.root","recreate");
  myFile->cd();
  modelConfig->Write();
  wspace->Write();
  myFile->Close();

  ///////////////////////
  // generate my own toys

  sigXSec.setVal(0);
  std::cout << "Cross section value: " << sigXSec.getVal() << std::endl;
  // check to make sure I have expectedEvents()
  std::cout << "Expected events (ias): " << model.expectedEvents(ias) << std::endl;

  // test stat - SB
  ProfileLikelihoodTestStat profll(model);
  profll.SetOneSided(true);
  //profll.SetStrategy(2);
  //profll.SetTolerance(1e-10);
  //profll.SetPrintLevel(3);
  // hists of test stat value
  TH1F* bgDataTestStatHist = new TH1F("bgDataTestStat","BG-only data test stat",500,0,5);
  TH1F* sbDataTestStatHist = new TH1F("sbDataTestStat","SB data test stat",500,0,5);
  double sigma1 = 0.00384;
  double sigma2 = 0.00390571;
  double sigma3 = 1e-3;

  //TH1* modelHistogram = model.createHistogram("ias",100);
  //modelHistogram->Scale(model.expectedEvents(ias));


  //// test stat on data
  //sigXSec.setVal(sigma3);
  //RooMsgService::instance().getStream(1).addTopic(Generation);
  //RooAbsData* testData = model.generate(RooArgSet(ias,Lumi),Extended());
  //sigXSec.setVal(sigma1);
  //profll.SetPrintLevel(3);
  //double testStatOnData = profll.Evaluate(*testData,RooArgSet(sigXSec));
  double testStatOnData = 3.15963; // we know this for sigma=0.00384
  //double testStatOnData = 3.21374; // we know this for sigma=0.00390571
  TH1F* dataTestStatHist = new TH1F("dataTestStat","SB data test stat",500,0,5);
  dataTestStatHist->Fill(testStatOnData);
  std::cout << "testStatOnData: " << testStatOnData << std::endl;

  profll.SetPrintLevel(0);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  int numToys = 2000;
  for(int i=0; i < numToys; ++i)
  {
    //
    sigXSec.setVal(sigma1);
    RooAbsData* toySBData = model.generate(RooArgSet(ias,Lumi),Extended());
    //TH1* toyHist = toyData->createHistogram("ias",50);
    sigXSec.setVal(sigma1);
    double testStatValue = profll.Evaluate(*toySBData,RooArgSet(sigXSec));
    sbDataTestStatHist->Fill(testStatValue);
    // BG data
    sigXSec.setVal(0);
    RooAbsData* toyBData = model.generate(RooArgSet(ias,Lumi),Extended());
    sigXSec.setVal(sigma1);
    testStatValue = profll.Evaluate(*toyBData,RooArgSet(sigXSec));
    bgDataTestStatHist->Fill(testStatValue);

    delete toySBData;
    delete toyBData;
  }
  //modelHistogram->Draw();
  //toyHist->SetMarkerColor(kRed);
  //toyHist->SetLineColor(kRed);
  //toyHist->Draw("sames");

  bgDataTestStatHist->SetLineColor(kBlue);
  bgDataTestStatHist->SetMarkerColor(kBlue);
  bgDataTestStatHist->Draw();
  sbDataTestStatHist->SetLineColor(kRed);
  sbDataTestStatHist->SetMarkerColor(kRed);
  sbDataTestStatHist->Draw("sames");
  dataTestStatHist->Draw("sames");

  double CLb = bgDataTestStatHist->Integral(bgDataTestStatHist->FindBin(testStatOnData),bgDataTestStatHist->GetNbinsX());
  double CLsb = sbDataTestStatHist->Integral(sbDataTestStatHist->FindBin(testStatOnData),sbDataTestStatHist->GetNbinsX());
  double CLs = CLsb/CLb;
  std::cout << "CLs: " << CLs << std::endl;
}

