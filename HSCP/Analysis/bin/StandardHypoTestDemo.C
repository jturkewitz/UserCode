// standard demo macro for just jypothesis test (no inversion !!!) using hybrid calculator 



#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;

#if ROOT_VERSION_CODE < ROOT_VERSION(5,32,02)
#define USE_OLD_ROOT
#endif

bool noSystematics = false;              // force all systematics to be off (i.e. set all nuisance parameters as constat
double nToysRatio = 5;                   // ratio Ntoys Null/ntoys ALT
//double poiValue = 1;                    // change poi snapshot value for S+B model
//double poiValue = 1e-2;

void StandardHypoTestDemo(const char* infile = "",
    const char* workspaceName = "combined",
    const char* modelSBName = "ModelConfig",
    const char* modelBName = "",
    const char* dataName = "obsData", 
    int calcType = 0, // 0 freq 1 hybrid, 2 asymptotic
    int testStatType = 3,   // 0 LEP, 1 TeV, 2 LHC, 3 LHC - one sided
    int ntoys = 25000,
    bool useNC = false, 
    const char * nuisPriorName = 0,
    const char* outputFile = "",
    double poiValue = 1e-2)
{

  //if (calcType == 1) 
  //  ToyMCSampler::SetAlwaysUseMultiGen(false);
  //else
  //  ToyMCSampler::SetAlwaysUseMultiGen(true);
  //SIC APR 20 -- never use multigen
  ToyMCSampler::SetAlwaysUseMultiGen(false);

  //SimpleLikelihoodRatioTestStat::SetAlwaysReuseNLL(true);
  RooRandom::randomGenerator()->SetSeed(0);

  /////////////////////////////////////////////////////////////
  // First part is just to access a user-defined file 
  // or create the standard example file if it doesn't exist
  ////////////////////////////////////////////////////////////
  const char* filename = "";
  if (!strcmp(infile,""))
    filename = "results/example_combined_GaussExample_model.root";
  else
    filename = infile;
  // Check if example input file exists
  TFile *file = TFile::Open(filename);

  // if input file was specified byt not found, quit
  if(!file && strcmp(infile,"")){
    cout <<"file not found" << endl;
    return;
  } 

  // if default file not found, try to create it
  if(!file ){
    // Normally this would be run on the command line
    cout <<"will run standard hist2workspace example"<<endl;
    gROOT->ProcessLine(".! prepareHistFactory .");
    gROOT->ProcessLine(".! hist2workspace config/example.xml");
    cout <<"\n\n---------------------"<<endl;
    cout <<"Done creating example input"<<endl;
    cout <<"---------------------\n\n"<<endl;
  }

  // now try to access the file again
  file = TFile::Open(filename);
  if(!file){
    // if it is still not there, then we can't continue
    cout << "Not able to run hist2workspace to create example input" <<endl;
    return;
  }


  /////////////////////////////////////////////////////////////
  // Tutorial starts here
  ////////////////////////////////////////////////////////////

  // get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get(workspaceName);
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }
  w->Print();

  // get the modelConfig out of the file
  ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);


  // get the modelConfig out of the file
  RooAbsData* data = w->data(dataName);

  // make sure ingredients are found
  if(!data || !sbModel){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }
  // make b model
  ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);


  // case of no systematics
  // remove nuisance parameters from model
  if (noSystematics) { 
    const RooArgSet * nuisPar = sbModel->GetNuisanceParameters();
    if (nuisPar && nuisPar->getSize() > 0) { 
      std::cout << "StandardHypoTestInvDemo" << "  -  Switch off all systematics by setting them constant to their initial values" << std::endl;
      RooStats::SetAllConstant(*nuisPar);
    }
    if (bModel) { 
      const RooArgSet * bnuisPar = bModel->GetNuisanceParameters();
      if (bnuisPar) 
        RooStats::SetAllConstant(*bnuisPar);
    }
  }


  if (!bModel ) {
    Info("StandardHypoTestInvDemo","The background model %s does not exist",modelBName);
    Info("StandardHypoTestInvDemo","Copy it from ModelConfig %s and set POI to zero",modelSBName);
    bModel = (ModelConfig*) sbModel->Clone();
    bModel->SetName(TString(modelSBName)+TString("B_only"));      
    RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
    if (!var) return;
    double oldval = var->getVal();
    var->setVal(0);
    //bModel->SetSnapshot( RooArgSet(*var, *w->var("lumi"))  );
    bModel->SetSnapshot( RooArgSet(*var)  );
    var->setVal(oldval);
  }

  if (!sbModel->GetSnapshot() || poiValue > 0) { 
    Info("StandardHypoTestDemo","Model %s has no snapshot  - make one using model poi",modelSBName);
    RooRealVar * var = dynamic_cast<RooRealVar*>(sbModel->GetParametersOfInterest()->first());
    if (!var) return;
    double oldval = var->getVal();
    if (poiValue > 0)  var->setVal(poiValue);
    //sbModel->SetSnapshot( RooArgSet(*var, *w->var("lumi") ) );
    sbModel->SetSnapshot( RooArgSet(*var) );
    if (poiValue > 0) var->setVal(oldval);
    //sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
  }


  ////XXX SIC mod
  //const RooArgSet* nuisParams = (bModel->GetNuisanceParameters() ) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
  //// set all gamma stat constraint terms constant
  //RooArgSet* statConstraints = nuisParams->selectByName("gamma_stat_*");
  //statConstraints->setAttribAll("Constant");


  // part 1, hypothesis testing 
  SimpleLikelihoodRatioTestStat * slrts = new SimpleLikelihoodRatioTestStat(*bModel->GetPdf(), *sbModel->GetPdf());
  // null parameters must includes snapshot of poi plus the nuisance values 
  RooArgSet nullParams(*bModel->GetSnapshot());
  if (bModel->GetNuisanceParameters()) nullParams.add(*bModel->GetNuisanceParameters());
  slrts->SetNullParameters(nullParams);
  RooArgSet altParams(*sbModel->GetSnapshot());
  if (sbModel->GetNuisanceParameters()) altParams.add(*sbModel->GetNuisanceParameters());
  slrts->SetAltParameters(altParams);


  ProfileLikelihoodTestStat * profll = new ProfileLikelihoodTestStat(*bModel->GetPdf());


  RatioOfProfiledLikelihoodsTestStat * 
    ropl = new RatioOfProfiledLikelihoodsTestStat(*bModel->GetPdf(), *sbModel->GetPdf(), sbModel->GetSnapshot());
  ropl->SetSubtractMLE(false);

  if (testStatType == 3) profll->SetOneSidedDiscovery(1);
  // profll.SetMinimizer(minimizerType.c_str());
  //profll.SetPrintLevel(2);

  // profll.SetReuseNLL(mOptimize);
  // slrts.SetReuseNLL(mOptimize);
  // ropl.SetReuseNLL(mOptimize);

  // SIC TEST OUTPUT
  //profll.SetPrintLevel(3);


  HypoTestCalculatorGeneric *  hypoCalc = 0;
  // note here Null is B and Alt is S+B
  if (calcType == 0) hypoCalc = new  FrequentistCalculator(*data, *sbModel, *bModel);
  else if (calcType == 1) hypoCalc= new  HybridCalculator(*data, *sbModel, *bModel);
  else if (calcType == 2) hypoCalc= new  AsymptoticCalculator(*data, *sbModel, *bModel);

  //if (calcType != 2)
  //  hypoCalc->SetToys(ntoys, ntoys/nToysRatio);
  // SIC
  if(calcType==0)
    ((FrequentistCalculator*) hypoCalc)->SetToys(ntoys,ntoys/nToysRatio);
  else if(calcType==1)
    ((HybridCalculator*) hypoCalc)->SetToys(ntoys,ntoys/nToysRatio);

  // check for nuisance prior pdf in case of nuisance parameters 
  if (calcType == 1 && (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters() )) {
    RooAbsPdf * nuisPdf = 0; 
    if (nuisPriorName) nuisPdf = w->pdf(nuisPriorName);
    // use prior defined first in bModel (then in SbModel)
    if (!nuisPdf)  { 
      Info("StandardHypoTestDemo","No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the   model");
#ifndef USE_OLD_ROOT
      if (bModel->GetPdf() && bModel->GetObservables() ) 
        nuisPdf = RooStats::MakeNuisancePdf(*bModel,"nuisancePdf_bmodel");
      else 
        nuisPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisancePdf_sbmodel");
#endif
    }   
    if (!nuisPdf ) {
      if (bModel->GetPriorPdf())  { 
        nuisPdf = bModel->GetPriorPdf();
        Info("StandardHypoTestDemo","No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf->GetName());            
      }
      else { 
        Error("StandardHypoTestDemo","Cannnot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived");
        return;
      }
    }
    assert(nuisPdf);
    Info("StandardHypoTestDemo","Using as nuisance Pdf ... " );
    nuisPdf->Print();

    const RooArgSet * nuisParams = (bModel->GetNuisanceParameters() ) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
    RooArgSet * np = nuisPdf->getObservables(*nuisParams);
    if (np->getSize() == 0) { 
      Warning("StandardHypoTestDemo","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
    }
    delete np;

    ((HybridCalculator*)hypoCalc)->ForcePriorNuisanceAlt(*nuisPdf);
    ((HybridCalculator*)hypoCalc)->ForcePriorNuisanceNull(*nuisPdf);

    //XXX SIC mod
    // set all gamma stat constraint terms constant
    RooArgSet* statConstraints = nuisParams->selectByName("gamma_stat_*");
    statConstraints->setAttribAll("Constant");
  }

  // hypoCalc->ForcePriorNuisanceAlt(*sbModel->GetPriorPdf());
  // hypoCalc->ForcePriorNuisanceNull(*bModel->GetPriorPdf());

  ToyMCSampler * sampler = (ToyMCSampler *)hypoCalc->GetTestStatSampler();
  if (useNC) sampler->SetNEventsPerToy(1); 
  // SIC APR 20 -- never use multigen
  sampler->SetUseMultiGen(false);


  if (testStatType == 0) sampler->SetTestStatistic(slrts); 
  if (testStatType == 1) sampler->SetTestStatistic(ropl); 
  if (testStatType == 2 || testStatType == 3) sampler->SetTestStatistic(profll); 


  HypoTestResult *  htr = hypoCalc->GetHypoTest();
  htr->SetPValueIsRightTail(true);
  htr->SetBackgroundAsAlt(false);
  htr->Print(); // how to get meaningfull CLs at this point?
  //cout << "Test statistic on data: " << htr->GetTestStatisticData() << endl;

  TFile* outFile = new TFile(outputFile,"recreate");
  outFile->cd();
  htr->Write();
  outFile->Close();

  delete sampler;
  delete slrts; 
  delete ropl; 
  delete profll;

  if (calcType != 2) {
    HypoTestPlot * plot = new HypoTestPlot(*htr,100);
    plot->SetLogYaxis(true);
    plot->Draw();
  }
  else { 
    std::cout << "Asymptotic results " << std::endl;

  }

  // look at expected significances 
  // found median of S+B distribution
  if (calcType != 2) { 

    SamplingDistribution * altDist = htr->GetAltDistribution();   
    HypoTestResult htExp("Expected Result");
    htExp.Append(htr);
    // find quantiles in alt (S+B) distribution 
    double p[5];
    double q[5];
    for (int i = 0; i < 5; ++i) { 
      double sig = -2  + i;
      p[i] = ROOT::Math::normal_cdf(sig,1);
    }
    std::vector<double> values = altDist->GetSamplingDistribution();
    TMath::Quantiles( values.size(), 5, &values[0], q, p, false);  

    for (int i = 0; i < 5; ++i) { 
      htExp.SetTestStatisticData( q[i] );
      double sig = -2  + i;      
      std::cout << " Expected p -value and significance at " << sig << " sigma = " 
        << htExp.NullPValue() << " significance " << htExp.Significance() << " sigma " << std::endl; 

    }
  }
  else { 
    // case of asymptotic calculator 
    for (int i = 0; i < 5; ++i) { 
      double sig = -2  + i;      
      // sigma is inverted here 
      double pval = AsymptoticCalculator::GetExpectedPValues( htr->NullPValue(), htr->AlternatePValue(), -sig, false);
      std::cout << " Expected p -value and significance at " << sig << " sigma = " 
        << pval << " significance " << ROOT::Math::normal_quantile_c(pval,1) << " sigma " << std::endl; 

    }
  }

}

