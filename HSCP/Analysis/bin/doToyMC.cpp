#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h" 
#include "TRandom3.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TArray.h"

//RooFit
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooMsgService.h"

#include "TSystem.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <vector>

namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle; class HSCPIsolation;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"

#include "PhysicsTools/FWLite/interface/CommandLineParser.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#endif

#include "commonFunctions.h"

// helper functions

// ****** main
int main(int argc, char ** argv)
{
  using namespace RooFit;

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // parse arguments
  if(argc < 2)
  {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  // get the python configuration
  const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  // now get each parameterset
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("doToyMC");
  // mass cut to use for the high-p high-Ih (search region) ias dist
  double massCut_ (ana.getParameter<double>("MassCut"));
  std::string backgroundPredictionRootFilename_ (ana.getParameter<std::string>("BackgroundPredictionInputRootFile"));
  std::string signalRootFilename_ (ana.getParameter<std::string>("SignalInputRootFile"));
  std::string outputRootFilename_ (ana.getParameter<std::string>("OutputRootFile"));
  // dE/dx calibration
  double dEdx_k (ana.getParameter<double>("dEdx_k"));
  double dEdx_c (ana.getParameter<double>("dEdx_c"));
  // definition of sidebands/search region
  double pSidebandThreshold (ana.getParameter<double>("PSidebandThreshold"));
  double ihSidebandThreshold (ana.getParameter<double>("IhSidebandThreshold"));
  // num background tracks D region after all skimming and preselections
  int numBackgroundTracksInDRegionInput (ana.getParameter<int>("NumberOfBackgroundTracksInDRegion"));
  // number of signal tracks expected in this int. lumi
  int numSignalTracksTotal (ana.getParameter<int>("NumberOfSignalTracks"));
  int numTrialsToDo (ana.getParameter<int>("NumberOfTrials"));
  int lowerNoM (ana.getParameter<int>("LowerNoMOfSlice"));
  double lowerEta (ana.getParameter<double>("LowerEtaOfSlice"));
  bool verbose_ (ana.getParameter<bool>("Verbose"));
  int reportEvery_ (ana.getParameter<int>("ReportEvery"));

  if(lowerNoM < 5 || lowerNoM > 21 || lowerEta < 0 || lowerEta > 2.2)
  {
    std::cout << "Eta/NoM slice specified (eta=" << lowerEta <<
      "-" << lowerEta+0.2 << ", nom=" << lowerNoM << "-" <<
      lowerNoM+1 << ") is out of range NoM [5-21], eta [0-2.2]. Exiting."
      << std::endl;
    return -1;
  }

  if(!verbose_) // kill all RooFit messages except warning and error
  {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }

  TFile* outputRootFile = new TFile(outputRootFilename_.c_str(),"recreate");
  TFile* backgroundPredictionRootFile = TFile::Open(backgroundPredictionRootFilename_.c_str());
  TFile* signalRootFile = TFile::Open(signalRootFilename_.c_str());

  // get this eta/nom hist
  string dirName="iasPredictionsVariableBins";
  //string dirName="iasPredictionsFixedBins";
  string getHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  string sigHistName = getHistName;
  sigHistName+="iasSignalVarBinHist";
  getHistName+="iasPredictionVarBinHist";
  //getHistName+="iasPredictionFixedHist";
  std::string iasSignalVarBinHistTitle = "Ias of signal";
  iasSignalVarBinHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);

  string fullPath = dirName;
  fullPath+="/";
  fullPath+=getHistName;

  TH1F* iasBackgroundPredictionMassCutNoMSliceHist =
    (TH1F*)backgroundPredictionRootFile->Get(fullPath.c_str());

  TH1F* iasSignalMassCutNoMSliceHist;
  if(iasBackgroundPredictionMassCutNoMSliceHist)
  {
    iasSignalMassCutNoMSliceHist =
      (TH1F*)iasBackgroundPredictionMassCutNoMSliceHist->Clone();
    iasSignalMassCutNoMSliceHist->Reset("ICESM");
    iasSignalMassCutNoMSliceHist->SetName(sigHistName.c_str());
    iasSignalMassCutNoMSliceHist->SetTitle(iasSignalVarBinHistTitle.c_str());
  }
  else
  {
    std::cout << "Cannot proceed further: no background prediction found (histogram " <<
      fullPath << ") in file: " << backgroundPredictionRootFilename_ << ")." << std::endl;
    return -2;
  }

  // get the background B region hist
  dirName="bRegionHistograms";
  getHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  getHistName+="bRegionHist";
  fullPath=dirName;
  fullPath+="/";
  fullPath+=getHistName;
  TH1F* bRegionBackgroundHist = (TH1F*)backgroundPredictionRootFile->Get(fullPath.c_str());
  if(!bRegionBackgroundHist)
  {
    std::cout << "Cannot proceed further: no b region hist found (histogram " <<
      fullPath << ") in file: " << backgroundPredictionRootFilename_ << ")." << std::endl;
    return -3;
  }

  // get the background entries in D region hist
  getHistName = "entriesInDRegion";
  TH2F* entriesInDRegionHist = (TH2F*)backgroundPredictionRootFile->Get(getHistName.c_str());
  if(!entriesInDRegionHist)
  {
    std::cout << "Cannot proceed further: no D region entries hist found (histogram " <<
      getHistName << ") in file: " << backgroundPredictionRootFilename_ << ")." << std::endl;
    return -4;
  }

  // define observables
  //TODO move to something like common functions
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",0,2.5);

  // get roodataset from signal file
  RooDataSet* rooDataSetAllSignal = (RooDataSet*)signalRootFile->Get("rooDataSetAll");
  if(!rooDataSetAllSignal)
  {
    std::cout << "Problem with RooDataSet named rooDataSetAll in signal file " <<
      signalRootFilename_.c_str() << std::endl;
    return -3;
  }
  
  // construct D region for signal
  std::string pSearchCutString = "rooVarP>";
  pSearchCutString+=floatToString(pSidebandThreshold);
  std::string ihSearchCutString = "rooVarIh>";
  ihSearchCutString+=floatToString(ihSidebandThreshold);
  RooDataSet* regionD1DataSetSignal = (RooDataSet*)rooDataSetAllSignal->reduce(Cut(pSearchCutString.c_str()));
  RooDataSet* regionDDataSetSignal = (RooDataSet*)regionD1DataSetSignal->reduce(Cut(ihSearchCutString.c_str()));
  int numSignalTracksInDRegion = regionDDataSetSignal->numEntries();
  // construct D region for signal in this eta/nom slice
  std::string nomCutString = "rooVarNoMias==";
  nomCutString+=intToString(lowerNoM);
  nomCutString+="||rooVarNoMias==";
  nomCutString+=intToString(lowerNoM+1);
  RooDataSet* nomCutDRegionDataSetSignal =
    (RooDataSet*)regionDDataSetSignal->reduce(Cut(nomCutString.c_str()));
  std::string etaCutString = "rooVarEta>";
  etaCutString+=floatToString(lowerEta);
  etaCutString+="&&rooVarEta<";
  etaCutString+=floatToString(lowerEta+0.2);
  RooDataSet* etaCutNomCutDRegionDataSetSignal =
    (RooDataSet*)nomCutDRegionDataSetSignal->reduce(Cut(etaCutString.c_str()));
  int numSignalTracksInDRegionThisSlice = etaCutNomCutDRegionDataSetSignal->numEntries();
  // construct signal prediction hist (with mass cut)
  const RooArgSet* argSet = etaCutNomCutDRegionDataSetSignal->get();
  RooRealVar* iasData = (RooRealVar*)argSet->find(rooVarIas.GetName());
  RooRealVar* ihData = (RooRealVar*)argSet->find(rooVarIh.GetName());
  RooRealVar* pData = (RooRealVar*)argSet->find(rooVarP.GetName());
  for(int evt=0; evt < etaCutNomCutDRegionDataSetSignal->numEntries(); ++evt)
  {
    etaCutNomCutDRegionDataSetSignal->get(evt);
    // apply mass cut
    float massSqr = (ihData->getVal()-dEdx_c)*pow(pData->getVal(),2)/dEdx_k;
    if(massSqr < 0)
      continue;
    else if(sqrt(massSqr) >= massCut_)
      iasSignalMassCutNoMSliceHist->Fill(iasData->getVal());
  }
  double numSignalTracksInDRegionMassCut = iasSignalMassCutNoMSliceHist->Integral();
  for(int bin=1; bin<=iasSignalMassCutNoMSliceHist->GetNbinsX(); ++bin)
  {
    double binWidthRatio =
      iasSignalMassCutNoMSliceHist->GetBinWidth(1)/
      iasSignalMassCutNoMSliceHist->GetBinWidth(bin);
    double binc = iasSignalMassCutNoMSliceHist->GetBinContent(bin);
    double bine = iasSignalMassCutNoMSliceHist->GetBinError(bin);
    iasSignalMassCutNoMSliceHist->SetBinContent(bin,binWidthRatio*binc);
    iasSignalMassCutNoMSliceHist->SetBinError(bin,binWidthRatio*bine);
  }

  RooDataHist* iasBackgroundDataHist;
  RooHistPdf* iasBackgroundHistPdf;
  RooDataHist* iasSignalDataHist;
  RooHistPdf* iasSignalHistPdf;
  // Make the dataHist/histPdf -- background
  std::string histName = iasBackgroundPredictionMassCutNoMSliceHist->GetName();
  std::string histTitle = iasBackgroundPredictionMassCutNoMSliceHist->GetTitle();
  std::string dataHistName=histName;
  dataHistName+="DataHistBG";
  std::string dataHistTitle=histTitle;
  dataHistTitle+=" dataHistBG";
  std::string histPdfName=histName;
  histPdfName+="HistPdfBG";
  std::string histPdfTitle=histTitle;
  histPdfTitle+=" histPdfBG";
  iasBackgroundDataHist = new RooDataHist(dataHistName.c_str(),dataHistTitle.c_str(),
      rooVarIas,iasBackgroundPredictionMassCutNoMSliceHist);
  iasBackgroundHistPdf = new RooHistPdf(histPdfName.c_str(),histPdfTitle.c_str(),
      rooVarIas, *iasBackgroundDataHist);
  // signal
  histName = iasSignalMassCutNoMSliceHist->GetName();
  histTitle = iasSignalMassCutNoMSliceHist->GetTitle();
  dataHistName=histName;
  dataHistName+="DataHistSig";
  dataHistTitle=histTitle;
  dataHistTitle+=" dataHistSig";
  histPdfName=histName;
  histPdfName+="HistPdfSig";
  histPdfTitle=histTitle;
  histPdfTitle+=" histPdfSig";
  iasSignalDataHist = new RooDataHist(dataHistName.c_str(),dataHistTitle.c_str(),
      rooVarIas,iasSignalMassCutNoMSliceHist);
  iasSignalHistPdf = new RooHistPdf(histPdfName.c_str(),histPdfTitle.c_str(),
      rooVarIas, *iasSignalDataHist);


  int nBins = iasBackgroundPredictionMassCutNoMSliceHist->GetNbinsX();
  RooBinning thisBinning(nBins,
      iasBackgroundPredictionMassCutNoMSliceHist->GetXaxis()->GetXbins()->GetArray());
  rooVarIas.setBinning(thisBinning);
  const int numTrials = numTrialsToDo;
  // figure out how many background tracks to generate in this eta/nom slice
  double fractionOfBGTracksPassingMassCutThisSlice =
    iasBackgroundPredictionMassCutNoMSliceHist->Integral()/bRegionBackgroundHist->Integral();
  double bgEntriesInDRegionThisSlice =
    entriesInDRegionHist->GetBinContent(entriesInDRegionHist->FindBin(lowerEta,lowerNoM));
  double fractionOfBGTracksInThisSlice = bgEntriesInDRegionThisSlice/entriesInDRegionHist->Integral();
  float numBackgroundTracksThisSlice =
    fractionOfBGTracksInThisSlice*fractionOfBGTracksPassingMassCutThisSlice*numBackgroundTracksInDRegionInput;
  // figure out how many signal tracks to generate in this eta/nom slice
  double fractionOfSigTracksPassingMassCutThisSlice =
    numSignalTracksInDRegionMassCut/numSignalTracksInDRegionThisSlice;
  double fractionOfSigTracksInThisSlice = numSignalTracksInDRegionThisSlice/(double)numSignalTracksInDRegion;
  float numSignalTracksThisSlice =
    fractionOfSigTracksInThisSlice*fractionOfSigTracksPassingMassCutThisSlice*numSignalTracksTotal;
  // sig + back PDF
  float totalNumTracksThisSlice = numSignalTracksThisSlice+numBackgroundTracksThisSlice;
  float sigFracThisSlice = numSignalTracksThisSlice/(double)totalNumTracksThisSlice;
  RooRealVar fsig("fsig","signal frac",sigFracThisSlice,0,1);
  RooHistPdf* iasBackgroundHistPdfForSBModel = new RooHistPdf("histPdfForSBModelOnly","Hist PDF for BG Model only",
      rooVarIas, *iasBackgroundDataHist);
  //RooAddPdf sigBackPdf("sigBackPdf","sigBackPdf",RooArgList(*iasSignalHistPdf,
  //      *iasBackgroundHistPdf),fsig);
  RooAddPdf sigBackPdf("sigBackPdf","sigBackPdf",RooArgList(*iasSignalHistPdf,
        *iasBackgroundHistPdfForSBModel),fsig);
  double totalSignalFraction =
    numSignalTracksTotal/(double)(numBackgroundTracksInDRegionInput+numSignalTracksTotal);
  bool generateBackgroundOnly = (totalSignalFraction > 0) ? false : true;

  // output
  //std::cout << "Found " << bRegionBackgroundHist->Integral() <<
  //  " tracks in bRegion background and " << iasBackgroundPredictionMassCutNoMSliceHist->Integral() <<
  //  " tracks in dRegion prediction after mass cut." << std::endl;
  std::cout << "Found " << fractionOfBGTracksPassingMassCutThisSlice*100 <<
    "% of background tracks passing mass cut in this slice. " << std::endl;
  std::cout << "Found " << fractionOfBGTracksInThisSlice*100 << "% of background tracks in this slice" <<
    std::endl;
  std::cout << "Generating " << numBackgroundTracksThisSlice << " background tracks" <<
    " from " << numBackgroundTracksInDRegionInput << " for entire D region." << std::endl;
  //std::cout << "Found " << bgEntriesInDRegionThisSlice << 
  //  " background tracks in D region this slice, calculated " << bgTracksInThisDRegionAfterMassCut <<
  //  " background tracks to generate in this slice (after mass cut)." << std::endl;
  std::cout << "Found " << fractionOfSigTracksPassingMassCutThisSlice*100 <<
    "% of signal tracks passing mass cut in this slice. " << std::endl;
  std::cout << "Found " << 100*fractionOfSigTracksInThisSlice << "% of signal tracks in this slice" <<
    std::endl;
  //std::cout << "Found " << numSignalTracksInDRegionThisSlice <<
  //  " signal tracks in D region for this slice and " <<
  //  numSignalTracksInDRegion << " signal tracks in D region. " <<
  //  std::endl;
  if(!generateBackgroundOnly)
    std::cout << "Generating " << numSignalTracksThisSlice << " signal tracks" <<
      " from " << numSignalTracksTotal << " for entire D region." << std::endl;
  else
    std::cout << "Generating 0 signal tracks" <<
      " from " << numSignalTracksTotal << " for entire D region." << std::endl;
  //std::cout << "Found " << numSignalTracksInDRegionThisSlice << 
  //  " signal tracks in D region this slice, calculated " << numSignalTracksThisSliceToGen <<
  //  " signal tracks to generate in this slice (after mass cut)." << std::endl;
  std::cout << "signal fraction: numSignalTracksThisSlice=" << numSignalTracksThisSlice << " / " <<
    " totalNumTracksThisSlice=" << totalNumTracksThisSlice <<
    " --> sigFracThisSlice =" << sigFracThisSlice << std::endl;

  // NLL B/sat hist
  std::string nllBSatHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllBSatHistName+="nllBSatHist";
  std::string nllBSatHistTitle = "-2*ln(B/sat)";
  nllBSatHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllBSatHist = new TH1F(nllBSatHistName.c_str(),nllBSatHistTitle.c_str(),1000,0,100);
  // NLL sat only hist
  std::string nllSatOnlyHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllSatOnlyHistName+="nllSatonlyHist";
  std::string nllSatOnlyHistTitle = "NLL sat only";
  nllSatOnlyHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllSatOnlyHist = new TH1F(nllSatOnlyHistName.c_str(),nllSatOnlyHistTitle.c_str(),30000,-30000,0);
  // NLL B only hist
  std::string nllBOnlyHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllBOnlyHistName+="nllBonlyHist";
  std::string nllBOnlyHistTitle = "NLL B only";
  nllBOnlyHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllBOnlyHist = new TH1F(nllBOnlyHistName.c_str(),nllBOnlyHistTitle.c_str(),30000,-30000,0);
  // NLL S only hist
  std::string nllSonlyHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllSonlyHistName+="nllSonlyHist";
  std::string nllSonlyHistTitle = "NLL S only";
  nllSonlyHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllSOnlyHist = new TH1F(nllSonlyHistName.c_str(),nllSonlyHistTitle.c_str(),30000,-30000,0);
  // NLL SB hist
  std::string nllSBHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllSBHistName+="nllSBHist";
  std::string nllSBHistTitle = "NLL SB";
  nllSBHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllSBHist = new TH1F(nllSBHistName.c_str(),nllSBHistTitle.c_str(),30000,-30000,0);
  // NLL S+B/B
  std::string nllSBOverBHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllSBOverBHistName+="nllSBOverBHist";
  std::string nllSBOverBHistTitle = "NLL -2*ln(SB/B)";
  nllSBOverBHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* nllSBOverBHist = new TH1F(nllSBOverBHistName.c_str(),nllSBOverBHistTitle.c_str(),2000,-100,100);
  // Num tracks per trial
  std::string numTracksPerTrialHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  numTracksPerTrialHistName+="numTracksPerTrial";
  std::string numTracksPerTrialHistTitle = "Tracks per trial";
  numTracksPerTrialHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH1F* numTracksPerTrialHist = new TH1F(numTracksPerTrialHistName.c_str(),numTracksPerTrialHistTitle.c_str(),
      500,totalNumTracksThisSlice-250,totalNumTracksThisSlice+250);
  // Num signal tracks vs. NLL S+B/B
  std::string numSigTracksVsNLLSBOverBHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  numSigTracksVsNLLSBOverBHistName+="numSigTracksVsNLLSBOverB";
  std::string numSigTracksVsNLLSBOverBHistTitle = "Number of Signal tracks vs. NLL -2*ln(SB/B)";
  numSigTracksVsNLLSBOverBHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH2F* numSigTracksVsNLLSBOverBHist = new TH2F(numSigTracksVsNLLSBOverBHistName.c_str(),numSigTracksVsNLLSBOverBHistTitle.c_str(),
      2000,-10,10,15,0,15);
  // Num tracks in last bin vs. NLL S+B/B
  std::string numTracksInLastBinVsSBOverBHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  numTracksInLastBinVsSBOverBHistName+="numTracksInLastBinVsNLLBSOverB";
  std::string numTracksInLastBinVsSBOverBHistTitle = "Number of tracks in last bin vs. NLL -2*ln(SB/B)";
  numTracksInLastBinVsSBOverBHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
  TH2F* numTracksInLastBinVsSBOverBHist = new TH2F(numTracksInLastBinVsSBOverBHistName.c_str(),
      numTracksInLastBinVsSBOverBHistTitle.c_str(),2000,-10,10,50,0,50);
  // background histPDF histogram
  TH1F* bgHistPdfHist = (TH1F*)iasBackgroundHistPdf->createHistogram("iasBackgroundPdfHist",rooVarIas);
  // signal histPDF histogram
  TH1F* sigHistPdfHist = (TH1F*)iasSignalHistPdf->createHistogram("iasSignalPdfHist",rooVarIas);
  // RooDataSet with nlls
  RooRealVar sigNLLRooVar("sigNLLRooVar","signal NLL",0,-100000,100000);
  RooRealVar backNLLRooVar("backNLLRooVar","background NLL",-100000,100000);
  RooRealVar sigBackNLLRooVar("sigBackNLLRooVar","signal+background NLL",-100000,100000);
  RooRealVar satModelNLLRooVar("satModelNLLRooVar","sat model NLL",-100000,100000);
  RooRealVar sigBackFitNLLRooVar("sigBackFitNLLRooVar","signal+background ML fit",-100000,100000);
  RooRealVar maxLSigFracRooVar("maxLSigFracRooVar","max L fit signal frac",0,1);
  RooRealVar assumedTotalSigTracksRooVar("assumedTotalSigTracksRooVar","assumed total signal tracks",0,1e4);
  RooRealVar assumedSigTracksThisSliceRooVar("assumedSigTracksThisSliceRooVar","assumed signal tracks this slice",0,1e4);
  RooRealVar trialIndexRooVar("trialIndex","trial index",0,1e8);
  std::string nllDatasetName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  nllDatasetName+="nllValues";
  RooDataSet* nllValues = new RooDataSet(nllDatasetName.c_str(),nllDatasetName.c_str(),
      RooArgSet(sigNLLRooVar,backNLLRooVar,sigBackNLLRooVar,satModelNLLRooVar,
        sigBackFitNLLRooVar,maxLSigFracRooVar,assumedTotalSigTracksRooVar,
        assumedSigTracksThisSliceRooVar,trialIndexRooVar));
  // actual events in last bin
  RooRealVar actualEventsInLastBinRooVar("actualEventsInLastBinRooVar","actual events in last bin",0,20000);
  RooRealVar signalEventsInLastBinRooVar("signalEventsInLastBinRooVar","signal events in last bin",0,0,20000);
  std::string datasetName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  datasetName+="observedEvents";
  RooDataSet* observedEvents = new RooDataSet(datasetName.c_str(),datasetName.c_str(),
      RooArgSet(actualEventsInLastBinRooVar,signalEventsInLastBinRooVar));
  // expected events this slice, total, expTotSigFrac
  RooRealVar expectedBackgroundEventsInLastBinRooVar("expectedBackgroundEventsInLastBinRooVar","expected bg events in last bin",0,20000);
  RooRealVar expectedSignalEventsInLastBinRooVar("expectedSignalEventsInLastBinRooVar","expected signal events in last bin",0,20000);
  RooRealVar expectedBackgroundEventsThisSliceRooVar("expectedBackgroundEventsThisSliceRooVar","expected bg events this slice",0,50000);
  RooRealVar expectedSignalEventsThisSliceRooVar("expectedSignalEventsThisSliceRooVar","expected signal events this slice",0,20000);
  RooRealVar expectedTotalSignalFractionRooVar("expectedTotalSignalFractionRooVar","expected total signal fraction",0,1);
  //RooRealVar numberOfTrialsRooVar("numberOfTrialsRooVar","number of trials",0,1e8);
  datasetName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
  datasetName+="expectedEvents";
  RooDataSet* expectedEvents = new RooDataSet(datasetName.c_str(),datasetName.c_str(),
      RooArgSet(expectedBackgroundEventsInLastBinRooVar,expectedSignalEventsInLastBinRooVar,
        expectedBackgroundEventsThisSliceRooVar,expectedSignalEventsThisSliceRooVar,
        expectedTotalSignalFractionRooVar));//,numberOfTrialsRooVar));

  // fill expected events RooVars
  expectedBackgroundEventsInLastBinRooVar = 
    numBackgroundTracksThisSlice*bgHistPdfHist->GetBinContent(bgHistPdfHist->GetNbinsX())
    *bgHistPdfHist->GetBinWidth(bgHistPdfHist->GetNbinsX())*bgHistPdfHist->GetNbinsX();
  //std::cout << "Expected bg events last bin: " << expectedEventsInLastBinRooVar.getVal() << std::endl;
  //std::cout << "binWidth last bin: " << bgHistPdfHist->GetBinWidth(bgHistPdfHist->GetNbinsX()) << std::endl;
  //std::cout << "bincontent last bin: " << bgHistPdfHist->GetBinContent(bgHistPdfHist->GetNbinsX()) << std::endl;
  //std::cout << "nbins: " << bgHistPdfHist->GetNbinsX() << std::endl;
  expectedSignalEventsInLastBinRooVar =
    numSignalTracksThisSlice*sigHistPdfHist->GetBinContent(sigHistPdfHist->GetNbinsX())
    *sigHistPdfHist->GetBinWidth(sigHistPdfHist->GetNbinsX())*sigHistPdfHist->GetNbinsX();
  expectedBackgroundEventsThisSliceRooVar = numBackgroundTracksThisSlice;
  expectedSignalEventsThisSliceRooVar = numSignalTracksThisSlice;
  expectedTotalSignalFractionRooVar = totalSignalFraction;
  //numberOfTrialsRooVar = numTrials;
  // add it into the dataset
  expectedEvents->add(RooArgSet(expectedBackgroundEventsInLastBinRooVar,expectedSignalEventsInLastBinRooVar,
                        expectedBackgroundEventsThisSliceRooVar,expectedSignalEventsThisSliceRooVar,
                        expectedTotalSignalFractionRooVar));//,numberOfTrialsRooVar));
  //// temp for display
  //outputRootFile->cd();
  ////int backgroundTracksThisSample = (int)numBackgroundTracksThisSlice;
  //int backgroundTracksThisSample = 1e6;
  //RooDataSet* displayDataSet = iasBackgroundHistPdf->generate(rooVarIas,backgroundTracksThisSample);
  //RooDataHist* displayDataHist = displayDataSet->binnedClone();
  //RooPlot* iasFrame = rooVarIas.frame();
  //displayDataHist->plotOn(iasFrame);
  //displayDataHist->statOn(iasFrame,Format("NEU",AutoPrecision(1)));
  ////sigBackPdf.fitTo(*sampleData);
  ////sigBackPdf.plotOn(iasFrame,Components(*iasSignalHistPdf),LineStyle(kDashed),
  ////      LineColor(kRed));
  ////sigBackPdf.plotOn(iasFrame,Components(*iasBackgroundHistPdf),LineStyle(kDashed),
  ////      LineColor(kBlue));
  //iasBackgroundHistPdf->plotOn(iasFrame,LineColor(kBlue));
  ////iasBackgroundHistPdf->paramOn(iasFrame,Label("fit result"),Format("NEU",AutoPrecision(1)));
  ////sigBackPdf.paramOn(iasFrame,Format("NEU",AutoPrecision(1)));
  //std::string name="Plot of sample dataset";
  //std::string title="Hist PDF for sample dataset";
  //iasFrame->SetName(name.c_str());
  //iasFrame->SetTitle(title.c_str());
  //iasFrame->Write();
  //// histogram
  //TH1F* histPdfHist = (TH1F*)iasBackgroundHistPdf->createHistogram("iasBackgroundPdfHist",rooVarIas);
  //histPdfHist->Write();
  //std::cout << "expected events in last bin: " <<
  //  backgroundTracksThisSample*histPdfHist->GetBinContent(histPdfHist->GetNbinsX())*0.7/0.0625 << std::endl;
  //// output
  //const RooArgSet* displayArgSet = displayDataSet->get();
  //RooRealVar* iasDisplayData = (RooRealVar*)displayArgSet->find(rooVarIas.GetName());
  //for(int evt=0; evt < displayDataSet->numEntries(); ++evt)
  //{
  //  displayDataSet->get(evt);
  //  std::cout << "ias generated: " << iasDisplayData->getVal() << std::endl;
  //}
  //// Examine the statistics of a binned dataset
  //cout << ">> number of bins in    : " << displayDataHist->numEntries() << endl ;
  //cout << ">> sum of weights in    : " << displayDataHist->sum(kFALSE)  << endl ;
  //cout << ">> integral over histogram: " << displayDataHist->sum(kTRUE)   << endl ; // accounts for bin volume
  //// Locate a bin from a set of coordinates and retrieve its properties
  //rooVarIas = 0.05 ;
  //cout << ">> retrieving the properties of the bin enclosing coordinate rooVarIas = (0.05) " << endl ;
  //cout << " bin center:" << endl ;
  //displayDataHist->get(RooArgSet(rooVarIas))->Print("v") ; // load bin center coordinates in internal buffer
  //cout << " weight = " << displayDataHist->weight() << endl ; // return weight of last loaded coordinates
  //cout << " bin volume = " << displayDataHist->binVolume(RooArgSet(rooVarIas)) << endl ; // return bin volume
  // Locate a bin from a set of coordinates and retrieve its properties
  //rooVarIas = 0.6 ;
  //cout << ">> retrieving the properties of the bin enclosing coordinate rooVarIas = (0.6) " << endl ;
  //cout << " bin center:" << endl ;
  //displayDataHist->get(RooArgSet(rooVarIas))->Print("v") ; // load bin center coordinates in internal buffer
  //cout << " weight = " << displayDataHist->weight() << endl ; // return weight of last loaded coordinates
  //cout << " bin volume = " << displayDataHist->binVolume(RooArgSet(rooVarIas)) << endl ; // return bin volume
  //delete iasFrame;
  //delete displayDataHist;
  //delete displayDataSet;


  TRandom3 myRandom;
  TDirectory* experimentsDir = outputRootFile->mkdir("experiments");

  // gen and fit for each trial
  for(int i=0; i<numTrials; ++i)
  {
    if(reportEvery_!=0 ? ((i+1)>0 && (i+1)%reportEvery_==0) : false)
      std::cout << "  generating/fitting for trial: " << i+1 << " of " << numTrials << std::endl;
    int backgroundTracksThisSample = (int)numBackgroundTracksThisSlice;
    RooDataSet* thisSampleDataSet = iasBackgroundHistPdf->generate(rooVarIas,backgroundTracksThisSample);
    RooDataSet* thisSampleSignalDataSet;
    if(!generateBackgroundOnly)
    {
      // fluctuate num signal tracks by poisson
      int signalTracksThisSample = myRandom.Poisson(numSignalTracksThisSlice);
      thisSampleSignalDataSet = iasSignalHistPdf->generate(rooVarIas,signalTracksThisSample);
      thisSampleDataSet->append(*thisSampleSignalDataSet);
      RooDataHist* signalDataHist = thisSampleSignalDataSet->binnedClone();
      rooVarIas = 0.99;
      signalDataHist->get(RooArgSet(rooVarIas));//->Print("v");
      signalEventsInLastBinRooVar = signalDataHist->weight();
    }
    RooDataHist* sampleData = thisSampleDataSet->binnedClone();
    // sic nov 24: keep nll values only, don't write datahists
    //experimentsDir->cd();
    //std::string sampleDataName = getHistNameBeg(lowerNoM,lowerEta);
    //sampleDataName+="sampleDataTrial";
    //sampleDataName+=intToString(i);
    //sampleData->SetName(sampleDataName.c_str());
    //sampleData->Write();
    trialIndexRooVar = i;

    //const RooDataSet* thisData = mcStudy.genData(i);
    //RooDataHist* sampleData = new RooDataHist("sampleData","sampleData",rooVarIas,*thisData);
    //TH1F* sampleDataHist = (TH1F*)sampleData->createHistogram("sampleDataHist",rooVarIas,Binning(thisBinning));
    //sampleDataHist->Write();
    //RooDataSet* sampleData = iasBackgroundHistPdf->generate(rooVarIas,2000);
    //RooDataHist* sampleData = iasBackgroundHistPdf->generateBinned(rooVarIas,
    //    10799);//,Extended(true));
    //sampleData->add(*iasSignalHistPdf->generateBinned(rooVarIas,
    //      11,Extended(true)));
    //numTracksPerTrialHist->Fill(thisData->sumEntries());
    numTracksPerTrialHist->Fill(sampleData->sumEntries());

    // loop over signal fractions
    //TODO make this configurable in some way
    int maxTestNumSigTracks = 25;
    for(int testNumTotalSigTracks = 0; testNumTotalSigTracks < maxTestNumSigTracks; 
        ++testNumTotalSigTracks)
    {
      //float testSigFracTot = testNumSigTracks/(double)(numBackgroundTracksInDRegionInput+testNumSigTracks);
      float testSignalTracksThisSlice =
        fractionOfSigTracksInThisSlice*fractionOfSigTracksPassingMassCutThisSlice*testNumTotalSigTracks;
      float testTotalNumTracksThisSlice = testSignalTracksThisSlice+numBackgroundTracksThisSlice;
      float testSigFracThisSlice = testSignalTracksThisSlice/(double)testTotalNumTracksThisSlice;
      // saturated model
      RooHistPdf saturatedModel("saturatedModel","saturatedModel",rooVarIas,*sampleData);
      // set signal frac
      fsig = testSigFracThisSlice;
      // NLLs
      //RooNLLVar nllVarSOnly("nllVarSOnly","NLL S only",*iasSignalHistPdf,*sampleData);
      RooNLLVar nllVarBOnly("nllVarBOnly","NLL B only",*iasBackgroundHistPdf,*sampleData);
      RooNLLVar nllVarSB("nllVarSB","NLL SB",sigBackPdf,*sampleData);
      RooNLLVar nllVarSatModel("nllVarSatModel","NLL Saturated Model",saturatedModel,*sampleData);

      //TODO ? hists
      //nllBSatHist->Fill(2*(nllVarBOnly.getVal()-nllVarSatModel.getVal()));
      //nllSatOnlyHist->Fill(nllVarSatModel.getVal());
      //nllBOnlyHist->Fill(nllVarBOnly.getVal());
      ////nllSOnlyHist->Fill(nllVarSOnly.getVal());
      //nllSBHist->Fill(nllVarSB.getVal());
      //nllSBOverBHist->Fill(2*(nllVarSB.getVal()-nllVarBOnly.getVal()));
      //numSigTracksVsNLLSBOverBHist->Fill(2*(nllVarSB.getVal()-nllVarBOnly.getVal()),signalTracksThisSample);

      // roo vars for dataset
      //sigNLLRooVar = nllVarSOnly.getVal();
      backNLLRooVar = nllVarBOnly.getVal();
      sigBackNLLRooVar = nllVarSB.getVal();
      satModelNLLRooVar = nllVarSatModel.getVal();
      assumedTotalSigTracksRooVar = testNumTotalSigTracks;
      assumedSigTracksThisSliceRooVar = testSignalTracksThisSlice;
      // fit the S+B PDF
      sigBackPdf.fitTo(*sampleData,PrintLevel(-1));
      if(verbose_)
        std::cout << "preferred fsig after fitting: " << fsig.getVal() << std::endl;
      maxLSigFracRooVar = fsig.getVal();
      RooNLLVar nllVarSBFit("nllVarSBFit","NLL SB Fit",sigBackPdf,*sampleData);
      sigBackFitNLLRooVar = nllVarSBFit.getVal();
      //numTracksInLastBinVsSBOverBHist->Fill(2*(nllVarSB.getVal()-nllVarBOnly.getVal()),
      //    actualEventsInLastBinRooVar.getVal());
      nllValues->add(RooArgSet(sigNLLRooVar,backNLLRooVar,sigBackNLLRooVar,
            satModelNLLRooVar,sigBackFitNLLRooVar,maxLSigFracRooVar,
            assumedTotalSigTracksRooVar,assumedSigTracksThisSliceRooVar,trialIndexRooVar));
    }

    // reset signal frac
    fsig = sigFracThisSlice;
    rooVarIas = 0.99;
    sampleData->get(RooArgSet(rooVarIas));
    actualEventsInLastBinRooVar = sampleData->weight();
    observedEvents->add(RooArgSet(actualEventsInLastBinRooVar,signalEventsInLastBinRooVar)); 

    if(verbose_)
    {
      std::cout << "Number of events in last bin: " << actualEventsInLastBinRooVar.getVal() << std::endl;
      std::cout << "Expected BG events in last bin: " << expectedBackgroundEventsInLastBinRooVar.getVal() << std::endl;
      std::cout << "Expected sig events in last bin: " << expectedSignalEventsInLastBinRooVar.getVal() << std::endl;
      ////std::cout << "NLL S ONLY: " << sigNLLRooVar.getVal() << std::endl;
      //std::cout << "NLL SB: " << sigBackNLLRooVar.getVal() << std::endl;
      //std::cout << "NLL SBfit: " << sigBackFitNLLRooVar.getVal() << std::endl;
      //std::cout << "NLL B ONLY: " << backNLLRooVar.getVal() << std::endl;
      //std::cout << "NLL saturated model: " << satModelNLLRooVar.getVal() << std::endl;
      //std::cout << "chi2 = 2*(NLL_B-Nll_sat) = " << 2*(backNLLRooVar.getVal()-satModelNLLRooVar.getVal()) << std::endl;
      //std::cout << "-2 ln SB/B = " << 2*(sigBackNLLRooVar.getVal()-backNLLRooVar.getVal()) << std::endl;
      //if(maxLSigFracRooVar.getVal() <= sigFracThisSlice)
      //  std::cout << "-2 ln L(mu)/L(mu_hat) = " << 2*(sigBackNLLRooVar.getVal()-sigBackFitNLLRooVar.getVal()) << std::endl;
      //else
      //  std::cout << "-2 ln L(mu)/L(mu_hat) = 0" << std::endl;
      std::cout << "sigFracThisSlice = " << sigFracThisSlice << std::endl;
      std::cout << "--------------------------------------------------------------------------" << std::endl;
    }

    if(i==0)
    {
      experimentsDir->cd();
      RooPlot* iasFrame = rooVarIas.frame();
      sampleData->plotOn(iasFrame);
      //sampleData->statOn(iasFrame,Label("data"),Format("NEU",AutoPrecision(1)));
      iasBackgroundHistPdf->plotOn(iasFrame,LineColor(kGreen));
      sigBackPdf.plotOn(iasFrame,Components(*iasSignalHistPdf),LineStyle(kDashed),
          LineColor(kRed));
      sigBackPdf.plotOn(iasFrame,Components(*iasBackgroundHistPdfForSBModel),LineStyle(kDashed),
          LineColor(kBlue));
      //iasBackgroundHistPdf->paramOn(iasFrame,Label("fit result"),Format("NEU",AutoPrecision(1)));
      //sigBackPdf.paramOn(iasFrame,Format("NEU",AutoPrecision(1)));
      std::string name = string(iasBackgroundHistPdf->GetName());
      std::string title = string(iasBackgroundHistPdf->GetTitle());
      name+="Plot";
      title+=" Hist PDF";
      iasFrame->SetName((name+intToString(i)).c_str());
      iasFrame->SetTitle((title+intToString(i)).c_str());
      iasFrame->Write();
      delete iasFrame;
    }

    delete sampleData;
    delete thisSampleDataSet;
    //delete thisData;
  }


  // Write the hists
  outputRootFile->cd();
  TDirectory* obsDatasetDir = outputRootFile->mkdir("observedEventsDatasets");
  obsDatasetDir->cd();
  observedEvents->Write();
  TDirectory* expDatasetDir = outputRootFile->mkdir("expectedEventsDatasets");
  expDatasetDir->cd();
  expectedEvents->Write();
  TDirectory* nllDatasetDir = outputRootFile->mkdir("nllDatasets");
  nllDatasetDir->cd();
  nllValues->Write();
  TDirectory* nllBSatDir = outputRootFile->mkdir("NllBOverSat");
  nllBSatDir->cd();
  nllBSatHist->Write();
  TDirectory* nllSatOnlyDir = outputRootFile->mkdir("NllSatOnly");
  nllSatOnlyDir->cd();
  nllSatOnlyHist->Write();
  TDirectory* nllBOnlyDir = outputRootFile->mkdir("NllBOnly");
  nllBOnlyDir->cd();
  nllBOnlyHist->Write();
  TDirectory* nllSOnlyDir = outputRootFile->mkdir("NllSOnly");
  nllSOnlyDir->cd();
  nllSOnlyHist->Write();
  TDirectory* nllSBDir = outputRootFile->mkdir("NllSB");
  nllSBDir->cd();
  nllSBHist->Write();
  TDirectory* nllSBOverBDir = outputRootFile->mkdir("NllSBOverB");
  nllSBOverBDir->cd();
  nllSBOverBHist->Write();
  TDirectory* iasBackgroundPredDir = outputRootFile->mkdir("iasBackgroundPredictions");
  iasBackgroundPredDir->cd();
  iasBackgroundPredictionMassCutNoMSliceHist->Write();
  iasBackgroundDataHist->Write();
  iasBackgroundHistPdf->Write();
  TDirectory* iasSignalDir = outputRootFile->mkdir("iasSignalPredictions");
  iasSignalDir->cd();
  iasSignalMassCutNoMSliceHist->Write();
  iasSignalDataHist->Write();
  iasSignalHistPdf->Write();
  TDirectory* tracksPerTrialDir = outputRootFile->mkdir("tracksPerTrial");
  tracksPerTrialDir->cd();
  numTracksPerTrialHist->Write();
  TDirectory* numSignalTracksVsNLLSBOverBDir = outputRootFile->mkdir("numSignalTracksVsNLLSBOverB");
  numSignalTracksVsNLLSBOverBDir->cd();
  numSigTracksVsNLLSBOverBHist->Write();
  TDirectory* numTracksLastBinVsNLLSBOverBDir = outputRootFile->mkdir("numTracksInLastBinVsNLLSBOverB");
  numTracksLastBinVsNLLSBOverBDir->cd();
  numTracksInLastBinVsSBOverBHist->Write();


  outputRootFile->Close();

  //std::cout << "Events with Ih >= 3.5 MeV/cm and P >= 50 GeV: " << eventsInHighIhHighPRegion
  //  << "; also passing mass >= " << massCutIasHighPHighIh_ << " : " << eventsInHighIhHighPRegionWithMassCut
  //  << std::endl;
}

