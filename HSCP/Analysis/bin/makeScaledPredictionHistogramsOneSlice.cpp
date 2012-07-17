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
#include <fstream>
#include <set>
#include <limits>

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

// globals
int minNoM = 0;
int maxNoM = 0;
float minEta = 0;
float maxEta = 0;
int numIasBins = 0;

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


// helper functions
int getEtaSliceFromLowerEta(float lowerEta)
{
  return lowerEta/0.2;
}
//
int getNoMSliceFromLowerNoM(int lowerNoM)
{
  return (lowerNoM-5)/2;
}
//
int getGlobalBinFromNomSliceEtaSliceIasBin(int nomSlice, int etaSlice, int iasBin)
{
  // ias bin will be from 1-d hist, so it goes from 1-400

  //int numNoMBins = getNoMSliceFromLowerNoM(maxNoM) + 1;
  int numEtaBins = getEtaSliceFromLowerEta(maxEta);

  int multFactorEtaSlice = numIasBins;
  int multFactorNoMSlice = multFactorEtaSlice*numEtaBins;

  return multFactorNoMSlice*nomSlice + multFactorEtaSlice*etaSlice + iasBin;
}
//
int getGlobalBinMax()
{
  int numNoMBins = getNoMSliceFromLowerNoM(maxNoM) + 1;
  int numEtaBins = getEtaSliceFromLowerEta(maxEta);

  int multFactorEtaSlice = numIasBins;
  int multFactorNoMSlice = multFactorEtaSlice*numEtaBins;

  return multFactorNoMSlice*(numNoMBins-1) + multFactorEtaSlice*(numEtaBins-1) + numIasBins;
}
// add the mass, p/pt, and ih SB thresholds to the output filename automatically
std::string generateFileNameEnd(double massCut, double pSideband, double ptSideband, bool usePt, double ihSideband,
    double iasSideband, bool useIas)
{
  std::string fileNameEnd = ".massCut";
  fileNameEnd+=intToString(((int)massCut));
  if(usePt)
  {
    fileNameEnd+=".pt";
    fileNameEnd+=intToString(((int)ptSideband));
  }
  else
  {
    fileNameEnd+=".p";
    fileNameEnd+=intToString(((int)pSideband));
  }
  if(useIas)
  {
    fileNameEnd+=".ias";
    // replace "." with "p" for the ih
    std::string iasSBStr = floatToString(iasSideband);
    iasSBStr.replace(iasSBStr.find("."),1,"p");
    fileNameEnd+=iasSBStr;
  }
  else
  {
    fileNameEnd+=".ih";
    // replace "." with "p" for the ih
    std::string ihSBStr = floatToString(ihSideband);
    ihSBStr.replace(ihSBStr.find("."),1,"p");
    fileNameEnd+=ihSBStr;
  }

  fileNameEnd+=".root";
  return fileNameEnd;
}



// ****** main
int main(int argc, char ** argv)
{
  using namespace RooFit;
  using namespace std;

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // parse arguments
  if(argc < 2)
  {
    cout << "Usage : " << argv[0] << " [parameters.py]" << endl;
    return -1;
  }

  // get the python configuration
  const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  // now get each parameterset
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("makeScaledPredictionHistograms");
  // mass cut to use for the high-p high-Ih (search region) ias dist
  double massCut_ (ana.getParameter<double>("MassCut"));
  string backgroundPredictionRootFilename_ (ana.getParameter<string>("BackgroundPredictionInputRootFile"));
  string signalTightRPCRootFilename_ (ana.getParameter<string>("SignalTightRPCInputRootFile"));
  string signalLooseRPCRootFilename_ (ana.getParameter<string>("SignalLooseRPCInputRootFile"));
  string outputRootFilename_ (ana.getParameter<string>("OutputRootFile"));
  // dE/dx calibration
  double dEdx_k (ana.getParameter<double>("dEdx_k"));
  double dEdx_c (ana.getParameter<double>("dEdx_c"));
  // definition of sidebands/search region
  double pSidebandThreshold (ana.getParameter<double>("PSidebandThreshold"));
  double ptSidebandThreshold (ana.getParameter<double>("PtSidebandThreshold"));
  bool usePtForSideband (ana.getParameter<bool>("UsePtForSideband"));
  double ihSidebandThreshold (ana.getParameter<double>("IhSidebandThreshold"));
  double integratedLumi (ana.getParameter<double>("IntegratedLuminosity")); // 1/pb
  double signalCrossSectionForEff (ana.getParameter<double>("SignalCrossSection")); // pb
  double iasCutForEffAcc (ana.getParameter<double>("IasCutForEfficiency"));
  bool useIasForSideband (ana.getParameter<bool>("UseIasForSideband"));

  // do mass?
  bool doMass = false;

  std::cout << "Using dE/dx k: " << dEdx_k << " c: " << dEdx_c << std::endl;

  // TODO configurable nom/eta limits
  // NB: always use upper edge of last bin for both
  minNoM = 5;
  maxNoM = 18;
  minEta = 0.0;
  maxEta = 1.6;


  // define observables
  //TODO move to something like common functions
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",0,2.5);
  RooRealVar rooVarRun("rooVarRun","run",0,4294967295);
  RooRealVar rooVarLumiSection("rooVarLumiSection","lumiSection",0,4294967295);
  RooRealVar rooVarEvent("rooVarEvent","event",0,4294967295);
  //TODO add mass in initial dataset generation?
  RooRealVar rooVarNumGenHSCPEvents("rooVarNumGenHSCPEvents","numGenHSCPEvents",0,1e10);
  RooRealVar rooVarNumGenHSCPTracks("rooVarNumGenHSCPTracks","numGenHSCPTracks",0,1e10);
  RooRealVar rooVarNumGenChargedHSCPTracks("rooVarNumGenChargedHSCPTracks","numGenChargedHSCPTracks",0,1e10);
  RooRealVar rooVarSignalEventCrossSection("rooVarSignalEventCrossSection","signalEventCrossSection",0,100); // pb
  RooRealVar rooVarNumGenHSCPEventsPUReweighted("rooVarNumGenHSCPEventsPUReweighted","numGenHSCPEventsPUReweighted",0,1e10);
  RooRealVar rooVarEventWeightSum("rooVarEventWeightSum","eventWeightSum",0,1e10);
  RooRealVar rooVarSampleWeight("rooVarSampleWeight","sampleWeight",0,1e10);
  //
  RooRealVar rooVarPUWeight("rooVarPUWeight","puWeight",-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
  RooRealVar rooVarPUSystFactor("rooVarPUSystFactor","puSystFactor",-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());

  //string fileNameEnd = generateFileNameEnd(massCut_,pSidebandThreshold,ptSidebandThreshold,usePtForSideband,ihSidebandThreshold);
  string fileNameEnd = "";
  TFile* outputRootFile = new TFile((outputRootFilename_+fileNameEnd).c_str(),"recreate");
  TFile* backgroundPredictionRootFile = TFile::Open(backgroundPredictionRootFilename_.c_str());
  TFile* signalTightRPCRootFile = TFile::Open(signalTightRPCRootFilename_.c_str());
  TFile* signalLooseRPCRootFile = TFile::Open(signalLooseRPCRootFilename_.c_str());

  // get roodataset from signal file -- tight RPC (period 0)
  //RooDataSet* rooDataSetAllSignal = (RooDataSet*)signalRootFile->Get("rooDataSetCandidates");
  RooDataSet* rooDataSetAllSignalTightRPC = (RooDataSet*)signalTightRPCRootFile->Get("rooDataSetOneCandidatePerEvent");
  bool countEvents  = true;
  if(!rooDataSetAllSignalTightRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetCandidates in signal file " <<
      signalTightRPCRootFilename_.c_str() << endl;
    return -3;
  }
  // get roodataset with gen information from signal file
  RooDataSet* rooDataSetGenSignalTightRPC = (RooDataSet*)signalTightRPCRootFile->Get("rooDataSetGenHSCPTracks");
  if(!rooDataSetGenSignalTightRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetGenSignal in signal file " <<
      signalTightRPCRootFilename_.c_str() << endl;
    return -3;
  }
  // get roodataset with pileup weight info
  RooDataSet* rooDataSetSignalPileupWeightsTightRPC = (RooDataSet*)signalTightRPCRootFile->Get("rooDataSetPileupWeights");
  if(!rooDataSetSignalPileupWeightsTightRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetPileupWeights in signal file " <<
      signalTightRPCRootFilename_.c_str() << endl;
    return -3;
  }
  // get roodataset from signal file -- loose RPC (period 1)
  RooDataSet* rooDataSetAllSignalLooseRPC = (RooDataSet*)signalLooseRPCRootFile->Get("rooDataSetOneCandidatePerEvent");
  if(!rooDataSetAllSignalLooseRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetCandidates in signal file " <<
      signalLooseRPCRootFilename_.c_str() << endl;
    return -3;
  }
  // get roodataset with gen information from signal file
  RooDataSet* rooDataSetGenSignalLooseRPC = (RooDataSet*)signalLooseRPCRootFile->Get("rooDataSetGenHSCPTracks");
  if(!rooDataSetGenSignalLooseRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetGenSignal in signal file " <<
      signalLooseRPCRootFilename_.c_str() << endl;
    return -3;
  }
  // get roodataset with pileup weight info
  RooDataSet* rooDataSetSignalPileupWeightsLooseRPC = (RooDataSet*)signalLooseRPCRootFile->Get("rooDataSetPileupWeights");
  if(!rooDataSetSignalPileupWeightsLooseRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetPileupWeights in signal file " <<
      signalLooseRPCRootFilename_.c_str() << endl;
    return -3;
  }

  //int numSignalTracksTotal = rooDataSetAllSignalLooseRPC->numEntries()+rooDataSetAllSignalTightRPC->numEntries();
  // construct D regions for signal
  string pSearchCutString = "rooVarP>";
  pSearchCutString+=floatToString(pSidebandThreshold);
  string ptSearchCutString = "rooVarPt>";
  ptSearchCutString+=floatToString(ptSidebandThreshold);
  string ihSearchCutString = "rooVarIh>";
  ihSearchCutString+=floatToString(ihSidebandThreshold);
  string iasSearchCutString = "rooVarIas>";
  iasSearchCutString+=floatToString(iasCutForEffAcc);
  RooDataSet* regionDDataSetSignalLooseRPC;
  RooDataSet* regionDDataSetSignalTightRPC;
  // SIC MAR 5 -- don't do P cut for signal D region -- will be taken care of by mass cut later
  if(useIasForSideband)
    regionDDataSetSignalLooseRPC = (RooDataSet*)rooDataSetAllSignalLooseRPC->reduce(Cut(iasSearchCutString.c_str()));
  else
    regionDDataSetSignalLooseRPC = (RooDataSet*)rooDataSetAllSignalLooseRPC->reduce(Cut(ihSearchCutString.c_str()));
  if(useIasForSideband)
    regionDDataSetSignalTightRPC = (RooDataSet*)rooDataSetAllSignalTightRPC->reduce(Cut(iasSearchCutString.c_str()));
  else
    regionDDataSetSignalTightRPC = (RooDataSet*)rooDataSetAllSignalTightRPC->reduce(Cut(ihSearchCutString.c_str()));
  //int numSignalTracksInDRegion = regionDDataSetSignalTightRPC->numEntries()+regionDDataSetSignalLooseRPC->numEntries();

  // get the background A region entries hist
  string fullPath = "entriesInARegion";
  TH2F* aRegionBackgroundEntriesHist = (TH2F*)backgroundPredictionRootFile->Get(fullPath.c_str());
  if(!aRegionBackgroundEntriesHist)
  {
    cout << "Cannot proceed further: no A region entries hist found (histogram " <<
      fullPath << ") in file: " << backgroundPredictionRootFilename_ << ")." << endl;
    return -3;
  }
  // get the background B region hist
  fullPath = "entriesInBRegion";
  TH2F* bRegionBackgroundEntriesHist = (TH2F*)backgroundPredictionRootFile->Get(fullPath.c_str());
  if(!bRegionBackgroundEntriesHist)
  {
    cout << "Cannot proceed further: no B region entries hist found (histogram " <<
      fullPath << ") in file: " << backgroundPredictionRootFilename_ << ")." << endl;
    return -3;
  }
  // get the background C region hist
  fullPath = "entriesInCRegion";
  TH2F* cRegionBackgroundEntriesHist = (TH2F*)backgroundPredictionRootFile->Get(fullPath.c_str());
  if(!cRegionBackgroundEntriesHist)
  {
    cout << "Cannot proceed further: no C region entries hist found (histogram " <<
      fullPath << ") in file: " << backgroundPredictionRootFilename_ << ")." << endl;
    return -3;
  }

  // can remove these and use only one loop
  vector<TH1D> bgLimitsHistsToUse;
  vector<TH1D> bgDiscoveryHistsToUse;
  vector<TH1D> dataHistsToUse;
  //vector<TH1F> bRegionHistsToUse;
  //vector<TH1F> cRegionHistsToUse;

  // for histograms already made
  //string dirName="iasPredictionsVariableBins";
  //string getHistName+="iasPredictionVarBinHist";
  string dirName="iasPredictionsFixedBins";
  string bgLimitsHistNameEnd="iasPredictionFixedLimitsHist";
  string bgDiscoveryHistNameEnd="iasPredictionFixedDiscoveryHist";
  string dataDirName="dRegionFixedBins";
  string dataHistNameEnd="dRegionFixedHist";

  // look at BG predictions to determine unrolled hist binning
  // if stepping made variable, must change getHistNameBeg calls
  for(int lowerNoM = minNoM; lowerNoM < maxNoM; lowerNoM+=2)
  {
    for(float lowerEta = minEta; lowerEta < maxEta; lowerEta+=0.2)
    {
      TH1D* iasBackgroundPredictionMassCutNoMSliceLimitsHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSliceDiscoveryHist = 0;
      //TH1F* bRegionHist = 0;
      //TH1F* cRegionHist = 0;
      // get this eta/nom hist
      string getLimitsHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
      string getDiscoveryHistName = getLimitsHistName;
      getLimitsHistName+=bgLimitsHistNameEnd;
      getDiscoveryHistName+=bgDiscoveryHistNameEnd;
      string fullPathLimits = dirName;
      fullPathLimits+="/";
      string fullPathDiscovery = fullPathLimits;
      fullPathLimits+=getLimitsHistName;
      fullPathDiscovery+=getDiscoveryHistName;

      TH1D* dRegionDataHist = 0;
      string dataHistGetName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
      dataHistGetName+=dataHistNameEnd;
      string fullDataHistPath = dataDirName;
      fullDataHistPath+="/";
      fullDataHistPath+=dataHistGetName;

      iasBackgroundPredictionMassCutNoMSliceLimitsHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathLimits.c_str());
      iasBackgroundPredictionMassCutNoMSliceDiscoveryHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathDiscovery.c_str());
        
      dRegionDataHist = (TH1D*)backgroundPredictionRootFile->Get(fullDataHistPath.c_str());

      // do the rest if the BG hist is found, which should always be the case (so we can remove these vectors)
      if(iasBackgroundPredictionMassCutNoMSliceLimitsHist && iasBackgroundPredictionMassCutNoMSliceDiscoveryHist)
      {
        // set ias bins or check hist consistency
        if(numIasBins < 1)
          numIasBins = iasBackgroundPredictionMassCutNoMSliceLimitsHist->GetNbinsX();
        else if(numIasBins != iasBackgroundPredictionMassCutNoMSliceLimitsHist->GetNbinsX())
        {
          cout << "ERROR: A histogram earlier had : " << numIasBins << " x bins and this one ("
            << iasBackgroundPredictionMassCutNoMSliceLimitsHist->GetName() << ") has : "
            << iasBackgroundPredictionMassCutNoMSliceLimitsHist->GetNbinsX() << " x bins...exiting."
            << endl;
          return -10;
        }
        //// debug
        //cout << "Found hist: " << iasBackgroundPredictionMassCutNoMSliceHist->GetName()
        //  << " with integral = " << iasBackgroundPredictionMassCutNoMSliceHist->Integral()
        //  << " pushing onto the vector." << endl;
        bgLimitsHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceLimitsHist);
        bgDiscoveryHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceDiscoveryHist);
        dataHistsToUse.push_back(*dRegionDataHist);
        // now get the b region hist
        //string getHistName = dirName;
        //getHistName+="/";
        //getHistName+=getHistNameBeg(lowerNoM,lowerEta);
        //getHistName+="bRegionHist";
        //bRegionHist = (TH1F*)backgroundPredictionRootFile->Get(getHistName.c_str());
        //if(bRegionHist && bRegionHist->Integral() > 0)
        //  bRegionHistsToUse.push_back(*bRegionHist);
        //// now get the c region hist
        //getHistName = dirName;
        //getHistName+="/";
        //getHistName+=getHistNameBeg(lowerNoM,lowerEta);
        //getHistName+="cRegionHist";
        //cRegionHist = (TH1F*)backgroundPredictionRootFile->Get(getHistName.c_str());
        //if(cRegionHist && cRegionHist->Integral() > 0)
        //  cRegionHistsToUse.push_back(*cRegionHist);
      }
      //else
      //{
      //  cout << "WARNING: no input found for this slice : lim or disc hist missing: " <<
      //    getHistNameBeg(lowerNoM,lowerEta) << " in file: " << backgroundPredictionRootFilename_ << ")." << 
      //    " Excluded -- should not happen. " << endl;
      //  return -11;
      //}
      delete iasBackgroundPredictionMassCutNoMSliceLimitsHist;
      delete iasBackgroundPredictionMassCutNoMSliceDiscoveryHist;
      //delete bRegionHist;
      //delete cRegionHist;
    }
  }

  assert(bgLimitsHistsToUse.size() == bgDiscoveryHistsToUse.size());

  // initialize the global unrolled hists
  int binsToExclude = ceil(numIasBins*iasCutForEffAcc);
  int numBinsEachSlice = numIasBins-binsToExclude;
  int numGlobalBins = numBinsEachSlice*bgLimitsHistsToUse.size();
  int numMassBins = 1000;
  int numGlobalMassBins = numMassBins*bgLimitsHistsToUse.size();
  cout << "INFO: numIasBins = " << numIasBins << endl;
  cout << "INFO: numBGLimitsHistsToUse = " << bgLimitsHistsToUse.size() << endl;
  cout << "INFO: numBGDiscoveryHistsToUse = " << bgDiscoveryHistsToUse.size() << endl;
  cout << "INFO: initializing unrolled histogram: bins=" << numGlobalBins << endl << endl;
  TH1D* backgroundAllNoMAllEtaUnrolledLimitsHist = new TH1D("backgroundAllNoMAllEtaUnrolledLimitsHist",
      "unrolled background hist (limits)",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledDiscoveryHist = new TH1D("backgroundAllNoMAllEtaUnrolledDiscoveryHist",
      "unrolled background hist (discovery)",numGlobalBins,1,numGlobalBins+1);
  //TH1D* backgroundAllNoMAllEtaUnrolledPlusOneSigmaHist = new TH1D("backgroundAllNoMAllEtaUnrolledPlusOneSigmaHist",
  //    "unrolled background hist (shape +1 sigma)",numGlobalBins,1,numGlobalBins+1);
  //TH1D* backgroundAllNoMAllEtaUnrolledMinusOneSigmaHist = new TH1D("backgroundAllNoMAllEtaUnrolledMinusOneSigmaHist",
  //    "unrolled background hist (shape -1 sigma)",numGlobalBins,1,numGlobalBins+1);
  TH1D* signalAllNoMAllEtaUnrolledHist = new TH1D("signalAllNoMAllEtaUnrolledHist",
      "unrolled signal hist",numGlobalBins,1,numGlobalBins+1);
  TH1D* dataAllNoMAllEtaUnrolledHist = new TH1D("dataAllNoMAllEtaUnrolledHist",
      "unrolled data hist, D region",numGlobalBins,1,numGlobalBins+1);
  // mass
  TH1D* backgroundAllNoMAllEtaUnrolledMassHist = new TH1D("backgroundAllNoMAllEtaUnrolledMassHist",
      "unrolled background hist (mass)",numGlobalMassBins,1,numGlobalMassBins+1);
  TH1D* signalAllNoMAllEtaUnrolledMassHist = new TH1D("signalAllNoMAllEtaUnrolledMassHist",
      "unrolled signal hist (mass)",numGlobalMassBins,1,numGlobalMassBins+1);
  TH2F* sigEffOverIasCutVsSliceHist = new TH2F("sigEffOverIasCutVsSlice","Sig eff. over ias cut;#eta;nom",12,0,2.4,9,5,23);
  sigEffOverIasCutVsSliceHist->GetYaxis()->SetNdivisions(509,false);
  TH1F* sigEffOverIasCutHist = new TH1F("sigEffOverIasCut","Signal eff. over ias cut",100,0,1);
  TH2F* backExpOverIasCutVsSliceHist = new TH2F("backExpOverIasCutVsSlice","Exp. background over ias cut;#eta;nom",12,0,2.4,9,5,23);
  backExpOverIasCutVsSliceHist->GetYaxis()->SetNdivisions(509,false);
  TH2F* backExpVsSliceHist = new TH2F("backExpVsSlice","Exp. background;#eta;nom",12,0,2.4,9,5,23);
  backExpVsSliceHist->GetYaxis()->SetNdivisions(509,false);
  TH1F* backExpOverIasCutHist = new TH1F("backExpOverIasCut","Exp. background over ias cut",100,0,1);
  TH2F* entriesSignalHist = new TH2F("entriesSignal","Signal tracks over mass cut in slice;#eta;nom",12,0,2.4,9,5,23);
  entriesSignalHist->GetYaxis()->SetNdivisions(509,false);
  TH2F* entriesSignalIasHist = new TH2F("entriesSignalIas","Signal tracks over mass/ias cuts in slice;#eta;nom",12,0,2.4,9,5,23);
  entriesSignalIasHist->GetYaxis()->SetNdivisions(509,false);

  double backgroundTracksOverIasCut = 0;
  double backgroundTracksOverIasCutErrorSqr = 0;
  double backgroundTracksOverIasCutNoApprox = 0;
  double backgroundTracksOverIasCutNoApproxErrorSqr = 0;
  int numBkOverIasCut = 0;
  double signalTracksOverIasCut = 0;
  double signalTracksLooseRPCTotal = 0;
  double signalTracksTightRPCTotal = 0;
  double signalEventsOverIasCutLooseRPC = 0;
  double signalEventsOverIasCutTightRPC = 0;
  double numSignalTracksInDRegionPassingMassCut = 0;
  set<EventInfo> selectedSignalEventsTightRPCSet;
  set<EventInfo> selectedSignalEventsLooseRPCSet;
  const RooArgSet* genArgSetLooseRPC = rooDataSetGenSignalLooseRPC->get();
  RooRealVar* numGenHSCPTracksLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarNumGenHSCPTracks.GetName());
  //RooRealVar* numGenHSCPEventsLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarNumGenHSCPEvents.GetName());
  //RooRealVar* signalCrossSectionLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarSignalEventCrossSection.GetName());
  //RooRealVar* numGenHSCPEventsPUReweightedLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarNumGenHSCPEventsPUReweighted.GetName());
  RooRealVar* eventWeightSumLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarEventWeightSum.GetName());
  RooRealVar* sampleWeightLooseRPCRooVar = (RooRealVar*)genArgSetLooseRPC->find(rooVarSampleWeight.GetName());
  rooDataSetGenSignalLooseRPC->get(0);
  signalTracksLooseRPCTotal = numGenHSCPTracksLooseRPCRooVar->getVal();

  const RooArgSet* genArgSetTightRPC = rooDataSetGenSignalTightRPC->get();
  RooRealVar* numGenHSCPTracksTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarNumGenHSCPTracks.GetName());
  //RooRealVar* numGenHSCPEventsTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarNumGenHSCPEvents.GetName());
  //RooRealVar* signalCrossSectionTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarSignalEventCrossSection.GetName());
  //RooRealVar* numGenHSCPEventsPUReweightedTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarNumGenHSCPEventsPUReweighted.GetName());
  RooRealVar* eventWeightSumTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarEventWeightSum.GetName());
  RooRealVar* sampleWeightTightRPCRooVar = (RooRealVar*)genArgSetTightRPC->find(rooVarSampleWeight.GetName());
  rooDataSetGenSignalTightRPC->get(0);
  signalTracksTightRPCTotal = numGenHSCPTracksTightRPCRooVar->getVal();
  //double totalGenHSCPTracks = signalTracksLooseRPCTotal+signalTracksTightRPCTotal;
  //double totalGenHSCPEvents = numGenHSCPEventsLooseRPCRooVar->getVal()+numGenHSCPEventsTightRPCRooVar->getVal();
  double totalGenHSCPEvents = eventWeightSumLooseRPCRooVar->getVal()+eventWeightSumTightRPCRooVar->getVal();
  cout << "event weight sum: " << totalGenHSCPEvents << endl;
  cout << "event weight sum loose: " << eventWeightSumLooseRPCRooVar->getVal() << " tight: " << eventWeightSumTightRPCRooVar->getVal() << endl;
  // get pileup-reweighted number of MC events -- loose RPC
  const RooArgSet* argSetPileupWeightsLooseRPC = rooDataSetSignalPileupWeightsLooseRPC->get();
  RooRealVar* eventNumPUDataLooseRPC = (RooRealVar*)argSetPileupWeightsLooseRPC->find(rooVarEvent.GetName());
  RooRealVar* lumiSecPUDataLooseRPC = (RooRealVar*)argSetPileupWeightsLooseRPC->find(rooVarLumiSection.GetName());
  RooRealVar* runNumPUDataLooseRPC = (RooRealVar*)argSetPileupWeightsLooseRPC->find(rooVarRun.GetName());
  RooRealVar* pileupWeightPUDataLooseRPC = (RooRealVar*)argSetPileupWeightsLooseRPC->find(rooVarPUWeight.GetName());
  //RooRealVar* pileupSystFactorPUDataLooseRPC = (RooRealVar*)argSetPileupWeightsLooseRPC->find(rooVarPUSystFactor.GetName());
  double sampleWeightLooseRPC = sampleWeightLooseRPCRooVar->getVal();
  // make map of pileup weights
  std::map<EventInfo,double> pileupWeightLooseRPCMap;
  std::map<EventInfo,double>::iterator pileupWeightLooseRPCMapItr;
  for(int idx=0; idx < rooDataSetSignalPileupWeightsLooseRPC->numEntries(); ++idx)
  {
    rooDataSetSignalPileupWeightsLooseRPC->get(idx);
    EventInfo evtInfo(runNumPUDataLooseRPC->getVal(),lumiSecPUDataLooseRPC->getVal(),eventNumPUDataLooseRPC->getVal());
    pileupWeightLooseRPCMap[evtInfo] = pileupWeightPUDataLooseRPC->getVal();
  }
  // get pileup-reweighted number of MC events -- tight RPC
  const RooArgSet* argSetPileupWeightsTightRPC = rooDataSetSignalPileupWeightsTightRPC->get();
  RooRealVar* eventNumPUDataTightRPC = (RooRealVar*)argSetPileupWeightsTightRPC->find(rooVarEvent.GetName());
  RooRealVar* lumiSecPUDataTightRPC = (RooRealVar*)argSetPileupWeightsTightRPC->find(rooVarLumiSection.GetName());
  RooRealVar* runNumPUDataTightRPC = (RooRealVar*)argSetPileupWeightsTightRPC->find(rooVarRun.GetName());
  RooRealVar* pileupWeightPUDataTightRPC = (RooRealVar*)argSetPileupWeightsTightRPC->find(rooVarPUWeight.GetName());
  //RooRealVar* pileupSystFactorPUDataTightRPC = (RooRealVar*)argSetPileupWeightsTightRPC->find(rooVarPUSystFactor.GetName());
  double sampleWeightTightRPC = sampleWeightTightRPCRooVar->getVal();
  // make map of pileup weights
  std::map<EventInfo,double> pileupWeightTightRPCMap;
  std::map<EventInfo,double>::iterator pileupWeightTightRPCMapItr;
  for(int idx=0; idx < rooDataSetSignalPileupWeightsTightRPC->numEntries(); ++idx)
  {
    rooDataSetSignalPileupWeightsTightRPC->get(idx);
    EventInfo evtInfo(runNumPUDataTightRPC->getVal(),lumiSecPUDataTightRPC->getVal(),eventNumPUDataTightRPC->getVal());
    pileupWeightTightRPCMap[evtInfo] = pileupWeightPUDataTightRPC->getVal();
  }


  // loop over limits hists to use
  int iteratorPos = 0;
  for(vector<TH1D>::iterator histItr = bgLimitsHistsToUse.begin();
      histItr != bgLimitsHistsToUse.end(); ++histItr)
  {
    // delete this and the matching delete at the bottom las
    //TH1F* iasBackgroundPredictionMassCutNoMSliceHist = 0;
    TH1D* iasSignalMassCutNoMSliceHist = 0;
    TH1D* iasSignalMassCutNoMSliceForEffHist = 0;
    string bgHistName = string(histItr->GetName());
    string bgHistNameBeg = bgHistName.substr(0,bgHistName.size()-bgLimitsHistNameEnd.size());
    int lowerNoM = getLowerNoMFromHistName(bgHistName,bgLimitsHistNameEnd);
    float lowerEta = getLowerEtaFromHistName(bgHistName,bgLimitsHistNameEnd);
    string sigHistName = bgHistNameBeg;
    //sigHistName+="iasSignalVarBinHist";
    sigHistName+="iasSignalFixedHist";
    string iasSignalHistTitle = "Ias of signal";
    iasSignalHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
    // get b region ias histogram
    string iasHistName = "bRegionHistograms/";
    iasHistName+=getHistNameBeg(lowerNoM,lowerEta,2,0.2);
    iasHistName+="iasBRegionHist";
    TH1D* iasHist = (TH1D*)backgroundPredictionRootFile->Get(iasHistName.c_str());
    // get c region tracks passing mass cut in ias bins
    string cRegionTracksOverMassCutProfileName="cRegionHistograms/";
    cRegionTracksOverMassCutProfileName+=getHistNameBeg(lowerNoM,lowerEta,2,0.2);
    cRegionTracksOverMassCutProfileName+="tracksInCeffOverMassCutProfile";
    TProfile* cRegionTracksOverMassCutProfile = (TProfile*)backgroundPredictionRootFile->Get(cRegionTracksOverMassCutProfileName.c_str());
    // debug
    //cout << " Looking at hist: " << bgHistName << " lowerNoM = " << lowerNoM
    //  << " lowerEta = " << lowerEta << endl;

    // copy binning and limits from background pred hist
    iasSignalMassCutNoMSliceHist =
      (TH1D*)histItr->Clone();
    iasSignalMassCutNoMSliceHist->Reset();
    iasSignalMassCutNoMSliceHist->SetName(sigHistName.c_str());
    iasSignalMassCutNoMSliceHist->SetTitle(iasSignalHistTitle.c_str());
    // mass
    TH1D* massBackgroundPredictionNoMSliceHist = 0;
    TH1D* massSignalNoMSliceHist = 0;
    if(doMass)
    {
      // bg mass hist
      string dirName="massPredictionsFixedBins";
      string bgHistNameEnd="massPredictionFixedHist";
      string getHistName = getHistNameBeg(lowerNoM,lowerEta,2,0.2);
      getHistName+=bgHistNameEnd;
      string fullPath = dirName;
      fullPath+="/";
      fullPath+=getHistName;
      massBackgroundPredictionNoMSliceHist = (TH1D*)backgroundPredictionRootFile->Get(fullPath.c_str());
      string sigMassHistName = bgHistNameBeg;
      sigMassHistName+="massSignalFixedHist";
      string signalMassHistTitle = "Mass of signal";
      signalMassHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,0);
      massSignalNoMSliceHist = new TH1D(sigMassHistName.c_str(),signalMassHistTitle.c_str(),1000,0,5000);
    }

    // construct D region for signal in this eta/nom slice
    //string nomCutString = "rooVarNoMias==";
    //nomCutString+=intToString(lowerNoM);
    //nomCutString+="||rooVarNoMias==";
    //nomCutString+=intToString(lowerNoM+1);
    //RooDataSet* nomCutDRegionDataSetSignal =
    //  (RooDataSet*)regionDDataSetSignal->reduce(Cut(nomCutString.c_str()));
    //std::string etaCutString = "(rooVarEta>";
    //etaCutString+=floatToString(lowerEta);
    //etaCutString+="&&rooVarEta<";
    //std::string upperEtaLimit;
    ////XXX TESTING eta 1.5 max
    //if(10*lowerEta==14)
    //  upperEtaLimit=floatToString(1.5);
    //else
    //  upperEtaLimit=floatToString(lowerEta+0.2);
    //etaCutString+=upperEtaLimit;
    //etaCutString+=")||(rooVarEta>-";
    //etaCutString+=upperEtaLimit;
    //etaCutString+="&&rooVarEta<-";
    //etaCutString+=floatToString(lowerEta);
    //etaCutString+=")";
    RooDataSet* etaCutNomCutDRegionDataSetSignalTightRPC =
      //(RooDataSet*)nomCutDRegionDataSetSignal->reduce(Cut(etaCutString.c_str()));
      (RooDataSet*)regionDDataSetSignalTightRPC->Clone();
    RooDataSet* etaCutNomCutDRegionDataSetSignalLooseRPC =
      //(RooDataSet*)nomCutDRegionDataSetSignal->reduce(Cut(etaCutString.c_str()));
      (RooDataSet*)regionDDataSetSignalLooseRPC->Clone();
    //int numSignalTracksInDRegionThisSlice = etaCutNomCutDRegionDataSetSignal->numEntries();

    // construct signal prediction hist (with mass cut) -- loose RPC
    const RooArgSet* argSetLooseRPC = etaCutNomCutDRegionDataSetSignalLooseRPC->get();
    RooRealVar* iasDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarIas.GetName());
    RooRealVar* nomIasDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarNoMias.GetName());
    RooRealVar* etaDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarEta.GetName());
    RooRealVar* ihDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarIh.GetName());
    RooRealVar* pDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarP.GetName());
    RooRealVar* ptDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarPt.GetName());
    RooRealVar* eventNumDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarEvent.GetName());
    RooRealVar* lumiSecDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarLumiSection.GetName());
    RooRealVar* runNumDataLooseRPC = (RooRealVar*)argSetLooseRPC->find(rooVarRun.GetName());
    for(int evt=0; evt < etaCutNomCutDRegionDataSetSignalLooseRPC->numEntries(); ++evt)
    {
      etaCutNomCutDRegionDataSetSignalLooseRPC->get(evt);
      if(fabs(etaDataLooseRPC->getVal()) > 1.5) continue;
      // find pileupWeight
      EventInfo evtInfo(runNumDataLooseRPC->getVal(),lumiSecDataLooseRPC->getVal(),eventNumDataLooseRPC->getVal());
      pileupWeightLooseRPCMapItr = pileupWeightLooseRPCMap.find(evtInfo);
      if(pileupWeightLooseRPCMapItr == pileupWeightLooseRPCMap.end())
      {
        std::cout << "ERROR: Could not find pileup weight for this event! Quitting." << std::endl;
        return -5;
      }
      double pileupWeight = pileupWeightLooseRPCMapItr->second;
      double eventWeight = sampleWeightLooseRPC*pileupWeight;

      // apply mass cut
      float massSqr = (ihDataLooseRPC->getVal()-dEdx_c)*pow(pDataLooseRPC->getVal(),2)/dEdx_k;
      if(massSqr < 0)
        continue;
      if(sqrt(massSqr) >= massCut_)
      {
        iasSignalMassCutNoMSliceHist->Fill(iasDataLooseRPC->getVal(),eventWeight);
        numSignalTracksInDRegionPassingMassCut+=eventWeight;
        entriesSignalHist->Fill(fabs(etaDataLooseRPC->getVal())+0.1,nomIasDataLooseRPC->getVal(),eventWeight);
        if(iasDataLooseRPC->getVal() > iasCutForEffAcc)
        {
          // if track over ias cut, put the event in the set (one track per event kept)
          // there should only be one track per event at this stage anyway
          pair<set<EventInfo>::iterator,bool> ret;
          ret = selectedSignalEventsLooseRPCSet.insert(evtInfo);
          if(ret.second)
          {
            //cout << "INFO: insert event - run: " << runNumDataLooseRPC->getVal() << " lumiSec: " << lumiSecDataLooseRPC->getVal()
            //  << " eventNum: " << eventNumDataLooseRPC->getVal() << endl;
            //cout << "add puweight=" << pileupWeight << "eventWeight=" << eventWeight << " to signalEventsOverIasCutLooseRPC; total=" << signalEventsOverIasCutLooseRPC << endl;
            signalEventsOverIasCutLooseRPC+=eventWeight;
            entriesSignalIasHist->Fill(fabs(etaDataLooseRPC->getVal())+0.1,nomIasDataLooseRPC->getVal(),eventWeight);
          }
        }
      }
      if(usePtForSideband)
      {
        if(ptDataLooseRPC->getVal() < ptSidebandThreshold)
          continue;
      }
      else
      {
        if(pDataLooseRPC->getVal() < pSidebandThreshold)
          continue;
      }
      // ih/ias cut already applied
      if(doMass)
        massSignalNoMSliceHist->Fill(sqrt(massSqr),eventWeight);
    }

    // construct signal prediction hist (with mass cut) -- tight RPC
    const RooArgSet* argSetTightRPC = etaCutNomCutDRegionDataSetSignalTightRPC->get();
    RooRealVar* iasDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarIas.GetName());
    RooRealVar* nomIasDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarNoMias.GetName());
    RooRealVar* etaDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarEta.GetName());
    RooRealVar* ihDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarIh.GetName());
    RooRealVar* pDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarP.GetName());
    RooRealVar* ptDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarPt.GetName());
    RooRealVar* eventNumDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarEvent.GetName());
    RooRealVar* lumiSecDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarLumiSection.GetName());
    RooRealVar* runNumDataTightRPC = (RooRealVar*)argSetTightRPC->find(rooVarRun.GetName());
    for(int evt=0; evt < etaCutNomCutDRegionDataSetSignalTightRPC->numEntries(); ++evt)
    {
      etaCutNomCutDRegionDataSetSignalTightRPC->get(evt);
      if(fabs(etaDataTightRPC->getVal()) > 1.5) continue;
      // find pileupWeight
      EventInfo evtInfo(runNumDataTightRPC->getVal(),lumiSecDataTightRPC->getVal(),eventNumDataTightRPC->getVal());
      pileupWeightTightRPCMapItr = pileupWeightTightRPCMap.find(evtInfo);
      if(pileupWeightTightRPCMapItr == pileupWeightTightRPCMap.end())
      {
        std::cout << "ERROR: Could not find pileup weight for this event! Quitting." << std::endl;
        return -5;
      }
      double pileupWeight = pileupWeightTightRPCMapItr->second;
      double eventWeight = sampleWeightTightRPC*pileupWeight;

      // apply mass cut
      float massSqr = (ihDataTightRPC->getVal()-dEdx_c)*pow(pDataTightRPC->getVal(),2)/dEdx_k;
      if(massSqr < 0)
        continue;
      if(sqrt(massSqr) >= massCut_)
      {
        iasSignalMassCutNoMSliceHist->Fill(iasDataTightRPC->getVal(),eventWeight);
        numSignalTracksInDRegionPassingMassCut+=eventWeight;
        entriesSignalHist->Fill(fabs(etaDataTightRPC->getVal())+0.1,nomIasDataTightRPC->getVal(),eventWeight);
        if(iasDataTightRPC->getVal() > iasCutForEffAcc)
        {
          // if track over ias cut, put the event in the set (one track per event kept)
          // there should only be one track per event at this stage anyway
          pair<set<EventInfo>::iterator,bool> ret;
          ret = selectedSignalEventsTightRPCSet.insert(evtInfo);
          if(ret.second)
          {
            //cout << "INFO: insert event - run: " << runNumDataTightRPC->getVal() << " lumiSec: " << lumiSecDataTightRPC->getVal()
            //  << " eventNum: " << eventNumDataTightRPC->getVal() << endl;
            //cout << "add puweight =" << pileupWeight << " eventWeight=" << eventWeight << " to signalEventsOverIasCutTightRPC; total=" << signalEventsOverIasCutTightRPC << endl;
            signalEventsOverIasCutTightRPC+=eventWeight;
            entriesSignalIasHist->Fill(fabs(etaDataTightRPC->getVal())+0.1,nomIasDataTightRPC->getVal(),eventWeight);
          }
        }
      }
      if(usePtForSideband)
      {
        if(ptDataTightRPC->getVal() < ptSidebandThreshold)
          continue;
      }
      else
      {
        if(pDataTightRPC->getVal() < pSidebandThreshold)
          continue;
      }
      if(doMass)
        massSignalNoMSliceHist->Fill(sqrt(massSqr),eventWeight);
    }

    double numSignalTracksInDRegionMassCutThisSlice = iasSignalMassCutNoMSliceHist->Integral();
    // adjust for bin width ratio if variable bins used
    //for(int bin=1; bin<=iasSignalMassCutNoMSliceHist->GetNbinsX(); ++bin)
    //{
    //  double binWidthRatio =
    //    iasSignalMassCutNoMSliceHist->GetBinWidth(1)/
    //    iasSignalMassCutNoMSliceHist->GetBinWidth(bin);
    //  double binc = iasSignalMassCutNoMSliceHist->GetBinContent(bin);
    //  double bine = iasSignalMassCutNoMSliceHist->GetBinError(bin);
    //  iasSignalMassCutNoMSliceHist->SetBinContent(bin,binWidthRatio*binc);
    //  iasSignalMassCutNoMSliceHist->SetBinError(bin,binWidthRatio*bine);
    //}
    iasSignalMassCutNoMSliceForEffHist = (TH1D*) iasSignalMassCutNoMSliceHist->Clone();

    // figure out overall normalization this slice, background from ABCD
    // if using one track per event, this is the number of events in the slice passing mass in D region
    double bgEntriesInARegionThisSlice =
      aRegionBackgroundEntriesHist->GetBinContent(aRegionBackgroundEntriesHist->FindBin(lowerEta+0.1,lowerNoM+1));
    //double bgEntriesInBRegionThisSlice = 
    //  bRegionBackgroundEntriesHist->GetBinContent(bRegionBackgroundEntriesHist->FindBin(lowerEta+0.1,lowerNoM+1));
    double bgEntriesInCRegionThisSlice = 
      cRegionBackgroundEntriesHist->GetBinContent(cRegionBackgroundEntriesHist->FindBin(lowerEta+0.1,lowerNoM+1));
    //double fractionOfBGTracksPassingMassCutThisSlice = (bgEntriesInBRegionThisSlice > 0) ?
    //  histItr->Integral()/bgEntriesInBRegionThisSlice : 0;
    //double bgTracksInDThisSlice = (bgEntriesInARegionThisSlice > 0) ?
    //  bgEntriesInBRegionThisSlice*bgEntriesInCRegionThisSlice/bgEntriesInARegionThisSlice : 0;
    //double numBackgroundTracksInDRegionPassingMassCutThisSlice =
    //  bgTracksInDThisSlice*fractionOfBGTracksPassingMassCutThisSlice;
    int etaBinCRegionEntriesHist = cRegionBackgroundEntriesHist->GetXaxis()->FindBin(lowerEta+0.1);
    int cRegionNoMSummedEntries = cRegionBackgroundEntriesHist->Integral(
        etaBinCRegionEntriesHist,etaBinCRegionEntriesHist,1,cRegionBackgroundEntriesHist->GetNbinsY());
        
    //errors
    //double fractionOfBGTracksPassingMassCutThisSliceError =
    //  fractionOfBGTracksPassingMassCutThisSlice*
    //  sqrt(pow(sqrt(histItr->Integral())/histItr->Integral(),2)+
    //      pow(sqrt(bgEntriesInBRegionThisSlice)/bgEntriesInBRegionThisSlice,2));
    //double bgTracksInDThisSliceError = bgTracksInDThisSlice*
    //  sqrt(pow(sqrt(bgEntriesInARegionThisSlice)/bgEntriesInARegionThisSlice,2)+
    //      pow(sqrt(bgEntriesInBRegionThisSlice)/bgEntriesInBRegionThisSlice,2)+
    //      pow(sqrt(bgEntriesInCRegionThisSlice)/bgEntriesInCRegionThisSlice,2));
    //double numBackgroundTracksInDRegionPassingMassCutThisSliceError =
    //  (fractionOfBGTracksPassingMassCutThisSlice > 0) ?
    //  numBackgroundTracksInDRegionPassingMassCutThisSlice*
    //  sqrt(pow(bgTracksInDThisSliceError/bgTracksInDThisSlice,2)+
    //      pow(fractionOfBGTracksPassingMassCutThisSliceError/fractionOfBGTracksPassingMassCutThisSlice,2)): 1;
    //
    // get typical number of c region tracks above mass cut
    double averageCTracksAboveMassCut = 0;
    int numberOfBinsInCTracksAboveMassCutSum = 0;
    for(int bin=1; bin < cRegionTracksOverMassCutProfile->GetNbinsX(); ++bin)
    {
      double binc = cRegionTracksOverMassCutProfile->GetBinContent(bin);
      if(binc > 0)
      {
        averageCTracksAboveMassCut+=binc;
        numberOfBinsInCTracksAboveMassCutSum++;
      }
    }
    if(numberOfBinsInCTracksAboveMassCutSum > 0)
      averageCTracksAboveMassCut/=numberOfBinsInCTracksAboveMassCutSum;
    int startBin = iasHist->FindBin(iasCutForEffAcc);
    double totalBkOverIas = 0;
    for(int bin=startBin; bin <= iasHist->GetNbinsX(); ++bin)
      totalBkOverIas+=iasHist->GetBinContent(bin);
    // output
    cout << "Slice: lowerNoM = " << lowerNoM << " lowerEta = " << lowerEta << endl;
    //cout << "bg entries this slice: A = " << bgEntriesInARegionThisSlice << " B = " << 
    //  bgEntriesInBRegionThisSlice << " C = " << bgEntriesInCRegionThisSlice << " ==> D = " <<
    //  bgTracksInDThisSlice << endl;
    //cout << "INFO: SumBkOverIas = " << totalBkOverIas << endl;
    //cout << "INFO: average Cj over mass cut = " << averageCTracksAboveMassCut << endl;
    //cout << "INFO: total C over NoM slices = " << cRegionNoMSummedEntries << endl;
    // error output
    // have to include SUM[Cj > mass cut] error here
    //cout << "INFO: numBackgroundTracksInDRegionPassingMassCutThisSlice: " <<
    //  numBackgroundTracksInDRegionPassingMassCutThisSlice <<
    //  " +/- " << sqrt(relativeSliceErrorSqr)*numBackgroundTracksInDRegionPassingMassCutThisSlice <<
    //  endl;
    //backgroundTracksOverIasCutErrorSqr+=
    //  pow(numBackgroundTracksInDRegionPassingMassCutThisSlice*sqrt(relativeSliceErrorSqr),2);
    //double sliceError = sqrt(totalBkOverIas);
    //if(bgEntriesInARegionThisSlice > 0)
    //  backgroundTracksOverIasCutErrorSqr+=
    //    pow(sliceError*averageCTracksAboveMassCut*bgEntriesInCRegionThisSlice/
    //        (bgEntriesInARegionThisSlice*cRegionNoMSummedEntries),2);
    numBkOverIasCut+=totalBkOverIas;
    // include c tracks passing mass cut error
    double sliceError = 0;
    if(totalBkOverIas > 0)
      sliceError+= 1/totalBkOverIas;
    if(averageCTracksAboveMassCut > 0)
      sliceError+= 1/averageCTracksAboveMassCut;
    if(bgEntriesInCRegionThisSlice > 0)
      sliceError+= 1/bgEntriesInCRegionThisSlice;
    if(bgEntriesInARegionThisSlice > 0)
      sliceError+= 1/bgEntriesInARegionThisSlice;
    if(cRegionNoMSummedEntries > 0)
      sliceError+= 1/cRegionNoMSummedEntries;
    // take sqrt
    sliceError = sqrt(sliceError);
    if(bgEntriesInARegionThisSlice > 0)
      backgroundTracksOverIasCutErrorSqr+=
        pow(sliceError*totalBkOverIas*averageCTracksAboveMassCut*bgEntriesInCRegionThisSlice/
            (bgEntriesInARegionThisSlice*cRegionNoMSummedEntries),2);
    // central value
    if(bgEntriesInARegionThisSlice > 0)
      backgroundTracksOverIasCut+=
        totalBkOverIas*averageCTracksAboveMassCut*bgEntriesInCRegionThisSlice/
        (bgEntriesInARegionThisSlice*cRegionNoMSummedEntries);

    // figure out overall normalization for signal this slice
    // must normalize as if sigma were 1 to measure the cross section
    double signalCrossSection = 1;
    double totalSignalEvents = integratedLumi*signalCrossSection;
    //FIXME track count is sort of wrong since it has to be weighted event-wise
    //double totalGenHSCPTracks = numGenHSCPTracksLooseRPCRooVar->getVal()+numGenHSCPTracksTightRPCRooVar->getVal();
    double totalGenHSCPEvents = eventWeightSumLooseRPCRooVar->getVal()+eventWeightSumTightRPCRooVar->getVal();
    //double fractionOfSigTracksInDRegionPassingMassCutThisSlice =
    //  numSignalTracksInDRegionMassCutThisSlice/totalGenHSCPTracks; // includes trigger eff.
    double fractionOfSigTracksInDRegionPassingMassCutThisSlice = 0;
    if(countEvents)
      fractionOfSigTracksInDRegionPassingMassCutThisSlice =
        numSignalTracksInDRegionMassCutThisSlice/totalGenHSCPEvents;
    double numSignalTracksInDRegionPassingMassCutThisSlice =
      fractionOfSigTracksInDRegionPassingMassCutThisSlice*totalSignalEvents;

    signalTracksOverIasCut+=iasSignalMassCutNoMSliceHist->Integral(
        iasSignalMassCutNoMSliceHist->FindBin(iasCutForEffAcc),
        iasSignalMassCutNoMSliceHist->GetNbinsX());

    // mass
    if(doMass)
    {
      double numSigTracksOverMassThisSlice = massSignalNoMSliceHist->
        Integral(massSignalNoMSliceHist->FindBin(massCut_),massSignalNoMSliceHist->GetNbinsX());
      double fracInThisSlice = numSigTracksOverMassThisSlice/totalGenHSCPEvents;
      double numSigEventsOverMassThisSlice = fracInThisSlice*totalSignalEvents;
      massSignalNoMSliceHist->Sumw2();
      double sigNormFactorMass = (massSignalNoMSliceHist->Integral() > 0) ?
        numSigEventsOverMassThisSlice/massSignalNoMSliceHist->Integral() : 0;
      massSignalNoMSliceHist->Scale(sigNormFactorMass);
    }

    //// figure out nominal sig norm for efficiency value
    //double totalSignalTracksThisSliceEff = integratedLumi*signalCrossSectionForEff;
    //double fractionOfSigTracksInDRegionPassingMassCutThisSliceEff =
    //  numSignalTracksInDRegionMassCutThisSlice/(double)numSignalTracksTotal;
    //double numSignalTracksInDRegionPassingMassCutThisSliceEff =
    //  fractionOfSigTracksInDRegionPassingMassCutThisSliceEff*totalSignalTracksThisSliceEff;
    //double sigNormFactorEff =
    //  numSignalTracksInDRegionPassingMassCutThisSliceEff/
    //  iasSignalMassCutNoMSliceHist->Integral();
    //iasSignalMassCutNoMSliceForEffHist->Scale(sigNormFactorEff);
    // FIX below if wanted
    //signalAcceptanceIasCut+=iasSignalMassCutNoMSliceForEffHist->Integral(
    //    iasSignalMassCutNoMSliceForEffHist->FindBin(iasCutForEffAcc),
    //    iasSignalMassCutNoMSliceForEffHist->GetNbinsX())/iasSignalMassCutNoMSliceForEffHist->Integral();

    //cout << "INFO: fraction of bg tracks passing mass cut this slice = " <<
    //  fractionOfBGTracksPassingMassCutThisSlice << endl << 
    //  " ==> numBackgroundTracksInDRegionPassingMassCutThisSlice = " <<
    //  numBackgroundTracksInDRegionPassingMassCutThisSlice << 
    //  " and integral = " << histItr->Integral() << " so bgNormFactor = ";
    //if(histItr->Integral() > 0)
    //  cout << numBackgroundTracksInDRegionPassingMassCutThisSlice/histItr->Integral();
    //else
    //  cout << "0";
    //cout << endl;
    //cout << "INFO: sum Bk over iasCut = " <<
    //  histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX()) << endl;
    //cout << "INFO: ias cut frac = " <<
    //  histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX())/histItr->Integral() << endl;

    //cout << "first Ias bin this slice: " << iteratorPos*numIasBins << endl;
    double sigTracksOverIasCutThisSlice = 
      iasSignalMassCutNoMSliceHist->Integral(
          iasSignalMassCutNoMSliceHist->FindBin(iasCutForEffAcc),
          iasSignalMassCutNoMSliceHist->GetNbinsX());
    double sigTracksTotalThisSlice = iasSignalMassCutNoMSliceHist->Integral();
    //cout << "INFO: signal tracks total this slice: " <<
    //  sigTracksTotalThisSlice <<
    //  " fraction passing mass cut this slice: " <<
    //  fractionOfSigTracksInDRegionPassingMassCutThisSlice <<
    //  " signal tracks over ias cut this slice: " <<
    //  sigTracksOverIasCutThisSlice <<
    //  " (normed) signal tracks passing mass cut this slice: " << 
    //  numSignalTracksInDRegionPassingMassCutThisSlice <<
    //  " integral: " <<
    //  iasSignalMassCutNoMSliceHist->Integral() <<
    //  endl << endl;
    cout << "Signal tracks over Ias cut this slice: " << sigTracksOverIasCutThisSlice << endl;
    if(doMass)
      cout << "(MASS) Signal tracks over Ias cut this slice: " << massSignalNoMSliceHist->Integral() << endl;

    sigEffOverIasCutVsSliceHist->Fill(lowerEta+0.1,lowerNoM,sigTracksOverIasCutThisSlice/sigTracksTotalThisSlice);
    sigEffOverIasCutHist->Fill(sigTracksOverIasCutThisSlice/sigTracksTotalThisSlice);


    // normalize BG hist
    //double iasShapeIntegral = (histItr->Integral() > 0)  ?
    //  histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX())/histItr->Integral() : 0;
    //double bgTracksInDThisSliceError = (bgEntriesInARegionThisSlice > 0) ?
    //  sqrt(iasShapeIntegral*bgEntriesInBRegionThisSlice)*bgEntriesInCRegionThisSlice/bgEntriesInARegionThisSlice: 0; // ignore A,C errors
    //double numBackgroundTracksInDRegionPassingMassCutThisSliceError =
    //  fractionOfBGTracksPassingMassCutThisSlice*bgTracksInDThisSliceError;
    //backgroundTracksOverIasCutErrorSqr+=pow(numBackgroundTracksInDRegionPassingMassCutThisSliceError,2);
    
    //XXX mar 1 no longer needed
    //histItr->Sumw2();
    //double bgNormFactor = (histItr->Integral() > 0) ?
    //  numBackgroundTracksInDRegionPassingMassCutThisSlice/
    //  histItr->Integral(): 0;
    //histItr->Scale(bgNormFactor);
    //

    //TODO: handle the hists with zero background in systematic
    // (assume an entry in the lowest Ias bin)
  
    //if(lowerNoM==11 && (int)(lowerEta*10)==4)
    //{
    //  cout << "backgroundTracksOverIasCut for NoM 11 eta 0.4: " << 
    //    histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX()) << endl;
    //}
    //XXX SIC FEB 20 switching to new background prediction at least for testing
    cout << "Background tracks over Ias/mass cuts this slice (limits hist): " << 
      histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX()) << endl;
    backgroundTracksOverIasCutNoApprox+=histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX());
    // compute error this slice -- assuming error in each ias bin is uncorrelated (if based only on Bk then it is)
    double errorSqrThisSlice = 0;
    for(int bin=histItr->FindBin(iasCutForEffAcc); bin <= histItr->GetNbinsX(); ++bin)
    {
      errorSqrThisSlice+=pow(histItr->GetBinError(bin),2);
      //std::cout << "bin: " << bin << " error = " << histItr->GetBinError(bin) << " binContent = "
      //  << histItr->GetBinContent(bin) << std::endl;
    }
    backgroundTracksOverIasCutNoApproxErrorSqr+=errorSqrThisSlice;
    backExpOverIasCutVsSliceHist->Fill(
        lowerEta+0.1,lowerNoM,histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX()));
    backExpVsSliceHist->Fill(lowerEta+0.1,lowerNoM,histItr->Integral());
    backExpOverIasCutHist->Fill(histItr->Integral(histItr->FindBin(iasCutForEffAcc),histItr->GetNbinsX()));
    // normalize sig hist
    iasSignalMassCutNoMSliceHist->Sumw2();
    double sigNormFactor = (iasSignalMassCutNoMSliceHist->Integral() > 0) ?
      numSignalTracksInDRegionPassingMassCutThisSlice/iasSignalMassCutNoMSliceHist->Integral() : 0;
    iasSignalMassCutNoMSliceHist->Scale(sigNormFactor);

    // Fill unrolled histograms
    // loop over ias bins, fill unrolled hists
    int firstIasBin = binsToExclude+1;
    for(int iasBin=firstIasBin; iasBin < numIasBins+1; ++iasBin)
    {
      int globalBinIndex = iteratorPos*numBinsEachSlice+iasBin-firstIasBin+1;
      double binc = histItr->GetBinContent(iasBin);
      double bine = histItr->GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,bine);
      //XXX SIC TEST MAY 10 -- SET LARGE ERROR BINS TO ZERO ERROR
      //if(bine/binc >= 1)
      //  backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,0);

      //std::cout << "INFO: Filling background hist: bin=" << globalBinIndex << " content=" << binc << " error=" << bine << std::endl;
      //// plus one sigma
      //backgroundAllNoMAllEtaUnrolledPlusOneSigmaHist->SetBinContent(globalBinIndex,binc+bine);
      //backgroundAllNoMAllEtaUnrolledPlusOneSigmaHist->SetBinError(globalBinIndex,bine);
      //// minus one sigma
      //if(binc-bine < 0)
      //  bine = 0.3*binc; // if 1-sigma shift is negative, shift down to 30% of nominal
      //backgroundAllNoMAllEtaUnrolledMinusOneSigmaHist->SetBinContent(globalBinIndex,binc-bine);
      //backgroundAllNoMAllEtaUnrolledMinusOneSigmaHist->SetBinError(globalBinIndex,bine);

      // discovery
      binc = bgDiscoveryHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgDiscoveryHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinError(globalBinIndex,bine);
      //XXX SIC TEST MAY 10 -- SET LARGE ERROR BINS TO ZERO ERROR
      //if(bine/binc >= 1)
      //  backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinError(globalBinIndex,0);

      binc = iasSignalMassCutNoMSliceHist->GetBinContent(iasBin);
      bine = iasSignalMassCutNoMSliceHist->GetBinError(iasBin);
      signalAllNoMAllEtaUnrolledHist->SetBinContent(globalBinIndex,binc);
      signalAllNoMAllEtaUnrolledHist->SetBinError(globalBinIndex,bine);

      // data D region
      binc = dataHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = dataHistsToUse[iteratorPos].GetBinError(iasBin);
      dataAllNoMAllEtaUnrolledHist->SetBinContent(globalBinIndex,binc);
      dataAllNoMAllEtaUnrolledHist->SetBinError(globalBinIndex,bine);

      //std::cout << "INFO: Filling signal hist: bin=" << globalBinIndex << " content=" << binc << " error=" << bine << std::endl;
    }
    //mass
    if(doMass)
    {
      for(int massBin=1; massBin < numMassBins+1; ++massBin)
      {
        int globalBinIndex = iteratorPos*numMassBins+massBin;
        double binc = massBackgroundPredictionNoMSliceHist->GetBinContent(massBin);
        double bine = massBackgroundPredictionNoMSliceHist->GetBinError(massBin);
        backgroundAllNoMAllEtaUnrolledMassHist->SetBinContent(globalBinIndex,binc);
        backgroundAllNoMAllEtaUnrolledMassHist->SetBinError(globalBinIndex,bine);
        // signal
        binc = massSignalNoMSliceHist->GetBinContent(massBin);
        bine = massSignalNoMSliceHist->GetBinError(massBin);
        signalAllNoMAllEtaUnrolledMassHist->SetBinContent(globalBinIndex,binc);
        signalAllNoMAllEtaUnrolledMassHist->SetBinError(globalBinIndex,bine);
      }
    }

    // write the histograms
    outputRootFile->cd();
    //scaledBGPredDir->cd();
    histItr->Write();
    bgDiscoveryHistsToUse[iteratorPos].Write();
    dataHistsToUse[iteratorPos].Write();
    //normedSignalPredDir->cd();
    iasSignalMassCutNoMSliceHist->Write();
    if(doMass)
    {
      massSignalNoMSliceHist->Write();
      massBackgroundPredictionNoMSliceHist->Write();
    }

    iteratorPos++;
    // cleanup
    //delete nomCutDRegionDataSetSignalLooseRPC;
    //delete nomCutDRegionDataSetSignalTightRPC;
    delete etaCutNomCutDRegionDataSetSignalLooseRPC;
    delete etaCutNomCutDRegionDataSetSignalTightRPC;
    //delete iasBackgroundPredictionMassCutNoMSliceHist;
    delete iasSignalMassCutNoMSliceHist;
    delete iasSignalMassCutNoMSliceForEffHist;
    delete iasHist;
  }

  // now do the standard analysis style prediction
  TH1F* etaBRegionHist = (TH1F*) backgroundPredictionRootFile->Get("etaBRegionHist");
  TH1F* etaARegionHist = (TH1F*) backgroundPredictionRootFile->Get("etaARegionHist");
  TH2F* pEtaCRegionHist = (TH2F*) backgroundPredictionRootFile->Get("pEtaCRegionHist");
  TH2F* pEtaCRegionOrigHist = (TH2F*) pEtaCRegionHist->Clone();
  pEtaCRegionOrigHist->SetName("pEtaCRegionOrig");
  TH1F* ihBRegionTotalHist = (TH1F*) backgroundPredictionRootFile->Get("ihBRegionTotalHist");
  double aRegionInt = etaARegionHist->Integral();
  double aRegionIntError = sqrt(aRegionInt);
  double bRegionInt = ihBRegionTotalHist->Integral();
  double bRegionIntError = sqrt(bRegionInt);
  double cRegionInt = pEtaCRegionHist->Integral();
  double cRegionIntError = sqrt(cRegionInt);
  double dPredTotal = 0;
  double dPredTotalError = 1;
  if(aRegionInt <= 0)
    std::cout << "ERROR: cannot do standard analysis prediction." << endl;
  else
  {
    dPredTotal = bRegionInt*cRegionInt/aRegionInt;
    dPredTotalError = dPredTotal*
      sqrt(pow(aRegionIntError/aRegionInt,2)+pow(bRegionIntError/bRegionInt,2)+pow(cRegionIntError/cRegionInt,2));
  }

  std::cout << "A = " << aRegionInt << " +/- " << aRegionIntError << std::endl <<
               "B = " << bRegionInt << " +/- " << bRegionIntError << std::endl <<
               "C = " << cRegionInt << " +/- " << cRegionIntError << std::endl <<
               "dPredTotal = " << dPredTotal << " +/- " << dPredTotalError << std::endl;

  // normalize
  etaBRegionHist->Scale(1.0/etaBRegionHist->Integral());
  etaARegionHist->Scale(1.0/etaARegionHist->Integral());

  for(int bin=0; bin <= pEtaCRegionHist->GetXaxis()->GetNbins(); ++bin)
  {
    double weight = 0;
    if(etaARegionHist->GetBinContent(bin) > 0)
    {
      weight = etaBRegionHist->GetBinContent(bin)/etaARegionHist->GetBinContent(bin);
      std::cout << "eta = " << etaARegionHist->GetBinCenter(bin) << " weight = " << weight << std::endl;
    }

    for(int y=0; y <= pEtaCRegionHist->GetYaxis()->GetNbins(); ++y)
      pEtaCRegionHist->SetBinContent(bin,y,pEtaCRegionHist->GetBinContent(bin,y)*weight);
  }
  //std::cout << "NOT doing eta reweighting." << std::endl;
  TH1D* reweightedPInCRegionHist = pEtaCRegionHist->ProjectionY("reweightedPInCRegion");
  reweightedPInCRegionHist->Sumw2();
  ihBRegionTotalHist->Sumw2();
  reweightedPInCRegionHist->Scale(1.0/reweightedPInCRegionHist->Integral());
  ihBRegionTotalHist->Scale(1.0/ihBRegionTotalHist->Integral());
  //std::cout << "Entries in reweightedPInCRegionHist = " << reweightedPInCRegionHist->GetEntries() << std::endl;

  // construct mass dist
  TH1F* massHist = new TH1F("massHist","predicted mass",200,0,2000);
  TH1F* massErrorHist = new TH1F("massErrorHist","predicted mass error",200,0,2000);
  for(int x=0; x <= reweightedPInCRegionHist->GetNbinsX(); ++x)
  {
    //std::cout << "DEBUG: reweightedPInCRegionHist bin = " << x << " binContent = " << reweightedPInCRegionHist->GetBinContent(x)
    //  << std::endl;
    if(reweightedPInCRegionHist->GetBinContent(x) <= 0.0)
      continue;
    double thisP = reweightedPInCRegionHist->GetBinCenter(x);
    for(int y=0; y <= ihBRegionTotalHist->GetNbinsX(); ++y)
    {
      //std::cout << "DEBUG: ihBRegionTotalHist bin = " << y << " binContent = " << ihBRegionTotalHist->GetBinContent(y)
      //  << std::endl;
      if(ihBRegionTotalHist->GetBinContent(y) <= 0.0)
        continue;
      double thisIh = ihBRegionTotalHist->GetBinCenter(y);
      double prob = reweightedPInCRegionHist->GetBinContent(x) * ihBRegionTotalHist->GetBinContent(y);
      double probError = prob*sqrt(
          pow(reweightedPInCRegionHist->GetBinError(x)/reweightedPInCRegionHist->GetBinContent(x),2)+
          pow(ihBRegionTotalHist->GetBinError(y)/ihBRegionTotalHist->GetBinContent(y),2));
      //std::cout << "\tDEBUG: prob = " << prob << std::endl;
      if(prob <= 0 || std::isnan(prob))
        continue;
      double massSqr = (thisIh-dEdx_c)*pow(thisP,2)/dEdx_k;
      //std::cout << "\tDEBUG: mass^2 = " << massSqr << std::endl;
      if(massSqr >= 0)
      {
        massHist->Fill(sqrt(massSqr),prob);
        massErrorHist->Fill(sqrt(massSqr),probError);
        //std::cout << "Fill mass = " << sqrt(massSqr) << " prob: " << prob << std::endl;
        //std::cout << "Fill mass = " << sqrt(massSqr) << " probError: " << probError << std::endl;
      }
    }
  }

  for(int x=0; x <= massHist->GetNbinsX(); ++x)
  {
    double scaledMass = massHist->GetBinContent(x) * dPredTotal;
    massHist->SetBinContent(x,scaledMass);
    if(massHist->GetBinContent(x) > 0)
      massErrorHist->SetBinContent(x,
          scaledMass*sqrt(pow(massErrorHist->GetBinContent(x)/massHist->GetBinContent(x),2)+
            pow(dPredTotalError/dPredTotal,2)));
      //massErrorHist->SetBinContent(x,massHist->GetBinContent(x)*dPredTotalError);
    else
      massErrorHist->SetBinContent(x,0);
    //cout << "DEBUG: massError this bin (" << x << ") = " << 
    //    scaledMass*sqrt(pow(massErrorHist->GetBinContent(x)/massHist->GetBinContent(x),2)+
    //      pow(dPredTotalError/dPredTotal,2)) << endl;
    //if(std::isnan(tmpH_Mass    ->GetBinContent(x) * PE_P)){printf("%f x %f\n",tmpH_Mass    ->GetBinContent(x),PE_P); fflush(stdout);exit(0);}
  }

  cout << endl << endl
    << "Number of background tracks in mass hist over ias cut (std ana method): "
    << massHist->Integral()
    << endl
    << "Number of background tracks expected over ias and mass cut (std ana method): "
    << massHist->Integral(massHist->FindBin(massCut_),massHist->GetNbinsX());
  double errorSqrSum = 0;
  for(int i=massHist->FindBin(massCut_); i <=massHist->GetNbinsX(); ++i)
    errorSqrSum+=pow(massErrorHist->GetBinContent(i),2);
  double massPredError = sqrt(errorSqrSum);
  cout << " +/- " << massPredError
       << endl << endl;
  // end of standard analysis style prediction

  double signalEventsOverIasCut = signalEventsOverIasCutLooseRPC+signalEventsOverIasCutTightRPC;
  cout << "Signal events > ias cut loose RPC: " << signalEventsOverIasCutLooseRPC <<
    " tight RPC: " << signalEventsOverIasCutTightRPC << endl;

  cout << endl << endl << "Ias cut = " << iasCutForEffAcc << endl << "\tfound " << backgroundTracksOverIasCutNoApprox
    << " +/- " << sqrt(backgroundTracksOverIasCutNoApproxErrorSqr)
    << " background tracks over ias cut from ias pred hist (limits)" << endl 
    << backgroundTracksOverIasCut
    << " +/- " << sqrt(backgroundTracksOverIasCutErrorSqr)
    << " tracks over ias cut from approximation" << endl
    << numBkOverIasCut << " Bk tracks over ias cut"
    << endl << endl
    << signalEventsOverIasCut << " signal events over ias cut" << endl
    << signalEventsOverIasCut/totalGenHSCPEvents << " signal efficiency (event level) or" << endl  << endl;
    //<< signalTracksOverIasCut/totalGenHSCPTracks << " signal efficiency (track level) with this ias cut. " << endl << endl;

  cout << "event weight sum: " << totalGenHSCPEvents << endl;
  cout << "event weight sum loose: " << eventWeightSumLooseRPCRooVar->getVal() << " tight: " << eventWeightSumTightRPCRooVar->getVal() << endl;
  //cout << "original generated events: " << totalGenHSCPEvents << " events over ias cut: " << signalEventsOverIasCut
  //  << endl;
  //cout << "original generated tracks: " << totalGenHSCPTracks << " tracks over ias cut: " << signalTracksOverIasCut
  //  << endl;
  //cout << "orig signal tracks passing trig/presel: " << numSignalTracksTotal << endl << " orig signal tracks in D region: "
  //  << numSignalTracksInDRegion << endl << " original signal tracks in D region passing mass cut: "
  //  << numSignalTracksInDRegionPassingMassCut
  //  << endl;
  //cout << "(eff) signal tracks passing trig/presel: " << numSignalTracksTotal/totalGenHSCPTracks << endl
  //  << " (eff) signal tracks in D region: " << numSignalTracksInDRegion/totalGenHSCPTracks << endl
  //  << " (eff) signal tracks in D region passing mass cut: " << numSignalTracksInDRegionPassingMassCut/totalGenHSCPTracks
  //  << endl;

//  // make toy data
//  TRandom3* rand = new TRandom3();
//  // background
//  TH1D* backgroundToyDataHist = (TH1D*) backgroundAllNoMAllEtaUnrolledLimitsHist->Clone();
//  backgroundToyDataHist->Reset();
//  backgroundToyDataHist->SetNameTitle("backgroundToyData","background toy data");
//  for(int bin=1; bin <= backgroundToyDataHist->GetNbinsX(); ++bin)
//    backgroundToyDataHist->SetBinContent(bin,rand->Poisson(backgroundAllNoMAllEtaUnrolledHist->GetBinContent(bin)));
//  // signal
//  TH1D* signalToyDataHist = (TH1D*) signalAllNoMAllEtaUnrolledHist->Clone();
//  signalToyDataHist->Reset();
//  signalToyDataHist->SetNameTitle("signalToyData","signal toy data");
//  for(int bin=1; bin <= signalToyDataHist->GetNbinsX(); ++bin)
//  {
//    // set signal zero bins to small value
//    if(signalAllNoMAllEtaUnrolledHist->GetBinContent(bin)==0)
//      signalAllNoMAllEtaUnrolledHist->SetBinContent(bin,1e-15);
//    signalToyDataHist->SetBinContent(bin,rand->Poisson(signalAllNoMAllEtaUnrolledHist->GetBinContent(bin)));
//  }
//  signalToyDataHist->Scale(signalCrossSectionForEff);  // scale by cross section
//  //signalToyDataHist->Scale(1e-4);
//  // signal and background
//  TH1D* signalAndBackgroundToyDataHist = (TH1D*) backgroundAllNoMAllEtaUnrolledHist->Clone();
//  signalAndBackgroundToyDataHist->Reset();
//  signalAndBackgroundToyDataHist->SetNameTitle("signalAndBackgroundToyDataHist","signal and background toy data");
//  signalAndBackgroundToyDataHist->Add(backgroundToyDataHist,signalToyDataHist);

  //check for zeros
  for(int bin=1; bin <= backgroundAllNoMAllEtaUnrolledLimitsHist->GetNbinsX(); ++bin)
  {
    double binc = backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinContent(bin);
    //double bincToy = signalAndBackgroundToyDataHist->GetBinContent(bin);
    //if(bincToy > 0)
    //{
    //  std::cout << "bincToy=" << bincToy << " binc=" << binc << std::endl;
    //  if(binc==0)
    //    std::cout << "\tDetected bin where binc==0!" << std::endl;
    //}
    if(binc == 0)
    {
      std::cout << "WARNING: limits hist: binc==0 : bin " << bin << std::endl;
      outputRootFile->cd();
      backgroundAllNoMAllEtaUnrolledLimitsHist->Write();
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->Write();
      outputRootFile->Close();
      return -11;
    }

    binc = backgroundAllNoMAllEtaUnrolledDiscoveryHist->GetBinContent(bin);
    if(binc == 0)
    {
      std::cout << "WARNING: discovery hist: binc==0 : bin " << bin << std::endl;
      outputRootFile->cd();
      backgroundAllNoMAllEtaUnrolledLimitsHist->Write();
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->Write();
      outputRootFile->Close();
      return -11;
    }
  }

  outputRootFile->cd();
  etaBRegionHist->Write();
  etaARegionHist->Write();
  pEtaCRegionHist->Write();
  pEtaCRegionOrigHist->Write();
  ihBRegionTotalHist->Write();
  reweightedPInCRegionHist->Write();
  massHist->Write();
  backgroundAllNoMAllEtaUnrolledLimitsHist->Write();
  backgroundAllNoMAllEtaUnrolledDiscoveryHist->Write();
  //backgroundAllNoMAllEtaUnrolledPlusOneSigmaHist->Write();
  //backgroundAllNoMAllEtaUnrolledMinusOneSigmaHist->Write();
  signalAllNoMAllEtaUnrolledHist->Write();
  dataAllNoMAllEtaUnrolledHist->Write();
  if(doMass)
  {
    backgroundAllNoMAllEtaUnrolledMassHist->Write();
    signalAllNoMAllEtaUnrolledMassHist->Write();
  }
  sigEffOverIasCutVsSliceHist->Write();
  sigEffOverIasCutHist->Write();
  backExpOverIasCutVsSliceHist->Write();
  backExpVsSliceHist->Write();
  backExpOverIasCutHist->Write();
  entriesSignalHist->Write();
  entriesSignalIasHist->Write();
  //// toy data
  //backgroundToyDataHist->Write();
  //signalToyDataHist->Write();
  //signalAndBackgroundToyDataHist->Write();
  outputRootFile->Close();

  //backgroundPredictionRootFile->Close();
  //signalRootFile->Close();

  if(backgroundAllNoMAllEtaUnrolledLimitsHist->Integral() <= 0
      || backgroundAllNoMAllEtaUnrolledDiscoveryHist->Integral() <= 0)
    return -1;
}

