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
#include <map>
#include <algorithm>

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
// NB: these must match what's used in makeIasPredictions
//TODO: make these configurable
int nomStep = 2;
float etaStep = 0.2;
int lastLowerNoM = 17;
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
// python config globals
double massCut_;
string backgroundPredictionRootFilename_;
string signalLooseRPCRootFilename_;
string outputRootFilename_;
double dEdx_k;
double dEdx_c;
double pSidebandThreshold;
double ptSidebandThreshold;
bool usePtForSideband;
double ihSidebandThreshold;
double integratedLumi;
double signalCrossSectionForEff;
double iasCutForEffAcc;
bool useIasForSideband;
bool doMass;
bool is8TeV = true;

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
  return lowerEta/etaStep;
}
//
int getNoMSliceFromLowerNoM(int lowerNoM)
{
  return (lowerNoM-5)/nomStep;
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
//
double getNoMWeightFromLowerNoM(int nom)
{
  //TODO update for 8 TeV
  // from study on pcminn26
  if(nom<=5)
    return 1.375364;
  else if(nom==6)
    return 1.756386;
  else if(nom==7)
    return 1.448682;
  else if(nom==8)
    return 1.350035;
  else if(nom==9)
    return 1.219157;
  else if(nom==10)
    return 1.127699;
  else if(nom==11)
    return 1.050045;
  else if(nom==12)
    return 1.011671;
  else if(nom==13)
    return 0.9925607;
  else if(nom==14)
    return 0.9688514;
  else if(nom==15)
    return 0.9585392;
  else if(nom==16)
    return 0.9565235;
  else if(nom==17)
    return 0.9436069;
  else if(nom==18)
    return 0.924295;
  else if(nom==19)
    return 0.8817048;
  else if(nom==20)
    return 0.8175668;
  else if(nom>=21)
    return 0.7557076;

    return 1; // nothing else to do
}

double getApproxMuValueFromMassCut(int massCut)
{
  if(massCut==180) return 0.003736;
  else if(massCut==230) return 0.002882;
  else if(massCut==290) return 0.002622;
  else if(massCut==360) return 0.002556;
  else if(massCut==410) return 0.002755;
  else if(massCut==460) return 0.002857;
  else if(massCut==500) return 0.002887;
  else if(massCut==540) return 0.003108;
  else if(massCut==570) return 0.003762;
  else if(massCut==600) return 0.004337;
  else if(massCut==70) return 0.005881;
  else if(massCut==130) return 0.003623;
  else if(massCut==190) return 0.003101;
  else if(massCut==250) return 0.002431;
  else if(massCut==310) return 0.002176;
  else if(massCut==360) return 0.001978;
  else if(massCut==410) return 0.00204;
  else if(massCut==450) return 0.001882;
  else if(massCut==20) return 0.00208;
  else if(massCut==50) return 0.002488;
  else if(massCut==80) return 0.001908;
  else if(massCut==120) return 0.00154;
  else if(massCut==150) return 0.001143;
  else if(massCut==190) return 0.001031;
  else if(massCut==230) return 0.00101;
  else if(massCut==260) return 0.000936;
  else if(massCut==300) return 0.000853;
  else return 0.001;
}



int doScaledPredictions(TFile* outputRootFile, TFile* backgroundPredictionRootFile, bool doIasShift)
{
  using namespace std;
  using namespace RooFit;
  TFile* signalLooseRPCRootFile = TFile::Open(signalLooseRPCRootFilename_.c_str());

  string signalDatasetName = doIasShift ? "rooDataSetIasShift" : "rooDataSetOneCandidatePerEvent";
  bool countEvents  = true;

  // get roodataset from signal file -- loose RPC (period 1)
  RooDataSet* rooDataSetAllSignalLooseRPC = (RooDataSet*)signalLooseRPCRootFile->Get(signalDatasetName.c_str());
  if(!rooDataSetAllSignalLooseRPC)
  {
    cout << "Problem with RooDataSet named rooDataSetCandidates in signal file " <<
      signalLooseRPCRootFilename_.c_str() << endl;
    return -3;
  }
  cout << "Testing 8 TeV " << endl;
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
  // SIC MAR 5 -- don't do P cut for signal D region -- will be taken care of by mass cut later
  if(useIasForSideband)
    regionDDataSetSignalLooseRPC = (RooDataSet*)rooDataSetAllSignalLooseRPC->reduce(Cut(iasSearchCutString.c_str()));
  else
    regionDDataSetSignalLooseRPC = (RooDataSet*)rooDataSetAllSignalLooseRPC->reduce(Cut(ihSearchCutString.c_str()));
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
  vector<TH1D> bgPlusOneSigmaLimitsHistsToUse;
  vector<TH1D> bgMinusOneSigmaLimitsHistsToUse;
  vector<TH1D> bgPlusOneSigmaDiscoveryHistsToUse;
  vector<TH1D> bgMinusOneSigmaDiscoveryHistsToUse;
  vector<TH1D> dataHistsToUse;
  //vector<TH1F> bRegionHistsToUse;
  //vector<TH1F> cRegionHistsToUse;

  // for histograms already made
  string dirName="iasPredictionsFixedBins";
  string bgLimitsHistNameEnd="iasPredictionFixedLimitsHist";
  string bgDiscoveryHistNameEnd="iasPredictionFixedDiscoveryHist";
  string bgPlusOneSigmaLimitsHistNameEnd="iasPredictionFixedSlopeAddOneSigmaLimitsHist";
  string bgMinusOneSigmaLimitsHistNameEnd="iasPredictionFixedSlopeMinusOneSigmaLimitsHist";
  string bgPlusOneSigmaDiscoveryHistNameEnd="iasPredictionFixedSlopeAddOneSigmaDiscoveryHist";
  string bgMinusOneSigmaDiscoveryHistNameEnd="iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist";
  string dataDirName="dRegionFixedBins";
  string dataHistNameEnd="dRegionFixedHist";

  // look at BG predictions to determine unrolled hist binning
  for(int lowerNoM = minNoM; lowerNoM < maxNoM; lowerNoM+=nomStep)
  {
    for(float lowerEta = minEta; lowerEta < maxEta; lowerEta+=etaStep)
    {
      TH1D* iasBackgroundPredictionMassCutNoMSliceLimitsHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSliceDiscoveryHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaLimitsHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaLimitsHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaDiscoveryHist = 0;
      TH1D* iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaDiscoveryHist = 0;
      //TH1F* bRegionHist = 0;
      //TH1F* cRegionHist = 0;
      // get this eta/nom hist
      string getLimitsHistName = getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep);
      string getDiscoveryHistName = getLimitsHistName;
      string getPlusOneSigmaLimitsHistName = getLimitsHistName;
      string getMinusOneSigmaLimitsHistName = getLimitsHistName;
      string getPlusOneSigmaDiscoveryHistName = getDiscoveryHistName;
      string getMinusOneSigmaDiscoveryHistName = getDiscoveryHistName;
      getLimitsHistName+=bgLimitsHistNameEnd;
      getDiscoveryHistName+=bgDiscoveryHistNameEnd;
      getPlusOneSigmaLimitsHistName+=bgPlusOneSigmaLimitsHistNameEnd;
      getMinusOneSigmaLimitsHistName+=bgMinusOneSigmaLimitsHistNameEnd;
      getPlusOneSigmaDiscoveryHistName+=bgPlusOneSigmaDiscoveryHistNameEnd;
      getMinusOneSigmaDiscoveryHistName+=bgMinusOneSigmaDiscoveryHistNameEnd;
      string fullPathLimits = dirName;
      fullPathLimits+="/";
      string fullPathDiscovery = fullPathLimits;
      string fullPathPlusOneSigmaLimits = fullPathLimits;
      string fullPathMinusOneSigmaLimits = fullPathLimits;
      string fullPathPlusOneSigmaDiscovery = fullPathLimits;
      string fullPathMinusOneSigmaDiscovery = fullPathLimits;
      fullPathLimits+=getLimitsHistName;
      fullPathDiscovery+=getDiscoveryHistName;
      fullPathPlusOneSigmaLimits+=getPlusOneSigmaLimitsHistName;
      fullPathMinusOneSigmaLimits+=getMinusOneSigmaLimitsHistName;
      fullPathPlusOneSigmaDiscovery+=getPlusOneSigmaDiscoveryHistName;
      fullPathMinusOneSigmaDiscovery+=getMinusOneSigmaDiscoveryHistName;

      TH1D* dRegionDataHist = 0;
      string dataHistGetName = getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep);
      dataHistGetName+=dataHistNameEnd;
      string fullDataHistPath = dataDirName;
      fullDataHistPath+="/";
      fullDataHistPath+=dataHistGetName;

      iasBackgroundPredictionMassCutNoMSliceLimitsHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathLimits.c_str());
      iasBackgroundPredictionMassCutNoMSliceDiscoveryHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathDiscovery.c_str());
      iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaLimitsHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathPlusOneSigmaLimits.c_str());
      iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaLimitsHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathMinusOneSigmaLimits.c_str());
      iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaDiscoveryHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathPlusOneSigmaDiscovery.c_str());
      iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaDiscoveryHist = (TH1D*)backgroundPredictionRootFile->Get(fullPathMinusOneSigmaDiscovery.c_str());
        
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
        //cout << "Found hist: " << iasBackgroundPredictionMassCutNoMSliceLimitsHist->GetName()
        //  << " with integral = " << iasBackgroundPredictionMassCutNoMSliceHist->Integral()
        //  << ", pushing onto the vector." << endl;
        bgLimitsHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceLimitsHist);
        bgDiscoveryHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceDiscoveryHist);
        //cout << "Pushing hist: " << fullPathPlusOneSigma << " onto the vector..." << endl;
        bgPlusOneSigmaLimitsHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaLimitsHist);
        bgMinusOneSigmaLimitsHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaLimitsHist);
        bgPlusOneSigmaDiscoveryHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaDiscoveryHist);
        bgMinusOneSigmaDiscoveryHistsToUse.push_back(*iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaDiscoveryHist);
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
      else
      {
        cout << "WARNING: no input found for this slice : lim or disc hist missing: " <<
          getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep) << " in file: " << backgroundPredictionRootFilename_ << ")." << 
          " Excluded -- should not happen. " << endl;
        return -11;
      }
      delete iasBackgroundPredictionMassCutNoMSliceLimitsHist;
      delete iasBackgroundPredictionMassCutNoMSliceDiscoveryHist;
      delete iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaLimitsHist;
      delete iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaLimitsHist;
      delete iasBackgroundPredictionMassCutNoMSlicePlusOneSigmaDiscoveryHist;
      delete iasBackgroundPredictionMassCutNoMSliceMinusOneSigmaDiscoveryHist;
      //delete bRegionHist;
      //delete cRegionHist;
    }
    if(lowerNoM==5) // step over 5 and the next slice as well
      lowerNoM+=nomStep;
  }

  assert(bgLimitsHistsToUse.size() == bgDiscoveryHistsToUse.size());

  string histNameSuffix = doIasShift ? "IasShiftHist" : "Hist";
  string histTitleSuffix = doIasShift ? " ias shifted" : "";

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
  // background and data hists get overwritten a few times, but I don't care
  TH1D* backgroundAllNoMAllEtaUnrolledLimitsHist = new TH1D("backgroundAllNoMAllEtaUnrolledLimitsHist",
      "unrolled background hist (limits)",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledDiscoveryHist = new TH1D("backgroundAllNoMAllEtaUnrolledDiscoveryHist",
      "unrolled background hist (discovery)",numGlobalBins,1,numGlobalBins+1);
  TH1D* dataAllNoMAllEtaUnrolledHist = new TH1D("dataAllNoMAllEtaUnrolledHist",
      "unrolled data hist, D region",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist = new TH1D("backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist",
      "unrolled background hist lims (shape +1 sigma)",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist = new TH1D("backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist",
      "unrolled background hist lims (shape -1 sigma)",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist = new TH1D("backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist",
      "unrolled background hist disc (shape +1 sigma)",numGlobalBins,1,numGlobalBins+1);
  TH1D* backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist = new TH1D("backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist",
      "unrolled background hist disc (shape -1 sigma)",numGlobalBins,1,numGlobalBins+1);
  TH1D* signalAllNoMAllEtaUnrolledHist = new TH1D(("signalAllNoMAllEtaUnrolled"+histNameSuffix).c_str(),
      ("unrolled signal hist"+histTitleSuffix).c_str(),numGlobalBins,1,numGlobalBins+1);
  TH1D* signalAllNoMAllEtaUnrolledNomShiftedHist = new TH1D(("signalAllNoMAllEtaUnrolledNoMShifted"+histNameSuffix).c_str(),
      ("unrolled signal hist, NoM weighted"+histTitleSuffix).c_str(),numGlobalBins,1,numGlobalBins+1);
  // mass
  TH1D* backgroundAllNoMAllEtaUnrolledMassHist = new TH1D("backgroundAllNoMAllEtaUnrolledMassHist",
      "unrolled background hist (mass)",numGlobalMassBins,1,numGlobalMassBins+1);
  TH1D* signalAllNoMAllEtaUnrolledMassHist = new TH1D(("signalAllNoMAllEtaUnrolledMass"+histNameSuffix).c_str(),
      ("unrolled signal hist (mass)"+histTitleSuffix).c_str(),numGlobalMassBins,1,numGlobalMassBins+1);
  string sigEffOverIasCutVsSliceHistTitle = "Sig eff. over ias cut";
  sigEffOverIasCutVsSliceHistTitle+=histTitleSuffix;
  TH2F* sigEffOverIasCutVsSliceHist = new TH2F(("sigEffOverIasCutVsSlice"+histNameSuffix).c_str(),(sigEffOverIasCutVsSliceHistTitle+";#eta;nom").c_str(),
      24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  int divisions = 500+int(13/nomStep)+1;
  sigEffOverIasCutVsSliceHist->GetYaxis()->SetNdivisions(divisions,false);
  TH1F* sigEffOverIasCutHist = new TH1F(("sigEffOverIasCut"+histNameSuffix).c_str(),("Signal eff. over ias cut"+histTitleSuffix).c_str(),100,0,1);

  TH2F* backExpOverIasCutVsSliceHist = new TH2F("backExpOverIasCutVsSlice","Exp. background over ias cut;#eta;nom",24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  backExpOverIasCutVsSliceHist->GetYaxis()->SetNdivisions(divisions,false);
  TH2F* backExpVsSliceHist = new TH2F("backExpVsSlice","Exp. background;#eta;nom",24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  backExpVsSliceHist->GetYaxis()->SetNdivisions(divisions,false);
  TH1F* backExpOverIasCutHist = new TH1F("backExpOverIasCut","Exp. background over ias cut",100,0,1);

  string entriesSignalHistTitle = "Signal tracks over mass cut in slice";
  entriesSignalHistTitle+=histTitleSuffix;
  TH2F* entriesSignalHist = new TH2F(("entriesSignal"+histNameSuffix).c_str(),(entriesSignalHistTitle+";#eta;nom").c_str(),24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  entriesSignalHist->GetYaxis()->SetNdivisions(divisions,false);
  string entriesSignalIasHistTitle = "Signal tracks over mass/ias cuts in slice";
  entriesSignalIasHistTitle+=histTitleSuffix;
  TH2F* entriesSignalIasHist = new TH2F(("entriesSignalIas"+histNameSuffix).c_str(),(entriesSignalIasHistTitle+";#eta;nom").c_str(),24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  entriesSignalIasHist->GetYaxis()->SetNdivisions(divisions,false);
  // sorting  [(delta B)/S] / [(B/S + mu)]^2
  map<double,int> deltaBOverSAndBinNumbersMap;

  double backgroundTracksOverIasCut = 0;
  double backgroundTracksOverIasCutErrorSqr = 0;
  double backgroundTracksOverIasCutNoApprox = 0;
  double backgroundTracksOverIasCutNoApproxErrorSqr = 0;
  int numBkOverIasCut = 0;
  double signalTracksOverIasCut = 0;
  double signalTracksLooseRPCTotal = 0;
  double signalEventsOverIasCutLooseRPC = 0;
  double numSignalTracksInDRegionPassingMassCut = 0;
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
  // SIC testing
  double likelihoodMuHatZero = 1;
  double likelihoodMuHatPos = 1;
  double likelihoodMuHatNeg = 1;

  //testing for 2012 only 1 period but I had set them as both
  double totalGenHSCPEvents = eventWeightSumLooseRPCRooVar->getVal();
  cout << "event weight sum: " << totalGenHSCPEvents << endl;
  cout << "event weight sum loose: " << eventWeightSumLooseRPCRooVar->getVal() << endl;
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
  // loop over limits hists to use
  int iteratorPos = 0;
  for(vector<TH1D>::iterator histItr = bgLimitsHistsToUse.begin();
      histItr != bgLimitsHistsToUse.end(); ++histItr)
  {
    // delete this and the matching delete at the bottom las
    //TH1F* iasBackgroundPredictionMassCutNoMSliceHist = 0;
    TH1D* iasSignalMassCutNoMSliceHist = 0;
    TH1D* iasSignalMassCutNoMSliceForEffHist = 0;
    TH1D* iasSignalMassCutNoMSliceNoMShiftedHist = 0;
    string bgHistName = string(histItr->GetName());
    string bgHistNameBeg = bgHistName.substr(0,bgHistName.size()-bgLimitsHistNameEnd.size());
    int lowerNoM = getLowerNoMFromHistName(bgHistName,bgLimitsHistNameEnd);
    float lowerEta = getLowerEtaFromHistName(bgHistName,bgLimitsHistNameEnd);
    string sigHistName = bgHistNameBeg;
    //sigHistName+="iasSignalVarBinHist";
    sigHistName+="iasSignalFixed";
    sigHistName+=histNameSuffix;
    string iasSignalHistTitle = "Ias of signal";
    iasSignalHistTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
    iasSignalHistTitle+=histTitleSuffix;
    // get b region ias histogram
    string iasHistName = "bRegionHistograms/";
    iasHistName+=getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep);
    iasHistName+="iasBRegionHist";
    TH1D* iasHist = (TH1D*)backgroundPredictionRootFile->Get(iasHistName.c_str());
    // get c region tracks passing mass cut in ias bins
    string cRegionTracksOverMassCutProfileName="cRegionHistograms/";
    cRegionTracksOverMassCutProfileName+=getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep);
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
    // nom shifted hist
    string sigHistNoMShiftName = bgHistNameBeg;
    sigHistNoMShiftName+="iasSignalFixed";
    sigHistNoMShiftName+=histNameSuffix;
    sigHistNoMShiftName+="NoMShifted";
    string iasSignalHistNoMShiftTitle = "Ias of signal";
    iasSignalHistNoMShiftTitle+=getHistTitleEnd(lowerNoM,lowerEta,massCut_);
    iasSignalHistNoMShiftTitle+=histTitleSuffix;
    iasSignalHistNoMShiftTitle+=" NoM weighted";
    iasSignalMassCutNoMSliceNoMShiftedHist = (TH1D*) histItr->Clone();
    iasSignalMassCutNoMSliceNoMShiftedHist->Reset();
    iasSignalMassCutNoMSliceNoMShiftedHist->SetNameTitle(sigHistNoMShiftName.c_str(),iasSignalHistNoMShiftTitle.c_str());
    // mass
    TH1D* massBackgroundPredictionNoMSliceHist = 0;
    TH1D* massSignalNoMSliceHist = 0;
    if(doMass)
    {
      // bg mass hist
      string dirName="massPredictionsFixedBins";
      string bgHistNameEnd="massPredictionFixedHist";
      string getHistName = getHistNameBeg(lowerNoM,lowerEta,nomStep,etaStep);
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
    std::string nomCutString = "rooVarNoMias>=";
    nomCutString+=intToString(lowerNoM);
    nomCutString+="&&rooVarNoMias<=";
    nomCutString+=intToString(lowerNoM+nomStep-1);
    if(lowerNoM==lastLowerNoM) // do nom lastLowerNoM+ in one slice
    {
      nomCutString = "rooVarNoMias>=";
      nomCutString+=intToString(lastLowerNoM);
    }
    if(lowerNoM==5) // combine 5-6 with next slice
    {
      nomCutString = "rooVarNoMias>=";
      nomCutString+=intToString(lowerNoM);
      nomCutString+="&&rooVarNoMias<=";
      nomCutString+=intToString(lowerNoM+2*nomStep-1);
    }
    RooDataSet* nomCutDRegionDataSetSignalLooseRPC =
      (RooDataSet*)regionDDataSetSignalLooseRPC->reduce(Cut(nomCutString.c_str()));
    std::string etaCutString = "(rooVarEta>";
    etaCutString+=floatToString(lowerEta);
    etaCutString+="&&rooVarEta<";
    std::string upperEtaLimit;
    //XXX eta 1.5 max
    if(lowerEta+etaStep>1.5)
      upperEtaLimit=floatToString(1.5);
    else
      upperEtaLimit=floatToString(lowerEta+etaStep);
    etaCutString+=upperEtaLimit;
    etaCutString+=")||(rooVarEta>-";
    etaCutString+=upperEtaLimit;
    etaCutString+="&&rooVarEta<-";
    etaCutString+=floatToString(lowerEta);
    etaCutString+=")";
    RooDataSet* etaCutNomCutDRegionDataSetSignalLooseRPC =
      (RooDataSet*)nomCutDRegionDataSetSignalLooseRPC->reduce(Cut(etaCutString.c_str()));
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
      // hard eta cut
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
        iasSignalMassCutNoMSliceNoMShiftedHist->Fill(iasDataLooseRPC->getVal(),eventWeight*getNoMWeightFromLowerNoM(nomIasDataLooseRPC->getVal()));
        numSignalTracksInDRegionPassingMassCut+=eventWeight;
        entriesSignalHist->Fill(fabs(etaDataLooseRPC->getVal())+0.01,nomIasDataLooseRPC->getVal(),eventWeight);
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
            entriesSignalIasHist->Fill(fabs(etaDataLooseRPC->getVal())+0.01,nomIasDataLooseRPC->getVal(),eventWeight);
          }
        }
      }
      //TODO: figure out why is this here - is pt made redundant by mass cut? JT
      //seems to be because no p cut is needed as the p cut is effectively done by the mass cut
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

    double numSignalTracksInDRegionMassCutThisSlice = iasSignalMassCutNoMSliceHist->Integral();
    double numSignalTracksInDRegionMassCutThisSliceNoMWeighted = iasSignalMassCutNoMSliceNoMShiftedHist->Integral();
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
      aRegionBackgroundEntriesHist->GetBinContent(aRegionBackgroundEntriesHist->FindBin(lowerEta+0.01,lowerNoM+1));
    if(lowerNoM==5)
      bgEntriesInARegionThisSlice+=aRegionBackgroundEntriesHist->GetBinContent(aRegionBackgroundEntriesHist->FindBin(lowerEta+0.01,lowerNoM+2));
    //double bgEntriesInBRegionThisSlice = 
    //  bRegionBackgroundEntriesHist->GetBinContent(bRegionBackgroundEntriesHist->FindBin(lowerEta+0.1,lowerNoM+1));
    double bgEntriesInCRegionThisSlice = 
      cRegionBackgroundEntriesHist->GetBinContent(cRegionBackgroundEntriesHist->FindBin(lowerEta+0.01,lowerNoM+1));
    if(lowerNoM==5)
      bgEntriesInCRegionThisSlice+=cRegionBackgroundEntriesHist->GetBinContent(cRegionBackgroundEntriesHist->FindBin(lowerEta+0.01,lowerNoM+2));
    //double fractionOfBGTracksPassingMassCutThisSlice = (bgEntriesInBRegionThisSlice > 0) ?
    //  histItr->Integral()/bgEntriesInBRegionThisSlice : 0;
    //double bgTracksInDThisSlice = (bgEntriesInARegionThisSlice > 0) ?
    //  bgEntriesInBRegionThisSlice*bgEntriesInCRegionThisSlice/bgEntriesInARegionThisSlice : 0;
    //double numBackgroundTracksInDRegionPassingMassCutThisSlice =
    //  bgTracksInDThisSlice*fractionOfBGTracksPassingMassCutThisSlice;
    int etaBinCRegionEntriesHist = cRegionBackgroundEntriesHist->GetXaxis()->FindBin(lowerEta+0.01);
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
    //testing in 2012 quick fix for loose / tight
    double totalGenHSCPEvents = eventWeightSumLooseRPCRooVar->getVal();
    //double fractionOfSigTracksInDRegionPassingMassCutThisSlice =
    //  numSignalTracksInDRegionMassCutThisSlice/totalGenHSCPTracks; // includes trigger eff.
    double fractionOfSigTracksInDRegionPassingMassCutThisSlice = 0;
    double fractionOfSigTracksInDRegionPassingMassCutThisSliceNoMWeighted = 0;
    if(countEvents)
    {
      fractionOfSigTracksInDRegionPassingMassCutThisSlice =
        numSignalTracksInDRegionMassCutThisSlice/totalGenHSCPEvents;
      fractionOfSigTracksInDRegionPassingMassCutThisSliceNoMWeighted =
        numSignalTracksInDRegionMassCutThisSliceNoMWeighted/totalGenHSCPEvents;
    }
    double numSignalTracksInDRegionPassingMassCutThisSlice =
      fractionOfSigTracksInDRegionPassingMassCutThisSlice*totalSignalEvents;
    // don't use - renormalize NoM shift using the unshifted normalization
    //double numSignalTracksInDRegionPassingMassCutThisSliceNoMWeighted = 
    //  fractionOfSigTracksInDRegionPassingMassCutThisSliceNoMWeighted*totalSignalEvents;

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
    // normalize sig NoM weighted/shifted hist
    iasSignalMassCutNoMSliceNoMShiftedHist->Sumw2();
    double sigNormFactorNoMWeighted = (iasSignalMassCutNoMSliceNoMShiftedHist->Integral() > 0) ?
      numSignalTracksInDRegionPassingMassCutThisSlice/iasSignalMassCutNoMSliceNoMShiftedHist->Integral() : 0;
      //numSignalTracksInDRegionPassingMassCutThisSliceNoMWeighted/iasSignalMassCutNoMSliceNoMShiftedHist->Integral() : 0;
      // SIC: normalize to nominal number of signal tracks to keep total integral the same with or without NoM shift
    iasSignalMassCutNoMSliceNoMShiftedHist->Scale(sigNormFactorNoMWeighted);

    int firstIasBin = binsToExclude+1;

    // Fill unrolled histograms
    // loop over ias bins, fill unrolled hists
    for(int iasBin=firstIasBin; iasBin < numIasBins+1; ++iasBin)
    {
      int globalBinIndex = iteratorPos*numBinsEachSlice+iasBin-firstIasBin+1;
      double binc = histItr->GetBinContent(iasBin);
      double bine = histItr->GetBinError(iasBin);
      // limits
      backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,bine);
      // SIC TESTING: DECREASE BACKGROUND PREDICTION BY 2%
      //backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinContent(globalBinIndex,0.98*binc);
      //backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,0.98*bine);
      //// SIC TESTING: DECREASE BACKGROUND PREDICTION BY 25%
      //backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinContent(globalBinIndex,binc*0.75);
      //backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,bine*0.75);
      // reset hist bins
      //histItr->SetBinContent(iasBin,binc*0.75);
      //histItr->SetBinError(iasBin,bine*0.75);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(globalBinIndex,0);
        histItr->SetBinError(iasBin,0);
      }
      double backExpThisBin = binc;
      double backExpErrorThisBin = histItr->GetBinError(iasBin);

      // +1 sigma slope shift
      binc = bgPlusOneSigmaLimitsHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgPlusOneSigmaLimitsHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist->SetBinError(globalBinIndex,bine);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist->SetBinError(globalBinIndex,0);
        bgPlusOneSigmaLimitsHistsToUse[iteratorPos].SetBinError(iasBin,0);
      }
      // -1 sigma slope shift
      binc = bgMinusOneSigmaLimitsHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgMinusOneSigmaLimitsHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist->SetBinError(globalBinIndex,bine);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist->SetBinError(globalBinIndex,0);
        bgMinusOneSigmaLimitsHistsToUse[iteratorPos].SetBinError(iasBin,0);
      }

      // discovery
      binc = bgDiscoveryHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgDiscoveryHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinError(globalBinIndex,bine);
      //// SIC TESTING: DECREASE BACKGROUND PREDICTION BY 2%
      //backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinContent(globalBinIndex,binc*0.98);
      //backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinError(globalBinIndex,bine*0.98);
      //bgDiscoveryHistsToUse[iteratorPos].SetBinContent(iasBin,binc*0.98);
      //bgDiscoveryHistsToUse[iteratorPos].SetBinError(iasBin,bine*0.98);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledDiscoveryHist->SetBinError(globalBinIndex,0);
        bgDiscoveryHistsToUse[iteratorPos].SetBinError(iasBin,0);
      }

      // +1 sigma slope shift
      binc = bgPlusOneSigmaDiscoveryHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgPlusOneSigmaDiscoveryHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist->SetBinError(globalBinIndex,bine);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist->SetBinError(globalBinIndex,0);
        bgPlusOneSigmaDiscoveryHistsToUse[iteratorPos].SetBinError(iasBin,0);
      }
      // -1 sigma slope shift
      binc = bgMinusOneSigmaDiscoveryHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = bgMinusOneSigmaDiscoveryHistsToUse[iteratorPos].GetBinError(iasBin);
      backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist->SetBinContent(globalBinIndex,binc);
      backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist->SetBinError(globalBinIndex,bine);
      // SIC TEST -- SET LARGE ERROR BINS TO ZERO ERROR
      if(bine/binc >= 1 && binc < 1)
      {
        backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist->SetBinError(globalBinIndex,0);
        bgMinusOneSigmaDiscoveryHistsToUse[iteratorPos].SetBinError(iasBin,0);
      }

      // signal
      binc = iasSignalMassCutNoMSliceHist->GetBinContent(iasBin);
      bine = iasSignalMassCutNoMSliceHist->GetBinError(iasBin);
      signalAllNoMAllEtaUnrolledHist->SetBinContent(globalBinIndex,binc);
      signalAllNoMAllEtaUnrolledHist->SetBinError(globalBinIndex,bine);
      double sigExpThisBin = binc;
      // NoM shifted version
      binc = iasSignalMassCutNoMSliceNoMShiftedHist->GetBinContent(iasBin);
      bine = iasSignalMassCutNoMSliceNoMShiftedHist->GetBinError(iasBin);
      signalAllNoMAllEtaUnrolledNomShiftedHist->SetBinContent(globalBinIndex,binc);
      signalAllNoMAllEtaUnrolledNomShiftedHist->SetBinError(globalBinIndex,bine);

      // push onto map
      double approxMuValue = getApproxMuValueFromMassCut(massCut_);
      //if(bine/binc < 1 || binc >= 1)
      if(backExpErrorThisBin > 0)
        deltaBOverSAndBinNumbersMap[(backExpErrorThisBin/sigExpThisBin)/pow((backExpThisBin/sigExpThisBin+approxMuValue),2)] = globalBinIndex;
      // data D region
      binc = dataHistsToUse[iteratorPos].GetBinContent(iasBin);
      bine = dataHistsToUse[iteratorPos].GetBinError(iasBin);
      dataAllNoMAllEtaUnrolledHist->SetBinContent(globalBinIndex,binc);
      dataAllNoMAllEtaUnrolledHist->SetBinError(globalBinIndex,bine);
      //// SIC TESTING -- putting one fake event in at high Ias
      //if(globalBinIndex==606)
      //{
      //  int eventsToAdd = 1;
      //  std::cout << "Data in bin 606=" << binc << std::endl;
      //  std::cout << "BG pred in bin 606=" << backgroundAllNoMAllEtaUnrolledDiscoveryHist->GetBinContent(606) << std::endl;
      //  std::cout << "sig pred in bin 606=" << signalAllNoMAllEtaUnrolledHist->GetBinContent(606) << std::endl;
      //  std::cout << "Put " << eventsToAdd << " fake event(s) into the data" << std::endl;
      //  dataAllNoMAllEtaUnrolledHist->SetBinContent(globalBinIndex,eventsToAdd);
      //  dataAllNoMAllEtaUnrolledHist->SetBinError(globalBinIndex,eventsToAdd);
      //  dataHistsToUse[iteratorPos].SetBinContent(iasBin,eventsToAdd);
      //  dataHistsToUse[iteratorPos].SetBinError(iasBin,eventsToAdd);
      //}
      double dataThisBin = dataHistsToUse[iteratorPos].GetBinContent(iasBin);

      likelihoodMuHatZero*=TMath::Poisson(dataThisBin,backExpThisBin);
      likelihoodMuHatPos*=TMath::Poisson(dataThisBin,backExpThisBin+0.1*sigExpThisBin);
      likelihoodMuHatNeg*=TMath::Poisson(dataThisBin,backExpThisBin-0.1*sigExpThisBin);
      //std::cout << "Likelihood 'by hand' this bin: data: "
      //  << dataThisBin << " backExp: " << backExpThisBin << " 0.1*sigExpThisBin " << 0.1*sigExpThisBin << " mu-hat = 0: " << TMath::Poisson(dataThisBin,backExpThisBin)
      //  << " mu-hat = -0.1: " << TMath::Poisson(dataThisBin,backExpThisBin-0.1*sigExpThisBin) << " mu-hat = 0.1: " << TMath::Poisson(dataThisBin,backExpThisBin+0.1*sigExpThisBin) << std::endl << std::endl;

      //std::cout << "INFO: Filling signal hist: bin=" << globalBinIndex << " content=" << binc << " error=" << bine << std::endl;
    }
    //std::cout << "Likelihood 'by hand' so far: mu-hat 0: " << likelihoodMuHatZero << " mu-hat -0.1: " << likelihoodMuHatNeg << " mu-hat 0.1: " << likelihoodMuHatPos << std::endl << std::endl;
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
    bgPlusOneSigmaLimitsHistsToUse[iteratorPos].Write();
    bgMinusOneSigmaLimitsHistsToUse[iteratorPos].Write();
    bgPlusOneSigmaDiscoveryHistsToUse[iteratorPos].Write();
    bgMinusOneSigmaDiscoveryHistsToUse[iteratorPos].Write();
    //normedSignalPredDir->cd();
    iasSignalMassCutNoMSliceHist->Write();
    iasSignalMassCutNoMSliceNoMShiftedHist->Write();
    if(doMass)
    {
      massSignalNoMSliceHist->Write();
      massBackgroundPredictionNoMSliceHist->Write();
    }

    iteratorPos++;
    // cleanup
    delete nomCutDRegionDataSetSignalLooseRPC;
    delete etaCutNomCutDRegionDataSetSignalLooseRPC;
    //delete iasBackgroundPredictionMassCutNoMSliceHist;
    delete iasSignalMassCutNoMSliceHist;
    delete iasSignalMassCutNoMSliceNoMShiftedHist;
    delete iasSignalMassCutNoMSliceForEffHist;
    delete iasHist;
  }

  double signalEventsOverIasCut = signalEventsOverIasCutLooseRPC;
  cout << "Signal events > ias cut loose RPC: " << signalEventsOverIasCutLooseRPC << endl;
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
  cout << "event weight sum loose: " << eventWeightSumLooseRPCRooVar->getVal() << endl;
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

  std::cout << "Size of map: " << deltaBOverSAndBinNumbersMap.size() << std::endl;
  // loop map of s/b values and bin numbers
  unsigned int maxBinsWithErrors = 50;
  if(deltaBOverSAndBinNumbersMap.size() > maxBinsWithErrors)
  {
    unsigned int count = 0;
    map<double,int>::reverse_iterator itr;
    for(itr = deltaBOverSAndBinNumbersMap.rbegin();
        itr != deltaBOverSAndBinNumbersMap.rend(); itr++)
    {

      //if(count >= maxBinsWithErrors) // keep error on first xxx bins with highest deltaB/S value
      //{
      //  if(backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinError(itr->second) > 0)
      //  {
      //    std::cout << "INFO: count = " << count << "; Setting error on background bin " << itr->second << " to zero." << std::endl;
      //    backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(itr->second,0);
      //  }
      //}
      //else
      //{
      //  std::cout << "INFO: count = " << count << "; Keeping error of " << backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinError(itr->second) <<
      //    " (deltaB/S(mu)) = " << itr->first << " on bin " << itr->second << std::endl;
      //}

      if(count < maxBinsWithErrors) // set first xxx largest to zero error
      {
        if(backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinError(itr->second) > 0)
        {
          std::cout << "INFO: count = " << count << "; Setting error on background bin " << itr->second << " to zero." << std::endl;
          backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(itr->second,0);
        }
      }
      else if(count >= 2*maxBinsWithErrors) // AND set all bins with deltaB/S numbers xxx and higher to zero error
      {
        if(backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinError(itr->second) > 0)
        {
          std::cout << "INFO: count = " << count << "; Setting error on background bin " << itr->second << " to zero." << std::endl;
          backgroundAllNoMAllEtaUnrolledLimitsHist->SetBinError(itr->second,0);
        }
      }
      else
      {
        std::cout << "INFO: count = " << count << "; Keeping error of " << backgroundAllNoMAllEtaUnrolledLimitsHist->GetBinError(itr->second) <<
          " (deltaB/S(mu)) = " << itr->first << " on bin " << itr->second << std::endl;
      }

      ++count;
    }
  }

  outputRootFile->cd();
  backgroundAllNoMAllEtaUnrolledLimitsHist->Write();
  backgroundAllNoMAllEtaUnrolledDiscoveryHist->Write();
  backgroundAllNoMAllEtaUnrolledPlusOneSigmaLimitsHist->Write();
  backgroundAllNoMAllEtaUnrolledMinusOneSigmaLimitsHist->Write();
  backgroundAllNoMAllEtaUnrolledPlusOneSigmaDiscoveryHist->Write();
  backgroundAllNoMAllEtaUnrolledMinusOneSigmaDiscoveryHist->Write();
  signalAllNoMAllEtaUnrolledHist->Write();
  signalAllNoMAllEtaUnrolledNomShiftedHist->Write();
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

  if(backgroundAllNoMAllEtaUnrolledLimitsHist->Integral() <= 0
      || backgroundAllNoMAllEtaUnrolledDiscoveryHist->Integral() <= 0)
    return -1;

  return 0;
}

void doStandardAnalysisPrediction(TFile* backgroundPredictionRootFile, TFile* outputRootFile)
{
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

  //std::cout << "Likelihood 'by hand': mu-hat = 0: " << likelihoodMuHatZero << " mu-hat = -0.1: " << likelihoodMuHatNeg << " mu-hat = 0.1: " << likelihoodMuHatPos << std::endl << std::endl;

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
      if(prob <= 0 || isnan((float)prob))
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
    //if(isnan((float)tmpH_Mass    ->GetBinContent(x) * PE_P)){printf("%f x %f\n",tmpH_Mass    ->GetBinContent(x),PE_P); fflush(stdout);exit(0);}
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
  outputRootFile->cd();
  etaBRegionHist->Write();
  etaARegionHist->Write();
  pEtaCRegionHist->Write();
  pEtaCRegionOrigHist->Write();
  ihBRegionTotalHist->Write();
  reweightedPInCRegionHist->Write();
  massHist->Write();
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
  massCut_ = ana.getParameter<double>("MassCut");
  backgroundPredictionRootFilename_ = ana.getParameter<string>("BackgroundPredictionInputRootFile");
  signalLooseRPCRootFilename_ = ana.getParameter<string>("SignalLooseRPCInputRootFile");
  outputRootFilename_ = ana.getParameter<string>("OutputRootFile");
  // dE/dx calibration
  dEdx_k = ana.getParameter<double>("dEdx_k");
  dEdx_c = ana.getParameter<double>("dEdx_c");
  // definition of sidebands/search region
  pSidebandThreshold = ana.getParameter<double>("PSidebandThreshold");
  ptSidebandThreshold = ana.getParameter<double>("PtSidebandThreshold");
  usePtForSideband = ana.getParameter<bool>("UsePtForSideband");
  ihSidebandThreshold = ana.getParameter<double>("IhSidebandThreshold");
  integratedLumi = ana.getParameter<double>("IntegratedLuminosity"); // 1/pb
  signalCrossSectionForEff = ana.getParameter<double>("SignalCrossSection"); // pb
  iasCutForEffAcc = ana.getParameter<double>("IasCutForEfficiency");
  useIasForSideband = ana.getParameter<bool>("UseIasForSideband");

  // do mass?
  doMass = false;

  std::cout << "Using dE/dx k: " << dEdx_k << " c: " << dEdx_c << std::endl;

  // TODO configurable nom/eta limits
  // NB: always use upper edge of last bin for both
  minNoM = 5;
  maxNoM = 18;
  minEta = 0.0;
  maxEta = 1.6;

  //string fileNameEnd = generateFileNameEnd(massCut_,pSidebandThreshold,ptSidebandThreshold,usePtForSideband,ihSidebandThreshold);
  string fileNameEnd = "";
  TFile* outputRootFile = new TFile((outputRootFilename_+fileNameEnd).c_str(),"recreate");
  TFile* backgroundPredictionRootFile = TFile::Open(backgroundPredictionRootFilename_.c_str());

  doScaledPredictions(outputRootFile,backgroundPredictionRootFile,false);
  doScaledPredictions(outputRootFile,backgroundPredictionRootFile,true);

  outputRootFile->cd();
  // get ahold of some needed hists
  TH1D* iasShiftedSignalHist = (TH1D*) outputRootFile->Get("signalAllNoMAllEtaUnrolledIasShiftHist");
  TH1D* nominalSignalHist = (TH1D*) outputRootFile->Get("signalAllNoMAllEtaUnrolledHist");
  TH1D* nomShiftedSignalHist = (TH1D*) outputRootFile->Get("signalAllNoMAllEtaUnrolledNoMShiftedHist");
  // new nominal
  TH1D* newNominalSignalHist = (TH1D*) nominalSignalHist->Clone();
  newNominalSignalHist->SetNameTitle("signalAllNoMAllEtaUnrolledHalfIasShiftHalfNoMShiftHist","unrolled signal hist, 1/2 Ias shift, 1/2 NoM shift (new nominal)");
  // make 1/2 ias shifted
  TH1D* halfIasShiftSignalHist = (TH1D*) nominalSignalHist->Clone();
  halfIasShiftSignalHist->Add(iasShiftedSignalHist,-1);
  halfIasShiftSignalHist->Scale(-0.5); // to take half difference
  // add half-diff to nominal for new nominal
  newNominalSignalHist->Add(halfIasShiftSignalHist,1);
  halfIasShiftSignalHist->SetNameTitle("signalAllNoMAllEtaUnrolledHalfIasShiftHist","unrolled signal hist, 1/2 Ias shift");
  // make 1/2 nom shifted
  TH1D* halfNoMShiftSignalHist = (TH1D*) nominalSignalHist->Clone();
  halfNoMShiftSignalHist->Add(nomShiftedSignalHist,-1);
  halfNoMShiftSignalHist->Scale(-0.5); // to take half difference
  // add half-diff to nominal for new nominal
  newNominalSignalHist->Add(halfNoMShiftSignalHist,1);
  halfNoMShiftSignalHist->SetNameTitle("signalAllNoMAllEtaUnrolledHalfNoMShiftHist","unrolled signal hist, 1/2 NoM shift");

  // make Ias variational hist
  TH1D* signalIasPlusOneSigmaHist = (TH1D*) iasShiftedSignalHist->Clone();
  signalIasPlusOneSigmaHist->Add(halfNoMShiftSignalHist,1);
  signalIasPlusOneSigmaHist->SetNameTitle("signalAllNoMAllEtaUnrolledPlusOneSigmaIasHist","unrolled signal hist, 1/2 nom shift +1#sigma Ias shift (full ias shift)");
  // make NoM variational hist
  TH1D* signalNoMPlusOneSigmaHist = (TH1D*) iasShiftedSignalHist->Clone();
  signalNoMPlusOneSigmaHist->Add(halfIasShiftSignalHist,1);
  signalNoMPlusOneSigmaHist->SetNameTitle("signalAllNoMAllEtaUnrolledPlusOneSigmaNoMHist","unrolled signal hist, 1/2 ias shift +1#sigma NoM shift (full nom shift)");

  // finish half-shifted signal hists
  // add on originals
  halfNoMShiftSignalHist->Add(nominalSignalHist,1);
  halfIasShiftSignalHist->Add(nominalSignalHist,1);
  // finish variation hists
  TH1D* signalIasMinusOneSigmaHist = (TH1D*) halfNoMShiftSignalHist->Clone();
  signalIasMinusOneSigmaHist->SetNameTitle("signalAllNoMAllEtaUnrolledMinusOneSigmaIasHist","unrolled signal hist, 1/2 nom shift -1#sigma Ias shift (old nominal)");
  TH1D* signalNoMMinusOneSigmaHist = (TH1D*) halfIasShiftSignalHist->Clone();
  signalNoMMinusOneSigmaHist->SetNameTitle("signalAllNoMAllEtaUnrolledMinusOneSigmaNoMHist","unrolled signal hist, 1/2 ias shift -1#sigma NoM shift (old nominal)");

  // check for zeros/negatives
  for(int bin=1; bin <= signalIasPlusOneSigmaHist->GetNbinsX(); ++bin)
  {
    double binc = signalIasPlusOneSigmaHist->GetBinContent(bin);
    if(binc <= 0)
    {
      signalIasPlusOneSigmaHist->SetBinContent(bin,1e-25);
      signalIasPlusOneSigmaHist->SetBinError(bin,0);
    }
    binc = signalIasMinusOneSigmaHist->GetBinContent(bin);
    if(binc <= 0)
    {
      signalIasMinusOneSigmaHist->SetBinContent(bin,1e-25);
      signalIasMinusOneSigmaHist->SetBinError(bin,0);
    }
    binc = signalNoMPlusOneSigmaHist->GetBinContent(bin);
    if(binc <= 0)
    {
      signalNoMPlusOneSigmaHist->SetBinContent(bin,1e-25);
      signalNoMPlusOneSigmaHist->SetBinError(bin,0);
    }
    binc = signalNoMMinusOneSigmaHist->GetBinContent(bin);
    if(binc <= 0)
    {
      signalNoMMinusOneSigmaHist->SetBinContent(bin,1e-25);
      signalNoMMinusOneSigmaHist->SetBinError(bin,0);
    }
    binc = halfIasShiftSignalHist->GetBinContent(bin);
    if(binc <= 0)
    {
      halfIasShiftSignalHist->SetBinContent(bin,1e-25);
      halfIasShiftSignalHist->SetBinError(bin,0);
    }
    binc = halfNoMShiftSignalHist->GetBinContent(bin);
    if(binc <= 0)
    {
      halfNoMShiftSignalHist->SetBinContent(bin,1e-25);
      halfNoMShiftSignalHist->SetBinError(bin,0);
    }
    binc = newNominalSignalHist->GetBinContent(bin);
    if(binc <= 0)
    {
      newNominalSignalHist->SetBinContent(bin,1e-25);
      newNominalSignalHist->SetBinError(bin,0);
    }
  }
  // write
  signalIasMinusOneSigmaHist->Write();
  signalIasPlusOneSigmaHist->Write();
  signalNoMMinusOneSigmaHist->Write();
  signalNoMPlusOneSigmaHist->Write();
  halfIasShiftSignalHist->Write();
  halfNoMShiftSignalHist->Write();
  newNominalSignalHist->Write();

  //FIXME seg fault for some reason
  //doStandardAnalysisPrediction(outputRootFile,backgroundPredictionRootFile);
  outputRootFile->Close();
  //backgroundPredictionRootFile->Close();
  //signalRootFile->Close();
}


