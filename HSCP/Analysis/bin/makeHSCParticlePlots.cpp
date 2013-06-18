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

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooPlot.h"

#include "TSystem.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <iomanip>
#include <vector>
#include <string>
#include <limits>

namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle; class HSCPIsolation;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag; class LumiReWeighting;}
namespace reweight{class PoissonMeanShifter;}

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/FWLite/interface/InputSource.h"
#include "DataFormats/FWLite/interface/OutputFiles.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/RefCore.h"
#include "DataFormats/Common/interface/RefCoreWithIndex.h"

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
#include "DataFormats/Math/interface/deltaR.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "commonFunctions.h"

// MC-data Ih peak shift functions
TF1* mcDataIhPeakShiftFuncs[8];
// ih width -- averaged over NoM, study on laptop
TF1* ihWidthSmearFunc; // quad diff
//Testing for 8 TeV using similar method for Ih and Ias shift
// MC-data Ias peak shift functions
TF1* mcDataIasPeakShiftFunc;
// ias width -- averaged over NoM, study on laptop
TF1* iasWidthSmearFunc; // 
// MC-data Ias peak shift graphs
TGraph* mcDataIasPeakShiftGraphs[8];
// MC-data Ias width shift graphs
TGraph* mcDataIasWidthShiftGraphs[8]; // quad diff

// adapted from AnalysisStep234.C (GetGenHSCPBeta)
int getNumGenHSCP(const std::vector<reco::GenParticle>& genColl, bool onlyCharged)
{
  int numGenHscp = 0;
  for(unsigned int g=0;g<genColl.size();g++)
  {
    if(genColl[g].pt()<5)continue;
    if(genColl[g].status()!=1)continue;
    int AbsPdg=abs(genColl[g].pdgId());
    if(AbsPdg<1000000)continue;
    if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))continue; //Skip neutral gluino RHadrons
    if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))continue;  //skip neutral stop RHadrons
    //if(beta1<0){beta1=genColl[g].p()/genColl[g].energy();}else if(beta2<0){beta2=genColl[g].p()/genColl[g].energy();return;}
    numGenHscp++;
  }
  return numGenHscp;
}

double GetSampleWeight(const double& IntegratedLuminosityInPb, const double& IntegratedLuminosityInPbBeforeTriggerChange,
    const double& CrossSection, const double& MCEvents, int period, bool is8TeV)
{
  double Weight = 1.0;
  if(IntegratedLuminosityInPb>=IntegratedLuminosityInPbBeforeTriggerChange && IntegratedLuminosityInPb>0)
  {
    double NMCEvents = MCEvents;
    //if(MaxEntry>0)NMCEvents=std::min(MCEvents,(double)MaxEntry);
    if(period==0)
      Weight = (CrossSection * IntegratedLuminosityInPbBeforeTriggerChange) / NMCEvents;
    else if(period==1 && !is8TeV)
      Weight = (CrossSection * (IntegratedLuminosityInPb-IntegratedLuminosityInPbBeforeTriggerChange)) / NMCEvents;
    else if(period==1 && is8TeV)
      Weight = (CrossSection * (IntegratedLuminosityInPb)) / NMCEvents;
  }
  return Weight;
}

double GetPUWeight(const fwlite::Event& ev, const std::string& pileup, double &PUSystFactor, edm::LumiReWeighting& LumiWeightsMC, edm::LumiReWeighting& LumiWeightsMCSyst){
   fwlite::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   PupInfo.getByLabel(ev, "addPileupInfo");
   if(!PupInfo.isValid()){printf("PileupSummaryInfo Collection NotFound\n");return 1.0;}
   double PUWeight_thisevent=1;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int npv = -1; float Tnpv = -1;

   if(pileup=="S4"){
      float sum_nvtx = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         npv = PVI->getPU_NumInteractions();
         sum_nvtx += float(npv);
      }
      float ave_nvtx = sum_nvtx/3.;
      PUWeight_thisevent = LumiWeightsMC.weight( ave_nvtx );
      if(PUWeight_thisevent==0) PUSystFactor=1;
      else PUSystFactor = LumiWeightsMCSyst.weight( ave_nvtx ) / PUWeight_thisevent;
   }else if(pileup=="S3"){
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      PUWeight_thisevent = LumiWeightsMC.weight( npv );
      if(PUWeight_thisevent==0) PUSystFactor=1;
      else PUSystFactor = LumiWeightsMCSyst.weight( npv ) / PUWeight_thisevent;
   }else if(pileup=="S10"){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
         Tnpv = PVI->getTrueNumInteractions();
         continue;
       }
     }
     PUWeight_thisevent = LumiWeightsMC.weight( Tnpv );
     if(PUWeight_thisevent==0) PUSystFactor=1;
     else PUSystFactor = LumiWeightsMCSyst.weight( Tnpv ) / PUWeight_thisevent;
   }
   else {
     printf("Can not find pile up scenario");
   }
   return PUWeight_thisevent;
}

void initializeIhPeakDiffFunctions(TF1* funcArray[], int arraySize, bool is8TeV)
{
  //functions should take momentum and return difference in Ih between MC and data (data - MC) JT
  for(int i=0; i<arraySize; ++i)
  {
    std::string funcName = "myFitNoM";
    funcName+=intToString(i);
    funcArray[i] = new TF1(funcName.c_str(),"pol1(0)");
  }
  if(!is8TeV)
  {
  // from study on laptop taking MC/data peak difference SC
  // sometimes excluding lowest momentum points
    funcArray[0]->SetParameters(0.6301,-0.3019); // NoM 5-6
    funcArray[1]->SetParameters(0.5195,-0.2264); // 7-8
    funcArray[2]->SetParameters(0.3924,-0.1238); // 9-10 <-- fitRange 0.68-1.42
    funcArray[3]->SetParameters(0.4338,-0.1765); // 11-12 <-- fitRange 0.76-1.41
    funcArray[4]->SetParameters(0.3947,-0.1671); // 13-14 <-- fitRange 0.84-1.4
    funcArray[5]->SetParameters(0.4378,-0.2382); // 15-16
    funcArray[6]->SetParameters(0.567,-0.3441);  // 17-18
    funcArray[7]->SetParameters(0.6271,-0.3795); // 18+
  }
  else
  {
    //testing, for now just have all NoM the same
    funcArray[0]->SetParameters(-0.08639,0.1674); // NoM 5-6
    funcArray[1]->SetParameters(-0.08639,0.1674); // 7-8
    funcArray[2]->SetParameters(-0.08639,0.1674); // 9-10 <-- fitRange 0.68-1.42
    funcArray[3]->SetParameters(-0.08639,0.1674); // 11-12 <-- fitRange 0.76-1.41
    funcArray[4]->SetParameters(-0.08639,0.1674); // 13-14 <-- fitRange 0.84-1.4
    funcArray[5]->SetParameters(-0.08639,0.1674); // 15-16
    funcArray[6]->SetParameters(-0.08639,0.1674);  // 17-18
    funcArray[7]->SetParameters(-0.08639,0.1674); // 18+
  }
}

void initializeIhWidthDiffFunction(bool is8TeV)
{
  if(!is8TeV)
  {
    ihWidthSmearFunc = new TF1("ihWidthSmearFunc","pol2(0)");
    ihWidthSmearFunc->SetParameters(1.035,-1.495,0.6101);
  }
  else
  {
    //testing
    ihWidthSmearFunc = new TF1("ihWidthSmearFunc","pol2(0)");
    ihWidthSmearFunc->SetParameters(-0.02842,0.1487,-0.09215);
  }
}
//testing
void initializeIasPeakDiffFunction()
{
  //functions should take momentum and return difference in Ih between MC and data (data - MC) JT
  mcDataIasPeakShiftFunc = new TF1("mcDataIasPeakShiftFunc","pol1(0)");
  mcDataIasPeakShiftFunc->SetParameters(0.153287,-0.178339);
}
//testing
void initializeIasWidthDiffFunction()
{
  iasWidthSmearFunc = new TF1("iasWidthSmearFunc","pol1(0)");
  iasWidthSmearFunc->SetParameters(-0.04935,0.07484);
}

void initializeIasShiftGraphs(std::string rootFileName)
{
  TFile inputGraphFile(rootFileName.c_str(),"r");
  // peak
  mcDataIasPeakShiftGraphs[0] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM0");
  mcDataIasPeakShiftGraphs[1] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM1");
  mcDataIasPeakShiftGraphs[2] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM2");
  mcDataIasPeakShiftGraphs[3] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM3");
  mcDataIasPeakShiftGraphs[4] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM4");
  mcDataIasPeakShiftGraphs[5] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM5");
  mcDataIasPeakShiftGraphs[6] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM6");
  mcDataIasPeakShiftGraphs[7] = (TGraph*) inputGraphFile.Get("mcDataPeakDiffGraphNoM7");
  // width (quad diff)
  mcDataIasWidthShiftGraphs[0] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM0");
  mcDataIasWidthShiftGraphs[1] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM1");
  mcDataIasWidthShiftGraphs[2] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM2");
  mcDataIasWidthShiftGraphs[3] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM3");
  mcDataIasWidthShiftGraphs[4] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM4");
  mcDataIasWidthShiftGraphs[5] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM5");
  mcDataIasWidthShiftGraphs[6] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM6");
  mcDataIasWidthShiftGraphs[7] = (TGraph*) inputGraphFile.Get("mcDataWidthQuadDiffGraphNoM7");
  inputGraphFile.Close();
}

int getSliceFromNoM(int nom)
{
  if(nom==5||nom==6)
    return 0;
  else if(nom==7||nom==8)
    return 1;
  else if(nom==9||nom==10)
    return 2;
  else if(nom==11||nom==12)
    return 3;
  else if(nom==13||nom==14)
    return 4;
  else if(nom==15||nom==16)
    return 5;
  else if(nom==17||nom==18)
    return 6;
  else if(nom>18)
    return 7;
  else
    return -1;
}

double getEquivProtonP(double betaHSCP)
{
  // p = gamma*beta*m
  double gamma = 1/sqrt(1-pow(betaHSCP,2));
  return gamma*betaHSCP*0.938;
}

double getIhPeakCorrection(int nom, double hscpBeta)
{
  int nomSlice = getSliceFromNoM(nom);
  double protonP = getEquivProtonP(hscpBeta);
  // special handling of apparent saturation effects
  // in NoM around eta approx. 1
  if(nomSlice==4 && protonP < 0.84)
  {
    if(protonP >= 0.8)
      return 0.1982;
    else if(protonP >= 0.76)
      return 0.2527;
    else if(protonP >= 0.72)
      return 0.2330;
    else
      return 0.1409;
  }
  else if(nomSlice==3 && protonP < 0.76)
  {
    if(protonP >= 0.72)
      return 0.2253;
    else if(protonP >= 0.68)
      return 0.1852;
    else
      return 0.1929;
  }
  else if(nomSlice==2 && protonP < 0.68)
  {
    if(protonP >= 0.64)
      return 0.2662;
    else if(protonP >= 0.60)
      return 0.2309;
    else
      return 0.3097;
  }

  return mcDataIhPeakShiftFuncs[nomSlice]->Eval(protonP);
}

double getIhWidthCorrection(double hscpBeta)
{
  double protonP = getEquivProtonP(hscpBeta);
  return ihWidthSmearFunc->Eval(protonP);
}

double getIasPeakCorrection(int nom, double hscpBeta, bool is8TeV)
{
  int nomSlice = getSliceFromNoM(nom);
  double protonP = getEquivProtonP(hscpBeta);
  if(!is8TeV)
  {
    return mcDataIasPeakShiftGraphs[nomSlice]->Eval(protonP);
  }
  else
  {
    return mcDataIasPeakShiftFunc->Eval(protonP);
  }
}

double getIasWidthCorrection(int nom, double hscpBeta, bool is8TeV)
{
  int nomSlice = getSliceFromNoM(nom);
  double protonP = getEquivProtonP(hscpBeta);
  if(!is8TeV)
  {
    return mcDataIasWidthShiftGraphs[nomSlice]->Eval(protonP);
  }
  else
  {
    return iasWidthSmearFunc->Eval(protonP);
  }
}


// ****** main
int main(int argc, char ** argv)
{
  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  // parse arguments
  if(argc < 2)
  {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }

  using namespace RooFit;

  // get the python configuration
  const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  fwlite::InputSource inputHandler_(process);
  fwlite::OutputFiles outputHandler_(process);
  int maxEvents_(inputHandler_.maxEvents());
  // now get each parameterset
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("makeHSCParticlePlots");
  double scaleFactor_(ana.getParameter<double>("ScaleFactor"));
  // mass cut to use for the high-p high-Ih (search region) ias dist
  double massCutIasHighPHighIh_ (ana.getParameter<double>("MassCut"));
  // dE/dx calibration
  double dEdx_k_ (ana.getParameter<double>("dEdx_k"));
  double dEdx_c_ (ana.getParameter<double>("dEdx_c"));
  // definition of sidebands/search region
  double pSidebandThreshold_ (ana.getParameter<double>("PSidebandThreshold"));
  double ihSidebandThreshold_ (ana.getParameter<double>("IhSidebandThreshold"));
  bool matchToHSCP_ (ana.getParameter<bool>("MatchToHSCP"));
  bool is8TeV_ (ana.getParameter<bool>("Is8TeV"));
  bool isMC_ (ana.getParameter<bool>("IsMC"));
  double signalEventCrossSection_ (ana.getParameter<double>("SignalEventCrossSection")); // pb
  double integratedLumi_ (ana.getParameter<double>("IntegratedLumi"));
  double integratedLumiBeforeTriggerChange_ (ana.getParameter<double>("IntegratedLumiBeforeTriggerChange"));

  // fileService initialization
  fwlite::TFileService fs = fwlite::TFileService(outputHandler_.file().c_str());
  TH1F* pDistributionHist = fs.make<TH1F>("pDistribution","P;GeV",200,0,2000);
  pDistributionHist->Sumw2();
  TH1F* ptDistributionHist = fs.make<TH1F>("ptDistribution","Pt;GeV",200,0,2000);
  ptDistributionHist->Sumw2();
  // p vs Ias
  TH2F* pVsIasDistributionHist;
  std::string pVsIasName = "trackPvsIas";
  std::string pVsIasTitle="Track P vs ias";
  pVsIasTitle+=";;GeV";
  pVsIasDistributionHist = fs.make<TH2F>(pVsIasName.c_str(),pVsIasTitle.c_str(),400,0,1,100,0,1000);
  pVsIasDistributionHist->Sumw2();
  // p vs Ih
  TH2F* pVsIhDistributionHist;
  std::string pVsIhName = "trackPvsIh";
  std::string pVsIhTitle="Track P vs ih";
  pVsIhTitle+=";MeV/cm;GeV";
  pVsIhDistributionHist = fs.make<TH2F>(pVsIhName.c_str(),pVsIhTitle.c_str(),400,0,10,100,0,1000);
  pVsIhDistributionHist->Sumw2();
  // p vs Ias
  TH2F* pVsIasHist = fs.make<TH2F>("pVsIas","P vs Ias (SearchRegion+MassCut);;GeV",20,0,1,40,0,1000);
  pVsIasHist->Sumw2();
  TH2F* trackEtaVsPHist = fs.make<TH2F>("trackEtaVsP","Track #eta vs. p (SearchRegion+MassCut);GeV",4000,0,2000,24,0,2.4);
  trackEtaVsPHist->Sumw2();
  // search region (only filled for MC)
//  TH1F* pDistributionSearchRegionHist = fs.make<TH1F>("pDistributionSearchRegion","SearchRegion+MassCutP;GeV",200,0,2000);
//  pDistributionSearchRegionHist->Sumw2();
//  TH1F* ptDistributionSearchRegionHist = fs.make<TH1F>("ptDistributionSearchRegion","SearchRegion+MassCutPt;GeV",200,0,2000);
//  ptDistributionSearchRegionHist->Sumw2();
//  // p vs Ias
//  TH2F* pVsIasDistributionSearchRegionHist;
//  std::string pVsIasSearchRegionName = "trackPvsIasSearchRegion";
//  std::string pVsIasSearchRegionTitle="Track P vs ias (SearchRegion+MassCut)";
//  pVsIasSearchRegionTitle+=";;GeV";
//  pVsIasDistributionSearchRegionHist = fs.make<TH2F>(pVsIasSearchRegionName.c_str(),pVsIasSearchRegionTitle.c_str(),400,0,1,100,0,1000);
//  pVsIasDistributionSearchRegionHist->Sumw2();
  // p vs Ih - SR
  TH2F* pVsIhDistributionSearchRegionHist;
  std::string pVsIhSearchRegionName = "trackPvsIhSearchRegion";
  std::string pVsIhSearchRegionTitle="Track P vs ih (SearchRegion+MassCut)";
  pVsIhSearchRegionTitle+=";MeV/cm;GeV";
  pVsIhDistributionSearchRegionHist = fs.make<TH2F>(pVsIhSearchRegionName.c_str(),pVsIhSearchRegionTitle.c_str(),400,0,10,100,0,1000);
  pVsIhDistributionSearchRegionHist->Sumw2();
  // p vs Ih - B region
  TH2F* pVsIhDistributionBRegionHist;
  std::string pVsIhBRegionName = "trackPvsIhBRegion";
  std::string pVsIhBRegionTitle="Track P vs ih (low P Region)";
  pVsIhBRegionTitle+=";MeV/cm;GeV";
  pVsIhDistributionBRegionHist = fs.make<TH2F>(pVsIhBRegionName.c_str(),pVsIhBRegionTitle.c_str(),400,0,10,100,0,1000);
  pVsIhDistributionBRegionHist->Sumw2();
  // p vs Ias
  TH2F* pVsIasSearchRegionHist = fs.make<TH2F>("pVsIasSearchRegion","P vs Ias (SearchRegion+MassCut);;GeV",20,0,1,40,0,1000);
  pVsIasSearchRegionHist->Sumw2();
  TH2F* trackEtaVsPSearchRegionHist = fs.make<TH2F>("trackEtaVsPSearchRegion","Track #eta vs. p (SearchRegion+MassCut);GeV",4000,0,2000,24,0,2.4);
  trackEtaVsPSearchRegionHist->Sumw2();
  // ias NoM
  //TH1F* iasNoMHist = fs.make<TH1F>("iasNoM","NoM (ias)",50,0,50);
  // ias NoM -- proton
  //TH1F* iasNoMProtonHist = fs.make<TH1F>("iasNoMProton","NoM proton mass cut (ias)",50,0,50);
  // ias NoM - no preselection
  TH1F* iasNoMNoPreselHist = fs.make<TH1F>("iasNoMNoPresel","NoM no preselection (ias)",50,0,50);
  // ias NoM - no preselection, proton mass cut
  TH1F* iasNoMNoPreselProtonHist = fs.make<TH1F>("iasNoMNoPreselProton","NoM no preselection, proton mass cut (ias)",50,0,50);
  // ToF SB - p vs Ias
  TH2F* pVsIasToFSBHist;
  std::string pVsIasToFSBName = "trackPvsIasToFSB";
  std::string pVsIasToFSBTitle="Track P vs ias (ToF SB: #beta>1.075)";
  pVsIasToFSBTitle+=";;GeV";
  pVsIasToFSBHist = fs.make<TH2F>(pVsIasToFSBName.c_str(),pVsIasToFSBTitle.c_str(),400,0,1,100,0,1000);
  pVsIasToFSBHist->Sumw2();
  // ToF SB - p vs Ih
  TH2F* pVsIhToFSBHist;
  std::string pVsIhToFSBName = "trackPvsIhToFSB";
  std::string pVsIhToFSBTitle="Track P vs ih (ToF SB: #beta>1.075)";
  pVsIhToFSBTitle+=";MeV/cm;GeV";
  pVsIhToFSBHist = fs.make<TH2F>(pVsIhToFSBName.c_str(),pVsIhToFSBTitle.c_str(),400,0,10,100,0,1000);
  pVsIhToFSBHist->Sumw2();
  // Ih vs Ias
//  TH2F* ihVsIasHist = fs.make<TH2F>("ihVsIas","Ih vs. Ias;;MeV/cm",400,0,1,400,0,10);
//  ihVsIasHist->Sumw2();
//  TH2F* ihVsIasHist_shifted = fs.make<TH2F>("ihVsIas_shifted","Ih vs. Ias;;MeV/cm",400,0,1,400,0,10);
//  ihVsIasHist_shifted->Sumw2();
  // p vs NoM
  //TH2F* pVsNoMHist = fs.make<TH2F>("pVsNoM","Track P vs. NoM (Ias);;GeV",50,0,50,100,0,1000);
  // p vs NoM, central eta only
  //TH2F* pVsNoMCentralEtaHist = fs.make<TH2F>("pVsNoMCentralEta","Track P vs. NoM (Ias), |#eta| < 0.9;;GeV",50,0,50,100,0,1000);
  // p vs NoM, eta slices
  TH2F* pVsNoMEtaSlice1Hist = fs.make<TH2F>("pVsNoMEtaSlice1","Track P vs. NoM (Ias), |#eta| < 0.2;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice2Hist = fs.make<TH2F>("pVsNoMEtaSlice2","Track P vs. NoM (Ias), 0.2 < |#eta| < 0.4;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice3Hist = fs.make<TH2F>("pVsNoMEtaSlice3","Track P vs. NoM (Ias), 0.4 < |#eta| < 0.6;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice4Hist = fs.make<TH2F>("pVsNoMEtaSlice4","Track P vs. NoM (Ias), 0.6 < |#eta| < 0.8;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice5Hist = fs.make<TH2F>("pVsNoMEtaSlice5","Track P vs. NoM (Ias), 0.8 < |#eta| < 1.0;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice6Hist = fs.make<TH2F>("pVsNoMEtaSlice6","Track P vs. NoM (Ias), 1.0 < |#eta| < 1.2;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice7Hist = fs.make<TH2F>("pVsNoMEtaSlice7","Track P vs. NoM (Ias), 1.2 < |#eta| < 1.4;;GeV",50,0,50,100,0,1000);
  TH2F* pVsNoMEtaSlice8Hist = fs.make<TH2F>("pVsNoMEtaSlice8","Track P vs. NoM (Ias), 1.4 < |#eta| < 1.6;;GeV",50,0,50,100,0,1000);
  //
  // ih vs. p in eta slices
  TH2F* ihVsPEtaSlice1Hist = fs.make<TH2F>("ihVsPEtaSlice1","Ih vs. p, |#eta| < 0.2;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice2Hist = fs.make<TH2F>("ihVsPEtaSlice2","Ih vs. p, 0.2 < |#eta| < 0.4;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice3Hist = fs.make<TH2F>("ihVsPEtaSlice3","Ih vs. p, 0.4 < |#eta| < 0.6;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice4Hist = fs.make<TH2F>("ihVsPEtaSlice4","Ih vs. p, 0.6 < |#eta| < 0.8;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice5Hist = fs.make<TH2F>("ihVsPEtaSlice5","Ih vs. p, 0.8 < |#eta| < 1.0;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice6Hist = fs.make<TH2F>("ihVsPEtaSlice6","Ih vs. p, 1.0 < |#eta| < 1.2;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice7Hist = fs.make<TH2F>("ihVsPEtaSlice7","Ih vs. p, 1.2 < |#eta| < 1.4;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPEtaSlice8Hist = fs.make<TH2F>("ihVsPEtaSlice8","Ih vs. p, 1.4 < |#eta| < 1.6;GeV;MeV/cm",100,0,4,200,0,20);
  // ias  vs. p in eta slices
  TH2F* iasVsPEtaSlice1Hist = fs.make<TH2F>("iasVsPEtaSlice1","Ias vs. p, |#eta| < 0.2;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice2Hist = fs.make<TH2F>("iasVsPEtaSlice2","Ias vs. p, 0.2 < |#eta| < 0.4;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice3Hist = fs.make<TH2F>("iasVsPEtaSlice3","Ias vs. p, 0.4 < |#eta| < 0.6;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice4Hist = fs.make<TH2F>("iasVsPEtaSlice4","Ias vs. p, 0.6 < |#eta| < 0.8;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice5Hist = fs.make<TH2F>("iasVsPEtaSlice5","Ias vs. p, 0.8 < |#eta| < 1.0;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice6Hist = fs.make<TH2F>("iasVsPEtaSlice6","Ias vs. p, 1.0 < |#eta| < 1.2;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice7Hist = fs.make<TH2F>("iasVsPEtaSlice7","Ias vs. p, 1.2 < |#eta| < 1.4;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPEtaSlice8Hist = fs.make<TH2F>("iasVsPEtaSlice8","Ias vs. p, 1.4 < |#eta| < 1.6;GeV;MeV/cm",100,0,4,100,0,1);
  // ih vs. p in nom slices
  TH2F* ihVsPNoMSlice1Hist = fs.make<TH2F>("ihVsPNoMSlice1","Ih vs. p, NoM 5-6;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice2Hist = fs.make<TH2F>("ihVsPNoMSlice2","Ih vs. p, NoM 7-8;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice3Hist = fs.make<TH2F>("ihVsPNoMSlice3","Ih vs. p, NoM 9-10;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice4Hist = fs.make<TH2F>("ihVsPNoMSlice4","Ih vs. p, NoM 11-12;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice5Hist = fs.make<TH2F>("ihVsPNoMSlice5","Ih vs. p, NoM 13-14;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice6Hist = fs.make<TH2F>("ihVsPNoMSlice6","Ih vs. p, NoM 15-16;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice7Hist = fs.make<TH2F>("ihVsPNoMSlice7","Ih vs. p, NoM 17-18;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPNoMSlice8Hist = fs.make<TH2F>("ihVsPNoMSlice8","Ih vs. p, NoM 19-20;GeV;MeV/cm",100,0,4,200,0,20);
  // ias  vs. p in nom slices
  TH2F* iasVsPNoMSlice1Hist = fs.make<TH2F>("iasVsPNoMSlice1","Ias vs. p, NoM 5-6;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice2Hist = fs.make<TH2F>("iasVsPNoMSlice2","Ias vs. p, NoM 7-8;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice3Hist = fs.make<TH2F>("iasVsPNoMSlice3","Ias vs. p, NoM 9-10;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice4Hist = fs.make<TH2F>("iasVsPNoMSlice4","Ias vs. p, NoM 11-12;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice5Hist = fs.make<TH2F>("iasVsPNoMSlice5","Ias vs. p, NoM 13-14;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice6Hist = fs.make<TH2F>("iasVsPNoMSlice6","Ias vs. p, NoM 15-16;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice7Hist = fs.make<TH2F>("iasVsPNoMSlice7","Ias vs. p, NoM 17-18;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPNoMSlice8Hist = fs.make<TH2F>("iasVsPNoMSlice8","Ias vs. p, NoM 19-20;GeV;MeV/cm",100,0,4,100,0,1);
  // ih vs. p in eta slices -- proton mass only
  TH2F* ihVsPProtonEtaSlice1Hist = fs.make<TH2F>("ihVsPProtonEtaSlice1","Ih vs. p, proton mass cut, |#eta| < 0.2;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice2Hist = fs.make<TH2F>("ihVsPProtonEtaSlice2","Ih vs. p, proton mass cut, 0.2 < |#eta| < 0.4;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice3Hist = fs.make<TH2F>("ihVsPProtonEtaSlice3","Ih vs. p, proton mass cut, 0.4 < |#eta| < 0.6;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice4Hist = fs.make<TH2F>("ihVsPProtonEtaSlice4","Ih vs. p, proton mass cut, 0.6 < |#eta| < 0.8;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice5Hist = fs.make<TH2F>("ihVsPProtonEtaSlice5","Ih vs. p, proton mass cut, 0.8 < |#eta| < 1.0;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice6Hist = fs.make<TH2F>("ihVsPProtonEtaSlice6","Ih vs. p, proton mass cut, 1.0 < |#eta| < 1.2;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice7Hist = fs.make<TH2F>("ihVsPProtonEtaSlice7","Ih vs. p, proton mass cut, 1.2 < |#eta| < 1.4;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonEtaSlice8Hist = fs.make<TH2F>("ihVsPProtonEtaSlice8","Ih vs. p, proton mass cut, 1.4 < |#eta| < 1.6;GeV;MeV/cm",100,0,4,200,0,20);
  // ias  vs. p in eta slices -- proton mass only
  TH2F* iasVsPProtonEtaSlice1Hist = fs.make<TH2F>("iasVsPProtonEtaSlice1","Ias vs. p, proton mass cut, |#eta| < 0.2;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice2Hist = fs.make<TH2F>("iasVsPProtonEtaSlice2","Ias vs. p, proton mass cut, 0.2 < |#eta| < 0.4;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice3Hist = fs.make<TH2F>("iasVsPProtonEtaSlice3","Ias vs. p, proton mass cut, 0.4 < |#eta| < 0.6;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice4Hist = fs.make<TH2F>("iasVsPProtonEtaSlice4","Ias vs. p, proton mass cut, 0.6 < |#eta| < 0.8;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice5Hist = fs.make<TH2F>("iasVsPProtonEtaSlice5","Ias vs. p, proton mass cut, 0.8 < |#eta| < 1.0;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice6Hist = fs.make<TH2F>("iasVsPProtonEtaSlice6","Ias vs. p, proton mass cut, 1.0 < |#eta| < 1.2;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice7Hist = fs.make<TH2F>("iasVsPProtonEtaSlice7","Ias vs. p, proton mass cut, 1.2 < |#eta| < 1.4;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonEtaSlice8Hist = fs.make<TH2F>("iasVsPProtonEtaSlice8","Ias vs. p, proton mass cut, 1.4 < |#eta| < 1.6;GeV;MeV/cm",100,0,4,100,0,1);
  // ih vs. p in nom slices -- proton mass only
  TH2F* ihVsPProtonNoMSlice1Hist = fs.make<TH2F>("ihVsPProtonNoMSlice1","Ih vs. p, proton mass cut, NoM 5-6;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice2Hist = fs.make<TH2F>("ihVsPProtonNoMSlice2","Ih vs. p, proton mass cut, NoM 7-8;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice3Hist = fs.make<TH2F>("ihVsPProtonNoMSlice3","Ih vs. p, proton mass cut, NoM 9-10;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice4Hist = fs.make<TH2F>("ihVsPProtonNoMSlice4","Ih vs. p, proton mass cut, NoM 11-12;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice5Hist = fs.make<TH2F>("ihVsPProtonNoMSlice5","Ih vs. p, proton mass cut, NoM 13-14;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice6Hist = fs.make<TH2F>("ihVsPProtonNoMSlice6","Ih vs. p, proton mass cut, NoM 15-16;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice7Hist = fs.make<TH2F>("ihVsPProtonNoMSlice7","Ih vs. p, proton mass cut, NoM 17-18;GeV;MeV/cm",100,0,4,200,0,20);
  TH2F* ihVsPProtonNoMSlice8Hist = fs.make<TH2F>("ihVsPProtonNoMSlice8","Ih vs. p, proton mass cut, NoM 19-20;GeV;MeV/cm",100,0,4,200,0,20);
  // ias  vs. p in nom slices -- proton mass only
  TH2F* iasVsPProtonNoMSlice1Hist = fs.make<TH2F>("iasVsPProtonNoMSlice1","Ias vs. p, proton mass cut, NoM 5-6;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice2Hist = fs.make<TH2F>("iasVsPProtonNoMSlice2","Ias vs. p, proton mass cut, NoM 7-8;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice3Hist = fs.make<TH2F>("iasVsPProtonNoMSlice3","Ias vs. p, proton mass cut, NoM 9-10;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice4Hist = fs.make<TH2F>("iasVsPProtonNoMSlice4","Ias vs. p, proton mass cut, NoM 11-12;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice5Hist = fs.make<TH2F>("iasVsPProtonNoMSlice5","Ias vs. p, proton mass cut, NoM 13-14;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice6Hist = fs.make<TH2F>("iasVsPProtonNoMSlice6","Ias vs. p, proton mass cut, NoM 15-16;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice7Hist = fs.make<TH2F>("iasVsPProtonNoMSlice7","Ias vs. p, proton mass cut, NoM 17-18;GeV;MeV/cm",100,0,4,100,0,1);
  TH2F* iasVsPProtonNoMSlice8Hist = fs.make<TH2F>("iasVsPProtonNoMSlice8","Ias vs. p, proton mass cut, NoM 19-20;GeV;MeV/cm",100,0,4,100,0,1);
  // NoH vs NoM
  //TH2F* nohVsNoMHist = fs.make<TH2F>("nohVsNoM","Number of hits vs. NoM (Ias);NoM;NoH",50,0,50,50,0,50);
  // NoH vs NoM, central eta
  //TH2F* nohVsNoMCentralEtaHist = fs.make<TH2F>("nohVsNoMCentral","Number of hits vs. NoM (Ias), |#eta| < 0.9;NoM;NoH",50,0,50,50,0,50);
  // p vs NoH
  //TH2F* pVsNoHHist = fs.make<TH2F>("pVsNoH","Track P vs. NoH;;GeV",50,0,50,100,0,1000);
  // p vs NoH, central eta only
  //TH2F* pVsNoHCentralEtaHist = fs.make<TH2F>("pVsNoHCentralEta","Track P vs. NoH, |#eta| < 0.9;;GeV",50,0,50,100,0,1000);
  // p vs relPerr
  //TH2F* pVsRelPerrHist = fs.make<TH2F>("pVsRelPerr","Track P vs. #Deltap/p;;GeV",100,0,1,100,0,1000);
  // p vs relPerr, central eta only
  //TH2F* pVsRelPerrCentralEtaHist = fs.make<TH2F>("pVsRelPerrCentralEta","Track P vs. #Deltap/p, |#eta| < 0.9;;GeV",100,0,1,100,0,1000);
  // num hscp gen per event
  TH1F* numHSCPGenPerEventHist = fs.make<TH1F>("numHSCPGenPerEvent","Number of gen HSCPs per event",4,0,4);
  // num hscp seen per event
  TH1F* numHSCPSeenPerEventHist = fs.make<TH1F>("numHSCPSeenPerEvent","Number of seen HSCPs per event",4,0,4);
  // num tracks passing presel in A region per event
  TH1F* numTracksPassingPreselARegionPerEventHist = fs.make<TH1F>("numTracksPassingPreselARegionPerEventHist",
      "Tracks passing preselection in A region per event",10,0,10);
  // num tracks passing presel in B region per event
  TH1F* numTracksPassingPreselBRegionPerEventHist = fs.make<TH1F>("numTracksPassingPreselBRegionPerEventHist",
      "Tracks passing preselection in B region per event",10,0,10);
  // num tracks passing presel in C region per event
  TH1F* numTracksPassingPreselCRegionPerEventHist = fs.make<TH1F>("numTracksPassingPreselCRegionPerEventHist",
      "Tracks passing preselection in C region per event",10,0,10);
  // HSCP eta
  //TH1F* hscpEtaHist = fs.make<TH1F>("hscpEtaHist","#eta (matched to HSCP)",60,-3,3);
  // nom c region
  //TH1F* nomCRegionHist = fs.make<TH1F>("nomCRegionHist","NoM C region (ias < 0.1)",16,5,21);
  // eta c region
  //TH1F* etaCRegionHist = fs.make<TH1F>("etaCRegionHist","Eta C region (ias < 0.1);#eta",60,-3,3);
  // pt c region
  //TH1F* ptCRegionHist = fs.make<TH1F>("ptCRegionHist","Pt C region (ias < 0.1);GeV",200,0,2000);
  // search region (only filled for MC)
  TH1F* pDistributionSearchRegionHist = fs.make<TH1F>("pDistributionSearchRegion","SearchRegion+MassCutP;GeV",200,0,2000);
  pDistributionSearchRegionHist->Sumw2();
  TH1F* ptDistributionSearchRegionHist = fs.make<TH1F>("ptDistributionSearchRegion","SearchRegion+MassCutPt;GeV",200,0,2000);
  ptDistributionSearchRegionHist->Sumw2();
  // beta
  TH1F* hscpBetaDistribution = fs.make<TH1F>("hscpBetaDistribution","HSCP Gen #beta",20,0,1);
  // eta
  TH1F* hscpEtaDistribution = fs.make<TH1F>("hscpEtaDistribution","HSCP Gen #eta",300,-3,3);
  // pt
  TH1F* hscpPtDistribution = fs.make<TH1F>("hscpPtDistribution","HSCP Gen Pt",100,0,2000);
  // p
  TH1F* hscpPDistribution = fs.make<TH1F>("hscpPDistribution","HSCP Gen P",100,0,2000);
  // p vs Ias
  TH2F* pVsIasDistributionSearchRegionHist;
  std::string pVsIasSearchRegionName = "trackPvsIasSearchRegion";
  std::string pVsIasSearchRegionTitle="Track P vs ias (SearchRegion+MassCut)";
  pVsIasSearchRegionTitle+=";;GeV";
  pVsIasDistributionSearchRegionHist = fs.make<TH2F>(pVsIasSearchRegionName.c_str(),pVsIasSearchRegionTitle.c_str(),400,0,1,100,0,1000);
  pVsIasDistributionSearchRegionHist->Sumw2();
  // pt vs Ias
  TH2F* ptVsIasHist;
  std::string ptVsIasName = "trackPtvsIasSignal";
  std::string ptVsIasTitle="Track Pt vs Ias (signal, trig, presel)";
  ptVsIasTitle+=";Ias;GeV";
  ptVsIasHist = fs.make<TH2F>(ptVsIasName.c_str(),ptVsIasTitle.c_str(),400,0,1,200,0,1000);
  ptVsIasHist->Sumw2();
  // pt vs Ias -- mass cut
  TH2F* ptVsIasMassCutHist;
  std::string ptVsIasMassCutName = "trackPtvsIasMassCutSignal";
  std::string ptVsIasMassCutTitle="Track Pt vs Ias after mass cut (signal, trig, presel)";
  ptVsIasMassCutTitle+=";Ias;GeV";
  ptVsIasMassCutHist = fs.make<TH2F>(ptVsIasMassCutName.c_str(),ptVsIasMassCutTitle.c_str(),400,0,1,200,0,1000);
  ptVsIasMassCutHist->Sumw2();
  // p vs ih -- mass cut
  TH2F* pVsIhMassCutHist;
  std::string pVsIhMassCutName = "trackPvsIhMassCutSignal";
  std::string pVsIhMassCutTitle="Track P vs Ih after mass cut (signal, trig, presel)";
  pVsIhMassCutTitle+=";MeV/cm;GeV";
  pVsIhMassCutHist = fs.make<TH2F>(pVsIhMassCutName.c_str(),pVsIhMassCutTitle.c_str(),400,0,20,200,0,1000);
  pVsIhMassCutHist->Sumw2();

  PlotStruct afterPreselectionPlots(fs,"AfterPreselection","after preselection");
  PlotStruct beforePreselectionPlots(fs,"BeforePreselection","before preselection");

  // RooFit observables and dataset
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIp("rooVarIp","ip",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",-2.5,2.5);
  RooRealVar rooVarRun("rooVarRun","run",0,4294967295);
  RooRealVar rooVarLumiSection("rooVarLumiSection","lumiSection",0,4294967295);
  RooRealVar rooVarEvent("rooVarEvent","event",0,4294967295);
  RooRealVar rooVarPUWeight("rooVarPUWeight","puWeight",-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
  RooRealVar rooVarPUSystFactor("rooVarPUSystFactor","puSystFactor",-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());

//  RooDataSet* rooDataSetCandidates = fs.make<RooDataSet>("rooDataSetCandidates","rooDataSetCandidates",
//      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
//      //RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarIp));
//  // Ip no longer included
//  RooDataSet* rooDataSetOneCandidatePerEvent = fs.make<RooDataSet>("rooDataSetOneCandidatePerEvent",
//      "rooDataSetOneCandidatePerEvent",
//      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
//  // keep track of pileup weights
  TTree *rooDataSetCandidatesTree = fs.make<TTree>("rooDataSetCandidatesTree","rooDataSetCandidatesTree");
  TTree *rooDataSetOneCandidatePerEventTree = fs.make<TTree>("rooDataSetOneCandidatePerEventTree","rooDataSetOneCandidatePerEventTree");
  //testing
  float ias_;
  float ih_;
  float trackP_;
  float trackPt_;
  float iasNoM_;
  float trackEta_;
  float eventNumber_;
  float lumiSection_;
  float runNumber_;
  rooDataSetCandidatesTree->Branch("ias_",&ias_,"ias_");
  rooDataSetCandidatesTree->Branch("ih_",&ih_,"ih_");
  rooDataSetCandidatesTree->Branch("trackP_",&trackP_,"trackP_");
  rooDataSetCandidatesTree->Branch("trackPt_",&trackPt_,"trackPt_");
  rooDataSetCandidatesTree->Branch("iasNoM_",&iasNoM_,"iasNoM_");
  rooDataSetCandidatesTree->Branch("trackEta_",&trackEta_,"trackEta_");
  rooDataSetCandidatesTree->Branch("runNumber_",&runNumber_,"runNumber_");
  rooDataSetCandidatesTree->Branch("lumiSection_",&lumiSection_,"lumiSection_");
  rooDataSetCandidatesTree->Branch("eventNumber_",&eventNumber_,"eventNumber_");

  rooDataSetOneCandidatePerEventTree->Branch("ias_",&ias_,"ias_");
  rooDataSetOneCandidatePerEventTree->Branch("ih_",&ih_,"ih_");
  rooDataSetOneCandidatePerEventTree->Branch("trackP_",&trackP_,"trackP_");
  rooDataSetOneCandidatePerEventTree->Branch("trackPt_",&trackPt_,"trackPt_");
  rooDataSetOneCandidatePerEventTree->Branch("iasNoM_",&iasNoM_,"iasNoM_");
  rooDataSetOneCandidatePerEventTree->Branch("trackEta_",&trackEta_,"trackEta_");
  rooDataSetOneCandidatePerEventTree->Branch("runNumber_",&runNumber_,"runNumber_");
  rooDataSetOneCandidatePerEventTree->Branch("lumiSection_",&lumiSection_,"lumiSection_");
  rooDataSetOneCandidatePerEventTree->Branch("eventNumber_",&eventNumber_,"eventNumber_");
  RooDataSet* rooDataSetCandidates = new RooDataSet ("rooDataSetCandidates","rooDataSetCandidates",
      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
      //RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarIp));
  // Ip no longer included
  RooDataSet* rooDataSetOneCandidatePerEvent = fs.make<RooDataSet>("rooDataSetOneCandidatePerEvent",
      "rooDataSetOneCandidatePerEvent",
      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
//  RooDataSet* rooDataSetOneCandidatePerEvent = new RooDataSet ("rooDataSetOneCandidatePerEvent",
//      "rooDataSetOneCandidatePerEvent",
//      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
  // keep track of pileup weights
  RooDataSet* rooDataSetPileupWeights = fs.make<RooDataSet>("rooDataSetPileupWeights","rooDataSetPileupWeights",
      RooArgSet(rooVarRun,rooVarLumiSection,rooVarEvent,rooVarPUWeight,rooVarPUSystFactor));

  // Dataset for IasShift
  RooDataSet* rooDataSetIasShift = fs.make<RooDataSet>("rooDataSetIasShift","rooDataSetIasShift",
      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
  // Dataset for IhShift
  RooDataSet* rooDataSetIhShift = fs.make<RooDataSet>("rooDataSetIhShift","rooDataSetIhShift",
      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));
  // Dataset for PtShift
  RooDataSet* rooDataSetPtShift = fs.make<RooDataSet>("rooDataSetPtShift","rooDataSetPtShift",
      RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarRun,rooVarLumiSection,rooVarEvent));

  // RooDataSet for number of original tracks, etc.
  RooRealVar rooVarNumGenHSCPEvents("rooVarNumGenHSCPEvents","numGenHSCPEvents",0,1e10);
  RooRealVar rooVarNumGenHSCPTracks("rooVarNumGenHSCPTracks","numGenHSCPTracks",0,1e10);
  RooRealVar rooVarNumGenChargedHSCPTracks("rooVarNumGenChargedHSCPTracks","numGenChargedHSCPTracks",0,1e10);
  RooRealVar rooVarSignalEventCrossSection("rooVarSignalEventCrossSection","signalEventCrossSection",0,100); // pb
  RooRealVar rooVarNumGenHSCPEventsPUReweighted("rooVarNumGenHSCPEventsPUReweighted","numGenHSCPEventsPUReweighted",0,1e10);
  RooRealVar rooVarEventWeightSum("rooVarEventWeightSum","eventWeightSum",0,1e10);
  RooRealVar rooVarSampleWeight("rooVarSampleWeight","sampleWeight",0,1e10);
  RooDataSet* rooDataSetGenHSCPTracks = fs.make<RooDataSet>("rooDataSetGenHSCPTracks","rooDataSetGenHSCPTracks",
      RooArgSet(rooVarNumGenHSCPEvents,rooVarNumGenHSCPTracks,rooVarNumGenChargedHSCPTracks,
        rooVarSignalEventCrossSection,rooVarNumGenHSCPEventsPUReweighted,rooVarEventWeightSum,
        rooVarSampleWeight));

  //for initializing PileUpReweighting utility.
  const float TrueDist2011_f[35] = {0.00285942, 0.0125603, 0.0299631, 0.051313, 0.0709713, 0.0847864, 0.0914627, 0.0919255, 0.0879994, 0.0814127, 0.0733995, 0.0647191, 0.0558327, 0.0470663, 0.0386988, 0.0309811, 0.0241175, 0.018241, 0.0133997, 0.00956071, 0.00662814, 0.00446735, 0.00292946, 0.00187057, 0.00116414, 0.000706805, 0.000419059, 0.000242856, 0.0001377, 7.64582e-05, 4.16101e-05, 2.22135e-05, 1.16416e-05, 5.9937e-06, 5.95542e-06};//from 2011 Full dataset
  const float Pileup_MC[35]= {1.45346E-01, 6.42802E-02, 6.95255E-02, 6.96747E-02, 6.92955E-02, 6.84997E-02, 6.69528E-02, 6.45515E-02, 6.09865E-02, 5.63323E-02, 5.07322E-02, 4.44681E-02, 3.79205E-02, 3.15131E-02, 2.54220E-02, 2.00184E-02, 1.53776E-02, 1.15387E-02, 8.47608E-03, 6.08715E-03, 4.28255E-03, 2.97185E-03, 2.01918E-03, 1.34490E-03, 8.81587E-04, 5.69954E-04, 3.61493E-04, 2.28692E-04, 1.40791E-04, 8.44606E-05, 5.10204E-05, 3.07802E-05, 1.81401E-05, 1.00201E-05, 5.80004E-06};
  const float TrueDist2012_f[60] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09};  // MB xsec = 69.3mb
  const float TrueDist2012_XSecShiftUp_f[60] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09}; // MB xsec = 73.5mb; observed in Z-->MuMu see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Calculating_Your_Pileup_Distribu
//  const   float TrueDist2012_XSecShiftDown_f[60] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09}; // MB xsec = 65.835mb
  const float Pileup_MC_Summer2012[60] = { 2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06};

  edm::LumiReWeighting LumiWeightsMC;
  edm::LumiReWeighting LumiWeightsMCSyst;
  std::vector<float> BgLumiMC; //MC                                           
  std::vector<float> TrueDist;  
  std::vector< float > TrueDistSyst;
  //initialize LumiReWeighting
  if(!is8TeV_)
  {
    for(int i=0; i<35; ++i)
      BgLumiMC.push_back(Pileup_MC[i]);
  }
  else
  {
    for(int i=0; i<60; ++i) 
      BgLumiMC.push_back(Pileup_MC_Summer2012[i]);
  }
  if(!is8TeV_)
  {
    for(int i=0; i<35; ++i)
      TrueDist.push_back(TrueDist2011_f[i]);
  }
  else
  {
    for(int i=0; i<60; ++i) 
      TrueDist.push_back(TrueDist2012_f[i]);
    for(int i=0; i<60; ++i) 
      TrueDistSyst.push_back(TrueDist2012_XSecShiftUp_f[i]);
  }
  LumiWeightsMC = edm::LumiReWeighting(BgLumiMC, TrueDist);
  LumiWeightsMCSyst = edm::LumiReWeighting(BgLumiMC, TrueDistSyst);
  //testing - not sure what this s4 or s10 means
  bool IsS4pileup = false;
  bool IsS10pileup = false;
  if(!is8TeV_)
  {
    IsS4pileup = true; // seems to be true for all 2011 signal MC
  }
  if(is8TeV_)
  {
    IsS10pileup = true; // seems to be true for all 2012 signal MC
  }

  reweight::PoissonMeanShifter PShift_(0.6);//0.6 for upshift, -0.6 for downshift

  // initialize ih shifting functions
  initializeIhPeakDiffFunctions(mcDataIhPeakShiftFuncs,8,is8TeV_);
  initializeIhWidthDiffFunction(is8TeV_);

  // initialize ias shifting functions
  std::string mcDataIasDiffGraphsRootFile = "iasDataMCDiffsCombined.may24.root"; // produced with my (SC) own macro (on laptop dedxSystematics dir)
  initializeIasShiftGraphs(mcDataIasDiffGraphsRootFile);
  initializeIasWidthDiffFunction(); //testing for 8 TeV
  initializeIasPeakDiffFunction(); //testing for 8 TeV

  double sampleWeight = 1;
  int period = 0;

  // get PU-reweighted total events
  double numPUReweightedMCEvents = 0;
  int ievt = 0;
  if(isMC_)
  {
    for(unsigned int iFile=0; iFile<inputHandler_.files().size(); ++iFile)
    {
      // open input file (can be located on castor)
      TFile* inFile = TFile::Open(inputHandler_.files()[iFile].c_str());
      if( inFile )
      {
        if(inputHandler_.files()[iFile].find("BX1") != string::npos || is8TeV_)
        {
          std::cout << "Found BX1-->period = 1" << std::endl;
          period = 1;
        }
        else
          std::cout << "Did not find BX1-->period = 0" << std::endl;
        fwlite::Event ev(inFile);
        for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
        {
          // break loop if maximal number of events is reached 
          if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;

          // Get PU Weight
          double PUWeight_thisevent;
          double PUSystFactor;
          if(IsS4pileup)
          {
            PUWeight_thisevent = GetPUWeight(ev, "S4", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
          }
          if(IsS10pileup)
          {
            PUWeight_thisevent = GetPUWeight(ev, "S10", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
          }
          if(ievt % 100 == 0)
            std::cout << "ientry=" << ievt << ", puwt=" << PUWeight_thisevent << ", NMCevents="
              << numPUReweightedMCEvents << std::endl;
          if(ievt < 100)
            std::cout << "ientry=" << ievt << ", puwt=" << PUWeight_thisevent << ", NMCevents="
              << numPUReweightedMCEvents << std::endl;
          numPUReweightedMCEvents+=PUWeight_thisevent;
        } // loop over events
        // close input file
        inFile->Close();
      }
      // break loop if maximal number of events is reached:
      // this has to be done twice to stop the file loop as well
      if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
    } // end file loop
    sampleWeight = GetSampleWeight(integratedLumi_,integratedLumiBeforeTriggerChange_,
        signalEventCrossSection_,numPUReweightedMCEvents,period,is8TeV_);
    rooVarSampleWeight = sampleWeight;

  } // is mc
  std::cout << "NumPUReweightedMCEvents: " << numPUReweightedMCEvents << std::endl;

  std::cout << "numTreeEntries: " << ievt << std::endl;
  std::cout << "IsS4pileup: " << IsS4pileup << std::endl;
  std::cout << "Int lumi = " << integratedLumi_ << " int lumi bef. trig. change = " << 
    integratedLumiBeforeTriggerChange_ << std::endl;
  std::cout << "Cross section = " << signalEventCrossSection_ << std::endl;
  std::cout << "NumPUReweightedMCEvents: " << numPUReweightedMCEvents << std::endl;
  std::cout << "Period = " << period << std::endl;
  std::cout << "SampleWeight: " << sampleWeight << std::endl;


  int numEventsPassingTrigger = 0;
  int numEventsWithOneTrackPassingPreselection = 0;
  int numEventsWithOneTrackOverIhCut = 0;
  int numEventsWithOneTrackOverIasCut = 0;
  int numEventsWithOneTrackOverIhAndIasCuts = 0;

  int numGenHSCPTracks = 0;
  int numGenChargedHSCPTracks = 0;
  int numGenHSCPEvents = 0;
  int numTracksPassingTrigger = 0;
  int numTracksPassingPreselection = 0;
  int numTracksInSearchRegion = 0;

  int numEventsWithOneHSCPGen = 0;
  int numEventsWithTwoHSCPGen = 0;
  int numEventsWithOneHSCPSeenMuonTrigger = 0;
  int numEventsWithTwoHSCPSeenMuonTrigger = 0;
  int numEventsWithOneHSCPSeenMetTrigger = 0;
  int numEventsWithTwoHSCPSeenMetTrigger = 0;
  double eventWeightSum = 0;
  // main event loop
  // loop the events
  ievt = 0;
  for(unsigned int iFile=0; iFile<inputHandler_.files().size(); ++iFile)
  {
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputHandler_.files()[iFile].c_str());
    if( inFile )
    {
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt)
      {
        double tempIas = -1;
        double tempIh = 0;
        double tempP = 0;
        double tempPt = 0;
        int tempNoMias = 0;
        double tempEta = 0;
        double tempRun = 0;
        double tempLumiSection = 0;
        double tempEvent = 0;
        double tempSystPt = 0;
        double tempSystP = 0;
        double tempSystIh = 0;
        double eventWeight = 1;
        // systematics datasets
        double iasSystIas = -1;
        double iasSystIh = 0;
        double iasSystP = 0;
        double iasSystPt = 0;
        double iasSystNoMias = 0;
        double iasSystEta = 0;
        double iasSystLumiSection = 0;
        double iasSystRun = 0;
        double iasSystEvent = 0;
        double highestShiftedIas = -1;

        // break loop if maximal number of events is reached 
        if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
        // simple event counter
        if(inputHandler_.reportAfter()!=0 ? (ievt>0 && ievt%inputHandler_.reportAfter()==0) : false) 
          std::cout << "  processing event: " << ievt << std::endl;

        // MC Gen Info
        fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
        std::vector<reco::GenParticle> genColl;
        if(isMC_)
        {
          genCollHandle.getByLabel(ev, "genParticles");
          if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
          genColl = *genCollHandle;
            // Get PU Weight
            double PUWeight_thisevent;
            double PUSystFactor;
            if(IsS4pileup)
            {
              PUWeight_thisevent = GetPUWeight(ev, "S4", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
            }
            if(IsS10pileup)
            {
              PUWeight_thisevent = GetPUWeight(ev, "S10", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
            }
            eventWeight = sampleWeight * PUWeight_thisevent;
            if(ievt < 100)
              std::cout << "ientry=" << ievt << " eventWeight = " << eventWeight << " totalE="
                << eventWeightSum << std::endl;
            eventWeightSum+=eventWeight;
            rooVarPUWeight = PUWeight_thisevent;
            rooVarPUSystFactor = PUSystFactor;
        }
        else
        {
          genCollHandle.getByLabel(ev, "genParticles");
          if(genCollHandle.isValid())
          {
            std::cout << "ERROR: GenParticle Collection Found when MC not indicated" << std::endl;
            return -5;
          }
        }

        fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
        hscpCollHandle.getByLabel(ev,"HSCParticleProducer");
        if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
        const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

        fwlite::Handle<reco::DeDxDataValueMap> dEdxSCollH;
        dEdxSCollH.getByLabel(ev, "dedxASmi");
        if(!dEdxSCollH.isValid()){printf("Invalid dEdx Selection collection\n");continue;}

        fwlite::Handle<reco::DeDxDataValueMap> dEdxIpCollH;
        dEdxIpCollH.getByLabel(ev, "dedxProd");
        if(!dEdxIpCollH.isValid()){printf("Invalid dEdx Selection (Ip) collection\n");continue;}

        fwlite::Handle<reco::DeDxDataValueMap> dEdxMCollH;
        dEdxMCollH.getByLabel(ev, "dedxHarm2");
        if(!dEdxMCollH.isValid()){printf("Invalid dEdx Mass collection\n");continue;}

        // ToF
        fwlite::Handle<reco::MuonTimeExtraMap> TOFCollH;
        TOFCollH.getByLabel(ev, "muontiming","combined");
        if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");continue;}

        fwlite::Handle<reco::MuonTimeExtraMap> TOFDTCollH;
        TOFDTCollH.getByLabel(ev, "muontiming","dt");
        if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

        fwlite::Handle<reco::MuonTimeExtraMap> TOFCSCCollH;
        TOFCSCCollH.getByLabel(ev, "muontiming","csc");
        if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");continue;}

        double genMass = -1;

        // Look for gen HSCP tracks if MC
        if(isMC_)
        {
          int numGenHSCPThisEvent = getNumGenHSCP(genColl,false); // consider also neutrals
          int numGenChargedHSCPThisEvent = getNumGenHSCP(genColl,true); // consider only charged HSCP
          numGenHSCPTracks+=numGenHSCPThisEvent;
          numGenChargedHSCPTracks+=numGenChargedHSCPThisEvent;
          if(numGenHSCPThisEvent>=1)
            numEventsWithOneHSCPGen++;
          if(numGenHSCPThisEvent>=2)
            numEventsWithTwoHSCPGen++;
          numHSCPGenPerEventHist->Fill(numGenHSCPThisEvent);

          int ClosestGen = -1;
          int numSeenHSCPThisEvt = 0;
          // attempt match track to gen hscp, look at gen/detected hscps
          for(unsigned int c=0;c<hscpColl.size();c++)
          {
            susybsm::HSCParticle hscp  = hscpColl[c];
            reco::MuonRef  muon  = hscp.muonRef();
            reco::TrackRef track = hscp.trackRef();
            if(track.isNull())
              continue;

            // match to gen hscp
            if(DistToHSCP(hscp, genColl, ClosestGen)>0.03)
              continue;
            else
              numSeenHSCPThisEvt++;

          }
          numHSCPSeenPerEventHist->Fill(numSeenHSCPThisEvt);
          // loop over gen particles
          for(unsigned int g=0;g<genColl.size();g++)
          {
            if(genColl[g].pt()<5)continue;
            if(genColl[g].status()!=1)continue;
            int AbsPdg=abs(genColl[g].pdgId());
            if(AbsPdg<1000000)continue;
            //if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))continue; //Skip neutral gluino RHadrons
            //if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))continue;  //skip neutral stop RHadrons
            //if(beta1<0){beta1=genColl[g].p()/genColl[g].energy();}else if(beta2<0){beta2=genColl[g].p()/genColl[g].energy();return;}
            genMass = genColl[g].mass();
            double genP = sqrt(pow(genColl[g].px(),2) + pow(genColl[g].py(),2) + pow(genColl[g].pz(),2) );
            hscpBetaDistribution->Fill(genP/sqrt(pow(genP,2)+pow(genMass,2)));
            hscpEtaDistribution->Fill(genColl[g].eta());
            hscpPtDistribution->Fill(genColl[g].pt());
            hscpPDistribution->Fill(genP);
          }

          // check triggers
          if(passesTrigger(ev,true,false,is8TeV_)) // consider mu trig only
          {
            if(numSeenHSCPThisEvt==1)
              numEventsWithOneHSCPSeenMuonTrigger++;
            if(numSeenHSCPThisEvt==2)
              numEventsWithTwoHSCPSeenMuonTrigger++;
          }
          //TODO fix this for 8 TeV
          if(passesTrigger(ev,false,true,is8TeV_)) // consider MET trig only
          {
            if(numSeenHSCPThisEvt==1)
              numEventsWithOneHSCPSeenMetTrigger++;
            if(numSeenHSCPThisEvt==2)
              numEventsWithTwoHSCPSeenMetTrigger++;
          }
        }

        numGenHSCPEvents++;

        double lumiSection = ev.id().luminosityBlock();
        double runNumber = ev.id().run();
        double eventNumber = ev.id().event();
        // check trigger
//        if(!passesTrigger(ev,true,true,is8TeV_))
//          continue;

        numEventsPassingTrigger++;

        int numTracksPassingPreselectionThisEvent = 0;
        int numTracksPassingIhCutThisEvent = 0;
        int numTracksPassingIasCutThisEvent = 0;
        int numTracksPassingIhAndIasCutsThisEvent = 0;
        int numTracksPassingPreselectionARegionThisEvent = 0;
        int numTracksPassingPreselectionBRegionThisEvent = 0;
        int numTracksPassingPreselectionCRegionThisEvent = 0;
        // loop over HSCParticles in this event
        for(unsigned int c=0;c<hscpColl.size();c++)
        {
          susybsm::HSCParticle hscp  = hscpColl[c];
          reco::MuonRef  muon  = hscp.muonRef();
          reco::TrackRef track = hscp.trackRef();
          if(track.isNull())
            continue;

          numTracksPassingTrigger++;

          const reco::DeDxData& dedxSObj  = dEdxSCollH->get(track.key());
          const reco::DeDxData& dedxMObj  = dEdxMCollH->get(track.key());
          const reco::DeDxData& dedxIpObj  = dEdxIpCollH->get(track.key());
          const reco::MuonTimeExtra* tof  = &TOFCollH->get(hscp.muonRef().key());
          const reco::MuonTimeExtra* dttof  = &TOFDTCollH->get(hscp.muonRef().key());
          const reco::MuonTimeExtra* csctof  = &TOFCSCCollH->get(hscp.muonRef().key());

          // for signal MC
          int NChargedHSCP = 0;
          if(isMC_ && matchToHSCP_)
          {
            int ClosestGen;
            if(DistToHSCP(hscp, genColl, ClosestGen)>0.03)
              continue;
           NChargedHSCP = getNumGenHSCP(genColl,true);
          }

          float ih = dedxMObj.dEdx();
          float ias = dedxSObj.dEdx();
          float ip = dedxIpObj.dEdx();
          float trackP = track->p();
          float trackPt = track->pt();
          int iasNoM = dedxSObj.numberOfMeasurements();
          float trackEta = track->eta();
          int trackNoH = track->found();
          float trackPtErr = track->ptError();
          
          // fill ias before preselection
          iasNoMNoPreselHist->Fill(iasNoM);
          // mass
          float massSqr = (ih-dEdx_c_)*pow(trackP,2)/dEdx_k_;
          float mass = (massSqr >= 0) ? sqrt(massSqr) : 0;
          // proton mass
          bool isProtonMass = (mass < 1.5 && mass > 0.65) ? true : false;
          if(isProtonMass)
            iasNoMNoPreselProtonHist->Fill(iasNoM);

          // apply preselections, not considering ToF
          //if(!passesPreselection(hscp,dedxSObj,dedxMObj,tof,dttof,csctof,ev,false,is8TeV_,beforePreselectionPlots))
//          if(!passesPreselection(hscp,dedxSObj,dedxMObj,tof,dttof,csctof,ev,false,true,beforePreselectionPlots))
          if(!passesPreselection(hscp,dedxSObj,dedxMObj,tof,dttof,csctof,ev,false,is8TeV_,beforePreselectionPlots))
            continue;

          // systematics datasets for MC
          // can look after applying preselection only since preselection does not use Ias, Ih, or P/Pt directly
          double shiftedIas = ias;
          double shiftedIh = ih;
          double shiftedPt = trackPt;
          double shiftedP = trackP;
          // apply shift, then choose highest ias (same track for P/Pt/Ih shifts, could be different for Ias shift)
          if(isMC_)
          {
            TRandom3 myRandom;
            // include TOF at some point?
            //rooVarIas = ias;
            //rooVarIp = ip;
            //rooVarIh = ih;
            //rooVarP = trackP;
            //rooVarPt = trackPt;
            //rooVarNoMias = iasNoM;
            //rooVarEta = fabs(trackEta);
            double beta = trackP/sqrt(pow(trackP,2)+pow(genMass,2));
            // ias shift
            //shiftedIas = ias + myRandom.Gaus(0,0.083) + 0.015; // from YK results Nov 21 2011 hypernews thread
            double iasPeakShift = getIasPeakCorrection(iasNoM,beta,is8TeV_);
            double iasWidthShift = getIasWidthCorrection(iasNoM,beta,is8TeV_);
            shiftedIas = ias + iasPeakShift + myRandom.Gaus(0,iasWidthShift);
            // ih shift
            //shiftedIh = ih*1.036; // from SIC results, Nov 10 2011 HSCP meeting
            double ihPeakShift = getIhPeakCorrection(iasNoM,beta);
            double ihWidthShift = getIhWidthCorrection(beta);
            shiftedIh = ih + ihPeakShift + myRandom.Gaus(0,ihWidthShift);
            // pt shift (and p shift) -- from MU-10-004-001 -- SIC report Jun 7 2011 HSCP meeting
            double newInvPt = 1/trackPt+0.000236-0.000135*pow(trackEta,2)+track->charge()*0.000282*TMath::Sin(track->phi()-1.337);
            shiftedPt = 1.0/newInvPt;
            shiftedP = shiftedPt/TMath::Sin(track->theta());
          }

          if(isMC_)
          {
            ptVsIasHist->Fill(ias,trackPt);
//            float massSqr = (ih-dEdx_c_)*pow(trackP,2)/dEdx_k_;
            if(sqrt(massSqr) >= massCutIasHighPHighIh_)
            {
              ptVsIasMassCutHist->Fill(ias,trackPt);
              pVsIhMassCutHist->Fill(ih,trackP);
            }
          }



          numTracksPassingPreselection++;
          numTracksPassingPreselectionThisEvent++;
          if(ih < ihSidebandThreshold_ && trackP < pSidebandThreshold_)
            numTracksPassingPreselectionARegionThisEvent++;
          else if(trackP < pSidebandThreshold_)
            numTracksPassingPreselectionBRegionThisEvent++;
          else if(ih < ihSidebandThreshold_)
            numTracksPassingPreselectionCRegionThisEvent++;

          if(ih > 3.5)
            numTracksPassingIhCutThisEvent++;
          if(ias > 0.35)
            numTracksPassingIasCutThisEvent++;
          if(ih > 3.5 && ias > 0.35)
            numTracksPassingIhAndIasCutsThisEvent++;

          if(trackP > pSidebandThreshold_ && ih > ihSidebandThreshold_)
          {
            numTracksInSearchRegion++;
            if(isMC_)
            {
              if(massSqr < 0)
                continue;
              else if(sqrt(massSqr) >= massCutIasHighPHighIh_)
              {
                pDistributionSearchRegionHist->Fill(trackP);
                ptDistributionSearchRegionHist->Fill(trackPt);
                trackEtaVsPSearchRegionHist->Fill(trackP,fabs(trackEta));
                pVsIasSearchRegionHist->Fill(ias,trackP);
                pVsIasDistributionSearchRegionHist->Fill(ias,trackP);
                pVsIhDistributionSearchRegionHist->Fill(ih,trackP);
              }
            }
          }

          // A/B region
          if(trackP <= pSidebandThreshold_)
            pVsIhDistributionBRegionHist->Fill(ih,trackP);

          // fill some overall hists (already passed preselection)
          afterPreselectionPlots.pDistributionHist->Fill(trackP);
          afterPreselectionPlots.ptDistributionHist->Fill(trackPt);
          afterPreselectionPlots.trackEtaVsPHist->Fill(trackP,fabs(trackEta));
          afterPreselectionPlots.pVsIasHist->Fill(ias,trackP);
          afterPreselectionPlots.pVsIhHist->Fill(ih,trackP);
          afterPreselectionPlots.iasNoMHist->Fill(iasNoM);
          afterPreselectionPlots.ihVsIasHist->Fill(ias,ih);
          afterPreselectionPlots.pVsNoMHist->Fill(iasNoM,trackP);
          if(isMC_)
          {
            afterPreselectionPlots.pVsIasHist_shifted->Fill(shiftedIas,trackP);
            afterPreselectionPlots.pVsIhHist_shifted->Fill(shiftedIh,trackP);
            afterPreselectionPlots.ihVsIasHist_shifted->Fill(shiftedIas,shiftedIh);
          }
//          float massSqr = (ih-dEdx_c_)*pow(trackP,2)/dEdx_k_;
          if(massSqr >= 0)
            afterPreselectionPlots.massHist->Fill(sqrt(massSqr));
          // fill pt vs ias only for ABC regions
          if(!(ias > 0.1 && trackPt > 50))
            afterPreselectionPlots.ptVsIasHist->Fill(ias,trackPt);
          if(fabs(trackEta) < 0.9)
            afterPreselectionPlots.pVsNoMCentralEtaHist->Fill(iasNoM,trackP);
          if(fabs(trackEta) < 0.2)
          {
            pVsNoMEtaSlice1Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice1Hist->Fill(trackP,ih);
            iasVsPEtaSlice1Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice1Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice1Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 0.4)
          {
            pVsNoMEtaSlice2Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice2Hist->Fill(trackP,ih);
            iasVsPEtaSlice2Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice2Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice2Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 0.6)
          {
            pVsNoMEtaSlice3Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice3Hist->Fill(trackP,ih);
            iasVsPEtaSlice3Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice3Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice3Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 0.8)
          {
            pVsNoMEtaSlice4Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice4Hist->Fill(trackP,ih);
            iasVsPEtaSlice4Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice4Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice4Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 1.0)
          {
            pVsNoMEtaSlice5Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice5Hist->Fill(trackP,ih);
            iasVsPEtaSlice5Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice5Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice5Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 1.2)
          {
            pVsNoMEtaSlice6Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice6Hist->Fill(trackP,ih);
            iasVsPEtaSlice6Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice6Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice6Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 1.4)
          {
            pVsNoMEtaSlice7Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice7Hist->Fill(trackP,ih);
            iasVsPEtaSlice7Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice7Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice7Hist->Fill(trackP,ias);
            }
          }
          else if(fabs(trackEta) < 1.6)
          {
            pVsNoMEtaSlice8Hist->Fill(iasNoM,trackP);
            ihVsPEtaSlice8Hist->Fill(trackP,ih);
            iasVsPEtaSlice8Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonEtaSlice8Hist->Fill(trackP,ih);
              iasVsPProtonEtaSlice8Hist->Fill(trackP,ias);
            }
          }
          // NoM
          if(iasNoM==5||iasNoM==6)
          {
            ihVsPNoMSlice1Hist->Fill(trackP,ih);
            iasVsPNoMSlice1Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice1Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice1Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==7||iasNoM==8)
          {
            ihVsPNoMSlice2Hist->Fill(trackP,ih);
            iasVsPNoMSlice2Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice2Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice2Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==9||iasNoM==10)
          {
            ihVsPNoMSlice3Hist->Fill(trackP,ih);
            iasVsPNoMSlice3Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice3Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice3Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==11||iasNoM==12)
          {
            ihVsPNoMSlice4Hist->Fill(trackP,ih);
            iasVsPNoMSlice4Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice4Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice4Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==13||iasNoM==14)
          {
            ihVsPNoMSlice5Hist->Fill(trackP,ih);
            iasVsPNoMSlice5Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice5Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice5Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==15||iasNoM==16)
          {
            ihVsPNoMSlice6Hist->Fill(trackP,ih);
            iasVsPNoMSlice6Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice6Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice6Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==17||iasNoM==18)
          {
            ihVsPNoMSlice7Hist->Fill(trackP,ih);
            iasVsPNoMSlice7Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice7Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice7Hist->Fill(trackP,ias);
            }
          }
          else if(iasNoM==19||iasNoM==20)
          {
            ihVsPNoMSlice8Hist->Fill(trackP,ih);
            iasVsPNoMSlice8Hist->Fill(trackP,ias);
            if(isProtonMass)
            {
              ihVsPProtonNoMSlice8Hist->Fill(trackP,ih);
              iasVsPProtonNoMSlice8Hist->Fill(trackP,ias);
            }
          }
          afterPreselectionPlots.nohVsNoMHist->Fill(iasNoM,trackNoH);
          if(fabs(trackEta) < 0.9)
            afterPreselectionPlots.nohVsNoMCentralEtaHist->Fill(iasNoM,trackNoH);
          afterPreselectionPlots.pVsNoHHist->Fill(trackNoH,trackP);
          if(fabs(trackEta) < 0.9)
            afterPreselectionPlots.pVsNoHCentralEtaHist->Fill(trackNoH,trackP);
          afterPreselectionPlots.pVsRelPerrHist->Fill(trackPtErr/trackPt,trackP);
          if(fabs(trackEta) < 0.9)
            afterPreselectionPlots.pVsRelPerrCentralEtaHist->Fill(trackPtErr/trackP,trackP);


          // fill roodataset
          //TODO include TOF at some point?
          rooVarIas = ias;
          rooVarIp = ip;
          rooVarIh = ih;
          rooVarP = trackP;
          rooVarPt = trackPt;
          rooVarNoMias = iasNoM;
          rooVarEta = trackEta;
          rooVarLumiSection = lumiSection;
          rooVarRun = runNumber;
          rooVarEvent = eventNumber;
          rooDataSetCandidates->add(RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarLumiSection,
                rooVarRun,rooVarEvent));
//          std::cout << "Ias: " << ias << std::endl;
          ias_ = ias;
          ih_ = ih;
          trackP_ = trackP;
          trackPt_ = trackPt;
          iasNoM_ = iasNoM;
          trackEta_ = trackEta;
          lumiSection_ = lumiSection;
          runNumber_ = runNumber;
          eventNumber_ = eventNumber;
          rooDataSetCandidatesTree->Fill();
//          rooDataSetTree->Print();
          if(ias > tempIas)
          {
            tempIas = ias;
            tempIh = ih;
            tempP = trackP;
            tempPt = trackPt;
            tempNoMias = iasNoM;
            tempEta = trackEta;
            tempLumiSection = lumiSection;
            tempRun = runNumber;
            tempEvent = eventNumber;
            tempSystPt = shiftedPt;
            tempSystP = shiftedP;
            tempSystIh = shiftedIh;
          }
          // for ias shift syst., have to select highest Ias track separately
          // for other shifts, highest ias track will be the same (use above section)
          if(isMC_ && shiftedIas > highestShiftedIas)
          {
            iasSystIas = ias;
            iasSystIh = ih;
            iasSystP = trackP;
            iasSystPt = trackPt;
            iasSystNoMias = iasNoM;
            iasSystEta = trackEta;
            iasSystLumiSection = lumiSection;
            iasSystRun = runNumber;
            iasSystEvent = eventNumber;
          }

          //// now consider the ToF
          //if(!passesPreselection(hscp,dedxSObj,dedxMObj,tof,dttof,csctof,ev,true,is8TeV_))
          //  continue;

          // require muon ref when accessing ToF, or the collections will be null
          if(hscp.muonRef().isNull())
            continue;
          // ndof check
          if(tof->nDof() < 8 && (dttof->nDof() < 6 || csctof->nDof() < 6) )
            continue;
          // select 1/beta low --> beta high
          if(tof->inverseBeta() >= 0.93)
            continue;
          // Max ToF err
          if(tof->inverseBetaErr() > 0.07)
            continue;

          // now we should have a sample where 1/beta < 0.93 --> beta > 1.075
          afterPreselectionPlots.pVsIasToFSBHist->Fill(ias,trackP);
          afterPreselectionPlots.pVsIhToFSBHist->Fill(ih,trackP);

        } // done looking at HSCParticle collection

        if(numTracksPassingPreselectionThisEvent > 0)
        {
          rooVarIas = tempIas;
          rooVarIh = tempIh;
          rooVarP = tempP;
          rooVarPt = tempPt;
          rooVarNoMias = tempNoMias;
          rooVarEta = tempEta;
          rooVarLumiSection = tempLumiSection;
          rooVarRun = tempRun;
          rooVarEvent = tempEvent;
          rooDataSetOneCandidatePerEvent->add(RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarLumiSection,
                rooVarRun,rooVarEvent));
          ias_ = tempIas;
          ih_ = tempIh;
          trackP_ = tempP;
          trackPt_ = tempPt;
          iasNoM_ = tempNoMias;
          trackEta_ = tempEta;
          lumiSection_ = tempLumiSection;
          runNumber_ = tempRun;
          eventNumber_ = tempEvent;
          rooDataSetOneCandidatePerEventTree->Fill();
          if(isMC_)
          {
            //cout << "INFO: insert event - run: " << rooVarRun.getVal() << " lumiSec: " << rooVarLumiSection.getVal()
            //  << " eventNum: " << rooVarEvent.getVal() << endl;
            //cout << "add puweight=" << rooVarPUWeight.getVal() << " eventWeight=" << rooVarPUWeight.getVal() * sampleWeight << std::endl;
            rooDataSetPileupWeights->add(RooArgSet(rooVarEvent,rooVarLumiSection,rooVarRun,rooVarPUWeight,rooVarPUSystFactor));
            // add tracks to systematics datasets with shifts
            // ih shift
            rooVarIh = tempSystIh;
            rooDataSetIhShift->add(RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarLumiSection,rooVarRun,rooVarEvent));
            rooVarIh = tempIh; // reset to original value
            // pt/p shift dataset
            rooVarPt = tempSystPt;
            rooVarP = tempSystP;
            rooDataSetPtShift->add(RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarLumiSection,rooVarRun,rooVarEvent));
            rooVarPt = tempPt; // reset
            rooVarP = tempP; // reset
            // ias shift dataset
            rooVarIas = iasSystIas;
            rooVarIh = iasSystIh;
            rooVarP = iasSystP;
            rooVarPt = iasSystPt;
            rooVarNoMias = iasSystNoMias;
            rooVarEta = iasSystEta;
            rooVarLumiSection = iasSystLumiSection;
            rooVarRun = iasSystRun;
            rooVarEvent = iasSystEvent;
            rooDataSetIasShift->add(RooArgSet(rooVarIas,rooVarIh,rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarLumiSection,rooVarRun,rooVarEvent));
            rooVarIas = tempIas; // reset to original value
            rooVarIh = tempIh;
            rooVarP = tempP;
            rooVarPt = tempPt;
            rooVarNoMias = tempNoMias;
            rooVarEta = tempEta;
            rooVarLumiSection = tempLumiSection;
            rooVarRun = tempRun;
            rooVarEvent = tempEvent;
          }

        }

        if(numTracksPassingPreselectionThisEvent > 0)
          numEventsWithOneTrackPassingPreselection++;
        if(numTracksPassingIhCutThisEvent > 0)
          numEventsWithOneTrackOverIhCut++;
        if(numTracksPassingIasCutThisEvent > 0)
          numEventsWithOneTrackOverIasCut++;
        if(numTracksPassingIhAndIasCutsThisEvent > 0)
          numEventsWithOneTrackOverIhAndIasCuts++;

      numTracksPassingPreselARegionPerEventHist->Fill(numTracksPassingPreselectionARegionThisEvent);
      numTracksPassingPreselBRegionPerEventHist->Fill(numTracksPassingPreselectionBRegionThisEvent);
      numTracksPassingPreselCRegionPerEventHist->Fill(numTracksPassingPreselectionCRegionThisEvent);

      } // loop over events

      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  // end of file/event loop
  afterPreselectionPlots.pDistributionHist->Scale(scaleFactor_);

  rooVarNumGenHSCPEvents = numGenHSCPEvents;
  rooVarNumGenHSCPTracks = numGenHSCPTracks;
  rooVarNumGenChargedHSCPTracks = numGenChargedHSCPTracks;
  rooVarSignalEventCrossSection = signalEventCrossSection_;
  rooVarNumGenHSCPEventsPUReweighted = numPUReweightedMCEvents;
  rooVarEventWeightSum = eventWeightSum;
  rooDataSetGenHSCPTracks->add(
      RooArgSet(rooVarNumGenHSCPEvents,rooVarNumGenHSCPTracks,rooVarNumGenChargedHSCPTracks,
        rooVarSignalEventCrossSection,rooVarNumGenHSCPEventsPUReweighted,rooVarEventWeightSum,
        rooVarSampleWeight));

  std::cout << "event weight sum = " << eventWeightSum << std::endl;
  std::cout << "Found: " << numEventsPassingTrigger << " events passing trigger selection." <<
    std::endl;
  std::cout << "Found: " << numEventsWithOneTrackPassingPreselection <<
    " events with at least one track passing preselection." << std::endl;
  std::cout << "Found: " << numEventsWithOneTrackOverIhCut <<
    " events with at least one track passing ih > 3.5" << std::endl;
  std::cout << "Found: " << numEventsWithOneTrackOverIasCut <<
    " events with at least one track passing ias > 0.35" << std::endl;
  std::cout << "Found: " << numEventsWithOneTrackOverIhAndIasCuts <<
    " events with at least one track passing ih > 3.5 and ias > 0.35" << std::endl;
  std::cout << std::endl;
  std::cout << "Found: " << numTracksPassingTrigger << " tracks passing trigger." <<
    std::endl;
  std::cout << "Found: " << numTracksPassingPreselection << " tracks passing preselections." <<
    std::endl;
  std::cout << "Added " << rooDataSetCandidates->numEntries() << " tracks into the dataset." << std::endl;


  double totalTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(1);

  if(isMC_)
  {
    totalTracks = numGenHSCPTracks;
    std::cout << std::endl;
    std::cout << "Number of initial MC events: " << numGenHSCPEvents
      << std::endl;

    std::cout << "fraction of events passing trigger selection: " << numEventsPassingTrigger/(float)numGenHSCPEvents
      << std::endl
      << std::endl << "tracks passing preselection/track passing trigger: " <<
      numTracksPassingPreselection/(float)numTracksPassingTrigger << std::endl;

    std::cout << "Found " << numTracksInSearchRegion << " tracks in search region ( P > " <<
      pSidebandThreshold_ << ", ih > " << ihSidebandThreshold_ << " )." << std::endl;

    std::cout << std::endl;

    std::cout << "Num events with 1 HSCP gen: " << numEventsWithOneHSCPGen << std::endl <<
      "\tnum events with 1 HSCP seen (Mu trigger): " << numEventsWithOneHSCPSeenMuonTrigger <<
      " eff = " << numEventsWithOneHSCPSeenMuonTrigger/(float)numEventsWithOneHSCPGen << std::endl <<
      "\tnum events with 1 HSCP seen (MET trigger): " << numEventsWithOneHSCPSeenMetTrigger <<
      " eff = " << numEventsWithOneHSCPSeenMetTrigger/(float)numEventsWithOneHSCPGen <<
      std::endl << "Num events with 2 HSCP gen: " << numEventsWithTwoHSCPGen << std::endl <<
      "\tnum events with 2 HSCP seen (Mu Trigger) : " << numEventsWithTwoHSCPSeenMuonTrigger <<
      " eff = " << numEventsWithTwoHSCPSeenMuonTrigger/(float)numEventsWithTwoHSCPGen << std::endl <<
      "\tnum events with 2 HSCP seen (MET Trigger) : " << numEventsWithTwoHSCPSeenMetTrigger <<
      " eff = " << numEventsWithTwoHSCPSeenMetTrigger/(float)numEventsWithTwoHSCPGen << std::endl;
  }
  double triggeredTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(1);
  double nomIhTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(2);
  double nohTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(3);
  double validFracTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(4);
  double pixelHitTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(5);
  double nomIasTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(6);
  double qualityMaskTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(7);
  double chi2NdfTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(8);
  double etaTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(9);
  double ptErrTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(10);
  double v3dTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(11);
  double trackIsoTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(12);
  double caloIsoTracks = beforePreselectionPlots.tracksVsPreselectionsHist->GetBinContent(13);
  std::cout << endl << endl;
  // preselection table
  std::cout << "Table of preselections" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << setw(10) << "Cut " << setw(10) << "Number" << setw(15) << "Percentage" << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  if(isMC_)
    std::cout << setw(10) << "Gen " << setw(10) << totalTracks << setw(15) << "100%" << std::endl;
  std::cout << setw(10) << "Trigger " << setw(10) << triggeredTracks << setw(15) << 100*triggeredTracks/totalTracks << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << setw(10) << "NoM Ih " << setw(10) << nomIhTracks << setw(15) << 100*nomIhTracks/triggeredTracks << std::endl;
  std::cout << setw(10) << "NoH " << setw(10) << nohTracks << setw(15) << 100*nohTracks/nomIhTracks << std::endl;
  std::cout << setw(10) << "validFrac " << setw(10) << validFracTracks << setw(15) << 100*validFracTracks/nohTracks << std::endl;
  std::cout << setw(10) << "pixelHit " << setw(10) << pixelHitTracks << setw(15) << 100*pixelHitTracks/validFracTracks << std::endl;
  std::cout << setw(10) << "NoM Ias " << setw(10) << nomIasTracks << setw(15) << 100*nomIasTracks/pixelHitTracks << std::endl;
  std::cout << setw(10) << "QualMask " << setw(10) << qualityMaskTracks << setw(15) << 100*qualityMaskTracks/nomIasTracks << std::endl;
  std::cout << setw(10) << "Chi2Ndf " << setw(10) << chi2NdfTracks << setw(15) << 100*chi2NdfTracks/qualityMaskTracks << std::endl;
  std::cout << setw(10) << "Eta " << setw(10) << etaTracks << setw(15) << 100*etaTracks/chi2NdfTracks << std::endl;
  std::cout << setw(10) << "PtErr " << setw(10) << ptErrTracks << setw(15) << 100*ptErrTracks/etaTracks << std::endl;
  std::cout << setw(10) << "v3D " << setw(10) << v3dTracks << setw(15) << 100*v3dTracks/ptErrTracks << std::endl;
  std::cout << setw(10) << "trackIso " << setw(10) << trackIsoTracks << setw(15) << 100*trackIsoTracks/v3dTracks << std::endl;
  std::cout << setw(10) << "caloIso " << setw(10) << caloIsoTracks << setw(15) << 100*caloIsoTracks/trackIsoTracks << std::endl;
  std::cout << "-------------------------------------------" << std::endl;
  std::cout << setw(10) << "All " << setw(10) << caloIsoTracks << setw(15) << 100*caloIsoTracks/triggeredTracks << std::endl;


}

