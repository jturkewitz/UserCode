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
#include "TRandom3.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooPlot.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

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
#include "DataFormats/Math/interface/deltaR.h"

#include "commonFunctions.h"

double emptyBinVal = 1e-25;
int numIasBins = 50;
int nomStep = 2;
float etaStep = 0.2;
int lastLowerNoM = 17;


std::string getHistTitleBeg(int lowerNom, float lowerEta, float pSB, float ihSB, float ptSB, float iasSB,
    bool usePt, bool useIas)
{
  std::string histTitle = "NoM ";
  histTitle+=intToString(lowerNom);
  histTitle+="-";
  if(lowerNom==lastLowerNoM)
    histTitle+="end";
  else if(lowerNom==5)
    histTitle+=intToString(lowerNom+2*nomStep-1);
  else
    histTitle+=intToString(lowerNom+nomStep-1);
  histTitle+=", ";
  histTitle+=floatToString(lowerEta);
  histTitle+=" < #eta < ";
  histTitle+=floatToString(lowerEta+etaStep);
  if(usePt)
  {
    histTitle+=", (pt < ";
    histTitle+=floatToString(ptSB);
  }
  else
  {
    histTitle+=", (p < ";
    histTitle+=floatToString(pSB);
  }
  if(useIas)
  {
    histTitle+=", Ias > ";
    histTitle+=floatToString(iasSB);
  }
  else
  {
    histTitle+=", Ih > ";
    histTitle+=floatToString(ihSB);
  }
  histTitle+=")";
  return histTitle;
}


// add the mass, p/pt, and ih SB thresholds to the output filename automatically
std::string generateFileNameEnd(double massCut, double pSideband, double ptSideband, bool usePt, double ihSideband)
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
  fileNameEnd+=".ih";
  // replace "." with "p" for the ih
  std::string ihSBStr = floatToString(ihSideband);
  ihSBStr.replace(ihSBStr.find("."),1,"p");
  fileNameEnd+=ihSBStr;

  fileNameEnd+=".root";
  return fileNameEnd;
}

void setBinLabelsEntriesHist(TH2F* hist)
{
  hist->GetXaxis()->SetBinLabel(1,"0.0-0.2");
  hist->GetXaxis()->SetBinLabel(2,"0.2-0.4");
  hist->GetXaxis()->SetBinLabel(3,"0.4-0.6");
  hist->GetXaxis()->SetBinLabel(4,"0.6-0.8");
  hist->GetXaxis()->SetBinLabel(5,"0.8-1.0");
  hist->GetXaxis()->SetBinLabel(6,"1.0-1.2");
  hist->GetXaxis()->SetBinLabel(7,"1.2-1.4");
  hist->GetXaxis()->SetBinLabel(8,"1.4-1.6");
  hist->GetXaxis()->SetBinLabel(9,"1.6-1.8");
  hist->GetXaxis()->SetBinLabel(10,"1.8-2.0");
  hist->GetXaxis()->SetBinLabel(11,"2.0-2.2");
  hist->GetXaxis()->SetBinLabel(12,"2.0-2.4");
  hist->GetYaxis()->SetBinLabel(1,"5-6");
  hist->GetYaxis()->SetBinLabel(2,"7-8");
  hist->GetYaxis()->SetBinLabel(3,"9-10");
  hist->GetYaxis()->SetBinLabel(4,"11-12");
  hist->GetYaxis()->SetBinLabel(5,"13-14");
  hist->GetYaxis()->SetBinLabel(6,"15-16");
  hist->GetYaxis()->SetBinLabel(7,"17-18");
  hist->GetYaxis()->SetBinLabel(8,"19-20");
  hist->GetYaxis()->SetBinLabel(9,"21-22");
}

//std::vector<double> computeVariableBins(TH1* hist)
//{
//  double binWidth = 0.02;
//  std::vector<double> binArray;
//
//  int lastBinWithContent = hist->GetNbinsX();
//  while(hist->GetBinContent(lastBinWithContent)<=0 && lastBinWithContent>0)
//  {
//    //std::cout << "BinContent for bin: " << lastBinWithContent << " = " <<
//    //  hist->GetBinContent(lastBinWithContent) << std::endl;
//    lastBinWithContent--;
//  }
//
//  if(lastBinWithContent==0)
//    lastBinWithContent=hist->GetNbinsX();
//
//  double xCoordOfLastBinWithContent = hist->GetBinLowEdge(lastBinWithContent);
//  for(double x = hist->GetXaxis()->GetXmin(); x < xCoordOfLastBinWithContent; x+=binWidth)
//    binArray.push_back(x);
//
//  if(xCoordOfLastBinWithContent-binArray.back() >= binWidth)
//    binArray.push_back(xCoordOfLastBinWithContent);
//  binArray.push_back(hist->GetXaxis()->GetXmax());
//
//  return binArray;
//}

// compute variable bin array from fixed bin hist
std::vector<double> computeVariableBins(TH1* hist)
{
  std::vector<double> binArray;
  binArray.push_back(hist->GetXaxis()->GetXmin());

  int binsExamined = 0;
  int binsToCombine = 2;
  while(binsExamined < hist->GetNbinsX())
  {
    int lastBin = hist->FindBin(binArray.back());
    int nextBin = lastBin+binsToCombine-1;

    double proposedBinContent = hist->Integral(lastBin,nextBin);
    // don't include that bin in the integral, but use its lower edge for the bin array
    nextBin++;

    ////debug
    //std::cout << "binsExamined: " << binsExamined << " lastBin: " << lastBin
    //  << " nextBin: " << nextBin
    //  << " lastBinLEdge: " << hist->GetBinLowEdge(lastBin)
    //  << " nextBin LEdge: " << hist->GetBinLowEdge(nextBin)
    //  << " proposedBinContent: " << proposedBinContent << " binsToCombine: " << binsToCombine
    //  << std::endl;
     
    binsExamined=nextBin;

    if(proposedBinContent > 0)
      binArray.push_back(hist->GetBinLowEdge(nextBin));
    else
      binsToCombine+=2;
  }

  //// debug
  //for(unsigned int i=0; i<binArray.size(); ++i)
  //  std::cout << "initial binArray[" << i << "]: " << binArray[i] << std::endl;

  while(hist->Integral(hist->FindBin(binArray.back()),hist->FindBin(hist->GetXaxis()->GetXmax())) <= 0
      && binArray.size() > 1)
  {
    binArray.erase(binArray.end()-1);
  }

  binArray.push_back(hist->GetXaxis()->GetXmax());

  //// debug
  //for(unsigned int i=0; i<binArray.size(); ++i)
  //  std::cout << "final binArray[" << i << "]: " << binArray[i] << std::endl;

  return binArray;
}

//// compute variable bin array from fixed bin hist
//std::vector<double> computeVariableBins(TH1* hist)
//{
//  std::vector<double> binArray;
//  binArray.push_back(hist->GetXaxis()->GetXmax());
//
//  // loop over x bins from right to left
//  for(int bin=hist->GetNbinsX(); bin > 0; --bin)
//  {
//    if(hist->GetBinContent(bin) > 0)
//      binArray.push_back(hist->GetBinLowEdge(bin));
//  }
//
//  binArray.erase(binArray.end()-1);
//  binArray.push_back(hist->GetXaxis()->GetXmin());
//
//  std::reverse(binArray.begin(),binArray.end());
//  return binArray;
//}

// ****** main
int main(int argc, char ** argv)
{

  // parse arguments
  if(argc < 2)
  {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return -1;
  }

  using namespace RooFit;

  // get the python configuration
  const edm::ParameterSet& process = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process");
  fwlite::InputSource inputHandler_(process);
  fwlite::OutputFiles outputHandler_(process);
  //int maxEvents_(inputHandler_.maxEvents());
  // now get each parameterset
  const edm::ParameterSet& ana = process.getParameter<edm::ParameterSet>("makeIasPredictions");
  // mass cut to use for the high-p high-Ih (search region) ias dist
  double massCutIasHighPHighIh_ (ana.getParameter<double>("MassCut"));
  // dE/dx calibration
  double dEdx_k (ana.getParameter<double>("dEdx_k"));
  double dEdx_c (ana.getParameter<double>("dEdx_c"));
  // definition of sidebands/search region
  double pSidebandThreshold (ana.getParameter<double>("PSidebandThreshold"));
  double ptSidebandThreshold (ana.getParameter<double>("PtSidebandThreshold"));
  bool usePtForSideband (ana.getParameter<bool>("UsePtForSideband"));
  double ihSidebandThreshold (ana.getParameter<double>("IhSidebandThreshold"));
  double iasSidebandThreshold (ana.getParameter<double>("IasSidebandThreshold"));
  bool useIasForSideband (ana.getParameter<bool>("UseIasForSideband"));
  double etaMin_ (ana.getParameter<double>("EtaMin"));
  double etaMax_ (ana.getParameter<double>("EtaMax"));
  int nomMin_ (ana.getParameter<int>("NoMMin"));
  int nomMax_ (ana.getParameter<int>("NoMMax"));

  // <ih> as a function of ias from MinBias data
  // TODO update for 8 TeV data
  TF1* avgIhFromIasFunc = new TF1("avgIhFromIasFunc","pol3(0)",0,1);
  avgIhFromIasFunc->SetParameters(2.782,7.069,-11.93,12.83);

  // fileService initialization
  //string fileNameEnd = generateFileNameEnd(massCutIasHighPHighIh_,pSidebandThreshold,ptSidebandThreshold,usePtForSideband,ihSidebandThreshold);
  string fileNameEnd = "";
  fwlite::TFileService fs = fwlite::TFileService((outputHandler_.file()+fileNameEnd).c_str());
  TFileDirectory successRateDir = fs.mkdir("successRates");
  TFileDirectory iasPredictionFixedBinsDir = fs.mkdir("iasPredictionsFixedBins");
  TFileDirectory dRegionFixedBinsDir = fs.mkdir("dRegionFixedBins");
  TFileDirectory massPredictionFixedBinsDir = fs.mkdir("massPredictionsFixedBins");
  TFileDirectory iasPredictionVarBinsDir = fs.mkdir("iasPredictionsVariableBins");
  TFileDirectory aRegionDir = fs.mkdir("aRegionHistograms");
  TFileDirectory bRegionDir = fs.mkdir("bRegionHistograms");
  TFileDirectory cRegionDir = fs.mkdir("cRegionHistograms");
  TFileDirectory ihMeanDir = fs.mkdir("ihMeanInIasBins");
  TFileDirectory minPCutMeanDir = fs.mkdir("minPCutInIasBins");
  TFileDirectory cRegionTracksOverMassCutDir = fs.mkdir("cRegionTracksOverIasCutInIasBins");
  // b region hist
  //TH1F* bRegionHist = fs.make<TH1F>("bRegionHist","bRegionHist",100,0,10);
  //bRegionHist->SetName("ihLowPsb_B");
  //bRegionHist->SetTitle("Ih in low P SB (B region);MeV/cm");
  // c region hist
  //TH1F* cRegionHist = fs.make<TH1F>("cRegionHist","cRegionHist",100,0,1000);
  //cRegionHist->SetName("pLowIhSB_C");
  //cRegionHist->SetTitle("P in low Ih SB (C region);GeV");
  // number of entries in datasets histos
  TH2F*  entriesInARegionHist = fs.make<TH2F>("entriesInARegion","Entries in A region (low P, low Ih);#eta;nom",24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  int divisions = 500+int(13/nomStep)+1;
  entriesInARegionHist->GetYaxis()->SetNdivisions(divisions,false);
  TH2F*  entriesInBRegionHist = fs.make<TH2F>("entriesInBRegion","Entries in B region (Ih in low P SB);#eta;nom",24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  entriesInBRegionHist->GetYaxis()->SetNdivisions(divisions,false);
  TH2F*  entriesInCRegionHist = fs.make<TH2F>("entriesInCRegion","Entries in C region (P in low Ih SB);#eta;nom",24/int(etaStep*10),0,2.4,int(13/nomStep)+1,5,18);
  entriesInCRegionHist->GetYaxis()->SetNdivisions(divisions,false);
  //TH2F*  entriesInDRegionHist = fs.make<TH2F>("entriesInDRegion","Entries in D region (search--> high Ih, high P);#eta;nom",12,0,2.4,9,5,23);
  //entriesInDRegionHist->GetYaxis()->SetNdivisions(509,false);
  // ias total distribution from B region
  std::string iasBRegionTotalHistName ="iasBRegionTotalHist";
  std::string iasBRegionTotalHistTitle = "Ias in B region";
  TH1F* iasBRegionTotalHist = fs.make<TH1F>(iasBRegionTotalHistName.c_str(),iasBRegionTotalHistTitle.c_str(),100,0,1);
  // ih total distribution from B region
  std::string ihBRegionTotalHistName ="ihBRegionTotalHist";
  std::string ihBRegionTotalHistTitle = "Ih in B region";
  TH1F* ihBRegionTotalHist = fs.make<TH1F>(ihBRegionTotalHistName.c_str(),ihBRegionTotalHistTitle.c_str(),200,3,15);
  // eta distribution from B region
  std::string etaBRegionHistName ="etaBRegionHist";
  std::string etaBRegionHistTitle = "Eta in B region";
  TH1F* etaBRegionHist = fs.make<TH1F>(etaBRegionHistName.c_str(),etaBRegionHistTitle.c_str(),50,-3,3);
  // eta distribution from A region
  std::string etaARegionHistName ="etaARegionHist";
  std::string etaARegionHistTitle = "Eta in A region";
  TH1F* etaARegionHist = fs.make<TH1F>(etaARegionHistName.c_str(),etaARegionHistTitle.c_str(),50,-3,3);
  // c region Pt vs. eta
  std::string pEtaCRegionHistName = "pEtaCRegionHist";
  std::string pEtaCRegionHistTitle = "P vs. eta in C Region";
  TH2F* pEtaCRegionHist = fs.make<TH2F>(pEtaCRegionHistName.c_str(),pEtaCRegionHistTitle.c_str(),50,-3,3,200,45,1200);

  //TODO move to common file
  // RooFit observables and dataset
//  RooRealVar rooVarIas("rooVarIas","ias",0,1);
//  RooRealVar rooVarIh("rooVarIh","ih",0,15);
//  RooRealVar rooVarP("rooVarP","p",0,5000);
//  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
//  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
//  RooRealVar rooVarEta("rooVarEta","eta",-2.5,2.5);
  RooRealVar ias_("ias_","ias",0,1);
  RooRealVar ih_("ih_","ih",0,15);
  RooRealVar trackP_("trackP_","p",0,5000);
  RooRealVar trackPt_("trackPt_","pt",0,5000);
  RooRealVar iasNoM_("iasNoM_","nom",0,30);
  RooRealVar trackEta_("trackEta_","eta",-2.5,2.5);
  TFile* inFile = TFile::Open(inputHandler_.files()[0].c_str());
  if(!inFile)
  {
    std::cout << "ERROR: input file " << inputHandler_.files()[0].c_str() << " not found." << std::endl;
    return -2;
  }

//  RooDataSet* rooDataSetInput = (RooDataSet*)inFile->Get("rooDataSetCandidates");
  TTree *rooDataSetCandidatesTree = (TTree*)inFile->Get("rooDataSetCandidatesTree");
//  rooDataSetCandidatesTree->SetBranchStatus("*",0);
//  rooDataSetCandidatesTree->SetBranchStatus("ias_",1);
//  rooDataSetCandidatesTree->SetBranchStatus("ih_",1);
//  rooDataSetCandidatesTree->SetBranchStatus("trackP_",1);
//  rooDataSetCandidatesTree->SetBranchStatus("trackPt_",1);
//  rooDataSetCandidatesTree->SetBranchStatus("iasNoM_",1);
//  rooDataSetCandidatesTree->SetBranchStatus("trackEta_",1);
//  TFile *newfile = new TFile("small.root","recreate");
//  TTree *rooDataSetCandidatesTreeReduced = rooDataSetCandidatesTree->CloneTree();
  RooDataSet* rooDataSetInput = new RooDataSet("rooDataSetCandidates","rooDataSetCandidates",rooDataSetCandidatesTree,
            RooArgSet(ias_,ih_,trackP_,trackPt_,iasNoM_,trackEta_));

  //RooDataSet* rooDataSetInput = (RooDataSet*)inFile->Get("rooDataSetOneCandidatePerEvent");
  if(rooDataSetInput->numEntries() < 1)
  {
    std::cout << "Problem with RooDataSet named rooDataSetCandidates in file " << inputHandler_.files()[0].c_str() << std::endl;
    return -3;
  }

  std::cout << "Doing ias predictions..." << std::endl;
  std::string pSBcutString = "trackP_<";
  pSBcutString+=floatToString(pSidebandThreshold);
  std::string pSearchCutString = "trackP_>";
  pSearchCutString+=floatToString(pSidebandThreshold);
  std::string ptSBcutString = "trackPt_<";
  ptSBcutString+=floatToString(ptSidebandThreshold);
  std::string ptSearchCutString = "trackPt_>";
  ptSearchCutString+=floatToString(ptSidebandThreshold);
  std::string ihSBcutString = "ih_<";
  ihSBcutString+=floatToString(ihSidebandThreshold);
  std::string ihSearchCutString = "ih_>";
  ihSearchCutString+=floatToString(ihSidebandThreshold);
  std::string iasSBcutString = "ias_<";
  iasSBcutString+=floatToString(iasSidebandThreshold);
  std::string iasSearchCutString = "ias_>";
  iasSearchCutString+=floatToString(iasSidebandThreshold);
  // TESTING CUT OFF PT < 50
  //RooDataSet* rooDataSetPtCut = (RooDataSet*) rooDataSetInput->reduce(Cut("trackPt_>50"));
  //std::cout << "INFO: Applying minimum Pt cut of 50" << std::endl;
  RooDataSet* rooDataSetPtCut = (RooDataSet*) rooDataSetInput->Clone();
  //XXX TURN OFF ETA > 1.5
  RooDataSet* rooDataSetAll = (RooDataSet*) rooDataSetPtCut->reduce(Cut("trackEta_<1.5&&trackEta_>-1.5"));
  std::cout << "INFO: Applying upper eta cut of 1.5" << std::endl;
  delete rooDataSetInput;
  delete rooDataSetPtCut;
  // A ==> low P/Pt, low Ih
  RooDataSet* regionA1DataSet;
  RooDataSet* regionADataSet;
  if(usePtForSideband)
    regionA1DataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ptSBcutString.c_str()));
  else
    regionA1DataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(pSBcutString.c_str()));
  if(useIasForSideband)
    regionADataSet = (RooDataSet*)regionA1DataSet->reduce(Cut(iasSBcutString.c_str()),
        SelectVars(RooArgSet(trackP_,trackPt_,trackEta_,ih_,iasNoM_)));
  else
    regionADataSet = (RooDataSet*)regionA1DataSet->reduce(Cut(ihSBcutString.c_str()),
        SelectVars(RooArgSet(trackP_,trackPt_,trackEta_,ih_,iasNoM_)));
  delete regionA1DataSet;
  // B ==> low P/Pt, high Ih
  RooDataSet* regionB1DataSet;
  RooDataSet* regionBDataSet;
  if(usePtForSideband)
    regionB1DataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ptSBcutString.c_str()));
  else
    regionB1DataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(pSBcutString.c_str()));
  if(useIasForSideband)
    regionBDataSet = (RooDataSet*)regionB1DataSet->reduce(Cut(iasSearchCutString.c_str()),
          SelectVars(RooArgSet(trackPt_,ih_,ias_,trackEta_,iasNoM_)));
  else
    regionBDataSet = (RooDataSet*)regionB1DataSet->reduce(Cut(ihSearchCutString.c_str()),
          SelectVars(RooArgSet(trackPt_,ih_,ias_,trackEta_,iasNoM_)));
  delete regionB1DataSet;
  // C ==> high P/Pt, low Ih
  RooDataSet* dedxSBDataSet;
  RooDataSet* regionCDataSet;
  if(useIasForSideband)
    dedxSBDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(iasSBcutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackPt_,iasNoM_,trackEta_,ih_)));
  else
    dedxSBDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ihSBcutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackPt_,iasNoM_,trackEta_,ih_)));
  //SIC FEB 10 -- cut off p > 1200
  //RooDataSet* pCutCDataSet = (RooDataSet*)regionC1DataSet->reduce(Cut("trackP_<1200"));
  //std::cout << "INFO: Applying max p of 1200 GeV cut in C region." << std::endl;
  std::cout << "INFO: NOT applying max p of 1200 GeV cut in C region." << std::endl;
  if(usePtForSideband)
   regionCDataSet = (RooDataSet*)dedxSBDataSet->reduce(Cut(ptSearchCutString.c_str()));
    //regionCDataSet = (RooDataSet*)pCutCDataSet->reduce(Cut(iasSBcutString.c_str()));
  else
   regionCDataSet = (RooDataSet*)dedxSBDataSet->reduce(Cut(pSearchCutString.c_str()));
    //regionCDataSet = (RooDataSet*)pCutCDataSet->reduce(Cut(ihSBcutString.c_str()));

  std::cout << "INFO: " << regionCDataSet->numEntries() << " entries in total C region." << std::endl;
  // SIC APR 14 -- define D region using mass and ias cuts only
  // no p cut needed since eff p cut is done by mass cut
  // D ==> search region, high dE/dx; apply mass cut later
  RooDataSet* highDeDxDataSet;
  if(useIasForSideband)
    highDeDxDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(iasSearchCutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackEta_,iasNoM_,ias_,ih_)));
  else
    highDeDxDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ihSearchCutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackEta_,iasNoM_,ias_,ih_)));
  delete rooDataSetAll;

  // broken
  int slicesToDo = (int)((10*etaMax_-10*etaMin_)/2.0)*(nomMax_-nomMin_+1)/2.0;

  int nomSlice = 0;
  // loop over given nom/eta slices
  for(int nom = nomMin_; nom < nomMax_; nom+=nomStep)
  {
    nomSlice++;
    int etaSlice = 0;
    for(float lowerEta = etaMin_; lowerEta < etaMax_; lowerEta+=etaStep)
    {
      etaSlice++;

      // book histograms -- b region
      std::string bRegionHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      bRegionHistName+="bRegionHist";
      std::string bRegionHistTitle = "Ih for ";
      bRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      bRegionHistTitle+="; MeV/cm";
      TH1F* bRegionHist = bRegionDir.make<TH1F>(bRegionHistName.c_str(),bRegionHistTitle.c_str(),100,0,10);
      // ias distribution from B region
      std::string iasBRegionHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasBRegionHistName+="iasBRegionHist";
      std::string iasBRegionHistTitle = "Ias in B region for ";
      iasBRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      TH1F* iasBRegionHist = bRegionDir.make<TH1F>(iasBRegionHistName.c_str(),iasBRegionHistTitle.c_str(),numIasBins,0,1);
      // c region hist
      std::string cRegionHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      cRegionHistName+="cRegionHist";
      std::string cRegionHistTitle = "P for ";
      cRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      cRegionHistTitle+="; GeV";
      TH1F* cRegionHist = cRegionDir.make<TH1F>(cRegionHistName.c_str(),cRegionHistTitle.c_str(),100,0,1000);
      // c region hist - cumulative
      std::string cRegionCumuHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      cRegionCumuHistName+="cRegionCumuHist";
      std::string cRegionCumuHistTitle = "P for ";
      cRegionCumuHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      cRegionCumuHistTitle+=" cumulative; GeV";
      TH1F* cRegionCumuHist = cRegionDir.make<TH1F>(cRegionCumuHistName.c_str(),cRegionCumuHistTitle.c_str(),100,0,1000);
      // ceff region hist
      std::string ceffRegionHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      ceffRegionHistName+="ceffRegionHist";
      std::string ceffRegionHistTitle = "P for ";
      ceffRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      ceffRegionHistTitle+="; GeV";
      TH1F* ceffRegionHist = cRegionDir.make<TH1F>(ceffRegionHistName.c_str(),ceffRegionHistTitle.c_str(),100,0,1000);
      // ceff region over mass cut in ias bins hist
      std::string ceffRegionTracksOverMassCutProfileName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      ceffRegionTracksOverMassCutProfileName+="tracksInCeffOverMassCutProfile";
      std::string ceffRegionTracksOverMassCutProfileTitle = "Tracks in Ceff Over Mass Cut in Ias bins for ";
      ceffRegionTracksOverMassCutProfileTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      ceffRegionTracksOverMassCutProfileTitle+="; Ias";
      TProfile* ceffRegionTracksOverMassCutProfile =
        cRegionDir.make<TProfile>(ceffRegionTracksOverMassCutProfileName.c_str(),
            ceffRegionTracksOverMassCutProfileTitle.c_str(),numIasBins,0,1);
      // ceff region over mass cut in ias bins hist -- using avg Ih from minBias data
      std::string ceffRegionTracksOverMassCutAvgIhProfileName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      ceffRegionTracksOverMassCutAvgIhProfileName+="tracksInCeffOverMassCutAvgIhProfile";
      std::string ceffRegionTracksOverMassCutAvgIhProfileTitle = "Tracks in Ceff Over Mass Cut (AvgIh from MinBias data) in Ias bins for ";
      ceffRegionTracksOverMassCutAvgIhProfileTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      ceffRegionTracksOverMassCutAvgIhProfileTitle+="; Ias";
      TProfile* ceffRegionTracksOverMassCutAvgIhProfile =
        cRegionDir.make<TProfile>(ceffRegionTracksOverMassCutAvgIhProfileName.c_str(),
            ceffRegionTracksOverMassCutAvgIhProfileTitle.c_str(),numIasBins,0,1);

      // B region dataset
      std::string nomCutString = "iasNoM_>=";
      nomCutString+=intToString(nom);
      nomCutString+="&&iasNoM_<=";
      nomCutString+=intToString(nom+nomStep-1);
      if(nom==lastLowerNoM) // do nom lastLowerNoM+ in one slice
      {
        nomCutString = "iasNoM_>=";
        nomCutString+=intToString(lastLowerNoM);
      }
      if(nom==5) // combine 5-6 with next slice
      {
        nomCutString = "iasNoM_>=";
        nomCutString+=intToString(nom);
        nomCutString+="&&iasNoM_<=";
        nomCutString+=intToString(nom+2*nomStep-1);
      }
      std::string etaCutString = "(trackEta_>";
      etaCutString+=floatToString(lowerEta);
      etaCutString+="&&trackEta_<";
      std::string upperEtaLimit;
      upperEtaLimit=floatToString(lowerEta+etaStep);
      etaCutString+=upperEtaLimit;
      etaCutString+=")||(trackEta_>-";
      etaCutString+=upperEtaLimit;
      etaCutString+="&&trackEta_<-";
      etaCutString+=floatToString(lowerEta);
      etaCutString+=")";
      RooDataSet* nomCutBRegionDataSet = (RooDataSet*)regionBDataSet->reduce(Cut(nomCutString.c_str()),
          SelectVars(RooArgSet(trackPt_,ih_,ias_,trackEta_)));
      RooDataSet* etaCutNomCutBRegionDataSet = (RooDataSet*)nomCutBRegionDataSet->reduce(Cut(etaCutString.c_str()));
      // C region dataset
      // SIC Dec. 23: study shows we can use all NoM slices
      //RooDataSet* nomCutCRegionDataSet = (RooDataSet*)regionCDataSet->reduce(Cut(nomCutString.c_str()));
      RooDataSet* etaCutCRegionDataSet = (RooDataSet*)regionCDataSet->reduce(Cut(etaCutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackPt_,iasNoM_,trackEta_,ih_)));
      RooDataSet* etaCutNomCutCRegionDataSet = (RooDataSet*)etaCutCRegionDataSet->reduce(Cut(nomCutString.c_str()));
      // dE/dx SB data set cut in eta --> used for Ceff region later
      RooDataSet* etaCutDeDxSBDataSet = (RooDataSet*)dedxSBDataSet->reduce(Cut(etaCutString.c_str()),
          SelectVars(RooArgSet(trackP_)));
      // High dE/dx dataset
      RooDataSet* nomCutHighDeDxDataSet = (RooDataSet*)highDeDxDataSet->reduce(Cut(nomCutString.c_str()),
          SelectVars(RooArgSet(trackP_,trackEta_,ih_,ias_)));
      RooDataSet* etaCutNomCutHighDeDxDataSet = (RooDataSet*)nomCutHighDeDxDataSet->reduce(Cut(etaCutString.c_str()));
      // A region dataset
      RooDataSet* etaCutARegionDataSet = (RooDataSet*)regionADataSet->reduce(Cut(etaCutString.c_str()),
          SelectVars(RooArgSet(trackPt_,trackEta_,ih_,iasNoM_)));
      RooDataSet* etaCutNomCutARegionDataSet = (RooDataSet*)etaCutARegionDataSet->reduce(Cut(nomCutString.c_str()));

      entriesInARegionHist->Fill(lowerEta+0.01,nom,etaCutNomCutARegionDataSet->numEntries());
      entriesInBRegionHist->Fill(lowerEta+0.01,nom,etaCutNomCutBRegionDataSet->numEntries());
      entriesInCRegionHist->Fill(lowerEta+0.01,nom,etaCutNomCutCRegionDataSet->numEntries());
      //entriesInDRegionHist->Fill(lowerEta+0.1,nom,etaCutNomCutDRegionDataSet->numEntries());
      //debug/output
      std::cout << "eta=" <<
        lowerEta << "-";
      std::cout << lowerEta+etaStep;
      if(nom==lastLowerNoM)
        std::cout << " nom=" << nom << "+" << std::endl;
      else if(nom==5)
        std::cout << " nom=5" << "-" << nom+2*nomStep-1 << std::endl;
      else
        std::cout << " nom=" << nom << "-" << nom+nomStep-1 << std::endl;
      std::cout << "A = " << etaCutNomCutARegionDataSet->numEntries() <<
        " B = " << etaCutNomCutBRegionDataSet->numEntries() <<
        " C = " << etaCutCRegionDataSet->numEntries() << std::endl;
       // " D = " << etaCutNomCutDRegionDataSet->numEntries() << std::endl;

      const RooArgSet* argSet_A = etaCutNomCutARegionDataSet->get();
      //RooRealVar* iasData_A = (RooRealVar*)argSet_A->find(ias_.GetName());
      RooRealVar* ihData_A = (RooRealVar*)argSet_A->find(ih_.GetName());
      RooRealVar* etaData_A = (RooRealVar*)argSet_A->find(trackEta_.GetName());
      RooRealVar* ptData_A = (RooRealVar*)argSet_A->find(trackPt_.GetName());
      for(int index=0; index < etaCutNomCutARegionDataSet->numEntries(); ++index)
      {
        etaCutNomCutARegionDataSet->get(index);
        if(ihData_A->getVal() > 3 && ptData_A->getVal() > 45)
          etaARegionHist->Fill(etaData_A->getVal());
      }
      //delete argSet_A;
      //delete iasData_A;
      //delete ihData_A;
      //delete etaData_A;
      //delete ptData_A;
      const RooArgSet* argSet_C1 = etaCutNomCutCRegionDataSet->get();
      RooRealVar* pData_C1 = (RooRealVar*)argSet_C1->find(trackP_.GetName());
      RooRealVar* ptData_C1 = (RooRealVar*)argSet_C1->find(trackPt_.GetName());
      RooRealVar* ihData_C1 = (RooRealVar*)argSet_C1->find(ih_.GetName());
      RooRealVar* etaData_C1 = (RooRealVar*)argSet_C1->find(trackEta_.GetName());
      for(int index=0; index < etaCutNomCutCRegionDataSet->numEntries(); ++index)
      {
        etaCutNomCutCRegionDataSet->get(index);
        //std::cout << "Fill hist: " << " eta = " << etaData_C1->getVal() << std::endl;
        //std::cout << " pt = " << ptData_C1->getVal() << std::endl;
        //std::cout << "Entries in hist: " << pEtaCRegionHist->GetEntries() << std::endl;
        if(ihData_C1->getVal() > 3 && ptData_C1->getVal() > 45)
          pEtaCRegionHist->Fill(etaData_C1->getVal(),pData_C1->getVal());
      }
      //delete argSet_C1;
      //delete pData_C1;
      //delete ptData_C1;
      //delete ihData_C1;
      //delete etaData_C1;

      // immediately delete no-longer-needed datasets
      //delete nomCutDRegionDataSet;
      //delete etaCutNomCutDRegionDataSet;
      delete etaCutNomCutARegionDataSet;
      delete etaCutNomCutCRegionDataSet;
      delete nomCutBRegionDataSet;

      const RooArgSet* argSet_B = etaCutNomCutBRegionDataSet->get();
      RooRealVar* ihData_B = (RooRealVar*)argSet_B->find(ih_.GetName());
      RooRealVar* iasData_B = (RooRealVar*)argSet_B->find(ias_.GetName());
      RooRealVar* etaData_B = (RooRealVar*)argSet_B->find(trackEta_.GetName());
      RooRealVar* ptData_B = (RooRealVar*)argSet_B->find(trackPt_.GetName());
      //remove eta/nom double checking
      //RooRealVar* nomData_B = (RooRealVar*)argSet_B->find(iasNoM_.GetName());
      // fill b region hist
      for(int index=0; index < etaCutNomCutBRegionDataSet->numEntries(); ++index)
      {
        etaCutNomCutBRegionDataSet->get(index);
        if(useIasForSideband)
        {
          if(iasData_B->getVal() < iasSidebandThreshold)
            std::cout << "ERROR: (in B region) ias=" << iasData_B->getVal() << std::endl;
        }
        else
        {
          if(ihData_B->getVal() < ihSidebandThreshold)
            std::cout << "ERROR: (in B region) ih=" << ihData_B->getVal() << std::endl;
        }
        //if(nomData_B->getVal() < nom || (nomData_B->getVal() > nom+1 && nom<21))
        //  std::cout << "ERROR: nom=" << nom << "-" << nom+1 << " and (in B region) data point nom="
        //    << nomData_B->getVal() << std::endl;
        //if(etaData_B->getVal() < lowerEta || etaData_B->getVal() > lowerEta+0.2)
        //  std::cout << "ERROR: eta=" << lowerEta << "-" << lowerEta+0.2 << " and (in B region) data point eta="
        //    << etaData_B->getVal() << std::endl;

        bRegionHist->Fill(ihData_B->getVal());
        iasBRegionHist->Fill(iasData_B->getVal());
        if(ihData_B->getVal() > 3 && ptData_B->getVal() > 45)
        {
          etaBRegionHist->Fill(etaData_B->getVal());
          iasBRegionTotalHist->Fill(iasData_B->getVal());
          ihBRegionTotalHist->Fill(ihData_B->getVal());
        }
        ////debug
        //std::cout << "Filled B region hist: " << bRegionHist->GetName()
        //  << " index = " << index << std::endl;

      }
      const RooArgSet* argSet_C = etaCutCRegionDataSet->get();
      RooRealVar* pData_C = (RooRealVar*)argSet_C->find(trackP_.GetName());
      RooRealVar* ptData_C = (RooRealVar*)argSet_C->find(trackPt_.GetName());
      //RooRealVar* etaData_C = (RooRealVar*)argSet_C->find(trackEta_.GetName());
      //RooRealVar* nomData_C = (RooRealVar*)argSet_C->find(iasNoM_.GetName());
      // fill c region hist
      for(int index=0; index < etaCutCRegionDataSet->numEntries(); ++index)
      {
        etaCutCRegionDataSet->get(index);
        cRegionHist->Fill(pData_C->getVal());

        if(usePtForSideband)
        {
          if(ptData_C->getVal() < ptSidebandThreshold)
            std::cout << "ERROR: (in C region) pt=" << ptData_C->getVal() << std::endl;
        }
        else
        {
          if(pData_C->getVal() < pSidebandThreshold)
            std::cout << "ERROR: (in C region) p=" << pData_C->getVal() << std::endl;
        }
        // no NoM checking in c region; using eta cut only
        //if(nomData_C->getVal() < nom || nomData_C->getVal() > nom+1)
        //  std::cout << "ERROR: nom=" << nom << "-" << nom+1 << " and (in C region) data point nom="
        //    << nomData_C->getVal() << std::endl;
        //if(etaData_C->getVal() < lowerEta || etaData_C->getVal() > lowerEta+0.2)
        //  std::cout << "ERROR: eta=" << lowerEta << "-" << lowerEta+0.2 << " and (in C region) data point eta="
        //    << etaData_C->getVal() << std::endl;
      }
      // fill ceff region hist
      const RooArgSet* argSet_Ceff = etaCutDeDxSBDataSet->get();
      RooRealVar* pData_Ceff = (RooRealVar*)argSet_Ceff->find(trackP_.GetName());
      for(int index=0; index < etaCutDeDxSBDataSet->numEntries(); ++index)
      {
        etaCutDeDxSBDataSet->get(index);
        ceffRegionHist->Fill(pData_Ceff->getVal());
      }

      // c region cumu hist
      cRegionHist->ComputeIntegral();
      Double_t *cRegionIntegral = cRegionHist->GetIntegral();
      cRegionCumuHist->SetContent(cRegionIntegral);

      // ias prediction histogram in this NoM/eta bin for limits -- underestimation of background
      std::string iasPredictionFixedLimitsHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedLimitsHistName+="iasPredictionFixedLimitsHist";
      std::string iasPredictionFixedLimitsHistTitle = "Ias prediction (limits) for ";
      iasPredictionFixedLimitsHistTitle+=intToString(nom);
      iasPredictionFixedLimitsHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedLimitsHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedLimitsHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedLimitsHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedLimitsHistTitle+=", ";
      iasPredictionFixedLimitsHistTitle+=floatToString(lowerEta);
      iasPredictionFixedLimitsHistTitle+=" < |#eta| < ";
      iasPredictionFixedLimitsHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedLimitsHistTitle+=", mass > ";
      iasPredictionFixedLimitsHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedLimitsHistTitle+=" GeV";
      TH1D* iasPredictionFixedLimitsHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedLimitsHistName.c_str(),iasPredictionFixedLimitsHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedLimitsHist->Sumw2();
      // ias prediction histogram for discovery  -- overestimation of background
      std::string iasPredictionFixedDiscoveryHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedDiscoveryHistName+="iasPredictionFixedDiscoveryHist";
      std::string iasPredictionFixedDiscoveryHistTitle = "Ias prediction (discovery) for ";
      iasPredictionFixedDiscoveryHistTitle+=intToString(nom);
      iasPredictionFixedDiscoveryHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedDiscoveryHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedDiscoveryHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedDiscoveryHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedDiscoveryHistTitle+=", ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(lowerEta);
      iasPredictionFixedDiscoveryHistTitle+=" < |#eta| < ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedDiscoveryHistTitle+=", mass > ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedDiscoveryHistTitle+=" GeV";
      TH1D* iasPredictionFixedDiscoveryHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedDiscoveryHistName.c_str(),iasPredictionFixedDiscoveryHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedDiscoveryHist->Sumw2();
      // ias prediction histogram in this NoM/eta bin -- add one sigma of slope in exp fit
      std::string iasPredictionFixedSlopeAddOneSigmaLimitsHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistName+="iasPredictionFixedSlopeAddOneSigmaLimitsHist";
      std::string iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle = "Ias prediction limits (+1 #sigma slope) for ";
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=intToString(nom);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=", ";
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=floatToString(lowerEta);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=" < |#eta| < ";
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=", mass > ";
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle+=" GeV";
      TH1D* iasPredictionFixedSlopeAddOneSigmaLimitsHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedSlopeAddOneSigmaLimitsHistName.c_str(),iasPredictionFixedSlopeAddOneSigmaLimitsHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedSlopeAddOneSigmaLimitsHist->Sumw2();
      // ias prediction histogram in this NoM/eta bin -- add one sigma of slope in exp fit
      std::string iasPredictionFixedSlopeMinusOneSigmaLimitsHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistName+="iasPredictionFixedSlopeMinusOneSigmaLimitsHist";
      std::string iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle = "Ias prediction limits (-1 #sigma slope) for ";
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=intToString(nom);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=", ";
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=floatToString(lowerEta);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=" < |#eta| < ";
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=", mass > ";
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle+=" GeV";
      TH1D* iasPredictionFixedSlopeMinusOneSigmaLimitsHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedSlopeMinusOneSigmaLimitsHistName.c_str(),iasPredictionFixedSlopeMinusOneSigmaLimitsHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedSlopeMinusOneSigmaLimitsHist->Sumw2();
      // ias prediction histogram in this NoM/eta bin -- add one sigma of slope in exp fit
      std::string iasPredictionFixedSlopeAddOneSigmaDiscoveryHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistName+="iasPredictionFixedSlopeAddOneSigmaDiscoveryHist";
      std::string iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle = "Ias prediction discovery (+1 #sigma slope) for ";
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=intToString(nom);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=", ";
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=floatToString(lowerEta);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=" < |#eta| < ";
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=", mass > ";
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle+=" GeV";
      TH1D* iasPredictionFixedSlopeAddOneSigmaDiscoveryHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedSlopeAddOneSigmaDiscoveryHistName.c_str(),iasPredictionFixedSlopeAddOneSigmaDiscoveryHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->Sumw2();
      // ias prediction histogram in this NoM/eta bin -- add one sigma of slope in exp fit
      std::string iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistName+="iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist";
      std::string iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle = "Ias prediction discovery (-1 #sigma slope) for ";
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=intToString(nom);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+="end";
      else if(nom==5)
        iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=intToString(nom+nomStep-1);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=", ";
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=floatToString(lowerEta);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=" < |#eta| < ";
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=floatToString(lowerEta+etaStep);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=", mass > ";
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle+=" GeV";
      TH1D* iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistName.c_str(),iasPredictionFixedSlopeMinusOneSigmaDiscoveryHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->Sumw2();
      // ias histogram points and fit in this NoM/eta bin
      std::string iasPointsAndFitHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      iasPointsAndFitHistName+="iasPointsAndFitHist";
      std::string iasPointsAndFitHistTitle = "Ias prediction for ";
      iasPointsAndFitHistTitle+=intToString(nom);
      iasPointsAndFitHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPointsAndFitHistTitle+="end";
      else if(nom==5)
        iasPointsAndFitHistTitle+=intToString(nom+2*nomStep-1);
      else
        iasPointsAndFitHistTitle+=intToString(nom+nomStep-1);
      iasPointsAndFitHistTitle+=", ";
      iasPointsAndFitHistTitle+=floatToString(lowerEta);
      iasPointsAndFitHistTitle+=" < |#eta| < ";
      iasPointsAndFitHistTitle+=floatToString(lowerEta+etaStep);
      iasPointsAndFitHistTitle+=", mass > ";
      iasPointsAndFitHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPointsAndFitHistTitle+=" GeV";
      TH1D* iasPointsAndFitHist = iasPredictionFixedBinsDir.make<TH1D>(iasPointsAndFitHistName.c_str(),iasPointsAndFitHistTitle.c_str(),numIasBins,0,1);
      iasPointsAndFitHist->Sumw2();
      // mass prediction histogram in this NoM/eta bin
      std::string massPredictionFixedHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      massPredictionFixedHistName+="massPredictionFixedHist";
      std::string massPredictionFixedHistTitle = "mass prediction for ";
      massPredictionFixedHistTitle+=intToString(nom);
      massPredictionFixedHistTitle+="-";
      if(nom==lastLowerNoM)
        massPredictionFixedHistTitle+="end";
      else if(nom==5)
        massPredictionFixedHistTitle+=intToString(nom+2*nomStep-1);
      else
        massPredictionFixedHistTitle+=intToString(nom+nomStep-1);
      massPredictionFixedHistTitle+=", ";
      massPredictionFixedHistTitle+=floatToString(lowerEta);
      massPredictionFixedHistTitle+=" < |#eta| < ";
      massPredictionFixedHistTitle+=floatToString(lowerEta+etaStep);
      massPredictionFixedHistTitle+=", mass > ";
      massPredictionFixedHistTitle+=floatToString(massCutIasHighPHighIh_);
      massPredictionFixedHistTitle+=" GeV";
      TH1D* massPredictionFixedHist = massPredictionFixedBinsDir.make<TH1D>(massPredictionFixedHistName.c_str(),massPredictionFixedHistTitle.c_str(),1000,0,5000);
      massPredictionFixedHist->Sumw2();
      // ias D region histogram in this NoM/eta bin
      std::string dRegionFixedHistName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      dRegionFixedHistName+="dRegionFixedHist";
      std::string dRegionFixedHistTitle = "Ias data in D region for ";
      dRegionFixedHistTitle+=intToString(nom);
      dRegionFixedHistTitle+="-";
      if(nom==lastLowerNoM)
        dRegionFixedHistTitle+="end";
      else if(nom==5)
        dRegionFixedHistTitle+=intToString(nom+2*nomStep-1);
      else
        dRegionFixedHistTitle+=intToString(nom+nomStep-1);
      dRegionFixedHistTitle+=", ";
      dRegionFixedHistTitle+=floatToString(lowerEta);
      dRegionFixedHistTitle+=" < |#eta| < ";
      dRegionFixedHistTitle+=floatToString(lowerEta+etaStep);
      dRegionFixedHistTitle+=", mass > ";
      dRegionFixedHistTitle+=floatToString(massCutIasHighPHighIh_);
      dRegionFixedHistTitle+=" GeV";
      TH1D* dRegionFixedHist = dRegionFixedBinsDir.make<TH1D>(dRegionFixedHistName.c_str(),dRegionFixedHistTitle.c_str(),numIasBins,0,1);
      dRegionFixedHist->Sumw2();

      std::cout << "Slice: eta=" << lowerEta << "-";
      // eta 1.5 max -- applied above
      if(10*lowerEta==14)
        std::cout << "1.5";
      else
        std::cout << lowerEta+etaStep;
      std::cout << " nom=" << nom << "-";
      if(nom==lastLowerNoM)
        std::cout << "end" << std::endl;
      else if(nom==5)
        std::cout << nom+2*nomStep-1 << std::endl;
      else
        std::cout << nom+nomStep-1 << std::endl;

      std::cout << "\t\t" << etaCutNomCutBRegionDataSet->numEntries() <<
        " entries in B region and " << etaCutCRegionDataSet->numEntries() <<
        " entries in C region (summed over NoM) and " << etaCutARegionDataSet->numEntries() <<
        " entries in A region (summed over NoM)" << endl <<
        " C/A = " << etaCutCRegionDataSet->numEntries()/(double)etaCutARegionDataSet->numEntries() << std::endl;

      // fill D region hist
      const RooArgSet* argSet_highDeDx = etaCutNomCutHighDeDxDataSet->get();
      RooRealVar* ihData_highDeDx = (RooRealVar*)argSet_highDeDx->find(ih_.GetName());
      RooRealVar* iasData_highDeDx = (RooRealVar*)argSet_highDeDx->find(ias_.GetName());
      RooRealVar* pData_highDeDx = (RooRealVar*)argSet_highDeDx->find(trackP_.GetName());
      for(int index=0; index < etaCutNomCutHighDeDxDataSet->numEntries(); ++index)
      {
        etaCutNomCutHighDeDxDataSet->get(index);
        double thisIh = ihData_highDeDx->getVal();
        double thisIas = iasData_highDeDx->getVal();
        double thisP = pData_highDeDx->getVal();
        if(thisIas < 0) continue; // require ias >= 0
        double massSqr = (thisIh-dEdx_c)*pow(thisP,2)/dEdx_k;
        if(sqrt(massSqr) > massCutIasHighPHighIh_)
          dRegionFixedHist->Fill(thisIas);
      }
      delete etaCutNomCutHighDeDxDataSet;
      delete nomCutHighDeDxDataSet;

      // require at least 1 entry in A/C datasets, 5 in B dataset to do prediction
      //if(etaCutNomCutBRegionDataSet->numEntries() < 1 || etaCutCRegionDataSet->numEntries() < 1)
      if(etaCutNomCutBRegionDataSet->numEntries() < 5 || etaCutCRegionDataSet->numEntries() < 1
          || etaCutARegionDataSet->numEntries() < 1)
      {
        std::cout << "WARNING: too few entries-not doing prediction for this slice" << std::endl;
        // fill with "empty" bins
        int firstIasBin = iasPredictionFixedLimitsHist->FindBin(iasSidebandThreshold);
        for(int bin=firstIasBin; bin < iasPredictionFixedLimitsHist->GetNbinsX()+1; ++bin)
        {
          if(iasPredictionFixedLimitsHist->GetBinContent(bin) == 0)
          {
            iasPredictionFixedLimitsHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPredictionFixedDiscoveryHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPointsAndFitHist->SetBinContent(bin,emptyBinVal);
            iasPointsAndFitHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
            iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
            iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
          }
        }
        continue;
      }

      // Ih mean in Ias bins histogram in this NoM/eta bin
      std::string ihMeanProfName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      ihMeanProfName+="ihMeanProf";
      std::string ihMeanProfTitle = "Ih mean for nom ";
      ihMeanProfTitle+=intToString(nom);
      ihMeanProfTitle+="-";
      if(nom==lastLowerNoM)
        ihMeanProfTitle+="end";
      else if(nom==5)
        ihMeanProfTitle+=intToString(nom+2*nomStep-1);
      else
        ihMeanProfTitle+=intToString(nom+nomStep-1);
      ihMeanProfTitle+=", ";
      ihMeanProfTitle+=floatToString(lowerEta);
      ihMeanProfTitle+=" < |#eta| < ";
      ihMeanProfTitle+=floatToString(lowerEta+etaStep);
      ihMeanProfTitle+=", mass > ";
      ihMeanProfTitle+=floatToString(massCutIasHighPHighIh_);
      ihMeanProfTitle+=" GeV; Ias";
      TProfile* ihMeanProf = ihMeanDir.make<TProfile>(ihMeanProfName.c_str(),ihMeanProfTitle.c_str(),numIasBins,0,1);
      // Min P cut mean in Ias bins histogram in this NoM/eta bin
      std::string minPCutMeanProfName = getHistNameBeg(nom,lowerEta,nomStep,etaStep);
      minPCutMeanProfName+="minPCutMeanProf";
      std::string minPCutMeanProfTitle = "Min P cut for nom ";
      minPCutMeanProfTitle+=intToString(nom);
      minPCutMeanProfTitle+="-";
      if(nom==lastLowerNoM)
        minPCutMeanProfTitle+="end";
      else if(nom==5)
        minPCutMeanProfTitle+=intToString(nom+2*nomStep-1);
      else
        minPCutMeanProfTitle+=intToString(nom+nomStep-1);
      minPCutMeanProfTitle+=", ";
      minPCutMeanProfTitle+=floatToString(lowerEta);
      minPCutMeanProfTitle+=" < |#eta| < ";
      minPCutMeanProfTitle+=floatToString(lowerEta+etaStep);
      minPCutMeanProfTitle+=", mass > ";
      minPCutMeanProfTitle+=floatToString(massCutIasHighPHighIh_);
      minPCutMeanProfTitle+=" GeV; Ias";
      TProfile* minPCutMeanProf = minPCutMeanDir.make<TProfile>(minPCutMeanProfName.c_str(),minPCutMeanProfTitle.c_str(),numIasBins,0,1);

      std::cout << "INFO doing slice ( " << (nomSlice-1)*(1+(nomMax_-nomMin_)/2)+etaSlice << " / " << slicesToDo <<
        " ): eta=" <<
        lowerEta << "-" << lowerEta+etaStep;
      if(nom==lastLowerNoM)
        std::cout << " nom=" << nom << "-end" << std::endl;
      else if(nom==5)
        std::cout << " nom=" << nom << "-" << nom+2*nomStep-1 << std::endl;
      else
        std::cout << " nom=" << nom << "-" << nom+nomStep-1 << std::endl;

// old method
      // make histogram
      TH1F* etaCutDeDxSBPHist = (TH1F*) etaCutDeDxSBDataSet->createHistogram("trackP_",10000);
      // loop over Ih/Ias and calculate average Cj over mass cut in each Ias bin
      for(int index=0; index < etaCutNomCutBRegionDataSet->numEntries(); ++index)
      {
        etaCutNomCutBRegionDataSet->get(index);
        double thisIh = ihData_B->getVal();
        double thisIas = iasData_B->getVal();
        if(thisIas < 0) continue; // require ias >= 0
        // compute minimum p to pass mass cut
        double momSqr = (pow(massCutIasHighPHighIh_,2)*dEdx_k)/(thisIh-dEdx_c);
        if(momSqr<0)
        {
          //std::cout << "ERROR: nom=" << nom << "-" << nom+1 << " and eta=" <<
          //  lowerEta << "-" << lowerEta+0.2 << "; momSqr for ih=" << thisIh <<
          //  " mass=" << massCutIasHighPHighIh_ << " = " << momSqr << std::endl;
          continue;
        }
        double minMomPassMass = sqrt(momSqr);
        int numMomPassingMass = etaCutDeDxSBPHist->Integral(etaCutDeDxSBPHist->FindBin(minMomPassMass),etaCutDeDxSBPHist->GetNbinsX()+1);
        ceffRegionTracksOverMassCutProfile->Fill(thisIas,numMomPassingMass);
        ihMeanProf->Fill(thisIas,thisIh);
        minPCutMeanProf->Fill(thisIas,minMomPassMass);
        //debug output
        //if(iasPredictionFixedLimitsHist->FindBin(thisIas)==12)
        //{
        //  std::cout << "For ias bin 12: thisIas = " << thisIas << std::endl <<
        //    "\tCtracksPassingMassCut: " << numMomPassingMass << " CtracksFailingMassCut: " << numMomFailingMass << 
        //    " CtracksTotal: " << etaCutCRegionDataSet->numEntries() << std::endl;
        //}
      }

      // for avg ih vs. ias
      for(int bin=1; bin <= ceffRegionTracksOverMassCutAvgIhProfile->GetNbinsX(); ++bin)
      {
        double ias = ceffRegionTracksOverMassCutAvgIhProfile->GetBinCenter(bin);
        double ih = avgIhFromIasFunc->Eval(ias);
        // compute minimum p to pass mass cut
        double momSqr = (pow(massCutIasHighPHighIh_,2)*dEdx_k)/(ih-dEdx_c);
        if(momSqr<0)
          continue;
        double minMomPassMass = sqrt(momSqr);
        int numMomPassingMass = etaCutDeDxSBPHist->Integral(etaCutDeDxSBPHist->FindBin(minMomPassMass),etaCutDeDxSBPHist->GetNbinsX()+1);
        ceffRegionTracksOverMassCutAvgIhProfile->Fill(ias,numMomPassingMass);
      }

// MASS
//      // mass prediction -- base on C region with lower Pt cut rather than Ceff without P cut
//      TH1F* etaCutNomCutBRegionIhHist = (TH1F*) etaCutNomCutBRegionDataSet->createHistogram("ih_",5000);
//      TH1F* etaCutCRegionPHist = (TH1F*) etaCutCRegionDataSet->createHistogram("trackP_",5000);
//      etaCutNomCutBRegionIhHist->Scale(1.0/etaCutNomCutBRegionIhHist->Integral());
//      etaCutCRegionPHist->Scale(1.0/etaCutCRegionPHist->Integral());
//      TH1F* massPredictionFixedTmpHist = (TH1F*) massPredictionFixedHist->Clone();
//      // loop over Ih values
//      for(int ihBin=1; ihBin < etaCutNomCutBRegionIhHist->GetNbinsX()+1; ++ihBin)
//      {
//        // loop over P values
//        for(int pBin=1; pBin < etaCutCRegionPHist->GetNbinsX()+1; ++pBin)
//        {
//          double prob = etaCutCRegionPHist->GetBinContent(pBin)*etaCutNomCutBRegionIhHist->GetBinContent(ihBin);
//          double thisIh = etaCutNomCutBRegionIhHist->GetBinCenter(ihBin);
//          double thisP = etaCutCRegionPHist->GetBinCenter(pBin);
//          double massSqr = (thisIh-dEdx_c)*pow(thisP,2)/dEdx_k;
//          if(massSqr >= 0)
//            massPredictionFixedTmpHist->Fill(sqrt(massSqr),prob);
//        }
//      }
//      double norm = etaCutNomCutBRegionDataSet->numEntries()
//        *etaCutCRegionDataSet->numEntries()/(double)etaCutARegionDataSet->numEntries();
//      for(int mBin=1; mBin < massPredictionFixedHist->GetNbinsX()+1; ++mBin)
//      {
//        massPredictionFixedHist->SetBinContent(mBin,massPredictionFixedTmpHist->GetBinContent(mBin)*norm);
//        // FIXME ERROR
//        massPredictionFixedHist->SetBinError(mBin,massPredictionFixedTmpHist->GetBinError(mBin)*norm);
//      }
//      // end mass


// do exp fit on the Bk's
      // find the last bin with more than 5 entries
      std::cout << "INFO: Doing exp fit" << std::endl;
      int lastDecentStatsBin = iasBRegionHist->FindLastBinAbove(5.0);
      lastDecentStatsBin-=2; // use at least 3 bins to do the fit
      int firstIasBin = iasBRegionHist->FindBin(iasSidebandThreshold);
      if(lastDecentStatsBin < firstIasBin)
        lastDecentStatsBin = firstIasBin;
      // fill bins, including empty bins
      int lastBinWithContent = iasBRegionHist->FindLastBinAbove(1);
      double lastBinError = iasBRegionHist->GetBinError(lastBinWithContent);
      for(int bin=firstIasBin; bin < iasPredictionFixedLimitsHist->GetNbinsX()+1; ++bin)
      {
        double Bk = iasBRegionHist->GetBinContent(bin);
        if(Bk != 0)
        {
          iasPredictionFixedLimitsHist->SetBinContent(bin,Bk);
          iasPredictionFixedLimitsHist->SetBinError(bin,sqrt(Bk));
          iasPredictionFixedDiscoveryHist->SetBinContent(bin,Bk);
          iasPredictionFixedDiscoveryHist->SetBinError(bin,sqrt(Bk));
          iasPointsAndFitHist->SetBinContent(bin,Bk);
          iasPointsAndFitHist->SetBinError(bin,sqrt(Bk));
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,Bk);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,sqrt(Bk));
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,Bk);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,sqrt(Bk));
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,Bk);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,sqrt(Bk));
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,Bk);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,sqrt(Bk));
        }
        else
        {
          //iasPredictionFixedLimitsHist->SetBinContent(bin,emptyBinVal);
          //iasPredictionFixedLimitsHist->SetBinError(bin,lastBinError);
          //iasPredictionFixedDiscoveryHist->SetBinContent(bin,emptyBinVal);
          //iasPredictionFixedDiscoveryHist->SetBinError(bin,lastBinError);
          //iasPointsAndFitHist->SetBinContent(bin,emptyBinVal);
          //iasPointsAndFitHist->SetBinError(bin,lastBinError);
          // set bin content to zero for likelihood fit
          iasPredictionFixedLimitsHist->SetBinContent(bin,0);
          iasPredictionFixedLimitsHist->SetBinError(bin,lastBinError);
          iasPredictionFixedDiscoveryHist->SetBinContent(bin,0);
          iasPredictionFixedDiscoveryHist->SetBinError(bin,lastBinError);
          iasPointsAndFitHist->SetBinContent(bin,0);
          iasPointsAndFitHist->SetBinError(bin,lastBinError);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,0);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,lastBinError);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,0);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,lastBinError);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,0);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,lastBinError);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,0);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,lastBinError);
        }
      }
      // now fit from here to the end with an exp and fill the empty bins
      double lowestFitIas = (iasPredictionFixedLimitsHist->GetBinCenter(lastDecentStatsBin) < iasSidebandThreshold) ?
        0.1 : iasPredictionFixedLimitsHist->GetBinCenter(lastDecentStatsBin);
      TF1* myExp = new TF1("myExp","expo(0)",lowestFitIas,1);
      myExp->SetLineColor(kRed);
      //TFitResultPtr r = iasPointsAndFitHist->Fit("myExp","RS");
      TFitResultPtr r = iasPointsAndFitHist->Fit("myExp","RSL");
      int fitStatus = r;
      if(fitStatus != 0)
      {
        std::cout << "ERROR: Fit status is " << fitStatus << " which is nonzero!  Not using the fit." << endl;
        continue;
      }
      // for p0, must constrain +/-1 sig exps to be equal to nominal where exp fit starts
      double fitStartPoint = iasPredictionFixedLimitsHist->GetBinCenter(lastDecentStatsBin+1);
      TF1* myExpPlusOneSigma = (TF1*)myExp->Clone(); // shallower slope
      myExpPlusOneSigma->SetParameters(myExp->GetParameter(0)-fitStartPoint*myExp->GetParError(1),myExp->GetParameter(1)+myExp->GetParError(1));
      TF1* myExpMinusOneSigma = (TF1*)myExp->Clone(); // steeper slope
      myExpMinusOneSigma->SetParameters(myExp->GetParameter(0)+fitStartPoint*myExp->GetParError(1),myExp->GetParameter(1)-myExp->GetParError(1));
      cout << "DEBUG: " << " original expVal at start point = " << myExp->Eval(fitStartPoint) << " +1sigma: " << myExpPlusOneSigma->Eval(iasPredictionFixedLimitsHist->GetBinCenter(lastDecentStatsBin+1))
        << " -1sigma: " << myExpMinusOneSigma->Eval(iasPredictionFixedLimitsHist->GetBinCenter(lastDecentStatsBin+1)) << endl;
      for(int bin=lastDecentStatsBin+1; bin <= iasPredictionFixedLimitsHist->GetNbinsX(); ++bin)
      {
        double expVal = myExp->Eval(iasPredictionFixedLimitsHist->GetBinCenter(bin));
        double absErr = sqrt(expVal);
        if(isnan((float)absErr) && expVal < emptyBinVal)
        {
          iasPredictionFixedLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
          iasPredictionFixedDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
          //std::cout << "WAS nan: bin = " << bin << " set bin content: " << expVal << " +/- " << "0" << std::endl;
        }
        else if(expVal < emptyBinVal)
        {
          iasPredictionFixedLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
          iasPredictionFixedDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
          //std::cout << "expVal of " << expVal << " is less than " << emptyBinVal << "; bin = " << bin << " set bin content: " << emptyBinVal << " +/- " << sqrt(emptyBinVal) << std::endl;
        }
        else
        {
          iasPredictionFixedLimitsHist->SetBinError(bin,absErr);
          iasPredictionFixedLimitsHist->SetBinContent(bin,expVal);
          iasPredictionFixedDiscoveryHist->SetBinError(bin,absErr);
          iasPredictionFixedDiscoveryHist->SetBinContent(bin,expVal);
        }
        //std::cout << "expVal=" << expVal << " par0Err=" << par0Err << " par1Err=" << par1Err << " cov=" << cov
        //  << " xVal=" << xVal << std::endl;
      }
      // now for +1 sigma exp slope
      for(int bin=lastDecentStatsBin+1; bin <= iasPredictionFixedLimitsHist->GetNbinsX(); ++bin)
      {
        double expVal = myExpPlusOneSigma->Eval(iasPredictionFixedLimitsHist->GetBinCenter(bin));
        double absErr = sqrt(expVal);
        if(isnan((float)absErr) && expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else if(expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else
        {
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,expVal);
          iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,absErr);
        }
      }
      // now for -1 sigma exp slope
      for(int bin=lastDecentStatsBin+1; bin <= iasPredictionFixedLimitsHist->GetNbinsX(); ++bin)
      {
        double expVal = myExpMinusOneSigma->Eval(iasPredictionFixedLimitsHist->GetBinCenter(bin));
        double absErr = sqrt(expVal);
        if(isnan((float)absErr) && expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else if(expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else
        {
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,expVal);
          iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,absErr);
        }
      }
      // now for +1 sigma exp slope
      for(int bin=lastDecentStatsBin+1; bin <= iasPredictionFixedDiscoveryHist->GetNbinsX(); ++bin)
      {
        double expVal = myExpPlusOneSigma->Eval(iasPredictionFixedDiscoveryHist->GetBinCenter(bin));
        double absErr = sqrt(expVal);
        if(isnan((float)absErr) && expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else if(expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else
        {
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,expVal);
          iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,absErr);
        }
      }
      // now for -1 sigma exp slope
      for(int bin=lastDecentStatsBin+1; bin <= iasPredictionFixedDiscoveryHist->GetNbinsX(); ++bin)
      {
        double expVal = myExpMinusOneSigma->Eval(iasPredictionFixedDiscoveryHist->GetBinCenter(bin));
        double absErr = sqrt(expVal);
        if(isnan((float)absErr) && expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else if(expVal < emptyBinVal)
        {
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,emptyBinVal);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,sqrt(emptyBinVal));
        }
        else
        {
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,expVal);
          iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,absErr);
        }
      }
      delete myExp;
      delete myExpPlusOneSigma;
      delete myExpMinusOneSigma;

      // now do the C/A scaling
      // for bins that had no Bk:
      //   1) For limits, underestimate background by using Ceff/A from last non-empty bin
      //   2) For discovery, overestimate background by using Ceff/A with Ceff defined based on avg ih in that ias bin from minBias data 2011
      double entriesInARegionNoM = etaCutARegionDataSet->numEntries();
      //double entriesInCRegionNoM = etaCutCRegionDataSet->numEntries();
      double cEffLastBin = ceffRegionTracksOverMassCutProfile->GetBinContent(
          ceffRegionTracksOverMassCutProfile->FindLastBinAbove(0)); // find last ias bin with nonzero entries
      for(int bin=firstIasBin; bin <= ceffRegionTracksOverMassCutProfile->GetNbinsX(); ++bin)
      {
        // Bk * < Cj over mass > / A
        double cEff = ceffRegionTracksOverMassCutProfile->GetBinContent(bin);
        double bkOriginal = iasPredictionFixedLimitsHist->GetBinContent(bin);
        double bk = bkOriginal;
        double cEffAvgIh = ceffRegionTracksOverMassCutAvgIhProfile->GetBinContent(bin);
        
        // (1) limits
        double binContent = emptyBinVal;
        if(cEff > 0)
            binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;

        double binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);

        iasPredictionFixedLimitsHist->SetBinContent(bin,binContent);
        iasPredictionFixedLimitsHist->SetBinError(bin,binError);
        std::cout << "(limits) Bin: " << bin << " Bk: " << bk << " Ceff: " << ceffRegionTracksOverMassCutProfile->GetBinContent(bin)
          << " A: " << entriesInARegionNoM << " cEffLastBin: " << cEffLastBin
          << " Bincontent: " << binContent << " +/- " << binError << std::endl;

        // (1') +1 sigma shifted hist --> use limits underestimation for C/A scaling
        bk = iasPredictionFixedSlopeAddOneSigmaLimitsHist->GetBinContent(bin);
        binContent = emptyBinVal;
        if(cEff > 0)
            binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;
        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);
        iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinContent(bin,binContent);
        iasPredictionFixedSlopeAddOneSigmaLimitsHist->SetBinError(bin,binError);

        // (1'') -1 sigma shifted hist --> use limits underestimation for C/A scaling
        bk = iasPredictionFixedSlopeMinusOneSigmaLimitsHist->GetBinContent(bin);
        binContent = emptyBinVal;
        if(cEff > 0)
            binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;
        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);
        iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinContent(bin,binContent);
        iasPredictionFixedSlopeMinusOneSigmaLimitsHist->SetBinError(bin,binError);

        // (2) discovery
        bk = bkOriginal;
        binContent = emptyBinVal;
        if(cEff > 0)
          binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffAvgIh > 0 && cEffAvgIh > cEffLastBin && cEffLastBin > 0)
          binContent = bk*cEffAvgIh / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;
        else
          binContent = bk / entriesInARegionNoM; // assume c = 1 to be conservative

        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffAvgIh > 0 && cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffAvgIh + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);
        else
          binError = binContent*sqrt(1.0/bk + 1.0/entriesInARegionNoM); //assume c = 1 to be conservative

        iasPredictionFixedDiscoveryHist->SetBinContent(bin,binContent);
        iasPredictionFixedDiscoveryHist->SetBinError(bin,binError);
        std::cout << "(disc) Bin: " << bin << " Bk: " << bk << " CeffAvgIh: " << ceffRegionTracksOverMassCutAvgIhProfile->GetBinContent(bin)
          << " A: " << entriesInARegionNoM << " cEffAvgIh: " << cEffAvgIh
          << " Bincontent: " << binContent << " +/- " << binError << std::endl;

        // (2') +1 sigma shifted hist --> use disc overestimation for C/A scaling
        bk = iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->GetBinContent(bin);
        binContent = emptyBinVal;
        if(cEff > 0)
          binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffAvgIh > 0 && cEffAvgIh > cEffLastBin && cEffLastBin > 0)
          binContent = bk*cEffAvgIh / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;
        else
          binContent = bk / entriesInARegionNoM; // assume c = 1 to be conservative
        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffAvgIh > 0 && cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffAvgIh + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);
        else
          binError = binContent*sqrt(1.0/bk + 1.0/entriesInARegionNoM); //assume c = 1 to be conservative
        iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinContent(bin,binContent);
        iasPredictionFixedSlopeAddOneSigmaDiscoveryHist->SetBinError(bin,binError);

        // (2'') -1 sigma shifted hist --> use disc overestimation for C/A scaling
        bk = iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->GetBinContent(bin);
        binContent = emptyBinVal;
        if(cEff > 0)
          binContent = bk*cEff / entriesInARegionNoM;
        else if(cEffAvgIh > 0 && cEffAvgIh > cEffLastBin && cEffLastBin > 0)
          binContent = bk*cEffAvgIh / entriesInARegionNoM;
        else if(cEffLastBin > 0)
          binContent = bk*cEffLastBin / entriesInARegionNoM;
        else
          binContent = bk / entriesInARegionNoM; // assume c = 1 to be conservative
        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(cEffAvgIh > 0 && cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffAvgIh + 1.0/entriesInARegionNoM);
        else if(cEffLastBin > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEffLastBin + 1.0/entriesInARegionNoM);
        else
          binError = binContent*sqrt(1.0/bk + 1.0/entriesInARegionNoM); //assume c = 1 to be conservative
        iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinContent(bin,binContent);
        iasPredictionFixedSlopeMinusOneSigmaDiscoveryHist->SetBinError(bin,binError);
      }

      // remove variable bin stuff - mar 1
      //// ias prediction histogram in this NoM/eta bin -- variable
      //std::string iasPredictionVarBinHistName = getHistNameBeg(nom,lowerEta);
      //iasPredictionVarBinHistName+="iasPredictionVarBinHist";
      //std::string iasPredictionVarBinHistTitle = "Ias prediction for nom ";
      //iasPredictionVarBinHistTitle+=intToString(nom);
      //iasPredictionVarBinHistTitle+="-";
      //if(nom==21)
      //  iasPredictionVarBinHistTitle+="end";
      //else
      //  iasPredictionVarBinHistTitle+=intToString(nom+1);
      //iasPredictionVarBinHistTitle+=", ";
      //iasPredictionVarBinHistTitle+=floatToString(lowerEta);
      //iasPredictionVarBinHistTitle+=" < |#eta| < ";
      //iasPredictionVarBinHistTitle+=floatToString(lowerEta+0.2);
      //iasPredictionVarBinHistTitle+=", mass > ";
      //iasPredictionVarBinHistTitle+=floatToString(massCutIasHighPHighIh_);
      //iasPredictionVarBinHistTitle+=" GeV";
      //// figure out variable bins
      //std::vector<double> binVec = computeVariableBins(iasPredictionFixedHist);
      //TH1F* iasPredictionVarBinHist = iasPredictionVarBinsDir.make<TH1F>(iasPredictionVarBinHistName.c_str(),
      //    iasPredictionVarBinHistTitle.c_str(),binVec.size()-1,&(*binVec.begin()));
      //iasPredictionVarBinHist->Sumw2();

      //// get bin contents/errors
      //const int numVarBins = iasPredictionVarBinHist->GetNbinsX();
      //std::vector<double> binContents;
      //std::vector<double> binErrors;
      //for(int i=0; i<numVarBins; ++i)
      //{
      //  binContents.push_back(0);
      //  binErrors.push_back(0);
      //}
      //for(int bin=1; bin<=iasPredictionFixedHist->GetNbinsX();++bin)
      //{
      //  double binc = iasPredictionFixedHist->GetBinContent(bin);
      //  double bine = iasPredictionFixedHist->GetBinError(bin);
      //  if(binc > 0)
      //  {
      //    binContents[iasPredictionVarBinHist->FindBin(iasPredictionFixedHist->GetBinCenter(bin))-1] += binc;
      //    binErrors[iasPredictionVarBinHist->FindBin(iasPredictionFixedHist->GetBinCenter(bin))-1] += (bine*bine);
      //  }
      //}
      //// fill it
      //for(int i=1; i<numVarBins+1; ++i)
      //{
      //  if(binContents[i-1] > 0)
      //  {
      //    // make bins contain number of events per horizontal unit
      //    double binWidthRatio = iasPredictionVarBinHist->GetBinWidth(1)/iasPredictionVarBinHist->GetBinWidth(i);
      //    iasPredictionVarBinHist->SetBinContent(i,binWidthRatio*binContents[i-1]);
      //    iasPredictionVarBinHist->SetBinError(i,binWidthRatio*sqrt(binErrors[i-1]));
      //  }
      //}

      //// if nonzero integral, check to make sure there are no empty bins in the prediction hist
      //if(iasPredictionVarBinHist->Integral() > 0)
      //{
      //  for(int bin=1; bin <= iasPredictionVarBinHist->GetNbinsX(); ++bin)
      //  {
      //    if(iasPredictionVarBinHist->GetBinContent(bin) <= 0)
      //    {
      //      std::cout << "ERROR: For this hist: " << iasPredictionVarBinHist->GetName()
      //        << " bin content for bin "
      //        << bin << " is " << iasPredictionVarBinHist->GetBinContent(bin)
      //        << " with binLowEdge= " << iasPredictionVarBinHist->GetBinLowEdge(bin)
      //        << " to binHighEdge="
      //        << iasPredictionVarBinHist->GetBinLowEdge(bin)+iasPredictionVarBinHist->GetBinWidth(bin)
      //        << ". Fix the binning.  Bailing out." << std::endl;
      //      return -9;
      //    }
      //  }
      //}

      delete etaCutNomCutBRegionDataSet;
      delete etaCutCRegionDataSet;
      delete etaCutDeDxSBDataSet;
      delete etaCutARegionDataSet;

    } // eta slices
    if(nom==lastLowerNoM)
      break; // don't do any more after lastLowerNoM, since we've done lastLowerNoM+ in that case
    if(nom==5)
      nom+=nomStep; // if nom==5, step over the next slice, since we've already done that one
  } // nom slices

}
