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
int lastLowerNoM = 17;


std::string getHistTitleBeg(int lowerNom, float lowerEta, float pSB, float ihSB, float ptSB, float iasSB,
    bool usePt, bool useIas)
{
  std::string histTitle = "NoM ";
  histTitle+=intToString(lowerNom);
  histTitle+="-";
  if(lowerNom==lastLowerNoM)
    histTitle+="end";
  else
    histTitle+=intToString(lowerNom+1);
  histTitle+=", ";
  histTitle+=floatToString(lowerEta);
  histTitle+=" < #eta < ";
  histTitle+=floatToString(lowerEta+0.2);
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
  //double etaMin_ (ana.getParameter<double>("EtaMin"));
  //double etaMax_ (ana.getParameter<double>("EtaMax"));
  //int nomMin_ (ana.getParameter<int>("NoMMin"));
  //int nomMax_ (ana.getParameter<int>("NoMMax"));

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
  TH2F*  entriesInARegionHist = fs.make<TH2F>("entriesInARegion","Entries in A region (low P, low Ih);#eta;nom",12,0,2.4,9,5,23);
  entriesInARegionHist->GetYaxis()->SetNdivisions(509,false);
  TH2F*  entriesInBRegionHist = fs.make<TH2F>("entriesInBRegion","Entries in B region (Ih in low P SB);#eta;nom",12,0,2.4,9,5,23);
  entriesInBRegionHist->GetYaxis()->SetNdivisions(509,false);
  TH2F*  entriesInCRegionHist = fs.make<TH2F>("entriesInCRegion","Entries in C region (P in low Ih SB);#eta;nom",12,0,2.4,9,5,23);
  entriesInCRegionHist->GetYaxis()->SetNdivisions(509,false);
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
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",-2.5,2.5);
  TFile* inFile = TFile::Open(inputHandler_.files()[0].c_str());
  if(!inFile)
  {
    std::cout << "ERROR: input file " << inputHandler_.files()[0].c_str() << " not found." << std::endl;
    return -2;
  }

  RooDataSet* rooDataSetInput = (RooDataSet*)inFile->Get("rooDataSetCandidates");
  //RooDataSet* rooDataSetInput = (RooDataSet*)inFile->Get("rooDataSetOneCandidatePerEvent");
  if(rooDataSetInput->numEntries() < 1)
  {
    std::cout << "Problem with RooDataSet named rooDataSetCandidates in file " << inputHandler_.files()[0].c_str() << std::endl;
    return -3;
  }

  std::cout << "Doing ias predictions..." << std::endl;
  std::string pSBcutString = "rooVarP<";
  pSBcutString+=floatToString(pSidebandThreshold);
  std::string pSearchCutString = "rooVarP>";
  pSearchCutString+=floatToString(pSidebandThreshold);
  std::string ptSBcutString = "rooVarPt<";
  ptSBcutString+=floatToString(ptSidebandThreshold);
  std::string ptSearchCutString = "rooVarPt>";
  ptSearchCutString+=floatToString(ptSidebandThreshold);
  std::string ihSBcutString = "rooVarIh<";
  ihSBcutString+=floatToString(ihSidebandThreshold);
  std::string ihSearchCutString = "rooVarIh>";
  ihSearchCutString+=floatToString(ihSidebandThreshold);
  std::string iasSBcutString = "rooVarIas<";
  iasSBcutString+=floatToString(iasSidebandThreshold);
  std::string iasSearchCutString = "rooVarIas>";
  iasSearchCutString+=floatToString(iasSidebandThreshold);
  //XXX TESTING CUT OFF PT < 50
  //RooDataSet* rooDataSetPtCut = (RooDataSet*) rooDataSetInput->reduce(Cut("rooVarPt>50"));
  //std::cout << "INFO: Applying minimum Pt cut of 50" << std::endl;
  RooDataSet* rooDataSetPtCut = (RooDataSet*) rooDataSetInput->Clone();
  //XXX TURN OFF ETA > 1.5
  RooDataSet* rooDataSetAll = (RooDataSet*) rooDataSetPtCut->reduce(Cut("rooVarEta<1.5&&rooVarEta>-1.5"));
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
        SelectVars(RooArgSet(rooVarP,rooVarPt,rooVarEta,rooVarIh,rooVarNoMias)));
  else
    regionADataSet = (RooDataSet*)regionA1DataSet->reduce(Cut(ihSBcutString.c_str()),
        SelectVars(RooArgSet(rooVarP,rooVarPt,rooVarEta,rooVarIh,rooVarNoMias)));
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
          SelectVars(RooArgSet(rooVarPt,rooVarIh,rooVarIas,rooVarEta,rooVarNoMias)));
  else
    regionBDataSet = (RooDataSet*)regionB1DataSet->reduce(Cut(ihSearchCutString.c_str()),
          SelectVars(RooArgSet(rooVarPt,rooVarIh,rooVarIas,rooVarEta,rooVarNoMias)));
  delete regionB1DataSet;
  // C ==> high P/Pt, low Ih
  RooDataSet* dedxSBDataSet;
  RooDataSet* regionCDataSet;
  if(useIasForSideband)
    dedxSBDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(iasSBcutString.c_str()),
          SelectVars(RooArgSet(rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarIh)));
  else
    dedxSBDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ihSBcutString.c_str()),
          SelectVars(RooArgSet(rooVarP,rooVarPt,rooVarNoMias,rooVarEta,rooVarIh)));
  //SIC FEB 10 -- cut off p > 1200
  //RooDataSet* pCutCDataSet = (RooDataSet*)regionC1DataSet->reduce(Cut("rooVarP<1200"));
  //std::cout << "INFO: Applying max p of 1200 GeV cut in C region." << std::endl;
  std::cout << "INFO: NOT applying max p of 1200 GeV cut in C region." << std::endl;
  if(usePtForSideband)
   regionCDataSet = (RooDataSet*)dedxSBDataSet->reduce(Cut(ptSearchCutString.c_str()));
    //regionCDataSet = (RooDataSet*)pCutCDataSet->reduce(Cut(iasSBcutString.c_str()));
  else
   regionCDataSet = (RooDataSet*)dedxSBDataSet->reduce(Cut(pSearchCutString.c_str()));
    //regionCDataSet = (RooDataSet*)pCutCDataSet->reduce(Cut(ihSBcutString.c_str()));

  std::cout << "INFO: " << regionCDataSet->numEntries() << " entries in total C region." << std::endl;
  //XXX SIC MAR 2 -- define D region later using mass and ias cuts only
  // no p cut needed since eff p cut is done by mass cut
  // D ==> search region, high dE/dx; apply mass cut later
  RooDataSet* highDeDxDataSet;
  if(useIasForSideband)
    highDeDxDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(iasSearchCutString.c_str()),
          SelectVars(RooArgSet(rooVarP,rooVarEta,rooVarNoMias,rooVarIas,rooVarIh)));
  else
    highDeDxDataSet = (RooDataSet*)rooDataSetAll->reduce(Cut(ihSearchCutString.c_str()),
          SelectVars(RooArgSet(rooVarP,rooVarEta,rooVarNoMias,rooVarIas,rooVarIh)));
  delete rooDataSetAll;

  int nom = 5;
  float lowerEta = 0.0;
      // book histograms -- b region
      std::string bRegionHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      bRegionHistName+="bRegionHist";
      std::string bRegionHistTitle = "Ih for ";
      bRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      bRegionHistTitle+="; MeV/cm";
      TH1F* bRegionHist = bRegionDir.make<TH1F>(bRegionHistName.c_str(),bRegionHistTitle.c_str(),100,0,10);
      // ias distribution from B region
      std::string iasBRegionHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      iasBRegionHistName+="iasBRegionHist";
      std::string iasBRegionHistTitle = "Ias in B region for ";
      iasBRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      TH1F* iasBRegionHist = bRegionDir.make<TH1F>(iasBRegionHistName.c_str(),iasBRegionHistTitle.c_str(),numIasBins,0,1);
      // c region hist
      std::string cRegionHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      cRegionHistName+="cRegionHist";
      std::string cRegionHistTitle = "P for ";
      cRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      cRegionHistTitle+="; GeV";
      TH1F* cRegionHist = cRegionDir.make<TH1F>(cRegionHistName.c_str(),cRegionHistTitle.c_str(),100,0,1000);
      // c region hist - cumulative
      std::string cRegionCumuHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      cRegionCumuHistName+="cRegionCumuHist";
      std::string cRegionCumuHistTitle = "P for ";
      cRegionCumuHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      cRegionCumuHistTitle+=" cumulative; GeV";
      TH1F* cRegionCumuHist = cRegionDir.make<TH1F>(cRegionCumuHistName.c_str(),cRegionCumuHistTitle.c_str(),100,0,1000);
      // ceff region hist
      std::string ceffRegionHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      ceffRegionHistName+="ceffRegionHist";
      std::string ceffRegionHistTitle = "P for ";
      ceffRegionHistTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      ceffRegionHistTitle+="; GeV";
      TH1F* ceffRegionHist = cRegionDir.make<TH1F>(ceffRegionHistName.c_str(),ceffRegionHistTitle.c_str(),100,0,1000);
      // ceff region over mass cut in ias bins hist
      std::string ceffRegionTracksOverMassCutProfileName = getHistNameBeg(nom,lowerEta,2,0.2);
      ceffRegionTracksOverMassCutProfileName+="tracksInCeffOverMassCutProfile";
      std::string ceffRegionTracksOverMassCutProfileTitle = "Tracks in Ceff Over Mass Cut in Ias bins for ";
      ceffRegionTracksOverMassCutProfileTitle+=getHistTitleBeg(nom,lowerEta,pSidebandThreshold,ihSidebandThreshold,ptSidebandThreshold,
          iasSidebandThreshold,usePtForSideband,useIasForSideband);
      ceffRegionTracksOverMassCutProfileTitle+="; Ias";
      TProfile* ceffRegionTracksOverMassCutProfile =
        cRegionDir.make<TProfile>(ceffRegionTracksOverMassCutProfileName.c_str(),
            ceffRegionTracksOverMassCutProfileTitle.c_str(),numIasBins,0,1);

      // B region dataset
      RooDataSet* nomCutBRegionDataSet = (RooDataSet*)regionBDataSet->Clone();
      RooDataSet* etaCutNomCutBRegionDataSet = (RooDataSet*)nomCutBRegionDataSet->Clone();
      // C region dataset
      // SIC Dec. 23: study shows we can use all NoM slices
      RooDataSet* etaCutCRegionDataSet = (RooDataSet*)regionCDataSet->Clone();
      RooDataSet* etaCutNomCutCRegionDataSet = (RooDataSet*)etaCutCRegionDataSet->Clone();
      // dE/dx SB data set cut in eta --> used for Ceff region later
      RooDataSet* etaCutDeDxSBDataSet = (RooDataSet*)dedxSBDataSet->Clone();
      // High dE/dx dataset
      RooDataSet* nomCutHighDeDxDataSet = (RooDataSet*)highDeDxDataSet->Clone();
      RooDataSet* etaCutNomCutHighDeDxDataSet = (RooDataSet*)nomCutHighDeDxDataSet->Clone();
      // A region dataset
      RooDataSet* etaCutARegionDataSet = (RooDataSet*)regionADataSet->Clone();
      RooDataSet* etaCutNomCutARegionDataSet = (RooDataSet*)etaCutARegionDataSet->Clone();

      entriesInARegionHist->Fill(lowerEta+0.1,nom,etaCutNomCutARegionDataSet->numEntries());
      entriesInBRegionHist->Fill(lowerEta+0.1,nom,etaCutNomCutBRegionDataSet->numEntries());
      entriesInCRegionHist->Fill(lowerEta+0.1,nom,etaCutNomCutCRegionDataSet->numEntries());
      //entriesInDRegionHist->Fill(lowerEta+0.1,nom,etaCutNomCutDRegionDataSet->numEntries());
      //debug/output
      std::cout << "eta=" <<
        lowerEta << "-";
      //XXX TESTING eta 1.5 max
      //if(10*lowerEta==14)
      //  std::cout << "1.5";
      //else
      //  std::cout << lowerEta+0.2;
      std::cout << lowerEta+0.2;
      std::cout << " nom=" << nom << "-" << nom+1 << std::endl;
      std::cout << "A = " << etaCutNomCutARegionDataSet->numEntries() <<
        " B = " << etaCutNomCutBRegionDataSet->numEntries() <<
        " C = " << etaCutCRegionDataSet->numEntries() << std::endl;
       // " D = " << etaCutNomCutDRegionDataSet->numEntries() << std::endl;

      const RooArgSet* argSet_A = etaCutNomCutARegionDataSet->get();
      //RooRealVar* iasData_A = (RooRealVar*)argSet_A->find(rooVarIas.GetName());
      RooRealVar* ihData_A = (RooRealVar*)argSet_A->find(rooVarIh.GetName());
      RooRealVar* etaData_A = (RooRealVar*)argSet_A->find(rooVarEta.GetName());
      RooRealVar* ptData_A = (RooRealVar*)argSet_A->find(rooVarPt.GetName());
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
      RooRealVar* pData_C1 = (RooRealVar*)argSet_C1->find(rooVarP.GetName());
      RooRealVar* ptData_C1 = (RooRealVar*)argSet_C1->find(rooVarPt.GetName());
      RooRealVar* ihData_C1 = (RooRealVar*)argSet_C1->find(rooVarIh.GetName());
      RooRealVar* etaData_C1 = (RooRealVar*)argSet_C1->find(rooVarEta.GetName());
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
      RooRealVar* ihData_B = (RooRealVar*)argSet_B->find(rooVarIh.GetName());
      RooRealVar* iasData_B = (RooRealVar*)argSet_B->find(rooVarIas.GetName());
      RooRealVar* etaData_B = (RooRealVar*)argSet_B->find(rooVarEta.GetName());
      RooRealVar* ptData_B = (RooRealVar*)argSet_B->find(rooVarPt.GetName());
      //remove eta/nom double checking
      //RooRealVar* nomData_B = (RooRealVar*)argSet_B->find(rooVarNoMias.GetName());
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
      RooRealVar* pData_C = (RooRealVar*)argSet_C->find(rooVarP.GetName());
      RooRealVar* ptData_C = (RooRealVar*)argSet_C->find(rooVarPt.GetName());
      //RooRealVar* etaData_C = (RooRealVar*)argSet_C->find(rooVarEta.GetName());
      //RooRealVar* nomData_C = (RooRealVar*)argSet_C->find(rooVarNoMias.GetName());
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
      RooRealVar* pData_Ceff = (RooRealVar*)argSet_Ceff->find(rooVarP.GetName());
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
      std::string iasPredictionFixedLimitsHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      iasPredictionFixedLimitsHistName+="iasPredictionFixedLimitsHist";
      std::string iasPredictionFixedLimitsHistTitle = "Ias prediction (limits) for ";
      iasPredictionFixedLimitsHistTitle+=intToString(nom);
      iasPredictionFixedLimitsHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedLimitsHistTitle+="end";
      else
        iasPredictionFixedLimitsHistTitle+=intToString(nom+1);
      iasPredictionFixedLimitsHistTitle+=", ";
      iasPredictionFixedLimitsHistTitle+=floatToString(lowerEta);
      iasPredictionFixedLimitsHistTitle+=" < |#eta| < ";
      iasPredictionFixedLimitsHistTitle+=floatToString(lowerEta+0.2);
      iasPredictionFixedLimitsHistTitle+=", mass > ";
      iasPredictionFixedLimitsHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedLimitsHistTitle+=" GeV";
      TH1D* iasPredictionFixedLimitsHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedLimitsHistName.c_str(),iasPredictionFixedLimitsHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedLimitsHist->Sumw2();
      // ias prediction histogram for discovery  -- overestimation of background
      std::string iasPredictionFixedDiscoveryHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      iasPredictionFixedDiscoveryHistName+="iasPredictionFixedDiscoveryHist";
      std::string iasPredictionFixedDiscoveryHistTitle = "Ias prediction (discovery) for ";
      iasPredictionFixedDiscoveryHistTitle+=intToString(nom);
      iasPredictionFixedDiscoveryHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPredictionFixedDiscoveryHistTitle+="end";
      else
        iasPredictionFixedDiscoveryHistTitle+=intToString(nom+1);
      iasPredictionFixedDiscoveryHistTitle+=", ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(lowerEta);
      iasPredictionFixedDiscoveryHistTitle+=" < |#eta| < ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(lowerEta+0.2);
      iasPredictionFixedDiscoveryHistTitle+=", mass > ";
      iasPredictionFixedDiscoveryHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPredictionFixedDiscoveryHistTitle+=" GeV";
      TH1D* iasPredictionFixedDiscoveryHist = iasPredictionFixedBinsDir.make<TH1D>(iasPredictionFixedDiscoveryHistName.c_str(),iasPredictionFixedDiscoveryHistTitle.c_str(),numIasBins,0,1);
      iasPredictionFixedDiscoveryHist->Sumw2();
      // ias histogram points and fit in this NoM/eta bin
      std::string iasPointsAndFitHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      iasPointsAndFitHistName+="iasPointsAndFitHist";
      std::string iasPointsAndFitHistTitle = "Ias prediction for ";
      iasPointsAndFitHistTitle+=intToString(nom);
      iasPointsAndFitHistTitle+="-";
      if(nom==lastLowerNoM)
        iasPointsAndFitHistTitle+="end";
      else
        iasPointsAndFitHistTitle+=intToString(nom+1);
      iasPointsAndFitHistTitle+=", ";
      iasPointsAndFitHistTitle+=floatToString(lowerEta);
      iasPointsAndFitHistTitle+=" < |#eta| < ";
      iasPointsAndFitHistTitle+=floatToString(lowerEta+0.2);
      iasPointsAndFitHistTitle+=", mass > ";
      iasPointsAndFitHistTitle+=floatToString(massCutIasHighPHighIh_);
      iasPointsAndFitHistTitle+=" GeV";
      TH1D* iasPointsAndFitHist = iasPredictionFixedBinsDir.make<TH1D>(iasPointsAndFitHistName.c_str(),iasPointsAndFitHistTitle.c_str(),numIasBins,0,1);
      iasPointsAndFitHist->Sumw2();
      // mass prediction histogram in this NoM/eta bin
      std::string massPredictionFixedHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      massPredictionFixedHistName+="massPredictionFixedHist";
      std::string massPredictionFixedHistTitle = "mass prediction for ";
      massPredictionFixedHistTitle+=intToString(nom);
      massPredictionFixedHistTitle+="-";
      if(nom==lastLowerNoM)
        massPredictionFixedHistTitle+="end";
      else
        massPredictionFixedHistTitle+=intToString(nom+1);
      massPredictionFixedHistTitle+=", ";
      massPredictionFixedHistTitle+=floatToString(lowerEta);
      massPredictionFixedHistTitle+=" < |#eta| < ";
      massPredictionFixedHistTitle+=floatToString(lowerEta+0.2);
      massPredictionFixedHistTitle+=", mass > ";
      massPredictionFixedHistTitle+=floatToString(massCutIasHighPHighIh_);
      massPredictionFixedHistTitle+=" GeV";
      TH1D* massPredictionFixedHist = massPredictionFixedBinsDir.make<TH1D>(massPredictionFixedHistName.c_str(),massPredictionFixedHistTitle.c_str(),1000,0,5000);
      massPredictionFixedHist->Sumw2();
      // success rate histogram in this NoM/eta bin -- removed APR 14
      // ias D region histogram in this NoM/eta bin
      std::string dRegionFixedHistName = getHistNameBeg(nom,lowerEta,2,0.2);
      dRegionFixedHistName+="dRegionFixedHist";
      std::string dRegionFixedHistTitle = "Ias data in D region for ";
      dRegionFixedHistTitle+=intToString(nom);
      dRegionFixedHistTitle+="-";
      if(nom==lastLowerNoM)
        dRegionFixedHistTitle+="end";
      else
        dRegionFixedHistTitle+=intToString(nom+1);
      dRegionFixedHistTitle+=", ";
      dRegionFixedHistTitle+=floatToString(lowerEta);
      dRegionFixedHistTitle+=" < |#eta| < ";
      dRegionFixedHistTitle+=floatToString(lowerEta+0.2);
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
        std::cout << lowerEta+0.2;
      std::cout << " nom=" << nom << "-";
      if(nom==lastLowerNoM)
        std::cout << "end" << std::endl;
      else
        std::cout << nom+1 << std::endl;

      std::cout << "\t\t" << etaCutNomCutBRegionDataSet->numEntries() <<
        " entries in B region and " << etaCutCRegionDataSet->numEntries() <<
        " entries in C region (summed over NoM) and " << etaCutARegionDataSet->numEntries() <<
        " entries in A region (summed over NoM)" << endl <<
        " C/A = " << etaCutCRegionDataSet->numEntries()/(double)etaCutARegionDataSet->numEntries() << std::endl;

      // fill D region hist
      const RooArgSet* argSet_highDeDx = etaCutNomCutHighDeDxDataSet->get();
      RooRealVar* ihData_highDeDx = (RooRealVar*)argSet_highDeDx->find(rooVarIh.GetName());
      RooRealVar* iasData_highDeDx = (RooRealVar*)argSet_highDeDx->find(rooVarIas.GetName());
      RooRealVar* pData_highDeDx = (RooRealVar*)argSet_highDeDx->find(rooVarP.GetName());
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
          }
        }
        return -1;
      }

      // Ih mean in Ias bins histogram in this NoM/eta bin
      std::string ihMeanProfName = getHistNameBeg(nom,lowerEta,2,0.2);
      ihMeanProfName+="ihMeanProf";
      std::string ihMeanProfTitle = "Ih mean for nom ";
      ihMeanProfTitle+=intToString(nom);
      ihMeanProfTitle+="-";
      if(nom==lastLowerNoM)
        ihMeanProfTitle+="end";
      else
        ihMeanProfTitle+=intToString(nom+1);
      ihMeanProfTitle+=", ";
      ihMeanProfTitle+=floatToString(lowerEta);
      ihMeanProfTitle+=" < |#eta| < ";
      ihMeanProfTitle+=floatToString(lowerEta+0.2);
      ihMeanProfTitle+=", mass > ";
      ihMeanProfTitle+=floatToString(massCutIasHighPHighIh_);
      ihMeanProfTitle+=" GeV; Ias";
      TProfile* ihMeanProf = ihMeanDir.make<TProfile>(ihMeanProfName.c_str(),ihMeanProfTitle.c_str(),numIasBins,0,1);
      // Min P cut mean in Ias bins histogram in this NoM/eta bin
      std::string minPCutMeanProfName = getHistNameBeg(nom,lowerEta,2,0.2);
      minPCutMeanProfName+="minPCutMeanProf";
      std::string minPCutMeanProfTitle = "Min P cut for nom ";
      minPCutMeanProfTitle+=intToString(nom);
      minPCutMeanProfTitle+="-";
      if(nom==17)
        minPCutMeanProfTitle+="end";
      else
        minPCutMeanProfTitle+=intToString(nom+1);
      minPCutMeanProfTitle+=", ";
      minPCutMeanProfTitle+=floatToString(lowerEta);
      minPCutMeanProfTitle+=" < |#eta| < ";
      minPCutMeanProfTitle+=floatToString(lowerEta+0.2);
      minPCutMeanProfTitle+=", mass > ";
      minPCutMeanProfTitle+=floatToString(massCutIasHighPHighIh_);
      minPCutMeanProfTitle+=" GeV; Ias";
      TProfile* minPCutMeanProf = minPCutMeanDir.make<TProfile>(minPCutMeanProfName.c_str(),minPCutMeanProfTitle.c_str(),numIasBins,0,1);

      //std::cout << "INFO doing slice ( " << (nomSlice-1)*(1+(nomMax_-nomMin_)/2)+etaSlice << " / " << slicesToDo <<
      //  " ): eta=" <<
      //  lowerEta << "-" << lowerEta+0.2;
      //if(nom==17)
      //  std::cout << " nom=" << nom << "-end" << std::endl;
      //else
      //  std::cout << " nom=" << nom << "-" << nom+1 << std::endl;

// old method
      // make histogram
      TH1F* etaCutDeDxSBPHist = (TH1F*) etaCutDeDxSBDataSet->createHistogram("rooVarP",10000);
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

// MASS
//      // mass prediction -- base on C region with lower Pt cut rather than Ceff without P cut
//      TH1F* etaCutNomCutBRegionIhHist = (TH1F*) etaCutNomCutBRegionDataSet->createHistogram("rooVarIh",5000);
//      TH1F* etaCutCRegionPHist = (TH1F*) etaCutCRegionDataSet->createHistogram("rooVarP",5000);
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

//XXX SIC MAR 1 new ias prediction
//
//      // create profile from etaCutNomCutBRegionDataSet, ias bins
//      for(int index=0; index < etaCutNomCutBRegionDataSet->numEntries(); ++index)
//      {
//        etaCutNomCutBRegionDataSet->get(index);
//        double thisIh = ihData_B->getVal();
//        double thisIas = iasData_B->getVal();
//        if(thisIas < 0) continue; // require ias >= 0
//        ihMeanProf->Fill(thisIas,thisIh);
//      }
//
//      // loop over the profile
//      for(int iasBin = 1; iasBin <= ihMeanProf->GetNbinsX(); ++iasBin)
//      {
//        double thisIas = ihMeanProf->GetBinCenter(iasBin);
//        double thisIh = ihMeanProf->GetBinContent(iasBin);
//        // compute minimum p to pass mass cut
//        double momSqr = (pow(massCutIasHighPHighIh_,2)*dEdx_k)/(thisIh-dEdx_c);
//        if(momSqr<0)
//        {
//          //std::cout << "ERROR: nom=" << nom << "-" << nom+1 << " and eta=" <<
//          //  lowerEta << "-" << lowerEta+0.2 << "; momSqr for ih=" << thisIh <<
//          //  " mass=" << massCutIasHighPHighIh_ << " = " << momSqr << std::endl;
//          continue;
//        }
//        double minMomPassMass = sqrt(momSqr);
//        RooDataSet* momPassDataSet = (RooDataSet*)etaCutDeDxSBDataSet->reduce(Cut(("rooVarP>"+floatToString(minMomPassMass)).c_str()),
//            SelectVars(rooVarP));
//        int numMomPassingMass = momPassDataSet->numEntries();
//        delete momPassDataSet;
//        cRegionTracksOverMassCutProfile->Fill(thisIas,numMomPassingMass);
//        minPCutMeanProf->Fill(thisIas,minMomPassMass);
//        
//        double bEntries = ihMeanProf->GetBinEntries(iasBin);
//        double aEntries = etaCutARegionDataSet->numEntries();
//        double binContent = bEntries*numMomPassingMass/aEntries;
//        double errSumSqr = 0;
//        if(bEntries > 0)
//          errSumSqr+=1/bEntries;
//        if(aEntries > 0)
//          errSumSqr+=1/aEntries;
//        if(numMomPassingMass > 0)
//          errSumSqr+=1/numMomPassingMass;
//        double binError = binContent*sqrt(errSumSqr);
//        std::cout << "Ias: " << thisIas << " binContent: " << binContent <<
//          " bEntries: " << ihMeanProf->GetBinEntries(iasBin) << 
//          " cEffEntries: " << numMomPassingMass <<
//          " aEntries: " << etaCutARegionDataSet->numEntries()
//          << std::endl;
//        iasPredictionFixedHist->SetBinContent(iasPredictionFixedHist->FindBin(thisIas),binContent);
//        iasPredictionFixedHist->SetBinError(iasPredictionFixedHist->FindBin(thisIas),binError);
//        //if(nom==11 && (int)(lowerEta*10)==4)
//        //{
//        //std::cout << "numMomPassingMass: " << numMomPassingMass << " numMomFailingMass: " <<
//        //  numMomFailingMass << " minMomPassMass:" << minMomPassMass << " thisIh: " << thisIh << std::endl;
//        //std::cout << "ias: " << thisIas << " successRate: " << successRate << " successRateError:"
//        //  << successRateError << " trials:"  << total << std::endl;
//        //}
//      }

// do exp fit on the Bk's
      // find the last bin with more than 5 entries
      std::cout << "INFO: Doing exp fit" << std::endl;
      int lastDecentStatsBin = iasBRegionHist->FindLastBinAbove(5.0);
      lastDecentStatsBin-=2; // use at least 3 bins to do the fit
      int firstIasBin = iasBRegionHist->FindBin(iasSidebandThreshold);
      if(lastDecentStatsBin < 1)
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
        return -1;
      }
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
      delete myExp;

      // now do the C/A scaling
      // for bins that had no Bk:
      //   1) For limits, underestimate background by using C/A from last non-empty bin
      //   2) For discovery, overestimate background by using C/A without mass cut on C
      double entriesInARegionNoM = etaCutARegionDataSet->numEntries();
      double entriesInCRegionNoM = etaCutCRegionDataSet->numEntries();
      double cEffLastBin = ceffRegionTracksOverMassCutProfile->GetBinContent(
          ceffRegionTracksOverMassCutProfile->FindLastBinAbove(10)); // find last ias bin with > 10 entries
      for(int bin=firstIasBin; bin <= ceffRegionTracksOverMassCutProfile->GetNbinsX(); ++bin)
      {
        // Bk * < Cj over mass > / A
        double cEff = ceffRegionTracksOverMassCutProfile->GetBinContent(bin);
        double bk = iasPredictionFixedLimitsHist->GetBinContent(bin);
        
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

        std::cout << "Bin: " << bin << " Bk: " << bk << " Ceff: " << ceffRegionTracksOverMassCutProfile->GetBinContent(bin)
          << " A: " << entriesInARegionNoM << " cEffLastBin: " << cEffLastBin
          << " Bincontent: " << binContent << " +/- " << binError << std::endl;

        binContent = emptyBinVal;
        if(cEff > 0)
          binContent = bk*cEff / entriesInARegionNoM;
        else if(entriesInCRegionNoM > 0)
          binContent = bk*entriesInCRegionNoM / entriesInARegionNoM;

        binError = sqrt(emptyBinVal);
        if(cEff > 0)
          binError = binContent*sqrt(1.0/bk + 1.0/cEff + 1.0/entriesInARegionNoM);
        else if(entriesInCRegionNoM > 0) // this check shouldn't be necessary
          binError = binContent*sqrt(1.0/bk + 1.0/entriesInCRegionNoM + 1.0/entriesInARegionNoM);

        iasPredictionFixedDiscoveryHist->SetBinContent(bin,binContent);
        iasPredictionFixedDiscoveryHist->SetBinError(bin,binError);
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

}
