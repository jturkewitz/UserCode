#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include <iostream>


void setBinLabelsPreselectionPlot(TH1F* hist)
{
  hist->GetXaxis()->SetBinLabel(1,"Trigger");
  hist->GetXaxis()->SetBinLabel(2,"NoM Ih");
  hist->GetXaxis()->SetBinLabel(3,"NoH");
  hist->GetXaxis()->SetBinLabel(4,"validFrac");
  hist->GetXaxis()->SetBinLabel(5,"PixelHits");
  hist->GetXaxis()->SetBinLabel(6,"NoM Ias");
  hist->GetXaxis()->SetBinLabel(7,"QualityMask");
  hist->GetXaxis()->SetBinLabel(8,"Chi2Ndf");
  hist->GetXaxis()->SetBinLabel(9,"Eta");
  hist->GetXaxis()->SetBinLabel(10,"PtErr");
  hist->GetXaxis()->SetBinLabel(11,"V3D");
  hist->GetXaxis()->SetBinLabel(12,"TrkIso");
  hist->GetXaxis()->SetBinLabel(13,"CaloIso");
}

// struct to hold before preselection plots
struct PlotStruct
{
  float pUpperLimit;
  float ihUpperLimit;
  TH1F* pDistributionHist;
  TH1F* ptDistributionHist;
  TH2F* pVsIhHist;
  TH2F* pVsIasHist;
  TH2F* ptVsIasHist;
  TH2F* trackEtaVsPHist;
  TH1F* iasNoMHist;
  TH2F* pVsIasToFSBHist;
  TH2F* pVsIhToFSBHist;
  TH2F* ihVsIasHist;
  TH2F* pVsNoMHist;
  TH2F* pVsNoMCentralEtaHist;
  TH2F* pVsNoMEtaSlice1Hist;
  TH2F* pVsNoMEtaSlice2Hist;
  TH2F* pVsNoMEtaSlice3Hist;
  TH2F* pVsNoMEtaSlice4Hist;
  TH2F* pVsNoMEtaSlice5Hist;
  TH2F* pVsNoMEtaSlice6Hist;
  TH2F* pVsNoMEtaSlice7Hist;
  TH2F* pVsNoMEtaSlice8Hist;
  TH2F* nohVsNoMHist;
  TH2F* nohVsNoMCentralEtaHist;
  TH2F* pVsNoHHist;
  TH2F* pVsNoHCentralEtaHist;
  TH2F* pVsRelPerrHist;
  TH2F* pVsRelPerrCentralEtaHist;
  TH1F* tracksVsPreselectionsHist;
  TH1F* massHist;

  PlotStruct(fwlite::TFileService& fs, std::string nameString, std::string titleString)
  {
    pUpperLimit = 2000.0; // GeV
    ihUpperLimit = 20; // MeV/cm  std::string dirName = "Plots";
    std::string dirName = "Plots";
    dirName+=nameString;
    TFileDirectory plotsDir = fs.mkdir(dirName.c_str());
    // p
    std::string pName = "pDistribution";
    pName+=nameString;
    std::string pTitle = "P ";
    pTitle+=titleString;
    pTitle+=";GeV";
    pDistributionHist = plotsDir.make<TH1F>(pName.c_str(),pTitle.c_str(),200,0,pUpperLimit);
    pDistributionHist->Sumw2();
    // pt
    std::string ptName = "ptDistribution";
    ptName+=nameString;
    std::string ptTitle = "Pt ";
    ptTitle+=titleString;
    ptTitle+=";GeV";
    ptDistributionHist = plotsDir.make<TH1F>(ptName.c_str(),ptTitle.c_str(),200,0,pUpperLimit);
    ptDistributionHist->Sumw2();
    // p vs Ih
    std::string pVsIhName = "trackPvsIh";
    pVsIhName+=nameString;
    std::string pVsIhTitle="Track P vs ih ";
    pVsIhTitle+=titleString;
    pVsIhTitle+=";MeV/cm;GeV";
    pVsIhHist = plotsDir.make<TH2F>(pVsIhName.c_str(),pVsIhTitle.c_str(),400,0,ihUpperLimit,200,0,pUpperLimit);
    pVsIhHist->Sumw2();
    // p vs Ias
    std::string pVsIasName = "pVsIas";
    pVsIasName+=nameString;
    std::string pVsIasTitle = "P vs Ias ";
    pVsIasTitle+=titleString;
    pVsIasTitle+=";;GeV";
    pVsIasHist = plotsDir.make<TH2F>(pVsIasName.c_str(),pVsIasTitle.c_str(),20,0,1,200,0,pUpperLimit);
    pVsIasHist->Sumw2();
    // pt vs Ias
    std::string ptVsIasName = "ptVsIas";
    ptVsIasName+=nameString;
    std::string ptVsIasTitle = "Pt vs Ias ";
    ptVsIasTitle+=titleString;
    ptVsIasTitle+=";;GeV";
    ptVsIasHist = plotsDir.make<TH2F>(ptVsIasName.c_str(),ptVsIasTitle.c_str(),20,0,1,200,0,pUpperLimit);
    ptVsIasHist->Sumw2();
    // eta vs p
    std::string etaVsPName = "trackEtaVsP";
    etaVsPName+=nameString;
    std::string etaVsPTitle = "Track #eta vs. p ";
    etaVsPTitle+=titleString;
    etaVsPTitle+=";GeV";
    trackEtaVsPHist = plotsDir.make<TH2F>(etaVsPName.c_str(),etaVsPTitle.c_str(),400,0,pUpperLimit,24,0,2.4);
    trackEtaVsPHist->Sumw2();
    // ias NoM
    std::string iasNoMName = "iasNoM";
    iasNoMName+=nameString;
    std::string iasNoMTitle = "NoM (ias) ";
    iasNoMTitle+=titleString;
    iasNoMHist = plotsDir.make<TH1F>(iasNoMName.c_str(),iasNoMTitle.c_str(),50,0,50);
    // ToF SB - p vs Ias
    std::string pVsIasToFSBName = "trackPvsIasToFSB";
    pVsIasToFSBName+=nameString;
    std::string pVsIasToFSBTitle="Track P vs ias (ToF SB: #beta>1.075) ";
    pVsIasToFSBTitle+=titleString;
    pVsIasToFSBTitle+=";;GeV";
    pVsIasToFSBHist = plotsDir.make<TH2F>(pVsIasToFSBName.c_str(),pVsIasToFSBTitle.c_str(),400,0,1,200,0,pUpperLimit);
    pVsIasToFSBHist->Sumw2();
    // ToF SB - p vs Ih
    std::string pVsIhToFSBName = "trackPvsIhToFSB";
    pVsIhToFSBName+=nameString;
    std::string pVsIhToFSBTitle="Track P vs ih (ToF SB: #beta>1.075) ";
    pVsIhToFSBTitle+=titleString;
    pVsIhToFSBTitle+=";MeV/cm;GeV";
    pVsIhToFSBHist = plotsDir.make<TH2F>(pVsIhToFSBName.c_str(),pVsIhToFSBTitle.c_str(),400,0,ihUpperLimit,200,0,pUpperLimit);
    pVsIhToFSBHist->Sumw2();
    // Ih vs Ias
    std::string ihVsIasName = "ihVsIas";
    ihVsIasName+=nameString;
    std::string ihVsIasTitle = "Ih vs. Ias ";
    ihVsIasTitle+=titleString;
    ihVsIasTitle+=";MeV/cm";
    ihVsIasHist = plotsDir.make<TH2F>(ihVsIasName.c_str(),ihVsIasTitle.c_str(),400,0,1,400,0,ihUpperLimit);
    ihVsIasHist->Sumw2();
    // p vs NoM
    std::string pVsNoMName = "pVsNoM";
    pVsNoMName+=nameString;
    std::string pVsNoMTitle = "Track P vs. NoM (Ias) ";
    pVsNoMTitle+=titleString;
    pVsNoMTitle+=";;GeV";
    pVsNoMHist = plotsDir.make<TH2F>(pVsNoMName.c_str(),pVsNoMTitle.c_str(),50,0,50,200,0,pUpperLimit);
    // p vs NoM, central eta only
    std::string pVsNoMCentralEtaName = "pVsNoMCentralEta";
    pVsNoMCentralEtaName+=nameString;
    std::string pVsNoMCentralEtaTitle = "Track P vs. NoM (Ias) ";
    pVsNoMCentralEtaTitle+=titleString;
    pVsNoMCentralEtaTitle+=";;GeV";
    pVsNoMCentralEtaHist = plotsDir.make<TH2F>(pVsNoMCentralEtaName.c_str(),pVsNoMCentralEtaTitle.c_str(),50,0,50,200,0,pUpperLimit);
    // p vs NoM, eta slices
    std::string pVsNoMEtaSliceBaseName = "pVsNoM";
    pVsNoMEtaSliceBaseName+=nameString;
    std::string pVsNoMEtaSliceBaseTitle = "Track P vs. NoM (Ias)";
    pVsNoMEtaSliceBaseTitle+=titleString;

    pVsNoMEtaSlice1Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice1").c_str(),(pVsNoMEtaSliceBaseTitle+" |#eta| < 0.2;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice2Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice2").c_str(),(pVsNoMEtaSliceBaseTitle+" 0.2 < |#eta| < 0.4;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice3Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice3").c_str(),(pVsNoMEtaSliceBaseTitle+" 0.4 < |#eta| < 0.6;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice4Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice4").c_str(),(pVsNoMEtaSliceBaseTitle+" 0.6 < |#eta| < 0.8;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice5Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice5").c_str(),(pVsNoMEtaSliceBaseTitle+" 0.8 < |#eta| < 1.0;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice6Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice6").c_str(),(pVsNoMEtaSliceBaseTitle+" 1.0 < |#eta| < 1.2;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice7Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice7").c_str(),(pVsNoMEtaSliceBaseTitle+" 1.2 < |#eta| < 1.4;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    pVsNoMEtaSlice8Hist = plotsDir.make<TH2F>((pVsNoMEtaSliceBaseName+"EtaSlice8").c_str(),(pVsNoMEtaSliceBaseTitle+" 1.4 < |#eta| < 1.6;;GeV").c_str(),50,0,50,100,0,pUpperLimit);
    // NoH vs NoM
    std::string nohVsNoMName = "nohVsNoM";
    nohVsNoMName+=nameString;
    std::string nohVsNoMTitle = "Number of hits vs. NoM (Ias) ";
    nohVsNoMTitle+=titleString;
    nohVsNoMTitle+=";NoM;NoH";
    nohVsNoMHist = plotsDir.make<TH2F>(nohVsNoMName.c_str(),nohVsNoMTitle.c_str(),50,0,50,50,0,50);
    // NoH vs NoM, central eta
    std::string nohVsNoMCentralName = "nohVsNoMCentral";
    nohVsNoMCentralName+=nameString;
    std::string nohVsNoMCentralTitle = "Number of hits vs. NoM (Ias) ";
    nohVsNoMCentralTitle+=titleString;
    nohVsNoMCentralTitle+=", |#eta| < 0.9;NoM;NoH";
    nohVsNoMCentralEtaHist = plotsDir.make<TH2F>(nohVsNoMCentralName.c_str(),nohVsNoMCentralTitle.c_str(),50,0,50,50,0,50);
    // p vs NoH
    std::string pVsNoHName = "pVsNoH";
    pVsNoHName+=nameString;
    std::string pVsNoHTitle = "Track P vs. NoH ";
    pVsNoHTitle+=titleString;
    pVsNoHTitle+=";;GeV";
    pVsNoHHist = plotsDir.make<TH2F>(pVsNoHName.c_str(),pVsNoHTitle.c_str(),50,0,50,100,0,pUpperLimit);
    // p vs NoH, central eta only
    std::string pVsNoHCentralEtaName = "pVsNoHCentralEta";
    pVsNoHCentralEtaName+=nameString;
    std::string pVsNoHCentralEtaTitle = "Track P vs. NoH ";
    pVsNoHCentralEtaTitle+=titleString;
    pVsNoHCentralEtaTitle+=", |#eta| < 0.9;;GeV";
    pVsNoHCentralEtaHist = plotsDir.make<TH2F>(pVsNoHCentralEtaName.c_str(),pVsNoHCentralEtaTitle.c_str(),50,0,50,100,0,pUpperLimit);
    // p vs relPerr
    std::string pVsRelPerrName = "pVsRelPerr";
    pVsRelPerrName+=nameString;
    std::string pVsRelPerrTitle = "Track P vs. #Deltap/p ";
    pVsRelPerrTitle+=titleString;
    pVsRelPerrTitle+=";;GeV";
    pVsRelPerrHist = plotsDir.make<TH2F>(pVsRelPerrName.c_str(),pVsRelPerrTitle.c_str(),100,0,1,100,0,pUpperLimit);
    // p vs relPerr, central eta only
    std::string pVsRelPerrCentralEtaName = "pVsRelPerrCentralEta";
    pVsRelPerrCentralEtaName+=nameString;
    std::string pVsRelPerrCentralEtaTitle = "Track P vs. #Deltap/p ";
    pVsRelPerrCentralEtaTitle+=titleString;
    pVsRelPerrCentralEtaTitle+=", |#eta| < 0.9;;GeV";
    pVsRelPerrCentralEtaHist = plotsDir.make<TH2F>(pVsRelPerrCentralEtaName.c_str(),pVsRelPerrCentralEtaTitle.c_str(),100,0,1,100,0,1000);
    // mass
    std::string massHistName = "mass";
    massHistName+=nameString;
    std::string massHistTitle = "Mass ";
    massHistTitle+=titleString;
    massHistTitle+=";GeV";
    massHist = plotsDir.make<TH1F>(massHistName.c_str(),massHistTitle.c_str(),2000,0,2000);
    // tracksVsPreselectionsHist
    tracksVsPreselectionsHist = plotsDir.make<TH1F>("tracksVsPreselections","Tracks at each preselection step",22,0,22);
    setBinLabelsPreselectionPlot(tracksVsPreselectionsHist);
  }
};


// preselection -- adapted from Analysis_Step234.C
bool passesPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj,
    const reco::DeDxData& dedxMObj, const reco::MuonTimeExtra* tof,
    const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof,
    const fwlite::Event& ev, bool considerToF, PlotStruct preselPlots)
{

  //TODO: use common py for these
  // preselections
  const int minTrackNoH = 11;
  const float minTrackValidFraction = 0.8;
  const int minValidPixelHits = 2;
  const unsigned int minIasNoM = 5;
  const unsigned int minIhNoM = 5;
  const int minTofNdof = 8;
  const int minTofNdofDt = 6;
  const int minTofNdofCsc = 6;
  const int minTrackQualityMask = 2;
  const float maxTrackChi2OverNdf = 5.0;
  //const float minTrackPt = 35.0;
  //const float minIas = 0.0;
  //const float minIh = 3.0;
  const float minTofInvBeta = 1.0;
  const float maxTofInvBetaErr = 0.07;
  const unsigned int minNumVertices = 1;
  const float maxV3D = 0.5;
  const float maxTrackEtIso = 50;
  const float maxCalEOverPIso = 0.3;
  const float maxTrackPtErr = 0.25;
  const float maxTrackEta = 2.5;

  reco::TrackRef track = hscp.trackRef();
  if(track.isNull())
    return false;

  preselPlots.tracksVsPreselectionsHist->Fill(0.5);
  // track
  if(dedxMObj.numberOfMeasurements() < minIhNoM)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(1.5);
  if(track->found() < minTrackNoH)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(2.5);
  if(track->validFraction() < minTrackValidFraction)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(3.5);
  if(track->hitPattern().numberOfValidPixelHits() < minValidPixelHits)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(4.5);
  if(dedxSObj.numberOfMeasurements() < minIasNoM)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(5.5);
  if(track->qualityMask() < minTrackQualityMask)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(6.5);
  if(track->chi2()/track->ndof() > maxTrackChi2OverNdf)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(7.5);
  if(fabs(track->eta()) > maxTrackEta)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(8.5);
  // ptError
  if(track->ptError()/track->pt() > maxTrackPtErr)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(9.5);

  // vertex
  fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
  vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
  if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
  const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
  if(vertexColl.size()<minNumVertices){printf("NO VERTEX\n"); return false;}
  double dz  = track->dz (vertexColl[0].position());
  double dxy = track->dxy(vertexColl[0].position());
  for(unsigned int i=1;i<vertexColl.size();i++){
    if(fabs(track->dz (vertexColl[i].position())) < fabs(dz) ){
      dz  = track->dz (vertexColl[i].position());
      dxy = track->dxy(vertexColl[i].position());
    }
  }
  double v3d = sqrt(dz*dz+dxy*dxy);
  if(v3d > maxV3D)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(10.5);

  // isolation
  fwlite::Handle<susybsm::HSCPIsolationValueMap> IsolationH;
  IsolationH.getByLabel(ev, "HSCPIsolation03");
  if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
  const edm::ValueMap<susybsm::HSCPIsolation>& IsolationMap = *IsolationH.product();
  susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
  if(hscpIso.Get_TK_SumEt() > maxTrackEtIso)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(11.5);
  double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
  if(EoP > maxCalEOverPIso)
    return false;
  preselPlots.tracksVsPreselectionsHist->Fill(12.5);

  // ToF
  // ndof check
  if(considerToF &&
      tof->nDof() < minTofNdof &&
      (dttof->nDof() < minTofNdofDt || csctof->nDof() < minTofNdofCsc) )
    return false;
  // Min ToF
  if(considerToF && tof->inverseBeta() < minTofInvBeta)
    return false;
  // Max ToF err
  if(considerToF && tof->inverseBetaErr() > maxTofInvBetaErr)
    return false;


  return true;
}
//
std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
std::string floatToString(float num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
std::string doubleToString(double num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}
//
int findFineNoMBin(int nom)
{
  if(nom >= 30) return 29;

  return nom-1;
}
//
std::string getHistNameBeg(int lowerNom, float lowerEta, int nomStep, float etaStep)
{
  int etaSlice = lowerEta*10;
  string histName="nom";
  histName+=intToString(lowerNom);
  histName+="to";
  if(lowerNom==5)
    histName+=intToString(lowerNom+2*nomStep-1);
  else
    histName+=intToString(lowerNom+nomStep-1);
  histName+="eta";
  histName+=intToString(etaSlice);
  histName+="to";
  histName+=intToString(etaSlice+etaStep*10);
  return histName;
}
//
int getLowerNoMFromHistName(std::string histName, std::string histNameEnd)
{
  std::string histNameBeg = histName.substr(0,histName.size()-histNameEnd.size());
  std::string nomString = histNameBeg.substr(3,histNameBeg.find("to")-3);
  return atoi(nomString.c_str());
}
//
float getLowerEtaFromHistName(std::string histName, std::string histNameEnd)
{
  std::string histNameBeg = histName.substr(0,histName.size()-histNameEnd.size());
  std::string etaPart = histNameBeg.substr(histNameBeg.find("eta"));
  std::string etaString = etaPart.substr(3,etaPart.find("to")-3);
  int etaSlice = atof(etaString.c_str());
  return etaSlice/10.0;
}
//
std::string getHistTitleEnd(int lowerNom, float lowerEta, float massCut)
{
  std::string title =" for ";
  title+=intToString(lowerNom);
  title+="-";
  if(lowerNom==5)
    title+=intToString(lowerNom+3);
  else
    title+=intToString(lowerNom+1);
  title+=", ";
  title+=floatToString(lowerEta);
  title+=" < #eta < ";
  title+=floatToString(lowerEta+0.2);
  title+=", mass > ";
  title+=floatToString(massCut);
  title+=" GeV";
  return title;
}

// taken from Analysis_Step234.C
double DistToHSCP(const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest){
  reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

  double RMin = 9999; IndexOfClosest=-1;
  for(unsigned int g=0;g<genColl.size();g++){
    if(genColl[g].pt()<5)continue;
    if(genColl[g].status()!=1)continue;
    int AbsPdg=abs(genColl[g].pdgId());
    if(AbsPdg<1000000)continue;    

    double dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
    if(dR<RMin){RMin=dR;IndexOfClosest=g;}
  }
  return RMin;
}

// modified from Analysis_Step234.C
int GetInitialNumberOfMCEvent(const vector<string>& fileNames)
{
  int Total = 0;
  fwlite::ChainEvent tree(fileNames);

  for(unsigned int f=0;f<fileNames.size();f++){
    TFile file(fileNames[f].c_str() );
    fwlite::LuminosityBlock ls( &file );
    for(ls.toBegin(); !ls.atEnd(); ++ls){
      fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
      nEventsTotalCounter.getByLabel(ls,"nEventsBefSkim");
      if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
      Total+= nEventsTotalCounter->value;
    }
  }
  return Total;
}

// taken from Analysis_Step234.C
bool passesTrigger(const fwlite::Event& ev, bool considerMuon = true, bool considerMET = true)
{
      edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");
      if(!tr.isValid())return false;

      if(considerMuon)
      {
        if(tr.accept(tr.triggerIndex("HscpPathSingleMu")))return true;
//      if(tr.accept(tr.triggerIndex("HscpPathDoubleMu")))return true;
      }
      if(considerMET)
      {
        if(tr.accept(tr.triggerIndex("HscpPathPFMet")))return true;
//      if(tr.accept(tr.triggerIndex("HscpPathCaloMet")))return true;
      }
      return false;
}

