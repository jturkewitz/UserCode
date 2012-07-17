#include "TH3.h"
#include "TH2.h"
#include "TProfile.h"
#include "TH1.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TDirectory.h"

#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>

// ****** helper functions
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
float getLastBinLowEdge(int nom)
{
  if(nom<7)
    return 0.9;
  else if(nom==7 || nom==8)
    return 0.65;
  else if(nom==9 || nom==10)
    return 0.5;
  else if(nom==11 || nom==12)
    return 0.4;
  else if(nom==13 || nom==14)
    return 0.4;
  else if(nom==15 || nom==16)
    return 0.3;
  else if(nom==17 || nom==18)
    return 0.25;
  else if(nom==19 || nom==20)
    return 0.2;
  else
    return 0.1;
}
//
float getSecondToLastBinLowEdge(int nom)
{
  if(nom<7)
    return 0.8;
  else if(nom==7 || nom==8)
    return 0.55;
  else if(nom==9 || nom==10)
    return 0.4;
  else if(nom==11 || nom==12)
    return 0.35;
  else if(nom==13 || nom==14)
    return 0.35;
  else if(nom==15 || nom==16)
    return 0.25;
  else if(nom==17 || nom==18)
    return 0.2;
  else if(nom==19 || nom==20)
    return 0.15;
  else
    return 0.05;
}
//
void moveBinsTH2F(TH2F* myprof)
{
  int nxbins = myprof->GetNbinsX();
  int nybins = myprof->GetNbinsY();

  for(int i=0; i<=(nxbins+2)*(nybins+2); i++ )
  {
    Double_t binents = myprof->GetBinContent(i);
    Double_t binerr = myprof->GetBinError(i);
    if(binents == 0 && binerr == 0)
    {
      myprof->SetBinContent(i,1);
      myprof->SetBinContent(i,-1000);
    }
  }
  return;
}



// start from 0 for slices
// 0-11 for NoM (slice 11--> NoM>=24)
// 0-11 for eta (slice 11--> 2.2 < |eta| < 2.4)
int main(int argc, char* argv[])
{
  using namespace std;
  if(argc < 5)
  {
    cout << "Usage error." << endl
      << "makeIasDists nomSlice etaSlice massCut rootfile" << endl
      << "nomSlice, etaSlice can be between 0 and 11" << endl;
    return -1;
  }

  float nomSlice = atof(argv[1]);
  float etaSlice = atof(argv[2]);
  float massCut = atof(argv[3]);
  char* rootFilename = argv[4];

  if(nomSlice < 0 || nomSlice > 11 || etaSlice < 0 || etaSlice > 11)
  {
    cout << "Usage error." << endl
      << "makeIasDists nomSlice etaSlice massCut rootfile" << endl
      << "nomSlice, etaSlice can be between 0 and 11" << endl;
    return -1;
  }

  string predHistName="trackIasFactorizedIhPNoMSlice";
  predHistName+=intToString(nomSlice);
  predHistName+="etaSlice";
  predHistName+=intToString(etaSlice);
  predHistName+="massCut";
  predHistName+=intToString(massCut);

  int nomLow = nomSlice*2+1;
  int nomHigh = nomSlice*2+2;
  float etaLow = etaSlice*0.2;
  float etaHigh = (etaSlice+1)*0.2;
  string predHistTitle="Track Ias from Ih,P SB hists with mass > ";
  predHistTitle+=floatToString(massCut);
  predHistTitle+=" for nom ";
  predHistTitle+=intToString(nomLow);
  predHistTitle+="-";
  predHistTitle+=intToString(nomHigh);
  predHistTitle+=", |#eta| ";
  predHistTitle+=floatToString(etaLow);
  predHistTitle+="-";
  predHistTitle+=floatToString(etaHigh);
  predHistTitle+=";Ias";

  string predHistFixedBinsName="fixedBinsIasFactorizedIhPNoMSlice";
  predHistFixedBinsName+=intToString(nomSlice);
  predHistFixedBinsName+="etaSlice";
  predHistFixedBinsName+=intToString(etaSlice);
  predHistFixedBinsName+="massCut";
  predHistFixedBinsName+=intToString(massCut);
  string predHistFixedBinsTitle="Track Ias from Ih,P SB hists with mass > ";
  predHistFixedBinsTitle+=floatToString(massCut);
  predHistFixedBinsTitle+=" for nom ";
  predHistFixedBinsTitle+=intToString(nomLow);
  predHistFixedBinsTitle+="-";
  predHistFixedBinsTitle+=intToString(nomHigh);
  predHistFixedBinsTitle+=", |#eta| ";
  predHistFixedBinsTitle+=floatToString(etaLow);
  predHistFixedBinsTitle+="-";
  predHistFixedBinsTitle+=floatToString(etaHigh);
  predHistFixedBinsTitle+=";Ias";

  string successRate2DHistName="successRate2DOfIasNoMSlice";
  successRate2DHistName+=intToString(nomSlice);
  successRate2DHistName+="etaSlice";
  successRate2DHistName+=intToString(etaSlice);
  successRate2DHistName+="massCut";
  successRate2DHistName+=intToString(massCut);
  string successRate2DHistTitle="Ias success rate with mass > ";
  successRate2DHistTitle+=floatToString(massCut);
  successRate2DHistTitle+=" for nom ";
  successRate2DHistTitle+=intToString(nomLow);
  successRate2DHistTitle+="-";
  successRate2DHistTitle+=intToString(nomHigh);
  successRate2DHistTitle+=", |#eta| ";
  successRate2DHistTitle+=floatToString(etaLow);
  successRate2DHistTitle+="-";
  successRate2DHistTitle+=floatToString(etaHigh);
  successRate2DHistTitle+=";Ias;Ih [MeV/cm]";

  string successRateHistName="successRateOfIasNoMSlice";
  successRateHistName+=intToString(nomSlice);
  successRateHistName+="etaSlice";
  successRateHistName+=intToString(etaSlice);
  successRateHistName+="massCut";
  successRateHistName+=intToString(massCut);
  string successRateHistTitle="Ias success rate with mass > ";
  successRateHistTitle+=floatToString(massCut);
  successRateHistTitle+=" for nom ";
  successRateHistTitle+=intToString(nomLow);
  successRateHistTitle+="-";
  successRateHistTitle+=intToString(nomHigh);
  successRateHistTitle+=", |#eta| ";
  successRateHistTitle+=floatToString(etaLow);
  successRateHistTitle+="-";
  successRateHistTitle+=floatToString(etaHigh);
  successRateHistTitle+=";Ias";

  string successIhHistName="successIhNoMSlice";
  successIhHistName+=intToString(nomSlice);
  successIhHistName+="etaSlice";
  successIhHistName+=intToString(etaSlice);
  successIhHistName+="massCut";
  successIhHistName+=intToString(massCut);
  string successIhHistTitle="Ih predicted with mass > ";
  successIhHistTitle+=floatToString(massCut);
  successIhHistTitle+=" for nom ";
  successIhHistTitle+=intToString(nomLow);
  successIhHistTitle+="-";
  successIhHistTitle+=intToString(nomHigh);
  successIhHistTitle+=", |#eta| ";
  successIhHistTitle+=floatToString(etaLow);
  successIhHistTitle+="-";
  successIhHistTitle+=floatToString(etaHigh);
  successIhHistTitle+=";Ih [MeV/cm]";

  string successRateIhHistName="successRateIhNoMSlice";
  successRateIhHistName+=intToString(nomSlice);
  successRateIhHistName+="etaSlice";
  successRateIhHistName+=intToString(etaSlice);
  successRateIhHistName+="massCut";
  successRateIhHistName+=intToString(massCut);
  string successRateIhHistTitle="Ih predicted with mass > ";
  successRateIhHistTitle+=floatToString(massCut);
  successRateIhHistTitle+=" for nom ";
  successRateIhHistTitle+=intToString(nomLow);
  successRateIhHistTitle+="-";
  successRateIhHistTitle+=intToString(nomHigh);
  successRateIhHistTitle+=", |#eta| ";
  successRateIhHistTitle+=floatToString(etaLow);
  successRateIhHistTitle+="-";
  successRateIhHistTitle+=floatToString(etaHigh);
  successRateIhHistTitle+=";Ih [MeV/cm]";

  string trialsHistName="trialsInIasNoMSlice";
  trialsHistName+=intToString(nomSlice);
  trialsHistName+="etaSlice";
  trialsHistName+=intToString(etaSlice);
  trialsHistName+="massCut";
  trialsHistName+=intToString(massCut);
  string trialsHistTitle="Trials with mass > ";
  trialsHistTitle+=floatToString(massCut);
  trialsHistTitle+=" for nom ";
  trialsHistTitle+=intToString(nomLow);
  trialsHistTitle+="-";
  trialsHistTitle+=intToString(nomHigh);
  trialsHistTitle+=", |#eta| ";
  trialsHistTitle+=floatToString(etaLow);
  trialsHistTitle+="-";
  trialsHistTitle+=floatToString(etaHigh);
  trialsHistTitle+=";Ias;Ih [MeV/cm]";

  string selectedIhvsPHistName="selectedIhvsPNoMSlice";
  selectedIhvsPHistName+=intToString(nomSlice);
  selectedIhvsPHistName+="etaSlice";
  selectedIhvsPHistName+=intToString(etaSlice);
  selectedIhvsPHistName+="massCut";
  selectedIhvsPHistName+=intToString(massCut);
  string selectedIhvsPHistTitle="Selected Ih vs. P for mass > ";
  selectedIhvsPHistTitle+=floatToString(massCut);
  selectedIhvsPHistTitle+=" for nom ";
  selectedIhvsPHistTitle+=intToString(nomLow);
  selectedIhvsPHistTitle+="-";
  selectedIhvsPHistTitle+=intToString(nomHigh);
  selectedIhvsPHistTitle+=", |#eta| ";
  selectedIhvsPHistTitle+=floatToString(etaLow);
  selectedIhvsPHistTitle+="-";
  selectedIhvsPHistTitle+=floatToString(etaHigh);
  selectedIhvsPHistTitle+=";GeV;MeV/cm";

  TFile* myOutputFile = new TFile((predHistName+".root").c_str(),"recreate");
  myOutputFile->cd();
  TH1F* massGeneratedHist = new TH1F("massGenerated","Masses generated;GeV",100,0,1000);
  TH1F* successRateHist = new TH1F(successRateHistName.c_str(),successRateHistTitle.c_str(),100,0,1);
  TH1D* successIhHist = new TH1D(successIhHistName.c_str(),successIhHistTitle.c_str(),400,0,10);
  TH1F* successRateIhHist = new TH1F(successRateIhHistName.c_str(),successRateIhHistTitle.c_str(),400,0,10);
  TH2F* successRate2DHist = new TH2F(successRate2DHistName.c_str(),successRate2DHistTitle.c_str(),100,0,1,400,0,10);
  TH2F* trialsHist = new TH2F(trialsHistName.c_str(),trialsHistTitle.c_str(),100,0,1,400,0,10);
  TH2F* selectedIhvsPHist = new TH2F(selectedIhvsPHistName.c_str(),selectedIhvsPHistTitle.c_str(),2000,0,2000,400,0,10);
  successRateHist->Sumw2();
  successRate2DHist->Sumw2();
  TFile* _file0 = new TFile(rootFilename);
  _file0->cd();
  const float k = 2.468;
  const float c = 2.679;

  string dedxHistName="trackDeDxE1VsDeDxD3LowPSBNoMSlice";
  dedxHistName+=intToString(nomSlice*2+1);
  TH3F* dedxHist = (TH3F*)_file0->Get(dedxHistName.c_str());
  string dedxHistName2="trackDeDxE1VsDeDxD3LowPSBNoMSlice";
  dedxHistName2+=intToString(nomSlice*2+2);
  TH3F* dedxHist2 = (TH3F*)_file0->Get(dedxHistName2.c_str());
  string pHistName="trackPMipSBNoMSlice";
  pHistName+=intToString(nomSlice*2+1);
  TH2F* pHist = (TH2F*)_file0->Get(pHistName.c_str());
  string pHistName2="trackPMipSBNoMSlice";
  pHistName2+=intToString(nomSlice*2+2);
  TH2F* pHist2 = (TH2F*)_file0->Get(pHistName2.c_str());
  ////XXX DEBUG
  //std::cout << dedxHistName<< "\t" << dedxHist->GetEntries() << std::endl;
  //std::cout << dedxHistName<< "\t" << dedxHist->Integral() << std::endl;
  //std::cout << dedxHistName2 << "\t" << dedxHist2->GetEntries() << std::endl;
  //std::cout << dedxHistName2 << "\t" << dedxHist2->Integral() << std::endl;
  //std::cout << pHistName<< "\t" << pHist->GetEntries() << std::endl;
  //std::cout << pHistName<< "\t" << pHist->Integral() << std::endl;
  //std::cout << pHistName2 << "\t" << pHist2->GetEntries() << std::endl;
  //std::cout << pHistName2 << "\t" << pHist2->Integral() << std::endl;
  //std::cout << "NOW ADDING THEM..." << std::endl;
  //std::cout << "NOW ADDING THEM..." << std::endl;
  ////XXX DEBUG
  dedxHist->Add(dedxHist2);
  pHist->Add(pHist2);
  // handle last NoM slice case (NoM>=23)
  if(nomSlice==11)
  {
    // in this case, we add NoM bins 23-30 together
    // already have 23-24, so add 25,26,27,28,29,30
    string dedxHistName25="trackDeDxE1VsDeDxD3LowPSBNoMSlice25";
    TH3F* dedxHist25 = (TH3F*)_file0->Get(dedxHistName25.c_str());
    string dedxHistName26="trackDeDxE1VsDeDxD3LowPSBNoMSlice26";
    TH3F* dedxHist26 = (TH3F*)_file0->Get(dedxHistName26.c_str());
    string dedxHistName27="trackDeDxE1VsDeDxD3LowPSBNoMSlice27";
    TH3F* dedxHist27 = (TH3F*)_file0->Get(dedxHistName27.c_str());
    string dedxHistName28="trackDeDxE1VsDeDxD3LowPSBNoMSlice28";
    TH3F* dedxHist28 = (TH3F*)_file0->Get(dedxHistName28.c_str());
    string dedxHistName29="trackDeDxE1VsDeDxD3LowPSBNoMSlice29";
    TH3F* dedxHist29 = (TH3F*)_file0->Get(dedxHistName29.c_str());
    string dedxHistName30="trackDeDxE1VsDeDxD3LowPSBNoMSlice30";
    TH3F* dedxHist30 = (TH3F*)_file0->Get(dedxHistName30.c_str());
    dedxHist->Add(dedxHist25);
    dedxHist->Add(dedxHist26);
    dedxHist->Add(dedxHist27);
    dedxHist->Add(dedxHist28);
    dedxHist->Add(dedxHist29);
    dedxHist->Add(dedxHist30);
    // now do the same for the P hist
    string pHistName25="trackPMipSBNoMSlice25";
    TH2F* pHist25 = (TH2F*)_file0->Get(pHistName25.c_str());
    string pHistName26="trackPMipSBNoMSlice26";
    TH2F* pHist26 = (TH2F*)_file0->Get(pHistName26.c_str());
    string pHistName27="trackPMipSBNoMSlice27";
    TH2F* pHist27 = (TH2F*)_file0->Get(pHistName27.c_str());
    string pHistName28="trackPMipSBNoMSlice28";
    TH2F* pHist28 = (TH2F*)_file0->Get(pHistName28.c_str());
    string pHistName29="trackPMipSBNoMSlice29";
    TH2F* pHist29 = (TH2F*)_file0->Get(pHistName29.c_str());
    string pHistName30="trackPMipSBNoMSlice30";
    TH2F* pHist30 = (TH2F*)_file0->Get(pHistName30.c_str());
    pHist->Add(pHist25);
    pHist->Add(pHist26);
    pHist->Add(pHist27);
    pHist->Add(pHist28);
    pHist->Add(pHist29);
    pHist->Add(pHist30);
  }
  ////XXX DEBUG
  //std::cout << dedxHistName<< "\t" << dedxHist->GetEntries() << std::endl;
  //std::cout << dedxHistName<< "\t" << dedxHist->Integral() << std::endl;
  //std::cout << dedxHistName2 << "\t" << dedxHist2->GetEntries() << std::endl;
  //std::cout << dedxHistName2 << "\t" << dedxHist2->Integral() << std::endl;
  //std::cout << pHistName<< "\t" << pHist->GetEntries() << std::endl;
  //std::cout << pHistName<< "\t" << pHist->Integral() << std::endl;
  //std::cout << pHistName2 << "\t" << pHist2->GetEntries() << std::endl;
  //std::cout << pHistName2 << "\t" << pHist2->Integral() << std::endl;
  ////XXX DEBUG

  myOutputFile->cd();
  string dirName="plotsNoMSlice";
  dirName+=intToString(nomSlice);
  dirName+="etaSlice";
  dirName+=intToString(etaSlice);
  dirName+="massCut";
  dirName+=intToString(massCut);
  TDirectory* thisEtaNoMSliceDir = myOutputFile->mkdir(dirName.c_str());
  thisEtaNoMSliceDir->cd();
  dedxHist->Write();
  _file0->cd();
  // set ranges to make the eta slice and project, dE/dx first
  dedxHist->GetZaxis()->SetRange(etaSlice*2+1,etaSlice*2+2);
  TH2F* etaSlicedNoMSlicedDeDxHist = (TH2F*) dedxHist->Project3D("yx");
  string histNameEtaSlicedDeDx="etaSlicedIhVsIasNoMSlice";
  histNameEtaSlicedDeDx+=intToString(nomSlice);
  histNameEtaSlicedDeDx+="etaSlice";
  histNameEtaSlicedDeDx+=intToString(etaSlice);
  string histTitleEtaSlicedDeDx="Ih vs. Ias from P SB for nom ";
  histTitleEtaSlicedDeDx+=intToString(nomLow);
  histTitleEtaSlicedDeDx+="-";
  histTitleEtaSlicedDeDx+=intToString(nomHigh);
  histTitleEtaSlicedDeDx+=", |#eta| ";
  histTitleEtaSlicedDeDx+=floatToString(etaLow);
  histTitleEtaSlicedDeDx+="-";
  histTitleEtaSlicedDeDx+=floatToString(etaHigh);
  etaSlicedNoMSlicedDeDxHist->SetNameTitle(histNameEtaSlicedDeDx.c_str(),histTitleEtaSlicedDeDx.c_str());
  // make the 1-D Ias for this eta/NoM (P SB)
  // first, only look at Ih > 3.5 MeV/cm, then project to Ias
  //XXX fixed possible bug here
  //etaSlicedNoMSlicedDeDxHist->GetYaxis()->SetRangeUser(3.5,10);
  etaSlicedNoMSlicedDeDxHist->GetYaxis()->SetRange(etaSlicedNoMSlicedDeDxHist->FindBin(3.5),etaSlicedNoMSlicedDeDxHist->GetNbinsY());
  TH1D* iasInThisEtaNoMBinBeforeMassCut = etaSlicedNoMSlicedDeDxHist->ProjectionX();
  string histNameIasEtaNoMBefMassCut="iasInNoMSliced";
  histNameIasEtaNoMBefMassCut+=intToString(nomSlice);
  histNameIasEtaNoMBefMassCut+="etaSlice";
  histNameIasEtaNoMBefMassCut+=intToString(etaSlice);
  histNameIasEtaNoMBefMassCut+="beforeMassCut";
  string histTitleIasEtaNoMBefMassCut="Ias from P SB for nom ";
  histTitleIasEtaNoMBefMassCut+=intToString(nomLow);
  histTitleIasEtaNoMBefMassCut+="-";
  histTitleIasEtaNoMBefMassCut+=intToString(nomHigh);
  histTitleIasEtaNoMBefMassCut+=", |#eta| ";
  histTitleIasEtaNoMBefMassCut+=floatToString(etaLow);
  histTitleIasEtaNoMBefMassCut+="-";
  histTitleIasEtaNoMBefMassCut+=floatToString(etaHigh);
  iasInThisEtaNoMBinBeforeMassCut->SetNameTitle(histNameIasEtaNoMBefMassCut.c_str(), histTitleIasEtaNoMBefMassCut.c_str());
  // same for ih
  TH1D* ihInThisEtaNoMBinBeforeMassCut = etaSlicedNoMSlicedDeDxHist->ProjectionY();
  string histNameIhEtaNoMBefMassCut="ihInNoMSliced";
  histNameIhEtaNoMBefMassCut+=intToString(nomSlice);
  histNameIhEtaNoMBefMassCut+="etaSlice";
  histNameIhEtaNoMBefMassCut+=intToString(etaSlice);
  histNameIhEtaNoMBefMassCut+="beforeMassCut";
  string histTitleIhEtaNoMBefMassCut="Ih ( > 3.5MeV/cm) from P SB for nom ";
  histTitleIhEtaNoMBefMassCut+=intToString(nomLow);
  histTitleIhEtaNoMBefMassCut+="-";
  histTitleIhEtaNoMBefMassCut+=intToString(nomHigh);
  histTitleIhEtaNoMBefMassCut+=", |#eta| ";
  histTitleIhEtaNoMBefMassCut+=floatToString(etaLow);
  histTitleIhEtaNoMBefMassCut+="-";
  histTitleIhEtaNoMBefMassCut+=floatToString(etaHigh);
  ihInThisEtaNoMBinBeforeMassCut->SetNameTitle(histNameIhEtaNoMBefMassCut.c_str(), histTitleIhEtaNoMBefMassCut.c_str());
  // p dist
  // set P range to be 100 GeV first
  pHist->GetXaxis()->SetRangeUser(100,2000);
  TH1D* etaSlicedNoMSlicedPHist = pHist->ProjectionX("etaSlicedNoMSlicedPHist_px",etaSlice*2+1,etaSlice*2+2);
  string histNameEtaSlicedP="pNoMSlice";
  histNameEtaSlicedP+=intToString(nomSlice);
  histNameEtaSlicedP+="etaSlice";
  histNameEtaSlicedP+=intToString(etaSlice);
  histNameEtaSlicedP+="massCut";
  histNameEtaSlicedP+=intToString(massCut);
  string histTitleEtaSlicedP="P from Ih SB with mass > ";
  histTitleEtaSlicedP+=floatToString(massCut);
  histTitleEtaSlicedP+=" for nom ";
  histTitleEtaSlicedP+=intToString(nomLow);
  histTitleEtaSlicedP+="-";
  histTitleEtaSlicedP+=intToString(nomHigh);
  histTitleEtaSlicedP+=", |#eta| ";
  histTitleEtaSlicedP+=floatToString(etaLow);
  histTitleEtaSlicedP+="-";
  histTitleEtaSlicedP+=floatToString(etaHigh);
  etaSlicedNoMSlicedPHist->SetNameTitle(histNameEtaSlicedP.c_str(),histTitleEtaSlicedP.c_str());

  // check to make sure everything went OK (not all hists checked for NoM>=23 case)
  if(etaSlicedNoMSlicedPHist->GetEntries() <= 0 || etaSlicedNoMSlicedDeDxHist->GetEntries() <= 0
      || pHist->GetEntries() <= 0 || dedxHist->GetEntries() <= 0
      || pHist2->GetEntries() <= 0 || dedxHist2->GetEntries() <= 0)
  {
    cout << "ERROR in histogram: " << endl
      << "\t etaSlicedNoMSlicedPHist\t\t" << etaSlicedNoMSlicedPHist->GetEntries() << endl
      << "\t etaSlicedNoMSlicedDeDxHist\t\t" << etaSlicedNoMSlicedDeDxHist->GetEntries() << endl
      << "\t " << pHistName << "\t\t\t" << pHist->GetEntries() << endl
      << "\t " << pHistName2 << "\t\t\t" << pHist2->GetEntries() << endl
      << "\t " << dedxHistName << "\t" << dedxHist->GetEntries() << endl
      << "\t " << dedxHistName2 << "\t" << dedxHist2->GetEntries() << endl;
    return -1;
  }

  TH1D* fluctuatedPHist = (TH1D*)etaSlicedNoMSlicedPHist->Clone();
  fluctuatedPHist->Reset();
  fluctuatedPHist->FillRandom(etaSlicedNoMSlicedPHist,(int)etaSlicedNoMSlicedPHist->Integral());
  //int numTrialPerIasBin = etaSlicedNoMSlicedPHist->GetEntries();
  //int numTrialPerIasBin = fluctuatedPHist->GetEntries();
  int numTrialPerIasBin = (int)fluctuatedPHist->Integral();
  // setup the arrays for keeping track of the pass/fail in each Ias bin
  const int numIasBins = 100;
  const int numIhBins = 250;
  int successCount[100][250];
  int trialsCount[100][250];
  for(int i=0; i<numIasBins;++i)
  {
    for(int j=0; j<numIhBins; ++j)
    {
      successCount[i][j] = 0;
      trialsCount[i][j] = 0;
    }
  }

  thisEtaNoMSliceDir->cd();
  TDirectory* ihDir = thisEtaNoMSliceDir->mkdir("ihDistsInIasSlices");
  ihDir->cd();

  for(int i=1; i<numIasBins+1; ++i)
  {
    // make 1-D Ih dist for this eta/NoM/Ias
    TH1D* ihForThisIasHist = etaSlicedNoMSlicedDeDxHist->ProjectionY(("ihForThisIasBin"+intToString(i)).c_str(),i,i);
    ihForThisIasHist->Write();
    if(ihForThisIasHist->GetEntries() < 1)
      continue;
    TH1D* fluctuatedIhHist = (TH1D*)ihForThisIasHist->Clone();
    fluctuatedIhHist->Reset();
    fluctuatedIhHist->FillRandom(ihForThisIasHist, (int)ihForThisIasHist->Integral());
    //XXX possible bug fix in range here
    //fluctuatedIhHist->GetXaxis()->SetRangeUser(3.5,25);
    fluctuatedIhHist->GetXaxis()->SetRange(fluctuatedIhHist->FindBin(3.5),fluctuatedIhHist->GetNbinsX());

    fluctuatedIhHist->SetName(("fluctuatedIhForIasBin"+intToString(i)).c_str());
    fluctuatedIhHist->Write();
    if(fluctuatedIhHist->GetEntries() < 1)
      continue;

    //debug
    cout << "Ias: bin=" << i << " entries=" << fluctuatedIhHist->GetEntries() <<
      " mean ih=" << fluctuatedIhHist->GetMean() << endl;

    int successesThisBin = 0;
    int trialsThisBin = 0;
    // generate the Ias dist, keeping track of pass/fail for each Ias bin
    for(int trial=1; trial < numTrialPerIasBin+1; ++trial)
    {
      //double ih = ihForThisIasHist->GetRandom();
      //double p = etaSlicedNoMSlicedPHist->GetRandom();
      //TODO for now, only using mean of the Ih
      double ih = fluctuatedIhHist->GetMean();
      double p = fluctuatedPHist->GetRandom();

      // require the candidate to be in the high Ih/high P region
      if(ih < 3.5 || p < 100)
        cout << "Ias: bin=" << i << " ERROR IN Ih/P range from hist: p=" << p << " mean ih=" << ih << endl;
      //if(ih < 3.5) continue;
      if(ih < 3.5) break;
      if(p < 100) continue;

      double massSqr = (ih-c)*pow(p,2)/k;
      trialsCount[i-1][ihForThisIasHist->FindBin(ih)-1]++;
      //debug
      //std::cout << "Generated points--  Ih=" << ih << "; p=" << p << "; mass=" << sqrt(massSqr) << std::endl;
      trialsThisBin++;

      if(massSqr < 0)
        continue;
      else if(sqrt(massSqr) < massCut)
        continue;
      
      successesThisBin++;
      massGeneratedHist->Fill(sqrt(massSqr));
      selectedIhvsPHist->Fill(p,ih);
      ////debug
      //std::cout << "Generated points--  Ih=" << ih << "; p=" << p << "; mass=" << sqrt(massSqr) << std::endl;
      //std::cout << "Ih bin: " << ihForThisIasHist->FindBin(ih) << std::endl;
      successCount[i-1][ihForThisIasHist->FindBin(ih)-1]++;
    }

    delete ihForThisIasHist;
    delete fluctuatedIhHist;
    double ias = successRateHist->GetBinCenter(i);
    double successRate = successesThisBin/(double)trialsThisBin;
    if(trialsThisBin > 0)
    {
      successRateHist->SetBinContent(i,successRate);
      successRateHist->SetBinError(i,sqrt(successRate*(1-successRate)/trialsThisBin));
    }
    //debug
    cout << "Ias: bin=" << i << " ias=" << ias
      << " successRate=" << successRate << " trials="
      << trialsThisBin << " error=" << sqrt(successRate*(1-successRate)/trialsThisBin)
      << endl;
  }

  for(int i=0; i<numIasBins;++i)
  {
    double weightSum = 0;
    double weightSqrSum = 0;
    double weightedAvg = 0;
    double errorAvg = 0;
    for(int j=0; j<numIhBins; ++j)
    {
      int thisBinTrialsCount = trialsCount[i][j];
      trialsHist->SetBinContent(i+1,j+1,thisBinTrialsCount);
      //double successRate = successCount[i][j]/(double)numTrialPerIasBin;
      // use binomial errors ( sigma^2 = n*p*q = n*p*(1-p) )
      //double binErr = sqrt(numTrialPerIasBin*successRate*(1-successRate))/(double)numTrialPerIasBin;
      if(thisBinTrialsCount>0)
      {
        double successRate = successCount[i][j]/(double)thisBinTrialsCount;
        successRate2DHist->SetBinContent(i+1,j+1,successRate);
        //double binErr = sqrt(thisBinTrialsCount*successRate*(1-successRate))/(double)thisBinTrialsCount;
        //double stdDev = sqrt(thisBinTrialsCount*successRate*(1-successRate));
        double stdDev = sqrt(successRate*(1-successRate));
        double binErr = stdDev/sqrt(thisBinTrialsCount);
        if(binErr==0)
          successRate2DHist->SetBinError(i+1,j+1,sqrt(thisBinTrialsCount));
        else 
          successRate2DHist->SetBinError(i+1,j+1,binErr);
      }
      ////debug
      ////if(thisBinTrialsCount>0)
      //if(i+1==12)
      //{
      //  cout << "binx: " << i+1 << " biny: " << j+1 << " successRate: " << successRate << " error: " << binErr
      //    << " successCount: " << successCount[i][j] << " thisBinTrialsCount: " << thisBinTrialsCount << endl;
      //}
      // create weighted average for this Ias slice
      //weightSum+=thisBinTrialsCount;
      //weightSqrSum+=thisBinTrialsCount*thisBinTrialsCount;
      //weightedAvg+=successRate*thisBinTrialsCount;
      //errorAvg+=binErr*binErr*thisBinTrialsCount*thisBinTrialsCount;
    }
    //if(weightSum>0)
    //{
    //  successRateHist->SetBinContent(i+1,weightedAvg/weightSum);
    //  double stdDev = sqrt(errorAvg/weightSqrSum);
    //  successRateHist->SetBinError(i+1,stdDev/sqrt(weightSum*weightSum/weightSqrSum));
    //  //successRateHist->SetBinContent(i+1,weightedAvg/weightSum);
    //  //TODO: calculate error
    //  //successRateHist->SetBinError(i+1,errorAvg/weightSum);
    //  //if(i+1==12)
    //  //  cout << "binx: " << i+1 << " weightedAvg= " << weightedAvg/weightSum << std::endl;
    //}
  }
  moveBinsTH2F(successRate2DHist);
  successRate2DHist->SetMinimum(-0.01);

  // do Ih with mass cut
  for(int bin=0; bin < numIhBins; ++bin)
  {
    int trialsThisBin = 0;
    int successesInThisBin = 0;
    if(ihInThisEtaNoMBinBeforeMassCut->GetBinContent(bin+1) < 1) continue;
    double ih = ihInThisEtaNoMBinBeforeMassCut->GetBinCenter(bin+1);
    if(ih < 3.5) continue;

    for(int trial=0; trial < (int)etaSlicedNoMSlicedPHist->Integral();
        ++trial)
    {
      //double ih = ihInThisEtaNoMBinBeforeMassCut->GetRandom();
      double p = etaSlicedNoMSlicedPHist->GetRandom();
      double massSqr = (ih-c)*pow(p,2)/k;
      // require the candidate to be in the high Ih/high P region
      if(ih < 3.5 || p < 100)
        cout << "Ih: ERROR IN Ih/P range from hist: p=" << p << " ih=" << ih << endl;
      if(ih < 3.5) continue;
      if(p < 100) continue;

      trialsThisBin++;

      if(massSqr < 0)
        continue;
      else if(sqrt(massSqr) < massCut)
        continue;

      successesInThisBin++;
    }
    double successRate = successesInThisBin/(double)trialsThisBin;
    if(trialsThisBin > 0)
    {
      successRateIhHist->SetBinContent(successRateIhHist->FindBin(ih),successRate);
      successRateIhHist->SetBinError(successRateIhHist->FindBin(ih),sqrt(successRate*(1-successRate)/trialsThisBin));
    }
    //debug
    //cout << "Ih: bin=" << bin << " ih=" << ihInThisEtaNoMBinBeforeMassCut->GetBinCenter(bin+1)
    //  << " successRate=" << successRate << " trials="
    //  << trialsThisBin << " error=" << sqrt(successRate*(1-successRate)/trialsThisBin)
    //  << endl;
  }

  successIhHist = ihInThisEtaNoMBinBeforeMassCut;
  successIhHist->Sumw2();
  successIhHist->Multiply(successRateIhHist);

  std::vector<double> binVec;
  //float ihBinArray[251]; // 250, 0 - 25
  //float etaBinArray[25]; // 24, 0 - 2.4
  //for(int i=0; i<251; ++i)
  //  ihBinArray[i] = i*0.1;
  //for(int i=0; i<25; ++i)
  //  etaBinArray[i] = i*0.1;
  //do variable bins
  double lastBinLowEdge = getLastBinLowEdge(nomLow);
  //float secondToLastBinLowEdge = getSecondToLastBinLowEdge(nomLow);
  //int numBins = floor(secondToLastBinLowEdge/0.05);
  int numBins = floor(lastBinLowEdge/0.05);
  for(int i=0; i<numBins; ++i)
    binVec.push_back(i*0.05);
  //binVec.push_back(secondToLastBinLowEdge);
  binVec.push_back(lastBinLowEdge);
  binVec.push_back(1);
  //std::cout << "numBins = " << numBins << " NoM = " << nom << std::endl;
  //for(unsigned int i=0; i<binVec.size(); ++i)
  //  std::cout << "i=" << i << " binVec[i] = " << binVec[i] << " binVec size = " << binVec.size() << std::endl;
  //trackDeDxE1VsDeDxD3LowPSBNoMSliceHists[nom-1] = fs.make<TH3F>(dedxName.c_str(),dedxTitle.c_str(),binVec.size()-1,&(*binVec.begin()),250,ihBinArray,24,etaBinArray);

  myOutputFile->cd();
  TH1D* tempHist = new TH1D(predHistFixedBinsName.c_str(),predHistFixedBinsTitle.c_str(),100,0,1);
  //TH1D* trackIasFactorizedIhPHist = new TH1D(predHistName.c_str(),predHistTitle.c_str(),100,0,1);
  //trackIasFactorizedIhPHist->Sumw2();
  tempHist->Sumw2();
  // the below uses uncorrelated errors between the two hists
  tempHist->Multiply(iasInThisEtaNoMBinBeforeMassCut,successRateHist);
  TH1D* trackIasFactorizedIhPHist;
  trackIasFactorizedIhPHist = (TH1D*)tempHist->Rebin(binVec.size()-1,predHistName.c_str(),&(*binVec.begin()));
  trackIasFactorizedIhPHist->SetTitle(predHistTitle.c_str());
  //debug
  std::cout << "iasInThisEtaNoMBinBeforeMassCut entries=" <<
    iasInThisEtaNoMBinBeforeMassCut->GetEntries() <<
    " binsX: " << iasInThisEtaNoMBinBeforeMassCut->GetNbinsX() <<
    " successRateHist entries=" <<
    successRateHist->GetEntries() <<
    " and final trackIasFactorizedIhPHist entries=" <<
    trackIasFactorizedIhPHist->GetEntries() <<
    " binsX: " << trackIasFactorizedIhPHist->GetNbinsX() <<
    std::endl;
  
  //std::cout << "Generated " << numEntries << " points using " << trialsCount << " trials: "
  //  << 100*(float)(trialsCount-numEntries)/trialsCount << "\% loss" << std::endl;

  //trackIasFactorizedIhPHist->Draw();
  //TCanvas t;
  //t.cd();
  //t.SetLogy();
  //trackIasFactorizedIhPHist->Draw();
  //t.Print((histName+".C").c_str());
  //t.Print((histName+".png").c_str());

  thisEtaNoMSliceDir->cd();
  ihInThisEtaNoMBinBeforeMassCut->Write();
  successIhHist->Write();
  successRateIhHist->Write();
  tempHist->Write();
  trialsHist->Write();
  trackIasFactorizedIhPHist->Write();
  pHist->Write();
  etaSlicedNoMSlicedPHist->Write();
  etaSlicedNoMSlicedDeDxHist->Write();
  massGeneratedHist->Write();
  iasInThisEtaNoMBinBeforeMassCut->Write();
  successRateHist->Write();
  successRate2DHist->Write();
  selectedIhvsPHist->Write();
  myOutputFile->Close();

}
