#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TF1.h"
#include "TProfile.h"
#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "HFitInterface.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <map>


//pasted in the example code from the combinedFit.C tutorial online
//IGNORE SLICE 8 FOR NOW (STRANGE SHAPE)
int iparEta1[5] = { 
  0, // x^0 common
  1, // x^1 common
  2, // x^2 common
  3, // x^3 common
  4  // x^4 common
};
int ipar8EtaSlices[7][6]; 
int ipar7EtaSlices[6][6]; 

struct GlobalChi2For8EtaSlices
{ 
  GlobalChi2For8EtaSlices(ROOT::Math::IMultiGenFunction& f1, ROOT::Math::IMultiGenFunction& f2,
      ROOT::Math::IMultiGenFunction& f3, ROOT::Math::IMultiGenFunction& f4,
      ROOT::Math::IMultiGenFunction& f5, ROOT::Math::IMultiGenFunction& f6,
      ROOT::Math::IMultiGenFunction& f7, ROOT::Math::IMultiGenFunction& f8
      ) : 
    fChi2_1(&f1),fChi2_2(&f2),fChi2_3(&f3),fChi2_4(&f4),fChi2_5(&f5),fChi2_6(&f6),
    fChi2_7(&f7),fChi2_8(&f8) {}

  double operator() (const double *par) const {
    double p1[5];
    double pOther[7][6];
    for(int i=0; i<5; ++i)
    {
      p1[i] = par[iparEta1[i]];
      for(int j=0; j<7; ++j)
        pOther[j][i] = par[ipar8EtaSlices[j][i]];
    }
    for(int i=0; i<7; ++i)
      pOther[i][5] = par[ipar8EtaSlices[i][5]];

    return (*fChi2_1)(p1) + (*fChi2_2)(pOther[0]) + (*fChi2_3)(pOther[1]) + (*fChi2_4)(pOther[2])
      + (*fChi2_5)(pOther[3]) + (*fChi2_6)(pOther[4]) + (*fChi2_7)(pOther[5]) + (*fChi2_8)(pOther[6]);
  } 

  const  ROOT::Math::IMultiGenFunction* fChi2_1;
  const  ROOT::Math::IMultiGenFunction* fChi2_2;
  const  ROOT::Math::IMultiGenFunction* fChi2_3;
  const  ROOT::Math::IMultiGenFunction* fChi2_4;
  const  ROOT::Math::IMultiGenFunction* fChi2_5;
  const  ROOT::Math::IMultiGenFunction* fChi2_6;
  const  ROOT::Math::IMultiGenFunction* fChi2_7;
  const  ROOT::Math::IMultiGenFunction* fChi2_8;
};

struct GlobalChi2For7EtaSlices
{ 
  GlobalChi2For7EtaSlices(ROOT::Math::IMultiGenFunction& f1, ROOT::Math::IMultiGenFunction& f2,
      ROOT::Math::IMultiGenFunction& f3, ROOT::Math::IMultiGenFunction& f4,
      ROOT::Math::IMultiGenFunction& f5, ROOT::Math::IMultiGenFunction& f6,
      ROOT::Math::IMultiGenFunction& f7
      ) : 
    fChi2_1(&f1),fChi2_2(&f2),fChi2_3(&f3),fChi2_4(&f4),fChi2_5(&f5),fChi2_6(&f6),
    fChi2_7(&f7){}

  double operator() (const double *par) const {
    double p1[5];
    double pOther[6][6];
    for(int i=0; i<5; ++i)
    {
      p1[i] = par[iparEta1[i]];
      for(int j=0; j<6; ++j)
        pOther[j][i] = par[ipar8EtaSlices[j][i]];
    }
    for(int i=0; i<6; ++i)
      pOther[i][5] = par[ipar8EtaSlices[i][5]];

    return (*fChi2_1)(p1) + (*fChi2_2)(pOther[0]) + (*fChi2_3)(pOther[1]) + (*fChi2_4)(pOther[2])
      + (*fChi2_5)(pOther[3]) + (*fChi2_6)(pOther[4]) + (*fChi2_7)(pOther[5]);
  } 

  const  ROOT::Math::IMultiGenFunction* fChi2_1;
  const  ROOT::Math::IMultiGenFunction* fChi2_2;
  const  ROOT::Math::IMultiGenFunction* fChi2_3;
  const  ROOT::Math::IMultiGenFunction* fChi2_4;
  const  ROOT::Math::IMultiGenFunction* fChi2_5;
  const  ROOT::Math::IMultiGenFunction* fChi2_6;
  const  ROOT::Math::IMultiGenFunction* fChi2_7;
};

std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}


// ********* MAIN
int main(int argc, char** argv)
{
  using namespace std;

  char* inputFile = argv[1];
  if(!inputFile)
  {
    cout << "Error: Missing input file." << endl;
    return -1;
  }

  gROOT->Macro("rootstyle.C");

  TFile* inputTFile = new TFile(inputFile);
  inputTFile->cd();

  TFile* outputTFile = new TFile("plotDeDxPdep.root","recreate");

  int nPBins = 200;
  int nNomBins = 30;
  int nEtaBins = 25;

  vector<double> p0s;
  vector<double> p1s;
  vector<double> p0Es;
  vector<double> p1Es;
  vector<double> etas;
  vector<double> etaEs;
  vector<double> p0sForSigma;
  vector<double> p1sForSigma;
  vector<double> p0EsForSigma;
  vector<double> p1EsForSigma;
  vector<double> etasForSigma;
  vector<double> etaEsForSigma;
  TGraphErrors* meanDeDxVsBetaTGraphsEta[25];
  TGraphErrors* sigmaDeDxVsBetaTGraphsEta[25];
  TGraphErrors* residualDeDxVsMeanTGraphsEta[25];
  TGraphErrors* secondResidualDeDxVsMeanTGraphsEta[25];
  TH2F* dedxHistsInNoMAndEtaSlices[30][25];
  for(int i=0; i<nEtaBins; ++i)
  {
    TGraphErrors* t = new TGraphErrors();
    meanDeDxVsBetaTGraphsEta[i] = t;
    sigmaDeDxVsBetaTGraphsEta[i] = t;
    residualDeDxVsMeanTGraphsEta[i] = t;
    secondResidualDeDxVsMeanTGraphsEta[i] = t;
    for(int j=0; j<nNomBins; ++j)
    {
      dedxHistsInNoMAndEtaSlices[j][i] = new TH2F();
    }
  }

  // fit func for beta*gamma dependence
  TF1* myFitFunc = new TF1("myFitFunc","[0]+[1]*((1-x^2)/x^2)");
  //TF1* myFitFunc = new TF1("myFitFunc","[0]*(1+((1-x^2)/x^2))"); // K/C == 1
  // fit for the residuals of the beta*gamma dependence of the mean in combined eta slices
  TF1* dedxFitMeanResidualsFuncEta1 = new TF1("dedxFitMeanResidualsFuncEta1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta2 = new TF1("dedxFitMeanResidualsFuncEta2","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta3 = new TF1("dedxFitMeanResidualsFuncEta3","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta4 = new TF1("dedxFitMeanResidualsFuncEta4","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta5 = new TF1("dedxFitMeanResidualsFuncEta5","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta6 = new TF1("dedxFitMeanResidualsFuncEta6","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta7 = new TF1("dedxFitMeanResidualsFuncEta7","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  //TF1* dedxFitMeanResidualsFuncEta8 = new TF1("dedxFitMeanResidualsFuncEta8","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*)x^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta9 = new TF1("dedxFitMeanResidualsFuncEta9","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta10 = new TF1("dedxFitMeanResidualsFuncEta10","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta11 = new TF1("dedxFitMeanResidualsFuncEta11","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta12 = new TF1("dedxFitMeanResidualsFuncEta12","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta13 = new TF1("dedxFitMeanResidualsFuncEta13","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta14 = new TF1("dedxFitMeanResidualsFuncEta14","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta15 = new TF1("dedxFitMeanResidualsFuncEta15","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta16 = new TF1("dedxFitMeanResidualsFuncEta16","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta17 = new TF1("dedxFitMeanResidualsFuncEta17","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta18 = new TF1("dedxFitMeanResidualsFuncEta18","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta19 = new TF1("dedxFitMeanResidualsFuncEta19","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta20 = new TF1("dedxFitMeanResidualsFuncEta20","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta21 = new TF1("dedxFitMeanResidualsFuncEta21","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta22 = new TF1("dedxFitMeanResidualsFuncEta22","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  TF1* dedxFitMeanResidualsFuncEta23 = new TF1("dedxFitMeanResidualsFuncEta23","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  //TF1* dedxFitMeanResidualsFuncEta23 = new TF1("dedxFitMeanResidualsFuncEta23","[0]+[1]*[5]*x+[2]*([5]*x)^2+[3]*([5]*x)^3+[4]*([5]*x)^4",3,25);
  ROOT::Math::WrappedMultiTF1 wmFuncEta1(*dedxFitMeanResidualsFuncEta1,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta2(*dedxFitMeanResidualsFuncEta2,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta3(*dedxFitMeanResidualsFuncEta3,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta4(*dedxFitMeanResidualsFuncEta4,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta5(*dedxFitMeanResidualsFuncEta5,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta6(*dedxFitMeanResidualsFuncEta6,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta7(*dedxFitMeanResidualsFuncEta7,1);
  //ROOT::Math::WrappedMultiTF1 wmFuncEta8(*dedxFitMeanResidualsFuncEta8,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta9(*dedxFitMeanResidualsFuncEta9,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta10(*dedxFitMeanResidualsFuncEta10,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta11(*dedxFitMeanResidualsFuncEta11,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta12(*dedxFitMeanResidualsFuncEta12,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta13(*dedxFitMeanResidualsFuncEta13,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta14(*dedxFitMeanResidualsFuncEta14,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta15(*dedxFitMeanResidualsFuncEta15,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta16(*dedxFitMeanResidualsFuncEta16,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta17(*dedxFitMeanResidualsFuncEta17,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta18(*dedxFitMeanResidualsFuncEta18,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta19(*dedxFitMeanResidualsFuncEta19,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta20(*dedxFitMeanResidualsFuncEta20,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta21(*dedxFitMeanResidualsFuncEta21,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta22(*dedxFitMeanResidualsFuncEta22,1);
  ROOT::Math::WrappedMultiTF1 wmFuncEta23(*dedxFitMeanResidualsFuncEta23,1);
  //init parameter index arrays
  for(int i=0; i<7; ++i)
  {
    for(int j=0; j<5; ++j)
      ipar8EtaSlices[i][j] = j;
  }
  for(int j=0; j<7; ++j)
    ipar8EtaSlices[j][5] = j+5;
  for(int i=0; i<6; ++i)
  {
    for(int j=0; j<5; ++j)
      ipar7EtaSlices[i][j] = j;
  }
  for(int j=0; j<6; ++j)
    ipar7EtaSlices[j][5] = j+5;
  //test
  //for(int i=0; i<5; ++i)
  //{
  //  cout << "p1[" << i << "] = " << iparEta1[i] << endl;
  //}
  //for(int i=0; i<6; ++i)
  //{
  //  for(int j=0; j<21; ++j)
  //    cout << "pOther[" << j << "][" << i<< "] = " << iparEtaSlices[j][i] << endl;
  //}

  // fitted beta*gamma dependence function in eta slices, with parameters
  // see plots in hscp/aug11 web directory; pcminn19:/data3/scooper/cmssw/426/HSCPStudiesAndToyMC/src/HSCP2011/ToyMC/bin/1.aug10.dedxPlots.aug11
  TF1* dedxMeanBetaDepFunc = new TF1("dedxMeanBetaDepFunc","[0]+[1]*((1-x^2)/x^2)",0.4,1);
  double dedxMeanBetaDepEtaParam0[25];
  double dedxMeanBetaDepEtaParam1[25];
  dedxMeanBetaDepEtaParam0[0] = 2.412;
  dedxMeanBetaDepEtaParam1[0] = 2.348;
  dedxMeanBetaDepEtaParam0[1] = 2.366; 
  dedxMeanBetaDepEtaParam1[1] = 2.393;
  dedxMeanBetaDepEtaParam0[2] = 2.406;
  dedxMeanBetaDepEtaParam1[2] = 2.32;
  dedxMeanBetaDepEtaParam0[3] = 2.46;
  dedxMeanBetaDepEtaParam1[3] = 2.23;
  dedxMeanBetaDepEtaParam0[4] = 2.351;
  dedxMeanBetaDepEtaParam1[4] = 2.354;
  dedxMeanBetaDepEtaParam0[5] = 2.291;
  dedxMeanBetaDepEtaParam1[5] = 2.498;
  dedxMeanBetaDepEtaParam0[6] = 2.215;
  dedxMeanBetaDepEtaParam1[6] = 2.645;
  dedxMeanBetaDepEtaParam0[7] = 2.776;
  dedxMeanBetaDepEtaParam1[7] = 1.46;
  dedxMeanBetaDepEtaParam0[8] = 2.292;
  dedxMeanBetaDepEtaParam1[8] = 2.518;
  dedxMeanBetaDepEtaParam0[9] = 2.574;
  dedxMeanBetaDepEtaParam1[9] = 1.913;
  dedxMeanBetaDepEtaParam0[10] = 2.672;
  dedxMeanBetaDepEtaParam1[10] = 1.701;
  dedxMeanBetaDepEtaParam0[11] = 2.483;
  dedxMeanBetaDepEtaParam1[11] = 2.043;
  dedxMeanBetaDepEtaParam0[12] = 2.54;
  dedxMeanBetaDepEtaParam1[12] = 1.816;
  dedxMeanBetaDepEtaParam0[13] = 2.472;
  dedxMeanBetaDepEtaParam1[13] = 1.974;
  dedxMeanBetaDepEtaParam0[14] = 2.469;
  dedxMeanBetaDepEtaParam1[14] = 1.921;
  dedxMeanBetaDepEtaParam0[15] = 2.382;
  dedxMeanBetaDepEtaParam1[15] = 2.118;
  dedxMeanBetaDepEtaParam0[16] = 2.402;
  dedxMeanBetaDepEtaParam1[16] = 2.028;
  dedxMeanBetaDepEtaParam0[17] = 2.481;
  dedxMeanBetaDepEtaParam1[17] = 1.828;
  dedxMeanBetaDepEtaParam0[18] = 2.427;
  dedxMeanBetaDepEtaParam1[18] = 1.936;
  dedxMeanBetaDepEtaParam0[19] = 2.541;
  dedxMeanBetaDepEtaParam1[19] = 1.76;
  dedxMeanBetaDepEtaParam0[20] = 2.593;
  dedxMeanBetaDepEtaParam1[20] = 1.742;
  dedxMeanBetaDepEtaParam0[21] = 2.604;
  dedxMeanBetaDepEtaParam1[21] = 1.815;
  dedxMeanBetaDepEtaParam0[22] = 2.683;
  dedxMeanBetaDepEtaParam1[22] = 1.695;
  dedxMeanBetaDepEtaParam0[23] = 2.789;
  dedxMeanBetaDepEtaParam1[23] = 1.455;
  dedxMeanBetaDepEtaParam0[24] = 2.664;
  dedxMeanBetaDepEtaParam1[24] = 1.902;

  TH1D* myHist;
  TH3F* rawHist;
  TH2F* dedxHistsInEtaBins[25]; // dE/dx vs. P


  string dir = "dedxHistsInNomBinsMatched";
  // loop over eta bins
  for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
  {
    inputTFile->cd();
    string name = "dedxE1MatchedNomBin";
    name+=intToString(1);

    string fullpath = dir;
    fullpath+="/";
    fullpath+=name;
    rawHist = (TH3F*) inputTFile->Get(fullpath.c_str());
    rawHist->GetZaxis()->SetRange(etaBin,etaBin);
    string projName = "yx";
    projName+="_etaBin";
    projName+=intToString(etaBin);
    TH2F* myHist2d = (TH2F*)rawHist->Project3D(projName.c_str());
    for(int nomBin=2; nomBin < nNomBins+1; ++nomBin)
    {
      string nameNom = "dedxE1MatchedNomBin";
      nameNom+=intToString(nomBin);
      string fullpathNom = dir;
      fullpathNom+="/";
      fullpathNom+=nameNom;
      rawHist = (TH3F*) inputTFile->Get(fullpathNom.c_str());
      rawHist->GetZaxis()->SetRange(etaBin,etaBin);
      string projNameNom = "yx";
      projNameNom+="_etaBin";
      projNameNom+=intToString(etaBin);
      projNameNom+="_nomBin";
      projNameNom+=intToString(nomBin);
      TH2F* thisHist = (TH2F*)rawHist->Project3D(projNameNom.c_str());
      //thisHist->SetName(nameNom.c_str());
      myHist2d->Add(thisHist);
      string histName = "dedxInNomBin";
      histName+=intToString(nomBin);
      histName+="EtaBin";
      histName+=intToString(etaBin);
      thisHist->SetName(histName.c_str());
      dedxHistsInNoMAndEtaSlices[nomBin-1][etaBin-1] = thisHist;
    }
    dedxHistsInEtaBins[etaBin-1] = myHist2d;
    // make output dir
    outputTFile->cd();
    string name2 = "dedxHistsIn";
    name2+="etaBin";
    name2+=intToString(etaBin);
    TDirectory* dedxE1MatchedEtaBinDir = outputTFile->mkdir(name2.c_str());
    dedxE1MatchedEtaBinDir->cd();
    dedxHistsInEtaBins[etaBin-1]->Write();

    vector<double> betas;
    vector<double> betaEs;
    vector<double> means;
    vector<double> sigmas;
    vector<double> meanEs;
    vector<double> sigmaEs;
    // Make graphs of dE/dx vs. beta in eta slices
    int minPbin = 40;
    //int minPbin = 135; // higher p (beta) only
    //if(etaBin < 6)
    //  minPbin = 100;
    for(int pBin=minPbin; pBin<nPBins; ++pBin)
    {
      string thisName2 = name2;
      thisName2+="pBin";
      thisName2+=intToString(pBin);
      //debug
      //cout << "Making projection, name=" << thisName2 << " pBin=" << pBin
      //  << " etaBin= " << etaBin << endl;
      myHist = dedxHistsInEtaBins[etaBin-1]->ProjectionX(thisName2.c_str(),pBin,pBin);

      if(myHist->GetEntries()>25)
      {
        TFitResultPtr r = myHist->Fit("gaus","SQ");
        //TFitResultPtr r = myHist->Fit("gaus","SQL");
        int ret = r;
        if(ret==0)
        {
          means.push_back(r->Parameter(1));
          sigmas.push_back(r->Parameter(2));
          meanEs.push_back(r->ParError(1));
          sigmaEs.push_back(r->ParError(2));
          float beta = (pBin*5+2.5)/sqrt(pow((pBin*5+2.5),2)+pow(500,2));
          betas.push_back(beta);
          float betaBinWidth = (pBin*5+5)/sqrt(pow((pBin*5+5),2)+pow(500,2))-(pBin*5)/sqrt(pow((pBin*5),2)+pow(500,2));
          betaEs.push_back(betaBinWidth/2);
          // write hist
          myHist->Write();
        }
      }
    }
    // make meanDeDx vs. beta graph for this eta, fit
    if(betas.size() <= 0)
      continue;
    TGraphErrors* t = new TGraphErrors(betas.size(),&(*betas.begin()),&(*means.begin()),&(*betaEs.begin()),&(*meanEs.begin()));
    TFitResultPtr r = t->Fit(myFitFunc,"SQ");
    int ret = r;
    if(ret!=0)
      cout << "Error in fitting 1-D hist from 2-d (no NoM slicing) mean; eta bin= " << etaBin << endl;
    else
    {
      etas.push_back(etaBin*0.1-0.05);
      etaEs.push_back(0.05);
      p0s.push_back(r->Parameter(0));
      p1s.push_back(r->Parameter(1));
      p0Es.push_back(r->ParError(0));
      p1Es.push_back(r->ParError(1));
    }
    // write hist
    string graphName = "meanDeDxVsBetaEtaBin";
    graphName+=intToString(etaBin);
    t->SetName(graphName.c_str());
    t->SetTitle(graphName.c_str());
    t->GetXaxis()->SetTitle("#beta");
    t->GetYaxis()->SetTitle("mean_Ih [MeV/cm]");
    t->Write();
    meanDeDxVsBetaTGraphsEta[etaBin-1] = t;

    // make residuals vs. expectedDeDx in this eta bin
    if(means.size() <= 0)
      continue;
    vector<double> residuals;
    vector<double> expDeDx;
    vector<double> expDeDxEs;
    vector<double>::const_iterator meanItr = means.begin();
    for(vector<double>::const_iterator betaItr = betas.begin();
        betaItr != betas.end(); ++betaItr)
    {
      dedxMeanBetaDepFunc->SetParameters(dedxMeanBetaDepEtaParam0[etaBin-1],dedxMeanBetaDepEtaParam1[etaBin-1]);
      double fitDeDx = dedxMeanBetaDepFunc->Eval(*betaItr);
      // SIC MOD over 5 only
      //if(fitDeDx >= 5)
      //{
        residuals.push_back(fitDeDx-*meanItr);
        expDeDx.push_back(fitDeDx);
        expDeDxEs.push_back(0);
      //}
      ////debug
      //if(etaBin==1)
      //cout << "For eta bin 1,  beta=" << *betaItr << " meanDeDx=" << *meanItr << " expectedDeDx=" << dedxMeanBetaDepFunc->Eval(*betaItr)
      //  << " difference = fit-mean = " << dedxMeanBetaDepFunc->Eval(*betaItr)-*meanItr
      //  << " param0=" << dedxMeanBetaDepFunc->GetParameter(0) << " param1=" << dedxMeanBetaDepFunc->GetParameter(1)
      //  << endl; 
      ++meanItr;
    }
    TCanvas tcan;
    tcan.cd();
    TGraphErrors* t3 = new TGraphErrors(expDeDx.size(),&(*expDeDx.begin()),&(*residuals.begin()),&(*expDeDxEs.begin()),&(*meanEs.begin()));
    // write hist
    string graphName3 = "residualDeDxVsMeanEtaBin";
    graphName3+=intToString(etaBin);
    t3->SetName(graphName3.c_str());
    t3->SetTitle(graphName3.c_str());
    t3->GetXaxis()->SetTitle("fit_Ih [MeV/cm]");
    t3->GetYaxis()->SetTitle("residual_Ih [MeV/cm]");
    t3->Draw("ap");
    //tcan.Print((graphName3+".C").c_str());
    t3->Write();
    residualDeDxVsMeanTGraphsEta[etaBin-1] = t3;

  }

  outputTFile->cd();

  // Make MultiGraphs of meanDeDx vs. beta (all eta slices together)
  TMultiGraph* meanDeDxVsBetaMG[6];
  TLegend* meanDeDxVsBetaLeg[6];
  for(int canv=0; canv<6; ++canv)
  {
    meanDeDxVsBetaMG[canv] = new TMultiGraph();
    meanDeDxVsBetaLeg[canv] = new TLegend(0.2,0.7,0.3,0.9);
    meanDeDxVsBetaLeg[canv]->SetHeader("#eta bin");
    int colorCounter = 1;
    for(int i=4*canv+1; i<4*canv+5; ++i)
    {
      if(meanDeDxVsBetaTGraphsEta[i-1]->GetN() > 0)
      {
        if(colorCounter==5)
          ++colorCounter;
        meanDeDxVsBetaTGraphsEta[i-1]->SetLineColor(colorCounter);
        meanDeDxVsBetaTGraphsEta[i-1]->SetMarkerColor(colorCounter);
        meanDeDxVsBetaMG[canv]->Add(meanDeDxVsBetaTGraphsEta[i-1],"p");
        meanDeDxVsBetaLeg[canv]->AddEntry(meanDeDxVsBetaTGraphsEta[i-1],intToString(i).c_str(),"l");
        ++colorCounter;
      }
    }
    if(canv==5)
    {
      meanDeDxVsBetaTGraphsEta[24]->SetLineColor(6);
      meanDeDxVsBetaTGraphsEta[24]->SetMarkerColor(6);
      meanDeDxVsBetaMG[canv]->Add(meanDeDxVsBetaTGraphsEta[24],"p");
      meanDeDxVsBetaLeg[canv]->AddEntry(meanDeDxVsBetaTGraphsEta[24],intToString(25).c_str(),"l");
    }
    TCanvas t;
    if(canv==0)
      t.SetName("meanDeDxVsBetaEtaSlices1_4");
    else if(canv==1)
      t.SetName("meanDeDxVsBetaEtaSlices5_8");
    else if(canv==2)
      t.SetName("meanDeDxVsBetaEtaSlices9_12");
    else if(canv==3)
      t.SetName("meanDeDxVsBetaEtaSlices13_16");
    else if(canv==4)
      t.SetName("meanDeDxVsBetaEtaSlices17_20");
    else if(canv==5)
      t.SetName("meanDeDxVsBetaEtaSlices21_25");
    meanDeDxVsBetaMG[canv]->Draw("a");
    meanDeDxVsBetaMG[canv]->GetXaxis()->SetTitle("#beta");
    meanDeDxVsBetaMG[canv]->GetYaxis()->SetTitle("mean_Ih [MeV/cm]");
    meanDeDxVsBetaLeg[canv]->Draw("same");
    t.Write();
  }

  TGraphErrors* p0VsEtaMeanTGraph = new TGraphErrors(etas.size(),&(*etas.begin()),&(*p0s.begin()),&(*etaEs.begin()),&(*p0Es.begin()));
  p0VsEtaMeanTGraph->SetName("p0VsEtaMean");
  p0VsEtaMeanTGraph->SetTitle("p0VsEta (mean)");
  p0VsEtaMeanTGraph->GetXaxis()->SetTitle("#eta");
  p0VsEtaMeanTGraph->GetYaxis()->SetTitle("p0");
  p0VsEtaMeanTGraph->Write();
  TGraphErrors* p1VsEtaMeanTGraph = new TGraphErrors(etas.size(),&(*etas.begin()),&(*p1s.begin()),&(*etaEs.begin()),&(*p1Es.begin()));
  p1VsEtaMeanTGraph->SetName("p1VsEtaMean");
  p1VsEtaMeanTGraph->SetTitle("p1VsEta (mean)");
  p1VsEtaMeanTGraph->GetXaxis()->SetTitle("#eta");
  p1VsEtaMeanTGraph->GetYaxis()->SetTitle("p1");
  p1VsEtaMeanTGraph->Write();

  // Set up fits for residualDeDx vs. meanDeDx
  ROOT::Fit::DataOptions dataOpt;
  ROOT::Fit::DataRange dataRange;
  dataRange.SetRange(3,25);
  ROOT::Fit::BinData* binData[23];
  for(int i=0; i<23; ++i)
  {
    if(i==7) continue;
    // set the data range
    binData[i] = new ROOT::Fit::BinData(dataOpt,dataRange);
    ROOT::Fit::FillData(*(binData[i]),residualDeDxVsMeanTGraphsEta[i]);
  }
  ROOT::Fit::Chi2Function chi2funcs_1(*(binData[0]),wmFuncEta1);
  ROOT::Fit::Chi2Function chi2funcs_2(*(binData[1]),wmFuncEta2);
  ROOT::Fit::Chi2Function chi2funcs_3(*(binData[2]),wmFuncEta3);
  ROOT::Fit::Chi2Function chi2funcs_4(*(binData[3]),wmFuncEta4);
  ROOT::Fit::Chi2Function chi2funcs_5(*(binData[4]),wmFuncEta5);
  ROOT::Fit::Chi2Function chi2funcs_6(*(binData[5]),wmFuncEta6);
  ROOT::Fit::Chi2Function chi2funcs_7(*(binData[6]),wmFuncEta7);
  //ROOT::Fit::Chi2Function chi2funcs_8(*(binData[7]),wmFuncEta8);
  ROOT::Fit::Chi2Function chi2funcs_9(*(binData[8]),wmFuncEta9);
  ROOT::Fit::Chi2Function chi2funcs_10(*(binData[9]),wmFuncEta10);
  ROOT::Fit::Chi2Function chi2funcs_11(*(binData[10]),wmFuncEta11);
  ROOT::Fit::Chi2Function chi2funcs_12(*(binData[11]),wmFuncEta12);
  ROOT::Fit::Chi2Function chi2funcs_13(*(binData[12]),wmFuncEta13);
  ROOT::Fit::Chi2Function chi2funcs_14(*(binData[13]),wmFuncEta14);
  ROOT::Fit::Chi2Function chi2funcs_15(*(binData[14]),wmFuncEta15);
  ROOT::Fit::Chi2Function chi2funcs_16(*(binData[15]),wmFuncEta16);
  ROOT::Fit::Chi2Function chi2funcs_17(*(binData[16]),wmFuncEta17);
  ROOT::Fit::Chi2Function chi2funcs_18(*(binData[17]),wmFuncEta18);
  ROOT::Fit::Chi2Function chi2funcs_19(*(binData[18]),wmFuncEta19);
  ROOT::Fit::Chi2Function chi2funcs_20(*(binData[19]),wmFuncEta20);
  ROOT::Fit::Chi2Function chi2funcs_21(*(binData[20]),wmFuncEta21);
  ROOT::Fit::Chi2Function chi2funcs_22(*(binData[21]),wmFuncEta22);
  ROOT::Fit::Chi2Function chi2funcs_23(*(binData[22]),wmFuncEta23);
  GlobalChi2For8EtaSlices globalChi2Eta1to9(chi2funcs_1,chi2funcs_2,chi2funcs_3,
      chi2funcs_4,chi2funcs_5,chi2funcs_6,chi2funcs_7,chi2funcs_9);
  GlobalChi2For7EtaSlices globalChi2Eta10to16(chi2funcs_10,chi2funcs_11,
      chi2funcs_12,chi2funcs_13,chi2funcs_14,chi2funcs_15,chi2funcs_16);
  GlobalChi2For7EtaSlices globalChi2Eta17to23(chi2funcs_17,chi2funcs_18,
      chi2funcs_19,chi2funcs_20,chi2funcs_21,chi2funcs_22,chi2funcs_23);

  ROOT::Fit::Fitter fitter;

  // First batch of eta slices
  const int Npar0 = 12;
  double par0[Npar0] = {3.1,-1.7,0.3,-0.02,0.0006,1,1,1,1,1,1,1};
  fitter.Config().SetParamsSettings(Npar0,par0);
  //fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(20000);
  fitter.Config().SetMinimizer("Minuit","Migrad"); 
  fitter.FitFCN(Npar0,globalChi2Eta1to9,par0);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  double residualChi2OverNdfs[25];
  dedxFitMeanResidualsFuncEta1->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4));
  dedxFitMeanResidualsFuncEta1->SetChisquare(residualDeDxVsMeanTGraphsEta[0]->Chisquare(dedxFitMeanResidualsFuncEta1));
  dedxFitMeanResidualsFuncEta1->SetNDF(residualDeDxVsMeanTGraphsEta[0]->GetN()-5);
  residualDeDxVsMeanTGraphsEta[0]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta1);
  residualChi2OverNdfs[0] = dedxFitMeanResidualsFuncEta1->GetChisquare()/dedxFitMeanResidualsFuncEta1->GetNDF();
  dedxFitMeanResidualsFuncEta2->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(5));
  dedxFitMeanResidualsFuncEta2->SetChisquare(residualDeDxVsMeanTGraphsEta[1]->Chisquare(dedxFitMeanResidualsFuncEta2));
  dedxFitMeanResidualsFuncEta2->SetNDF(residualDeDxVsMeanTGraphsEta[1]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[1]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta2);
  residualChi2OverNdfs[1] = dedxFitMeanResidualsFuncEta2->GetChisquare()/dedxFitMeanResidualsFuncEta2->GetNDF();
  dedxFitMeanResidualsFuncEta3->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(6));
  dedxFitMeanResidualsFuncEta3->SetChisquare(residualDeDxVsMeanTGraphsEta[2]->Chisquare(dedxFitMeanResidualsFuncEta3));
  dedxFitMeanResidualsFuncEta3->SetNDF(residualDeDxVsMeanTGraphsEta[2]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[2]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta3);
  residualChi2OverNdfs[2] = dedxFitMeanResidualsFuncEta3->GetChisquare()/dedxFitMeanResidualsFuncEta3->GetNDF();
  dedxFitMeanResidualsFuncEta4->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(7));
  dedxFitMeanResidualsFuncEta4->SetChisquare(residualDeDxVsMeanTGraphsEta[3]->Chisquare(dedxFitMeanResidualsFuncEta4));
  dedxFitMeanResidualsFuncEta4->SetNDF(residualDeDxVsMeanTGraphsEta[3]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[3]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta4);
  residualChi2OverNdfs[3] = dedxFitMeanResidualsFuncEta4->GetChisquare()/dedxFitMeanResidualsFuncEta4->GetNDF();
  dedxFitMeanResidualsFuncEta5->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(8));
  dedxFitMeanResidualsFuncEta5->SetChisquare(residualDeDxVsMeanTGraphsEta[4]->Chisquare(dedxFitMeanResidualsFuncEta5));
  dedxFitMeanResidualsFuncEta5->SetNDF(residualDeDxVsMeanTGraphsEta[4]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[4]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta5);
  residualChi2OverNdfs[4] = dedxFitMeanResidualsFuncEta5->GetChisquare()/dedxFitMeanResidualsFuncEta5->GetNDF();
  dedxFitMeanResidualsFuncEta6->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(9));
  dedxFitMeanResidualsFuncEta6->SetChisquare(residualDeDxVsMeanTGraphsEta[5]->Chisquare(dedxFitMeanResidualsFuncEta6));
  dedxFitMeanResidualsFuncEta6->SetNDF(residualDeDxVsMeanTGraphsEta[5]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[5]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta6);
  residualChi2OverNdfs[5] = dedxFitMeanResidualsFuncEta6->GetChisquare()/dedxFitMeanResidualsFuncEta6->GetNDF();
  dedxFitMeanResidualsFuncEta7->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(10));
  dedxFitMeanResidualsFuncEta7->SetChisquare(residualDeDxVsMeanTGraphsEta[6]->Chisquare(dedxFitMeanResidualsFuncEta7));
  dedxFitMeanResidualsFuncEta7->SetNDF(residualDeDxVsMeanTGraphsEta[6]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[6]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta7);
  residualChi2OverNdfs[6] = dedxFitMeanResidualsFuncEta7->GetChisquare()/dedxFitMeanResidualsFuncEta7->GetNDF();
  //residualDeDxVsMeanTGraphsEta[7]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta8);
  dedxFitMeanResidualsFuncEta9->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(11));
  dedxFitMeanResidualsFuncEta9->SetChisquare(residualDeDxVsMeanTGraphsEta[8]->Chisquare(dedxFitMeanResidualsFuncEta9));
  dedxFitMeanResidualsFuncEta9->SetNDF(residualDeDxVsMeanTGraphsEta[8]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[8]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta9);
  residualChi2OverNdfs[8] = dedxFitMeanResidualsFuncEta9->GetChisquare()/dedxFitMeanResidualsFuncEta9->GetNDF();

  // Second batch of eta slices
  const int Npar1 = 11;
  double par1[Npar1] = {-0.68,-0.23,0.081,-0.0035,9.6e-5,1,1,1,1,1,1};
  fitter.Config().SetParamsSettings(Npar1,par1);
  //fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(20000);
  fitter.Config().SetMinimizer("Minuit","Migrad"); 
  fitter.FitFCN(Npar1,globalChi2Eta10to16,par1);
  result = fitter.Result();
  result.Print(std::cout);

  dedxFitMeanResidualsFuncEta10->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4));
  dedxFitMeanResidualsFuncEta10->SetChisquare(residualDeDxVsMeanTGraphsEta[9]->Chisquare(dedxFitMeanResidualsFuncEta10));
  dedxFitMeanResidualsFuncEta10->SetNDF(residualDeDxVsMeanTGraphsEta[9]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[9]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta10);
  residualChi2OverNdfs[9] = dedxFitMeanResidualsFuncEta10->GetChisquare()/dedxFitMeanResidualsFuncEta10->GetNDF();
  dedxFitMeanResidualsFuncEta11->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(5));
  dedxFitMeanResidualsFuncEta11->SetChisquare(residualDeDxVsMeanTGraphsEta[10]->Chisquare(dedxFitMeanResidualsFuncEta11));
  dedxFitMeanResidualsFuncEta11->SetNDF(residualDeDxVsMeanTGraphsEta[10]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[10]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta11);
  residualChi2OverNdfs[10] = dedxFitMeanResidualsFuncEta11->GetChisquare()/dedxFitMeanResidualsFuncEta11->GetNDF();
  dedxFitMeanResidualsFuncEta12->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(6));
  dedxFitMeanResidualsFuncEta12->SetChisquare(residualDeDxVsMeanTGraphsEta[11]->Chisquare(dedxFitMeanResidualsFuncEta12));
  dedxFitMeanResidualsFuncEta12->SetNDF(residualDeDxVsMeanTGraphsEta[11]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[11]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta12);
  residualChi2OverNdfs[11] = dedxFitMeanResidualsFuncEta12->GetChisquare()/dedxFitMeanResidualsFuncEta12->GetNDF();
  dedxFitMeanResidualsFuncEta13->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(7));
  dedxFitMeanResidualsFuncEta13->SetChisquare(residualDeDxVsMeanTGraphsEta[12]->Chisquare(dedxFitMeanResidualsFuncEta13));
  dedxFitMeanResidualsFuncEta13->SetNDF(residualDeDxVsMeanTGraphsEta[12]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[12]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta13);
  residualChi2OverNdfs[12] = dedxFitMeanResidualsFuncEta13->GetChisquare()/dedxFitMeanResidualsFuncEta13->GetNDF();
  dedxFitMeanResidualsFuncEta14->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(8));
  dedxFitMeanResidualsFuncEta14->SetChisquare(residualDeDxVsMeanTGraphsEta[13]->Chisquare(dedxFitMeanResidualsFuncEta14));
  dedxFitMeanResidualsFuncEta14->SetNDF(residualDeDxVsMeanTGraphsEta[13]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[13]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta14);
  residualChi2OverNdfs[13] = dedxFitMeanResidualsFuncEta14->GetChisquare()/dedxFitMeanResidualsFuncEta14->GetNDF();
  dedxFitMeanResidualsFuncEta15->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(9));
  dedxFitMeanResidualsFuncEta15->SetChisquare(residualDeDxVsMeanTGraphsEta[14]->Chisquare(dedxFitMeanResidualsFuncEta15));
  dedxFitMeanResidualsFuncEta15->SetNDF(residualDeDxVsMeanTGraphsEta[14]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[14]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta15);
  residualChi2OverNdfs[14] = dedxFitMeanResidualsFuncEta15->GetChisquare()/dedxFitMeanResidualsFuncEta15->GetNDF();
  dedxFitMeanResidualsFuncEta16->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(10));
  dedxFitMeanResidualsFuncEta16->SetChisquare(residualDeDxVsMeanTGraphsEta[15]->Chisquare(dedxFitMeanResidualsFuncEta16));
  dedxFitMeanResidualsFuncEta16->SetNDF(residualDeDxVsMeanTGraphsEta[15]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[15]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta16);
  residualChi2OverNdfs[15] = dedxFitMeanResidualsFuncEta16->GetChisquare()/dedxFitMeanResidualsFuncEta16->GetNDF();


  // Last batch of eta slices
  const int Npar2 = 11;
  double par2[Npar2] = {-1.86,1.23,-0.34,0.038,-0.00012,1,1,1,1,1,1};
  fitter.Config().SetParamsSettings(Npar2,par2);
  //fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().MinimizerOptions().SetMaxFunctionCalls(20000);
  fitter.Config().SetMinimizer("Minuit","Migrad"); 
  fitter.FitFCN(Npar2,globalChi2Eta17to23,par2);
  result = fitter.Result();
  result.Print(std::cout);

  dedxFitMeanResidualsFuncEta17->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4));
  dedxFitMeanResidualsFuncEta17->SetChisquare(residualDeDxVsMeanTGraphsEta[16]->Chisquare(dedxFitMeanResidualsFuncEta17));
  dedxFitMeanResidualsFuncEta17->SetNDF(residualDeDxVsMeanTGraphsEta[16]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[16]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta17);
  residualChi2OverNdfs[16] = dedxFitMeanResidualsFuncEta17->GetChisquare()/dedxFitMeanResidualsFuncEta17->GetNDF();
  dedxFitMeanResidualsFuncEta18->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(5));
  dedxFitMeanResidualsFuncEta18->SetChisquare(residualDeDxVsMeanTGraphsEta[17]->Chisquare(dedxFitMeanResidualsFuncEta18));
  dedxFitMeanResidualsFuncEta18->SetNDF(residualDeDxVsMeanTGraphsEta[17]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[17]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta18);
  residualChi2OverNdfs[17] = dedxFitMeanResidualsFuncEta18->GetChisquare()/dedxFitMeanResidualsFuncEta18->GetNDF();
  dedxFitMeanResidualsFuncEta19->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(6));
  dedxFitMeanResidualsFuncEta19->SetChisquare(residualDeDxVsMeanTGraphsEta[18]->Chisquare(dedxFitMeanResidualsFuncEta19));
  dedxFitMeanResidualsFuncEta19->SetNDF(residualDeDxVsMeanTGraphsEta[18]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[18]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta19);
  residualChi2OverNdfs[18] = dedxFitMeanResidualsFuncEta19->GetChisquare()/dedxFitMeanResidualsFuncEta19->GetNDF();
  dedxFitMeanResidualsFuncEta20->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(7));
  dedxFitMeanResidualsFuncEta20->SetChisquare(residualDeDxVsMeanTGraphsEta[19]->Chisquare(dedxFitMeanResidualsFuncEta20));
  dedxFitMeanResidualsFuncEta20->SetNDF(residualDeDxVsMeanTGraphsEta[19]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[19]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta20);
  residualChi2OverNdfs[19] = dedxFitMeanResidualsFuncEta20->GetChisquare()/dedxFitMeanResidualsFuncEta20->GetNDF();
  dedxFitMeanResidualsFuncEta21->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(8));
  dedxFitMeanResidualsFuncEta21->SetChisquare(residualDeDxVsMeanTGraphsEta[20]->Chisquare(dedxFitMeanResidualsFuncEta21));
  dedxFitMeanResidualsFuncEta21->SetNDF(residualDeDxVsMeanTGraphsEta[20]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[20]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta21);
  residualChi2OverNdfs[20] = dedxFitMeanResidualsFuncEta21->GetChisquare()/dedxFitMeanResidualsFuncEta21->GetNDF();
  dedxFitMeanResidualsFuncEta22->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(9));
  dedxFitMeanResidualsFuncEta22->SetChisquare(residualDeDxVsMeanTGraphsEta[21]->Chisquare(dedxFitMeanResidualsFuncEta22));
  dedxFitMeanResidualsFuncEta22->SetNDF(residualDeDxVsMeanTGraphsEta[21]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[21]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta22);
  residualChi2OverNdfs[21] = dedxFitMeanResidualsFuncEta22->GetChisquare()/dedxFitMeanResidualsFuncEta22->GetNDF();
  dedxFitMeanResidualsFuncEta23->SetParameters(result.Parameter(0),result.Parameter(1),result.Parameter(2),result.Parameter(3),result.Parameter(4),result.Parameter(10));
  dedxFitMeanResidualsFuncEta23->SetChisquare(residualDeDxVsMeanTGraphsEta[22]->Chisquare(dedxFitMeanResidualsFuncEta23));
  dedxFitMeanResidualsFuncEta23->SetNDF(residualDeDxVsMeanTGraphsEta[22]->GetN()-6);
  residualDeDxVsMeanTGraphsEta[22]->GetListOfFunctions()->Add(dedxFitMeanResidualsFuncEta23);
  residualChi2OverNdfs[22] = dedxFitMeanResidualsFuncEta23->GetChisquare()/dedxFitMeanResidualsFuncEta23->GetNDF();
  for(int i=0; i<25; ++i)
    residualDeDxVsMeanTGraphsEta[i]->Write();
  for(int i=0; i<25; ++i)
  {
    cout << "Eta slice: " << i+1 << " Chi2/Ndf: " << residualChi2OverNdfs[i];
    if(i==7 || i==23 || i==24)
      cout << "\t<-- not included in fit";
    cout << endl;
  }

  // Make MultiGraphs of residualDeDx vs. meanDeDx (all eta slices together)
  TMultiGraph* residualDeDxVsMeanMG[6];
  TLegend* residualDeDxVsMeanLeg[6];
  for(int canv=0; canv<6; ++canv)
  {
    residualDeDxVsMeanMG[canv] = new TMultiGraph();
    residualDeDxVsMeanLeg[canv] = new TLegend(0.2,0.7,0.3,0.9);
    residualDeDxVsMeanLeg[canv]->SetHeader("#eta bin");
    int colorCounter = 1;
    for(int i=4*canv+1; i<4*canv+5; ++i)
    {
      if(residualDeDxVsMeanTGraphsEta[i-1]->GetN() > 0)
      {
        if(colorCounter==5)
          ++colorCounter;
        residualDeDxVsMeanTGraphsEta[i-1]->SetLineColor(colorCounter);
        residualDeDxVsMeanTGraphsEta[i-1]->SetMarkerColor(colorCounter);
        residualDeDxVsMeanMG[canv]->Add(residualDeDxVsMeanTGraphsEta[i-1],"p");
        residualDeDxVsMeanLeg[canv]->AddEntry(residualDeDxVsMeanTGraphsEta[i-1],intToString(i).c_str(),"l");
        ++colorCounter;
      }
    }
    if(canv==5)
    {
      residualDeDxVsMeanTGraphsEta[24]->SetLineColor(6);
      residualDeDxVsMeanTGraphsEta[24]->SetMarkerColor(6);
      residualDeDxVsMeanMG[canv]->Add(residualDeDxVsMeanTGraphsEta[24],"p");
      residualDeDxVsMeanLeg[canv]->AddEntry(residualDeDxVsMeanTGraphsEta[24],intToString(25).c_str(),"l");
    }
    TCanvas t;
    if(canv==0)
      t.SetName("residualDeDxVsMeanEtaSlices1_4");
    else if(canv==1)
      t.SetName("residualDeDxVsMeanEtaSlices5_8");
    else if(canv==2)
      t.SetName("residualDeDxVsMeanEtaSlices9_12");
    else if(canv==3)
      t.SetName("residualDeDxVsMeanEtaSlices13_16");
    else if(canv==4)
      t.SetName("residualDeDxVsMeanEtaSlices17_20");
    else if(canv==5)
      t.SetName("residualDeDxVsMeanEtaSlices21_25");
    residualDeDxVsMeanMG[canv]->Draw("a");
    residualDeDxVsMeanMG[canv]->GetYaxis()->SetTitle("#Delta=fit-measured");
    residualDeDxVsMeanMG[canv]->GetYaxis()->SetRangeUser(-1,10);
    residualDeDxVsMeanMG[canv]->GetXaxis()->SetTitle("fit Ih [MeV/cm]");
    residualDeDxVsMeanMG[canv]->GetXaxis()->SetRangeUser(2,18);
    residualDeDxVsMeanLeg[canv]->Draw("same");
    t.Write();
  }

}

