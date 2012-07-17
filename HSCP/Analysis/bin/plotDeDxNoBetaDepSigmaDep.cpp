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
#include <iomanip>


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

  gROOT->Macro("../rootstyle.C");

  TFile* inputTFile = new TFile(inputFile);
  inputTFile->cd();

  TFile* outputTFile = new TFile("plotDeDxNoBetaDepSigmaDep.root","recreate");

  int nPBins = 200;
  const int nNomBins = 30;
  const int nEtaBins = 25;

  double dedxSigmaBetaDepEtaParam0[nEtaBins];
  TH1F* dedxNoBetaDepHists[nEtaBins][nNomBins];
  TGraphErrors* meanVsNoMGraphsEtaSlices[nEtaBins];
  TGraphErrors* sigmaVsNoMGraphsEtaSlices[nEtaBins];
  TH2F* dedxNoBetaDepHistsInEtaBins[nEtaBins];
  TGraphErrors* meanDeDxVsBetaTGraphsEta[nEtaBins];
  TGraphErrors* sigmaDeDxVsBetaTGraphsEta[nEtaBins];
  TGraphErrors* residualDeDxVsSigmaTGraphsEta[nEtaBins];
  for(int i=0; i<nEtaBins; ++i)
  {
    sigmaDeDxVsBetaTGraphsEta[i] = new TGraphErrors();
    meanDeDxVsBetaTGraphsEta[i] = new TGraphErrors();
    residualDeDxVsSigmaTGraphsEta[i] = new TGraphErrors();
  }
  // parameters from 1-parameter fit to sigma(Ih) vs. beta with restricted beta range
  // see: pcminn19:/data3/scooper/cmssw/426/HSCPStudiesAndToyMC/src/HSCP2011/ToyMC/bin/1.aug25.dedxSigmaPlots.takeOutNoMDep.aug26.1parBetaFit
  // also posted on the web at: hscp/aug26/1.aug25.dedxSigmaPlots.takeOutNoMDep.aug26.1parBetaFit
  //   example: sigmaIhMeasMinusIhFitPlusResNoMOutVsPfit_etaBin1.png
  dedxSigmaBetaDepEtaParam0[0] = 0.1519;
  dedxSigmaBetaDepEtaParam0[1] = 0.1504;
  dedxSigmaBetaDepEtaParam0[2] = 0.1502;
  dedxSigmaBetaDepEtaParam0[3] = 0.147;
  dedxSigmaBetaDepEtaParam0[4] = 0.1451;
  dedxSigmaBetaDepEtaParam0[5] = 0.1497;
  dedxSigmaBetaDepEtaParam0[6] = 0.1208;
  dedxSigmaBetaDepEtaParam0[7] = 0;        // ignored for now
  dedxSigmaBetaDepEtaParam0[8] = 0.1489;
  dedxSigmaBetaDepEtaParam0[9] = 0.1347;
  dedxSigmaBetaDepEtaParam0[10] = 0.1392;
  dedxSigmaBetaDepEtaParam0[11] = 0.1359;
  dedxSigmaBetaDepEtaParam0[12] = 0.1428;
  dedxSigmaBetaDepEtaParam0[13] = 0.1555;
  dedxSigmaBetaDepEtaParam0[14] = 0.1433;
  dedxSigmaBetaDepEtaParam0[15] = 0.1456;
  dedxSigmaBetaDepEtaParam0[16] = 0.1519;
  dedxSigmaBetaDepEtaParam0[17] = 0.1585;
  dedxSigmaBetaDepEtaParam0[18] = 0.1532;
  dedxSigmaBetaDepEtaParam0[19] = 0.1556;
  dedxSigmaBetaDepEtaParam0[20] = 0.1575;
  dedxSigmaBetaDepEtaParam0[21] = 0.1646;
  dedxSigmaBetaDepEtaParam0[22] = 0.1659;
  dedxSigmaBetaDepEtaParam0[23] = 0;       // ignored
  dedxSigmaBetaDepEtaParam0[24] = 0;       // ignored

  // fit func for NoM dependence
  // normalize to nom=10 for eta 1-14; nom=12 for eta > 15
  TF1* myNomDepFunc = new TF1("myNomDepFunc","[0]/TMath::Sqrt(x/[1])",5,25);
  // fit func for gaussians
  TF1* myGaus = new TF1("myGaus","gaus(0)",-5,5);
  // fit func for beta-dependence of sigma
  //TF1* myFitFunc = new TF1("myFitFunc","[0]+[1]/x^2",0.5,1);
  TF1* myFitFunc = new TF1("myFitFunc","[0]/x^2",0.5,1);
  // fit func for residuals
  TF1* residualSigmaDeDxVsSigmaDeDxFitFunc = new TF1("residualSigmaDeDxVsSigmaDeDxFitFunc","pol1(0)",0.1,1);

  TH3F* rawHist;

  string dir = "dedxHistsInNomBinsMatched";
  // loop over nom bins
  for(int nomBin=1; nomBin < nNomBins+1; ++nomBin)
  {
    inputTFile->cd();
    string name = "dedxE1NoPBiasMatchedNomBin";
    name+=intToString(nomBin);
    string fullpath = dir;
    fullpath+="/";
    fullpath+=name;
    rawHist = (TH3F*) inputTFile->Get(fullpath.c_str());

    // loop over eta bins
    for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
    {
      rawHist->GetZaxis()->SetRange(etaBin,etaBin);
      // set beta/p range (limit saturation effects)
      int minPbin = 135; // higher p (beta) only
      if(etaBin < 6)
        minPbin = 100;
      rawHist->GetYaxis()->SetRange(minPbin,nPBins);

      string projName = "x";
      projName+="_EtaBin";
      projName+=intToString(etaBin);
      dedxNoBetaDepHists[etaBin-1][nomBin-1] = (TH1F*)rawHist->Project3D(projName.c_str());
      string histName = "dedxE1NoPBiasMatchedEtaBin";
      histName+=intToString(etaBin);
      histName+="NomBin";
      histName+=intToString(nomBin);
      dedxNoBetaDepHists[etaBin-1][nomBin-1]->SetName(histName.c_str());
    }
  }

  double sigmaChi2OverNdfs[25];

  cout << setiosflags(ios::left);
  cout << setw(8) << "etaBin" << setw(12) << "chi2" << setw(6) << "Ndf" << setw(6) << "chi2/Ndf" << endl;
  // loop over the newly-stored hists -- eta bins first
  for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
  {
    vector<double> means;
    vector<double> sigmas;
    vector<double> meanEs;
    vector<double> sigmaEs;
    vector<double> noms;
    vector<double> nomEs;

    // make output dir
    outputTFile->cd();
    string name2 = "dedxE1NoPBiasMatchedInEtaBin";
    name2+=intToString(etaBin);
    TDirectory* dedxE1MatchedEtaBinDir = outputTFile->mkdir(name2.c_str());
    dedxE1MatchedEtaBinDir->cd();

    // loop over nom bins, fit make graphs
    for(int nomBin=1; nomBin < nNomBins+1; ++nomBin)
    {
      if(dedxNoBetaDepHists[etaBin-1][nomBin-1]->GetEntries() > 25)
      {
        //myGaus->SetParameter(0,500);
        myGaus->SetParameter(1,0);
        myGaus->SetParameter(2,0.5);
        //XXX limit sigma to be positive?
        myGaus->SetParLimits(2,0,1);
        TFitResultPtr r = dedxNoBetaDepHists[etaBin-1][nomBin-1]->Fit("myGaus","SQR");
        int ret = r;
        if(ret==0)
        {
          means.push_back(r->Parameter(1));
          sigmas.push_back(r->Parameter(2));
          meanEs.push_back(r->ParError(1));
          sigmaEs.push_back(r->ParError(2));
          noms.push_back(nomBin);
          nomEs.push_back(0);
        }
        dedxNoBetaDepHists[etaBin-1][nomBin-1]->GetXaxis()->SetTitle("Ih_meas-Ih_fit+res [MeV/cm]");
        dedxNoBetaDepHists[etaBin-1][nomBin-1]->Write();
      }
    }
    // Make graphs, fit, write
    string meanGraphName="meanVsNoMEtaSlice";
    meanGraphName+=intToString(etaBin);
    meanVsNoMGraphsEtaSlices[etaBin-1] = new TGraphErrors(noms.size(),&(*noms.begin()),&(*means.begin()),&(*nomEs.begin()),&(*meanEs.begin()));
    meanVsNoMGraphsEtaSlices[etaBin-1]->SetName(meanGraphName.c_str());
    meanVsNoMGraphsEtaSlices[etaBin-1]->SetTitle(meanGraphName.c_str());
    meanVsNoMGraphsEtaSlices[etaBin-1]->GetXaxis()->SetTitle("NoM");
    meanVsNoMGraphsEtaSlices[etaBin-1]->GetYaxis()->SetTitle("Ih_meas-Ih_fit+res [MeV/cm]");
    meanVsNoMGraphsEtaSlices[etaBin-1]->Write();
    TCanvas tMean;
    tMean.cd();
    meanVsNoMGraphsEtaSlices[etaBin-1]->GetYaxis()->SetRangeUser(-0.75,0.75);
    meanVsNoMGraphsEtaSlices[etaBin-1]->Draw("ap");
    string canName = "meanIhMeasMinusIhFitPlusResVsNoMfit_etaBin";
    canName+=intToString(etaBin);
    tMean.Print((canName+".png").c_str());

    string sigmaGraphName="sigmaVsNoMEtaSlice";
    sigmaGraphName+=intToString(etaBin);
    sigmaVsNoMGraphsEtaSlices[etaBin-1] = new TGraphErrors(noms.size(),&(*noms.begin()),&(*sigmas.begin()),&(*nomEs.begin()),&(*sigmaEs.begin()));
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->SetName(sigmaGraphName.c_str());
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->SetTitle(sigmaGraphName.c_str());
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->GetXaxis()->SetTitle("NoM");
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->GetYaxis()->SetTitle("#sigma(Ih_meas-Ih_fit+res) [MeV/cm]");

    unsigned int nomNormIndex = 0;
    int nomToNormalizeTo = 10;
    if(etaBin>=15)
      nomToNormalizeTo = 12;
    while(nomNormIndex<noms.size() && noms[nomNormIndex] != nomToNormalizeTo)
      ++nomNormIndex;
    if(nomNormIndex>=noms.size())
    {
      cout << "ERROR: couldn't find nom==" << nomToNormalizeTo << " in this eta slice: " << etaBin << endl;
      continue;
    }
    myNomDepFunc->SetParameter(0,sigmas[nomNormIndex]);
    myNomDepFunc->SetParameter(1,nomToNormalizeTo);

    myNomDepFunc->SetChisquare(sigmaVsNoMGraphsEtaSlices[etaBin-1]->Chisquare(myNomDepFunc));
    myNomDepFunc->SetNDF(sigmaVsNoMGraphsEtaSlices[etaBin-1]->GetN()-1);
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->GetListOfFunctions()->Add(myNomDepFunc);
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->Write();

    TCanvas tSig;
    tSig.cd();
    sigmaVsNoMGraphsEtaSlices[etaBin-1]->Draw("ap");
    myNomDepFunc->Draw("same");
    canName = "sigmaIhMeasMinusIhFitPlusResVsNoMfit_etaBin";
    canName+=intToString(etaBin);
    tSig.Print((canName+".png").c_str());
    sigmaChi2OverNdfs[etaBin-1] = myNomDepFunc->GetChisquare()/myNomDepFunc->GetNDF();
    cout << setw(8) << etaBin << setw(12) << sigmaVsNoMGraphsEtaSlices[etaBin-1]->Chisquare(myNomDepFunc) << setw(6) << sigmaVsNoMGraphsEtaSlices[etaBin-1]->GetN()-1;
    cout << setw(6) << sigmaChi2OverNdfs[etaBin-1] << endl;

  }

  cout << endl;
  cout << setiosflags(ios::left);
  cout << setw(8) << "etaBin" << setw(12) << "chi2" << setw(6) << "Ndf" << setw(6) << "chi2/Ndf" << endl;
  TH1D* myHist;
  vector<double> etas;
  vector<double> etaEs;
  vector<double> p0s;
  vector<double> p0Es;
  vector<double> p1s;
  vector<double> p1Es;
  // now to get the beta dependence; could probably be combined with above loop
  // loop over eta bins
  for(int etaBin=1; etaBin<nEtaBins+1; ++etaBin)
  {
    inputTFile->cd();
    //string name = "dedxE1NoPBiasMatchedNomBin";
    string name = "dedxE1NoPBiasNoNomBiasMatchedNomBin";
    name+=intToString(1);

    string fullpath = dir;
    fullpath+="/";
    fullpath+=name;
    rawHist = (TH3F*) inputTFile->Get(fullpath.c_str());
    rawHist->GetZaxis()->SetRange(etaBin,etaBin);
    rawHist->GetYaxis()->SetRange(1,nPBins);
    string projName = "yx";
    projName+="_etaBin";
    projName+=intToString(etaBin);
    TH2F* myHist2d = (TH2F*)rawHist->Project3D(projName.c_str());
    outputTFile->cd();
    //myHist2d->Write();
    for(int nomBin=2; nomBin < nNomBins+1; ++nomBin)
    {
      inputTFile->cd();
      string nameNom = "dedxE1NoPBiasNoNomBiasMatchedNomBin";
      //string nameNom = "dedxE1NoPBiasMatchedNomBin";
      nameNom+=intToString(nomBin);
      string fullpathNom = dir;
      fullpathNom+="/";
      fullpathNom+=nameNom;
      rawHist = (TH3F*) inputTFile->Get(fullpathNom.c_str());
      rawHist->GetZaxis()->SetRange(etaBin,etaBin);
      rawHist->GetYaxis()->SetRange(1,nPBins);
      string projNameNom = "yx";
      projNameNom+="_etaBin";
      projNameNom+=intToString(etaBin);
      projNameNom+="_nomBin";
      projNameNom+=intToString(nomBin);
      TH2F* thisHist = (TH2F*)rawHist->Project3D(projNameNom.c_str());
      // thisHist is p (y) vs. dE/dx E1 (x)
      outputTFile->cd();
      //thisHist->Write();
      myHist2d->Add(thisHist);
      string histName = "dedxInNomBin";
      histName+=intToString(nomBin);
      histName+="EtaBin";
      histName+=intToString(etaBin);
      thisHist->SetName(histName.c_str());
    }
    dedxNoBetaDepHistsInEtaBins[etaBin-1] = myHist2d;
    // make output dir
    outputTFile->cd();
    string name2 = "dedxNoBetaDepHistsIn";
    name2+="etaBin";
    name2+=intToString(etaBin);
    TDirectory* dedxE1MatchedEtaBinDir = outputTFile->mkdir(name2.c_str());
    dedxE1MatchedEtaBinDir->cd();
    dedxNoBetaDepHistsInEtaBins[etaBin-1]->Write();

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
      myHist = dedxNoBetaDepHistsInEtaBins[etaBin-1]->ProjectionX(thisName2.c_str(),pBin,pBin);

      if(myHist->GetEntries()>10)
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
    // write hist
    string graphName = "meanDeDxVsBetaEtaBin";
    graphName+=intToString(etaBin);
    t->SetName(graphName.c_str());
    t->SetTitle(graphName.c_str());
    t->GetXaxis()->SetTitle("#beta");
    t->GetYaxis()->SetTitle("mean_Ih [MeV/cm]");
    t->Write();
    meanDeDxVsBetaTGraphsEta[etaBin-1] = t;

    t = new TGraphErrors(betas.size(),&(*betas.begin()),&(*sigmas.begin()),&(*betaEs.begin()),&(*sigmaEs.begin()));
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
    graphName = "sigmaDeDxVsBetaEtaBin";
    graphName+=intToString(etaBin);
    t->SetName(graphName.c_str());
    t->SetTitle(graphName.c_str());
    t->GetXaxis()->SetTitle("#beta");
    t->GetYaxis()->SetTitle("sigma_Ih [MeV/cm]");
    t->Write();
    sigmaDeDxVsBetaTGraphsEta[etaBin-1] = t;
    TCanvas tSig;
    tSig.cd();
    sigmaDeDxVsBetaTGraphsEta[etaBin-1]->GetYaxis()->SetRangeUser(0,1);
    sigmaDeDxVsBetaTGraphsEta[etaBin-1]->Draw("ap");
    myFitFunc->Draw("same");
    string canName = "sigmaIhMeasMinusIhFitPlusResNomOutVsPfit_etaBin";
    canName+=intToString(etaBin);
    tSig.Print((canName+".png").c_str());

    cout << setw(8) << etaBin << setw(12) << myFitFunc->GetChisquare() << setw(6) << myFitFunc->GetNDF();
    cout << setw(6) << myFitFunc->GetChisquare() / myFitFunc->GetNDF() << endl;

    // make residuals vs. expectedSigmaDeDx in this eta bin
    if(sigmas.size() <= 0)
      continue;
    vector<double> residuals;
    vector<double> expDeDx;
    vector<double> expDeDxEs;
    vector<double>::const_iterator sigmaItr = sigmas.begin();
    for(vector<double>::const_iterator betaItr = betas.begin();
        betaItr != betas.end(); ++betaItr)
    {
      myFitFunc->SetParameter(0,dedxSigmaBetaDepEtaParam0[etaBin-1]);
      double fitDeDx = myFitFunc->Eval(*betaItr);
      // SIC MOD over 5 only
      //if(fitDeDx >= 5)
      //{
        residuals.push_back(fitDeDx-*sigmaItr);
        expDeDx.push_back(fitDeDx);
        expDeDxEs.push_back(0);
      //}
      ////debug
      //if(etaBin==1)
      //cout << "For eta bin 1,  beta=" << *betaItr << " meanDeDx=" << *meanItr << " expectedDeDx=" << myFitFunc->Eval(*betaItr)
      //  << " difference = fit-mean = " << myFitFunc->Eval(*betaItr)-*meanItr
      //  << " param0=" << myFitFunc->GetParameter(0) << " param1=" << myFitFunc->GetParameter(1)
      //  << endl; 
      ++sigmaItr;
    }
    TCanvas tcan;
    tcan.cd();
    TGraphErrors* t3 = new TGraphErrors(expDeDx.size(),&(*expDeDx.begin()),&(*residuals.begin()),&(*expDeDxEs.begin()),&(*meanEs.begin()));
    t3->Fit("residualSigmaDeDxVsSigmaDeDxFitFunc","Q");
    cout << "residualFit:" << setw(8) << etaBin << setw(12) << residualSigmaDeDxVsSigmaDeDxFitFunc->GetChisquare() << setw(6) << residualSigmaDeDxVsSigmaDeDxFitFunc->GetNDF();
    cout << setw(6) << residualSigmaDeDxVsSigmaDeDxFitFunc->GetChisquare() / residualSigmaDeDxVsSigmaDeDxFitFunc->GetNDF() << endl;
    // write hist
    string graphName3 = "residualDeDxVsSigmaEtaBin";
    graphName3+=intToString(etaBin);
    t3->SetName(graphName3.c_str());
    t3->SetTitle(graphName3.c_str());
    t3->GetXaxis()->SetTitle("fit_#sigma(Ih) [MeV/cm]");
    t3->GetYaxis()->SetTitle("residual_#sigma(Ih) [MeV/cm]");
    t3->GetYaxis()->SetRangeUser(-1,1);
    t3->Draw("ap");
    tcan.Print((graphName3+".png").c_str());
    t3->Write();
    residualDeDxVsSigmaTGraphsEta[etaBin-1] = t3;

  }

  outputTFile->cd();

  // Make MultiGraphs of sigmaDeDx vs. beta (all eta slices together)
  TMultiGraph* sigmaDeDxVsBetaMG[6];
  TLegend* sigmaDeDxVsBetaLeg[6];
  for(int canv=0; canv<6; ++canv)
  {
    sigmaDeDxVsBetaMG[canv] = new TMultiGraph();
    sigmaDeDxVsBetaLeg[canv] = new TLegend(0.2,0.7,0.3,0.9);
    sigmaDeDxVsBetaLeg[canv]->SetHeader("#eta bin");
    int colorCounter = 1;
    for(int i=4*canv+1; i<4*canv+5; ++i)
    {
      if(sigmaDeDxVsBetaTGraphsEta[i-1]->GetN() > 0)
      {
        if(colorCounter==5)
          ++colorCounter;
        sigmaDeDxVsBetaTGraphsEta[i-1]->SetLineColor(colorCounter);
        sigmaDeDxVsBetaTGraphsEta[i-1]->SetMarkerColor(colorCounter);
        sigmaDeDxVsBetaMG[canv]->Add(sigmaDeDxVsBetaTGraphsEta[i-1],"p");
        sigmaDeDxVsBetaLeg[canv]->AddEntry(sigmaDeDxVsBetaTGraphsEta[i-1],intToString(i).c_str(),"l");
        ++colorCounter;
      }
    }
    if(canv==5)
    {
      sigmaDeDxVsBetaTGraphsEta[24]->SetLineColor(6);
      sigmaDeDxVsBetaTGraphsEta[24]->SetMarkerColor(6);
      sigmaDeDxVsBetaMG[canv]->Add(sigmaDeDxVsBetaTGraphsEta[24],"p");
      sigmaDeDxVsBetaLeg[canv]->AddEntry(sigmaDeDxVsBetaTGraphsEta[24],intToString(25).c_str(),"l");
    }
    TCanvas t;
    if(canv==0)
      t.SetName("sigmaDeDxVsBetaEtaSlices1_4");
    else if(canv==1)
      t.SetName("sigmaDeDxVsBetaEtaSlices5_8");
    else if(canv==2)
      t.SetName("sigmaDeDxVsBetaEtaSlices9_12");
    else if(canv==3)
      t.SetName("sigmaDeDxVsBetaEtaSlices13_16");
    else if(canv==4)
      t.SetName("sigmaDeDxVsBetaEtaSlices17_20");
    else if(canv==5)
      t.SetName("sigmaDeDxVsBetaEtaSlices21_25");
    sigmaDeDxVsBetaMG[canv]->Draw("a");
    sigmaDeDxVsBetaMG[canv]->GetXaxis()->SetTitle("#beta");
    sigmaDeDxVsBetaMG[canv]->GetYaxis()->SetTitle("sigma_Ih [MeV/cm]");
    sigmaDeDxVsBetaLeg[canv]->Draw("same");
    t.Write();
  }

  TGraphErrors* p0VsEtasigmaTGraph = new TGraphErrors(etas.size(),&(*etas.begin()),&(*p0s.begin()),&(*etaEs.begin()),&(*p0Es.begin()));
  p0VsEtasigmaTGraph->SetName("p0VsEtasigma");
  p0VsEtasigmaTGraph->SetTitle("p0VsEta (sigma)");
  p0VsEtasigmaTGraph->GetXaxis()->SetTitle("#eta");
  p0VsEtasigmaTGraph->GetYaxis()->SetTitle("p0");
  p0VsEtasigmaTGraph->Write();
  TGraphErrors* p1VsEtasigmaTGraph = new TGraphErrors(etas.size(),&(*etas.begin()),&(*p1s.begin()),&(*etaEs.begin()),&(*p1Es.begin()));
  p1VsEtasigmaTGraph->SetName("p1VsEtasigma");
  p1VsEtasigmaTGraph->SetTitle("p1VsEta (sigma)");
  p1VsEtasigmaTGraph->GetXaxis()->SetTitle("#eta");
  p1VsEtasigmaTGraph->GetYaxis()->SetTitle("p1");
  p1VsEtasigmaTGraph->Write();

  outputTFile->Close();

}

