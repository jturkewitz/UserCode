#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include <iostream>
#include <memory>
#include <string>
#include <sstream>

#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"

#include "RooRealVar.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

#include "dedxSingleHitDists.h"


// Created FEB 23 2011 by SIC
// Sample a Landau with a number of points and take the harmonic^2 mean, like in the HSCP note
// Modified SEP 8 2011
// Use the single-hit dE/dx distribution from 2011 data
using namespace RooFit;

//
std::string intToString(int num)
{
    using namespace std;
    ostringstream myStream;
    myStream << num << flush;
    return (myStream.str ());
}


int main(int argc, char* argv[])
{
    using namespace std;

    bool exampleFit = false;
    int numSamplePoints320 = 7; // number of measurements of the 320 um single-hit dist for harm^2 avg.
    int numSamplePoints500 = 7; // number of measurements of the 500 um single-hit dist for harm^2 avg.
    int numSamplePointsAll = 7; // number of measurements of the all hits single-hit dist for harm^2 avg.
    int numSamplePointsMixed320 = 3; //2; // nom for mixing 500 and 320 um hits
    int numSamplePointsMixed500 = 4; //3; // nom for mixing 500 and 320 um hits

    int numberOfTrials = 50000;  // number of tracks to simulate

    RooRealVar hitDeDx("hitDeDx","hitDeDx",0,50,"MeV/cm");
    //RooRealVar landauDeDxMuonMean("landauDeDxMuonMean","landauDeDxMuonMean",2.816);
    //RooRealVar landauDeDxMuonSigma("landauDeDxMuonSigma","landauDeDxMuonSigma",0.2896);
    //RooRealVar landauDeDxMuonMean("landauDeDxMuonMean","landauDeDxMuonMean",3);
    //RooRealVar landauDeDxMuonSigma("landauDeDxMuonSigma","landauDeDxMuonSigma",0.3);
    //RooLandau landauDeDxMuon("landauDeDxMuon","landauDeDxMuon",hitDeDx,landauDeDxMuonMean,landauDeDxMuonSigma);

    // init single hit dE/dx dists from data -- in include file above
    initHists();
    // dataHist/histPdf 320 um
    RooDataHist dedx320DataHist("dedx320DataHist","dedx320DataHist",hitDeDx,dedxSingleHit320etaBin10Hist);
    RooHistPdf dedx320HistPDF("dedx320HistPDF","dedx320HistPDF",hitDeDx,dedx320DataHist);
    // dataHist/histPdf 500 um
    RooDataHist dedx500DataHist("dedx500DataHist","dedx500DataHist",hitDeDx,dedxSingleHit500etaBin10Hist);
    RooHistPdf dedx500HistPDF("dedx500HistPDF","dedx500HistPDF",hitDeDx,dedx500DataHist);
    // dataHist/histPdf All um
    RooDataHist dedxAllDataHist("dedxAllDataHist","dedxAllDataHist",hitDeDx,dedxSingleHitAlletaBin10Hist);
    RooHistPdf dedxAllHistPDF("dedxAllHistPDF","dedxAllHistPDF",hitDeDx,dedxAllDataHist);

    TFile* outfile = new TFile("sampleLandauToyMC.root","recreate");
    outfile->cd();
    TH1F* harm2Avg320Hist = new TH1F("harm2Avg320Hist","Dist. of harmonic^2 mean (Ih) from 320 um hits; MeV/cm",200,0,20);
    TH1F* harm2Avg500Hist = new TH1F("harm2Avg500Hist","Dist. of harmonic^2 mean (Ih) from 500 um hits; MeV/cm",200,0,20);
    TH1F* harm2AvgAllHist = new TH1F("harm2AvgAllHist","Dist. of harmonic^2 mean (Ih) from all hits; MeV/cm",200,0,20);
    string mixedTitle="Dist. of harmonic^2 mean (Ih) from mixture of ";
    mixedTitle+=intToString(numSamplePointsMixed320);
    mixedTitle+=" hits from 320 um layers and ";
    mixedTitle+=intToString(numSamplePointsMixed500);
    mixedTitle+=" hits from 500 um layers; MeV/cm";
    TH1F* harm2AvgMixedHist = new TH1F("harm2AvgMixedHist",mixedTitle.c_str(),200,0,20);
    TCanvas t;
    t.SetName("exampleFit");
    t.SetTitle("exampleFit");

    for(int trial=0; trial < numberOfTrials; ++trial)
    {
      RooDataSet data320(*dedx320HistPDF.generate(hitDeDx,numSamplePoints320));
      RooDataSet data500(*dedx500HistPDF.generate(hitDeDx,numSamplePoints500));
      RooDataSet dataAll(*dedxAllHistPDF.generate(hitDeDx,numSamplePointsAll));

      RooDataSet dataMixed(*dedx320HistPDF.generate(hitDeDx,numSamplePointsMixed320));
      dataMixed.append(*dedx500HistPDF.generate(hitDeDx,numSamplePointsMixed500));

      //if(exampleFit)
      //{
      //  t.cd();
      //  landauDeDxMuon.fitTo(*data);
      //  RooPlot* dedxFrame = hitDeDx.frame();
      //  data->plotOn(dedxFrame);
      //  landauDeDxMuon.plotOn(dedxFrame);
      //  dedxFrame->Draw("e0");
      //  exampleFit = false;
      //}

      float harmMean320 = 0;
      // Loop over the points and calculate the harmonic^2 mean, 320 um
      for(int entry = 0; entry < data320.numEntries(); ++entry)
      {
        const RooArgSet temp = *(data320.get(entry));
        hitDeDx = temp.getRealValue("hitDeDx");
        harmMean320+=pow(hitDeDx.getVal(),-2);
        //std::cout << "hitDeDx: " << hitDeDx.getVal() << std::endl;
      }
      harmMean320/=data320.numEntries();
      harmMean320=pow(harmMean320,-0.5);
      harm2Avg320Hist->Fill(harmMean320);

      float harmMean500 = 0;
      // Loop over the points and calculate the harmonic^2 mean, 500 um
      for(int entry = 0; entry < data500.numEntries(); ++entry)
      {
        const RooArgSet temp = *(data500.get(entry));
        hitDeDx = temp.getRealValue("hitDeDx");
        harmMean500+=pow(hitDeDx.getVal(),-2);
        //std::cout << "hitDeDx: " << hitDeDx.getVal() << std::endl;
      }
      harmMean500/=data500.numEntries();
      harmMean500=pow(harmMean500,-0.5);
      harm2Avg500Hist->Fill(harmMean500);

      float harmMeanAll = 0;
      // Loop over the points and calculate the harmonic^2 mean, all hits
      for(int entry = 0; entry < dataAll.numEntries(); ++entry)
      {
        const RooArgSet temp = *(dataAll.get(entry));
        hitDeDx = temp.getRealValue("hitDeDx");
        harmMeanAll+=pow(hitDeDx.getVal(),-2);
        //std::cout << "hitDeDx: " << hitDeDx.getVal() << std::endl;
      }
      harmMeanAll/=dataAll.numEntries();
      harmMeanAll=pow(harmMeanAll,-0.5);
      harm2AvgAllHist->Fill(harmMeanAll);

      float harmMeanMixed = 0;
      // Loop over the points and calculate the harmonic^2 mean, all hits
      for(int entry = 0; entry < dataMixed.numEntries(); ++entry)
      {
        const RooArgSet temp = *(dataMixed.get(entry));
        hitDeDx = temp.getRealValue("hitDeDx");
        harmMeanMixed+=pow(hitDeDx.getVal(),-2);
        //std::cout << "hitDeDx: " << hitDeDx.getVal() << std::endl;
      }
      harmMeanMixed/=dataMixed.numEntries();
      harmMeanMixed=pow(harmMeanMixed,-0.5);
      harm2AvgMixedHist->Fill(harmMeanMixed);

    }

    harm2Avg320Hist->Write();
    harm2Avg500Hist->Write();
    harm2AvgAllHist->Write();
    harm2AvgMixedHist->Write();
    t.Write();
    outfile->Close();

}
