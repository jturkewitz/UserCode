#ifndef __CINT__
#include "RooGlobalFunc.h"
#else
// Refer to a class implemented in libRooFit to force its loading
// via the autoloader.
class Roo2DKeysPdf;
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/DistFunc.h"
#ifndef __CINT__
#include "RooCFunction1Binding.h" 
#include "RooCFunction3Binding.h"
#endif
#include "RooTFnBinding.h"

using namespace RooFit;
void plotPDistsInNoMSlices()
{
  TCanvas* t = new TCanvas();
  t->cd();
  t->SetLogy();

  // configure which eta slice to look at here
  int slice = 8;
  string fileNameBeg;
  TH2F* pVsNoMInEtaSlice;

  if(slice==1)
  {
    fileNameBeg = "eta0to0p2_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice1");
  }
  else if(slice==2)
  {
    fileNameBeg = "eta0p2to0p4_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice2");
  }
  else if(slice==3)
  {
    fileNameBeg = "eta0p4to0p6_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice3");
  }
  else if(slice==4)
  {
    fileNameBeg = "eta0p6to0p8_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice4");
  }
  else if(slice==5)
  {
    fileNameBeg = "eta0p8to1p0_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice5");
  }
  else if(slice==6)
  {
    fileNameBeg = "eta1p0to1p2_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice6");
  }
  else if(slice==7)
  {
    fileNameBeg = "eta1p2to1p4_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice7");
  }
  else if(slice==8)
  {
    fileNameBeg = "eta1p4to1p6_";
    pVsNoMInEtaSlice = (TH2F*) _file0->Get("PlotsAfterPreselection/pVsNoMAfterPreselectionEtaSlice8");
  }

  // 7-8
  TH1D* nom7to8 = pVsNoMInEtaSlice->ProjectionY("nom7to8",8,9);
  // 9-10
  TH1D* nom9to10 = pVsNoMInEtaSlice->ProjectionY("nom9to10",10,11);
  // 11-12
  TH1D* nom11to12 = pVsNoMInEtaSlice->ProjectionY("nom11to12",12,13);
  // 13-14
  TH1D* nom13to14 = pVsNoMInEtaSlice->ProjectionY("nom13to14",14,15);
  // 15-16
  TH1D* nom15to16 = pVsNoMInEtaSlice->ProjectionY("nom15to16",16,17);
  // 17-18
  TH1D* nom17to18 = pVsNoMInEtaSlice->ProjectionY("nom17to18",18,19);
  // all
  TH1D* nomAll = pVsNoMInEtaSlice->ProjectionY("nomAll",1,pVsNoMInEtaSlice->GetNbinsX());

  
  //TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  //leg->AddEntry("nom7to8","NoM 7-8","l");
  //leg->AddEntry("nom11to12","NoM 11-12","l");
  //leg->AddEntry("nom15to16","NoM 15-16","l");
  //nom7to8->DrawNormalized();
  //nom11to12->SetLineColor(2);
  //nom11to12->DrawNormalized("sames");
  //nom15to16->SetLineColor(4);
  //nom15to16->DrawNormalized("sames");
  //leg->Draw("same");

  double fitRangeMin = 300;
  double fitRangeMax = 1000;
  TF1* myExpo = new TF1("myExpo","expo(0)",fitRangeMin,fitRangeMax);
  nomAll->Fit("myExpo","RL");
  t->Print((fileNameBeg+"nomAll_plot.C").c_str());
  
  TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
  nom7to8->SetTitle("");
  nom7to8->SetMarkerStyle(20);
  nom7to8->GetYaxis()->SetRangeUser(1e-1,1);
  nom7to8->GetXaxis()->SetTitle("p [GeV]");
  nom7to8->DrawNormalized("p");
  nom9to10->SetLineColor(2);
  nom9to10->SetMarkerColor(2);
  nom9to10->SetMarkerSize(1.0);
  nom9to10->SetMarkerStyle(23);
  nom9to10->DrawNormalized("samep");
  nom11to12->SetLineColor(4);
  nom11to12->SetMarkerColor(4);
  nom11to12->SetMarkerSize(1.0);
  nom11to12->SetMarkerStyle(22);
  nom11to12->DrawNormalized("samep");
  nom15to16->SetLineColor(8);
  nom15to16->SetMarkerColor(8);
  //nom15to16->SetMarkerSize(1.0);
  nom15to16->SetMarkerStyle(21);
  nom15to16->DrawNormalized("samep");
  gStyle->SetOptStat(0);
  //gStyle->SetErrorX(0);
  leg->AddEntry("nom7to8","NoM 7-8","p");
  leg->AddEntry("nom9to10","NoM 9-10","p");
  leg->AddEntry("nom11to12","NoM 11-12","p");
  leg->AddEntry("nom15to16","NoM 15-16","p");
  leg->SetBorderSize(0);
  leg->Draw("same");
  t->SetLogy();
  t->Print((fileNameBeg+"selectedNoM_plot.png").c_str());
  return; // don't move on here

  RooRealVar rooVarP("rooVarP","p",fitRangeMin,fitRangeMax);
  RooDataHist* nom7to8DataHist = new RooDataHist("nom7to8dataHist","nom7to8dataHist",
      rooVarP, nom7to8);
  RooHistPdf* nom7to8HistPdf = new RooHistPdf("nom7to8histPdf","nom7to8histPdf",
      rooVarP, *nom7to8DataHist);
  RooDataHist* nom9to10DataHist = new RooDataHist("nom9to10dataHist","nom9to10dataHist",
      rooVarP, nom9to10);
  RooHistPdf* nom9to10HistPdf = new RooHistPdf("nom9to10histPdf","nom9to10histPdf",
      rooVarP, *nom9to10DataHist);
  RooDataHist* nom11to12DataHist = new RooDataHist("nom11to12dataHist","nom11to12dataHist",
      rooVarP, nom11to12);
  RooHistPdf* nom11to12HistPdf = new RooHistPdf("nom11to12histPdf","nom11to12histPdf",
      rooVarP, *nom11to12DataHist);
  RooDataHist* nom13to14DataHist = new RooDataHist("nom13to14dataHist","nom13to14dataHist",
      rooVarP, nom13to14);
  RooHistPdf* nom13to14HistPdf = new RooHistPdf("nom13to14histPdf","nom13to14histPdf",
      rooVarP, *nom13to14DataHist);
  RooDataHist* nom15to16DataHist = new RooDataHist("nom15to16dataHist","nom15to16dataHist",
      rooVarP, nom15to16);
  RooHistPdf* nom15to16HistPdf = new RooHistPdf("nom15to16histPdf","nom15to16histPdf",
      rooVarP, *nom15to16DataHist);
  RooDataHist* nom17to18DataHist = new RooDataHist("nom17to18dataHist","nom17to18dataHist",
      rooVarP, nom17to18);
  RooHistPdf* nom17to18HistPdf = new RooHistPdf("nom17to18histPdf","nom17to18histPdf",
      rooVarP, *nom17to18DataHist);

  RooAbsPdf* myExpPdf = bindPdf(myExpo,rooVarP);
  //
  RooHistPdf saturatedModelNoM7to8("saturatedModelNoM7to8","saturatedModel NoM7to8",rooVarP,*nom7to8DataHist);
  RooNLLVar nllVarExpoModelNoM7to8("nllVarExpoModelNoM7to8","NLL expo model NoM7to8",*myExpPdf,*nom7to8DataHist);//,Range("fitRange"));
  RooNLLVar nllVarSatModelNoM7to8("nllVarSatModelNoM7to8","NLL Saturated Model NoM7to8",saturatedModelNoM7to8,*nom7to8DataHist);//,Range("fitRange"));
  //
  RooHistPdf saturatedModelNoM9to10("saturatedModelNoM9to10","saturatedModel NoM9to10",rooVarP,*nom9to10DataHist);
  RooNLLVar nllVarExpoModelNoM9to10("nllVarExpoModelNoM9to10","NLL expo model NoM9to10",*myExpPdf,*nom9to10DataHist);//,Range("fitRange"));
  RooNLLVar nllVarSatModelNoM9to10("nllVarSatModelNoM9to10","NLL Saturated Model NoM9to10",saturatedModelNoM9to10,*nom9to10DataHist);//,Range("fitRange"));
  //
  RooHistPdf saturatedModelNoM11to12("saturatedModelNoM11to12","saturatedModel NoM11to12",rooVarP,*nom11to12DataHist);
  RooNLLVar nllVarExpoModelNoM11to12("nllVarExpoModelNoM11to12","NLL expo model NoM11to12",*myExpPdf,*nom11to12DataHist);//,Range("fitRange"));
  RooNLLVar nllVarSatModelNoM11to12("nllVarSatModelNoM11to12","NLL Saturated Model NoM11to12",saturatedModelNoM11to12,*nom11to12DataHist);//,Range("fitRange"));
  //
  RooHistPdf saturatedModelNoM13to14("saturatedModelNoM13to14","saturatedModel NoM13to14",rooVarP,*nom13to14DataHist);
  RooNLLVar nllVarExpoModelNoM13to14("nllVarExpoModelNoM13to14","NLL expo model NoM13to14",*myExpPdf,*nom13to14DataHist);//,Range("fitRange"));
  RooNLLVar nllVarSatModelNoM13to14("nllVarSatModelNoM13to14","NLL Saturated Model NoM13to14",saturatedModelNoM13to14,*nom13to14DataHist);//,Range("fitRange"));
  //
  RooHistPdf saturatedModelNoM15to16("saturatedModelNoM15to16","saturatedModel NoM15to16",rooVarP,*nom15to16DataHist);
  RooNLLVar nllVarExpoModelNoM15to16("nllVarExpoModelNoM15to16","NLL expo model NoM15to16",*myExpPdf,*nom15to16DataHist);//,Range("fitRange"));
  RooNLLVar nllVarSatModelNoM15to16("nllVarSatModelNoM15to16","NLL Saturated Model NoM15to16",saturatedModelNoM15to16,*nom15to16DataHist);//,Range("fitRange"));
  //
  //RooHistPdf saturatedModelNoM17to18("saturatedModelNoM17to18","saturatedModel NoM17to18",rooVarP,*nom17to18DataHist);
  //RooNLLVar nllVarExpoModelNoM17to18("nllVarExpoModelNoM17to18","NLL expo model NoM17to18",*myExpPdf,*nom17to18DataHist);//,Range("fitRange"));
  //RooNLLVar nllVarSatModelNoM17to18("nllVarSatModelNoM17to18","NLL Saturated Model NoM17to18",saturatedModelNoM17to18,*nom17to18DataHist);//,Range("fitRange"));

  double satValNoM7to8 = 2*(nllVarExpoModelNoM7to8.getVal()-nllVarSatModelNoM7to8.getVal());
  int numBinsNoM7to8 = nom7to8DataHist->numEntries();
  //std::cout << "NoM7-8: nllVarExpoModel = " << nllVarExpoModel.getVal() << " nllVarSatModel: " <<
  //  nllVarSatModel.getVal() << std::endl;
  std::cout << "NoM7to8: satVal: " << satValNoM7to8 << " numBins: " << numBinsNoM7to8
    << " satVal/numBins = " << satValNoM7to8/numBinsNoM7to8 << std::endl;

  RooPlot* momPlot = rooVarP.frame();
  momPlot->SetAxisRange(fitRangeMin,fitRangeMax);
  nom7to8DataHist->plotOn(momPlot);
  myExpPdf->plotOn(momPlot);
  //saturatedModel.plotOn(momPlot);
  momPlot->Draw();
  t->Print((fileNameBeg+"nom7to8_plot.C").c_str());
  //
  momPlot = rooVarP.frame();
  momPlot->SetAxisRange(fitRangeMin,fitRangeMax);
  nom9to10DataHist->plotOn(momPlot);
  myExpPdf->plotOn(momPlot);
  momPlot->Draw();
  t->Print((fileNameBeg+"nom9to10_plot.C").c_str());
  //
  momPlot = rooVarP.frame();
  momPlot->SetAxisRange(fitRangeMin,fitRangeMax);
  nom11to12DataHist->plotOn(momPlot);
  myExpPdf->plotOn(momPlot);
  momPlot->Draw();
  //t->Print((fileNameBeg+"nom11to12_plot.png").c_str());
  t->Print((fileNameBeg+"nom11to12_plot.C").c_str());
  //
  momPlot = rooVarP.frame();
  momPlot->SetAxisRange(fitRangeMin,fitRangeMax);
  nom13to14DataHist->plotOn(momPlot);
  myExpPdf->plotOn(momPlot);
  momPlot->Draw();
  //t->Print((fileNameBeg+"nom13to14_plot.png").c_str());
  t->Print((fileNameBeg+"nom13to14_plot.C").c_str());
  //
  momPlot = rooVarP.frame();
  momPlot->SetAxisRange(fitRangeMin,fitRangeMax);
  nom15to16DataHist->plotOn(momPlot);
  myExpPdf->plotOn(momPlot);
  momPlot->Draw();
  //t->Print((fileNameBeg+"nom15to16_plot.png").c_str());
  t->Print((fileNameBeg+"nom15to16_plot.C").c_str());

  double satValNoM9to10 = 2*(nllVarExpoModelNoM9to10.getVal()-nllVarSatModelNoM9to10.getVal());
  int numBinsNoM9to10 = nom9to10DataHist->numEntries();
  std::cout << "NoM9to10: satVal: " << satValNoM7to8 << " numBins: " << numBinsNoM7to8
    << " satVal/numBins = " << satValNoM7to8/numBinsNoM7to8 << std::endl;
  double satValNoM11to12 = 2*(nllVarExpoModelNoM11to12.getVal()-nllVarSatModelNoM11to12.getVal());
  int numBinsNoM11to12 = nom11to12DataHist->numEntries();
  std::cout << "NoM11to12: satVal: " << satValNoM7to8 << " numBins: " << numBinsNoM7to8
    << " satVal/numBins = " << satValNoM7to8/numBinsNoM7to8 << std::endl;
  double satValNoM13to14 = 2*(nllVarExpoModelNoM13to14.getVal()-nllVarSatModelNoM13to14.getVal());
  int numBinsNoM13to14 = nom13to14DataHist->numEntries();
  std::cout << "NoM13to14: satVal: " << satValNoM13to14 << " numBins: " << numBinsNoM13to14
    << " satVal/numBins = " << satValNoM13to14/numBinsNoM13to14 << std::endl;
  double satValNoM15to16 = 2*(nllVarExpoModelNoM15to16.getVal()-nllVarSatModelNoM15to16.getVal());
  int numBinsNoM15to16 = nom15to16DataHist->numEntries();
  std::cout << "NoM15to16: satVal: " << satValNoM15to16 << " numBins: " << numBinsNoM15to16
    << " satVal/numBinsNoM15to16 = " << satValNoM15to16/numBinsNoM15to16 << std::endl;
  //double satValNoM17to18 = 2*(nllVarExpoModelNoM17to18.getVal()-nllVarSatModelNoM17to18.getVal());
  //int numBinsNoM17to18 = nom17to18DataHist->numEntries();
  //std::cout << "NoM17to18: satVal: " << satValNoM17to18 << " numBins: " << numBinsNoM17to18
  //  << " satVal/numBins = " << satValNoM17to18/numBinsNoM17to18 << std::endl;

  // make graph of satVals
  float lowerNoMArr[5] = {7,9,11,13,15};
  float satValArr[5] = {satValNoM7to8,satValNoM9to10,satValNoM11to12,satValNoM13to14,satValNoM15to16};
  float satValOverNumBinsArr[5] = {satValArr[0]/(float)numBinsNoM7to8,
    satValArr[1]/(float)numBinsNoM9to10,satValArr[2]/(float)numBinsNoM11to12,
    satValArr[3]/(float)numBinsNoM13to14,satValArr[4]/(float)numBinsNoM15to16};

  t->Clear();
  TGraph* satValGraph = new TGraph(5,lowerNoMArr,satValArr);
  satValGraph->SetName("satValGraph");
  satValGraph->SetMarkerColor(kRed);
  satValGraph->SetLineColor(kRed);
  TGraph* satValNumBinsGraph = new TGraph(5,lowerNoMArr,satValOverNumBinsArr);
  satValNumBinsGraph->SetName("satValNumBinsGraph");
  satValNumBinsGraph->SetTitle("satVal/NumBins;lowerNoM");
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry("satValGraph","saturated value","l");
  leg->AddEntry("satValNumBinsGraph","saturated value/numBins","l");

  //TMultiGraph* mg = new TMultiGraph();
  //mg->SetTitle("Sat. val and sat. val divided by num bins;lowerNoM");
  //mg->Add(satValGraph);
  //mg->Add(satValNumBinsGraph);
  //mg->Draw("ap");
  //leg->Draw("same");
  t->SetLogy(0);
  satValNumBinsGraph->Draw("ap");

  //t->Print((fileNameBeg+"satVal_plot.png").c_str());
  t->Print((fileNameBeg+"satVal_plot.C").c_str());

}
