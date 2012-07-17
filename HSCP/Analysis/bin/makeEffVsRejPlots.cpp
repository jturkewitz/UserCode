#include "TH1F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMultiGraph.h"

#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <math.h>

std::string intToString(int num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}

std::string floatToString(float num)
{
  using namespace std;
  ostringstream myStream;
  myStream << num << flush;
  return(myStream.str()); //returns the string form of the stringstream object
}

int main (int argc, char** argv)
{
  using namespace std;
  using namespace RooFit;

  char* inFileSig = argv[1];
  char* inFileBack = argv[2];

  if(!inFileSig)
  {
    cout << " No input signal file specified !" << endl;
    cout << "Usage: makeEffVsRejPlots signalFile backgroundFile" << endl;
    return -1;
  }
  if(!inFileBack)
  {
    cout << " No input background file specified !" << endl;
    cout << "Usage: makeEffVsRejPlots signalFile backgroundFile" << endl;
    return -1;
  }

  // RooFit observables and dataset
  RooRealVar rooVarIas("rooVarIas","ias",0,1);
  RooRealVar rooVarIp("rooVarIp","ip",0,1);
  RooRealVar rooVarIh("rooVarIh","ih",0,15);
  RooRealVar rooVarP("rooVarP","p",0,5000);
  RooRealVar rooVarPt("rooVarPt","pt",0,5000);
  RooRealVar rooVarNoMias("rooVarNoMias","nom",0,30);
  RooRealVar rooVarEta("rooVarEta","eta",0,2.5);
  TFile* tFileBack = TFile::Open(inFileBack);
  if(!tFileBack)
  {
    std::cout << "ERROR: input file " << inFileBack << " not found." << std::endl;
    return -2;
  }
  RooDataSet* rooDataSetBack1 = (RooDataSet*)tFileBack->Get("rooDataSetAll");
  if(rooDataSetBack1->numEntries() < 1)
  {
    std::cout << "Problem with RooDataSet named rooDataSetAll in file " << inFileBack << std::endl;
    return -3;
  }
  TFile* tFileSig = TFile::Open(inFileSig);
  if(!tFileSig)
  {
    std::cout << "ERROR: input file " << inFileSig << " not found." << std::endl;
    return -2;
  }
  RooDataSet* rooDataSetSig1 = (RooDataSet*)tFileSig->Get("rooDataSetAll");
  if(rooDataSetSig1->numEntries() < 1)
  {
    std::cout << "Problem with RooDataSet named rooDataSetAll in file " << inFileSig << std::endl;
    return -3;
  }

  RooDataSet* bRegionBackgroundDataSet1 =
    (RooDataSet*) rooDataSetBack1->reduce(Cut("rooVarP<100"));
  RooDataSet* bRegionBackgroundDataSet2 =
    (RooDataSet*) bRegionBackgroundDataSet1->reduce(Cut("rooVarIh>3"));
  RooDataSet* bRegionBackgroundDataSet =
    (RooDataSet*) bRegionBackgroundDataSet2->reduce(Cut("rooVarEta<1.6"));

  RooDataSet* rooDataSetSig2 = (RooDataSet*) rooDataSetSig1->reduce(("rooVarEta<1.6"));
  RooDataSet* rooDataSetSig = (RooDataSet*) rooDataSetSig2->reduce(("rooVarIh>3"));
  //const RooArgSet* argSet_backB = bRegionBackgroundDataSet->get();
  //RooRealVar* ihData_backB = (RooRealVar*)argSet_backB->find(rooVarIh.GetName());
  //RooRealVar* iasData_backB = (RooRealVar*)argSet_backB->find(rooVarIas.GetName());
  ////RooRealVar* ipData_backB = (RooRealVar*)argSet_backB->find(rooVarIp.GetName());
  ////RooRealVar* nomData_backB = (RooRealVar*)argSet_backB->find(rooVarNoMias.GetName());
  ////RooRealVar* etaData_backB = (RooRealVar*)argSet_backB->find(rooVarEta.GetName());
  //const RooArgSet* argSet_sig = rooDataSetSig->get();
  //RooRealVar* ihData_sig = (RooRealVar*)argSet_sig->find(rooVarIh.GetName());
  //RooRealVar* iasData_sig = (RooRealVar*)argSet_sig->find(rooVarIas.GetName());
  ////RooRealVar* ipData_sig = (RooRealVar*)argSet_sig->find(rooVarIp.GetName());
  ////RooRealVar* nomData_sig = (RooRealVar*)argSet_sig->find(rooVarNoMias.GetName());
  ////RooRealVar* etaData_sig = (RooRealVar*)argSet_sig->find(rooVarEta.GetName());

  // output
  TFile* outputTFile = new TFile("effVsRejPlots.root","recreate");
  outputTFile->cd();
  TGraph* effVsCutValIasGraph;
  TGraph* effVsCutValIpGraph;
  TGraph* effVsCutValIhGraph;
  TGraph* rejVsCutValIasGraph;
  TGraph* rejVsCutValIpGraph;
  TGraph* rejVsCutValIhGraph;
  // Graphs for eff vs rej
  TGraph* efficiencyVsRejectionIasGraph;
  TGraph* efficiencyVsRejectionIpGraph;
  TGraph* efficiencyVsRejectionIhGraph;

  // calculate quantities
  vector<double> cutValsIas;
  vector<double> cutValsIp;
  vector<double> cutValsIh;
  vector<double> effByCutValIas;
  vector<double> effByCutValIp;
  vector<double> effByCutValIh;
  vector<double> rejByCutValIas;
  vector<double> rejByCutValIp;
  vector<double> rejByCutValIh;

  int numPointsInBRegionBackground = bRegionBackgroundDataSet->numEntries();
  int numPointsInSignal = rooDataSetSig->numEntries();
  cout << "num initial entries background b region: " << numPointsInBRegionBackground
    << " and signal (all): " << numPointsInSignal << endl;

  // Ias
  for(float iasCut=0; iasCut<1; iasCut+=0.01) // loop over cut values
  {
    string iasCutString = "rooVarIas>";
    iasCutString+=floatToString(iasCut);
    cutValsIas.push_back(iasCut);

    RooDataSet* thisCutDataSetBack = (RooDataSet*)(bRegionBackgroundDataSet->reduce(Cut(iasCutString.c_str())));
    if(thisCutDataSetBack->numEntries() > 0)
      rejByCutValIas.push_back(-1*log(thisCutDataSetBack->numEntries()/(float)numPointsInBRegionBackground));

    RooDataSet* thisCutDataSetSig = (RooDataSet*)(rooDataSetSig->reduce(Cut(iasCutString.c_str())));
    effByCutValIas.push_back(thisCutDataSetSig->numEntries()/(float)numPointsInSignal);

    delete thisCutDataSetBack;
    delete thisCutDataSetSig;
  }
  // Ip
  for(float ipCut=0; ipCut<1; ipCut+=0.01) // loop over cut values
  {
    string ipCutString = "rooVarIp>";
    ipCutString+=floatToString(ipCut);
    cutValsIp.push_back(ipCut);

    RooDataSet* thisCutDataSetBack = (RooDataSet*)(bRegionBackgroundDataSet->reduce(Cut(ipCutString.c_str())));
    if(thisCutDataSetBack->numEntries() > 0)
      rejByCutValIp.push_back(-1*log(thisCutDataSetBack->numEntries()/(float)numPointsInBRegionBackground));

    RooDataSet* thisCutDataSetSig = (RooDataSet*)(rooDataSetSig->reduce(Cut(ipCutString.c_str())));
    effByCutValIp.push_back(thisCutDataSetSig->numEntries()/(float)numPointsInSignal);

    delete thisCutDataSetBack;
    delete thisCutDataSetSig;
  }
  // Ih
  for(float ihCut=2.8; ihCut<10; ihCut+=0.05) // loop over cut values
  {
    string ihCutString = "rooVarIh>";
    ihCutString+=floatToString(ihCut);
    cutValsIh.push_back(ihCut);

    RooDataSet* thisCutDataSetBack = (RooDataSet*)(bRegionBackgroundDataSet->reduce(Cut(ihCutString.c_str())));
    if(thisCutDataSetBack->numEntries() > 0)
      rejByCutValIh.push_back(-1*log(thisCutDataSetBack->numEntries()/(float)numPointsInBRegionBackground));

    RooDataSet* thisCutDataSetSig = (RooDataSet*)(rooDataSetSig->reduce(Cut(ihCutString.c_str())));
    effByCutValIh.push_back(thisCutDataSetSig->numEntries()/(float)numPointsInSignal);

    delete thisCutDataSetBack;
    delete thisCutDataSetSig;
  }

  // Make the graphs
  effVsCutValIasGraph = new TGraph(effByCutValIas.size(),&(*cutValsIas.begin()),&(*effByCutValIas.begin()));
  effVsCutValIasGraph->SetName("effVsCutIasGraph");
  effVsCutValIasGraph->SetTitle("Signal efficiency vs. Ias cut");
  rejVsCutValIasGraph = new TGraph(rejByCutValIas.size(),&(*cutValsIas.begin()),&(*rejByCutValIas.begin()));
  rejVsCutValIasGraph->SetName("rejByCutIasGraph");
  rejVsCutValIasGraph->SetTitle("-Log(rej) Vs Cut Ias");
  efficiencyVsRejectionIasGraph = new TGraph(rejByCutValIas.size(),&(*rejByCutValIas.begin()),&(*effByCutValIas.begin()));
  efficiencyVsRejectionIasGraph->SetName("efficiencyVsRejectionIasGraph");
  efficiencyVsRejectionIasGraph->SetTitle("Signal efficiency vs. -log(rej)");
  // ip
  effVsCutValIpGraph = new TGraph(effByCutValIp.size(),&(*cutValsIp.begin()),&(*effByCutValIp.begin()));
  effVsCutValIpGraph->SetName("effVsCutIpGraph");
  effVsCutValIpGraph->SetTitle("Signal efficiency vs. Ip cut");
  rejVsCutValIpGraph = new TGraph(rejByCutValIp.size(),&(*cutValsIp.begin()),&(*rejByCutValIp.begin()));
  rejVsCutValIpGraph->SetName("rejByCutIpGraph");
  rejVsCutValIpGraph->SetTitle("-Log(rej) Vs Cut Ip");
  efficiencyVsRejectionIpGraph = new TGraph(rejByCutValIp.size(),&(*rejByCutValIp.begin()),&(*effByCutValIp.begin()));
  efficiencyVsRejectionIpGraph->SetName("efficiencyVsRejectionIpGraph");
  efficiencyVsRejectionIpGraph->SetTitle("Signal efficiency vs. -log(rej)");
  // ih
  effVsCutValIhGraph = new TGraph(effByCutValIh.size(),&(*cutValsIh.begin()),&(*effByCutValIh.begin()));
  effVsCutValIhGraph->SetName("effVsCutIhGraph");
  effVsCutValIhGraph->SetTitle("Signal efficiency vs. Ih cut");
  rejVsCutValIhGraph = new TGraph(rejByCutValIh.size(),&(*cutValsIh.begin()),&(*rejByCutValIh.begin()));
  rejVsCutValIhGraph->SetName("rejByCutIhGraph");
  rejVsCutValIhGraph->SetTitle("-Log(rej) Vs Cut Ih");
  efficiencyVsRejectionIhGraph = new TGraph(rejByCutValIh.size(),&(*rejByCutValIh.begin()),&(*effByCutValIh.begin()));
  efficiencyVsRejectionIhGraph->SetName("efficiencyVsRejectionIhGraph");
  efficiencyVsRejectionIhGraph->SetTitle("Signal efficiency vs. -log(rej)");

  // multigraph
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(efficiencyVsRejectionIasGraph,"lp+");
  efficiencyVsRejectionIhGraph->SetLineColor(4);
  efficiencyVsRejectionIhGraph->SetMarkerColor(4);
  mg->Add(efficiencyVsRejectionIhGraph,"lp+");
  efficiencyVsRejectionIpGraph->SetLineColor(2);
  efficiencyVsRejectionIpGraph->SetMarkerColor(2);
  mg->Add(efficiencyVsRejectionIpGraph,"lp+");
  mg->SetName("efficiencyVsRejectionOverlayGraph");

  // Write
  outputTFile->cd();
  effVsCutValIasGraph->Write();
  rejVsCutValIasGraph->Write();
  efficiencyVsRejectionIasGraph->Write();
  effVsCutValIhGraph->Write();
  rejVsCutValIhGraph->Write();
  efficiencyVsRejectionIhGraph->Write();
  effVsCutValIpGraph->Write();
  rejVsCutValIpGraph->Write();
  efficiencyVsRejectionIpGraph->Write();
  mg->Write();
}
