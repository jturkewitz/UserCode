#include <iostream>
#include <vector>
#include <math.h>

#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TFile.h"

int main (int argc, char** argv)
{
  std::vector<float> betas;
  std::vector<float> betaErrors;
  std::vector<float> means;
  std::vector<float> meanErrors;
  std::vector<float> sigmas;
  std::vector<float> sigmaErrors;

  for(float i=0.25; i < 0.9; i+=0.1)
  {
    betas.push_back(i);
    betaErrors.push_back(0);
  }
  means.push_back(11.87);
  means.push_back(9.73);
  means.push_back(8.327);
  means.push_back(6.899);
  means.push_back(5.554);
  means.push_back(4.374);
  means.push_back(3.505);
  meanErrors.push_back(0.1);
  meanErrors.push_back(0.05);
  meanErrors.push_back(0.03);
  meanErrors.push_back(0.025);
  meanErrors.push_back(0.016);
  meanErrors.push_back(0.015);
  meanErrors.push_back(0.018);

  sigmas.push_back(1.901);
  sigmas.push_back(1.298);
  sigmas.push_back(1.064);
  sigmas.push_back(0.7731);
  sigmas.push_back(0.5676);
  sigmas.push_back(0.4724);
  sigmas.push_back(0.3688);
  sigmaErrors.push_back(0.079);
  sigmaErrors.push_back(0.027);
  sigmaErrors.push_back(0.017);
  sigmaErrors.push_back(0.0151);
  sigmaErrors.push_back(0.0112);
  sigmaErrors.push_back(0.0117);
  sigmaErrors.push_back(0.0139);

  TFile* outFile = new TFile("plotDeDxInBetaSlices.root","recreate");
  outFile->cd();
  TGraphErrors* meanGraph = new TGraphErrors(betas.size(),&(*betas.begin()),&(*means.begin()),&(*betaErrors.begin()),&(*meanErrors.begin()));
  meanGraph->SetName("meanGraph");
  TGraphErrors* sigmaGraph = new TGraphErrors(betas.size(),&(*betas.begin()),&(*sigmas.begin()),&(*betaErrors.begin()),&(*sigmaErrors.begin()));
  sigmaGraph->SetName("sigmaGraph");
  meanGraph->Write();
  sigmaGraph->Write();
  outFile->Close();

}
