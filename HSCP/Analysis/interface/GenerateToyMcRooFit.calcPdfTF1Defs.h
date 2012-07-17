#include "TF1.h"
#include "TMath.h"

// **************************************************************************************
// Create my fit function -- Pt/Nom slice 5
Double_t fitTailSlice5(Double_t *x, Double_t *par)
{
  double splicePoint1 = 3.5;
  double splicePoint2 = 4.5;
  double fitVal = 0;

  if(x[0] < splicePoint1)
  {
    double gaus1 = par[3]*TMath::Gaus(splicePoint1,par[4],par[5]);
    double exp1 = TMath::Exp(par[2]+par[1]*splicePoint1+par[0]*splicePoint1*splicePoint1);
    double par7 = (TMath::Log(exp1-gaus1)-par[6])/splicePoint1;
    fitVal = par[3]*TMath::Gaus(x[0],par[4],par[5])+TMath::Exp(par[6]+par7*x[0]);
  }
  else if(x[0] < splicePoint2)
  {
    fitVal = TMath::Exp(par[2]+par[1]*x[0]+par[0]*x[0]*x[0]);
  }
  else
  {
    // constrained exp: exp^(alpha*x+beta)
    double alpha = 2*par[0]*splicePoint2+par[1];
    double beta = -1*par[0]*splicePoint2*splicePoint2+par[2];
    fitVal = TMath::Exp(alpha*x[0]+beta);
  }
  return fitVal;
}

// **************************************************************************************
// Create my fit function -- Pt/Nom slice 4
Double_t fitTailSlice4(Double_t *x, Double_t *par)
{
  double splicePoint1 = 3.5;
  double splicePoint2 = 5;
  double fitVal = 0;

  if(x[0] < splicePoint1)
  {
    double gaus1 = par[3]*TMath::Gaus(splicePoint1,par[4],par[5]);
    double exp1 = TMath::Exp(par[2]+par[1]*splicePoint1+par[0]*splicePoint1*splicePoint1);
    double par7 = (TMath::Log(exp1-gaus1)-par[6])/splicePoint1;
    fitVal = par[3]*TMath::Gaus(x[0],par[4],par[5])+TMath::Exp(par[6]+par7*x[0]);
  }
  else if(x[0] < splicePoint2)
  {
    fitVal = TMath::Exp(par[2]+par[1]*x[0]+par[0]*x[0]*x[0]);
  }
  else
  {
    // constrained exp: exp^(alpha*x+beta)
    double alpha = 2*par[0]*splicePoint2+par[1];
    double beta = -1*par[0]*splicePoint2*splicePoint2+par[2];
    fitVal = TMath::Exp(alpha*x[0]+beta);
  }
  return fitVal;
}

// **************************************************************************************
// Create my fit function -- Pt/Nom slice 3
Double_t fitTailSlice3(Double_t *x, Double_t *par)
{
  double splicePoint1 = 5.;
  double splicePoint2 = 7;
  double fitVal = 0;

  if(x[0] < splicePoint1)
  {
    double gaus1 = par[3]*TMath::Gaus(splicePoint1,par[4],par[5]);
    double exp1 = TMath::Exp(par[2]+par[1]*splicePoint1+par[0]*splicePoint1*splicePoint1);
    double par7 = (TMath::Log(exp1-gaus1)-par[6])/splicePoint1;
    fitVal = par[3]*TMath::Gaus(x[0],par[4],par[5])+TMath::Exp(par[6]+par7*x[0]);
  }
  else if(x[0] < splicePoint2)
  {
    fitVal = TMath::Exp(par[2]+par[1]*x[0]+par[0]*x[0]*x[0]);
  }
  else
  {
    // constrained exp: exp^(alpha*x+beta)
    double alpha = 2*par[0]*splicePoint2+par[1];
    double beta = -1*par[0]*splicePoint2*splicePoint2+par[2];
    fitVal = TMath::Exp(alpha*x[0]+beta);
  }
  return fitVal;
}

// **************************************************************************************
// Create my fit function -- Pt/Nom slice 2
Double_t fitTailSlice2(Double_t *x, Double_t *par)
{
  double splicePoint1 = 5.5;
  double splicePoint2 = 9;
  double fitVal = 0;

  if(x[0] < splicePoint1)
  {
    double gaus1 = par[3]*TMath::Gaus(splicePoint1,par[4],par[5]);
    double exp1 = TMath::Exp(par[2]+par[1]*splicePoint1+par[0]*splicePoint1*splicePoint1);
    double par7 = (TMath::Log(exp1-gaus1)-par[6])/splicePoint1;
    fitVal = par[3]*TMath::Gaus(x[0],par[4],par[5])+TMath::Exp(par[6]+par7*x[0]);
  }
  else if(x[0] < splicePoint2)
  {
    fitVal = TMath::Exp(par[2]+par[1]*x[0]+par[0]*x[0]*x[0]);
  }
  else
  {
    // constrained exp: exp^(alpha*x+beta)
    double alpha = 2*par[0]*splicePoint2+par[1];
    double beta = -1*par[0]*splicePoint2*splicePoint2+par[2];
    fitVal = TMath::Exp(alpha*x[0]+beta);
  }
  return fitVal;
}

