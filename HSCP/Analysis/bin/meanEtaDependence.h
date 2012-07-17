#include "TF1.h"
#include "TMath.h"

// Create my fit function
Double_t fitMeanEtaDependence(Double_t *x, Double_t *par)
{
  double splicePoint1 = 1.05;
  double splicePoint2 = 1.4;
  double fitVal = 0;

  // [0]+[1]*cos([2]*x); ap + par[3]*x; app + par[4]*x*par[5]*x^2+par[6]*x^3+par[7]*x^4

  double ap = par[0]+par[1]*TMath::Cos(par[2]*splicePoint1)-par[3]*splicePoint1;
  double app = ap+par[3]*splicePoint2-par[4]*splicePoint2-par[5]*pow(splicePoint2,2)-par[6]*pow(splicePoint2,3)-par[7]*pow(splicePoint2,4);

  if(x[0] < splicePoint1)
  {
    fitVal = par[0]+par[1]*TMath::Cos(par[2]*x[0]);
  }
  else if(x[0] < splicePoint2)
  {
    fitVal = ap+par[3]*x[0];
  }
  else
  {
    fitVal = app+par[4]*x[0]+par[5]*x[0]*x[0]+par[6]*x[0]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0]*x[0];
  }

  return fitVal;
}

