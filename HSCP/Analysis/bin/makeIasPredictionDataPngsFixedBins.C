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


makeIasPredictionDataPngsFixedBins()
{
  TCanvas t("myCanvas","myCanvas",100,100,800,600);
  t.cd();
  t.SetLogy();

  for(int nom=5; nom < 19; nom+=2)
  {
    for(int etaSlice=0; etaSlice < 16; etaSlice+=2)
    {
      string histName="nom";
      histName+=intToString(nom);
      histName+="to";
      if(nom==5)
        histName+=intToString(nom+3);
      else
        histName+=intToString(nom+1);
      histName+="eta";
      histName+=intToString(etaSlice);
      histName+="to";
      histName+=intToString(etaSlice+2);
      string histNameBeg = histName;
      histName+="iasPredictionFixedLimitsHist";
      //histName+="iasPredictionFixedDiscoveryHist";

      //cout << "Getting hist: " << histName << endl;
      TH1D* iasPredHistThisSlice = (TH1D*)_file0->Get(histName.c_str());

      histName="nom";
      histName+=intToString(nom);
      histName+="to";
      if(nom==5)
        histName+=intToString(nom+3);
      else
        histName+=intToString(nom+1);
      histName+="eta";
      histName+=intToString(etaSlice);
      histName+="to";
      histName+=intToString(etaSlice+2);
      histName+="dRegionFixedHist";

      //cout << "Getting hist: " << histName << endl;
      TH1D* iasDataHistThisSlice = (TH1D*)_file0->Get(histName.c_str());

      iasDataHistThisSlice->SetLineColor(kBlack);
      iasDataHistThisSlice->SetMarkerColor(kBlack);
      iasDataHistThisSlice->SetMarkerStyle(20);
      iasDataHistThisSlice->Draw();
      iasPredHistThisSlice->SetLineColor(kBlue);
      iasPredHistThisSlice->SetMarkerColor(kBlue);
      iasPredHistThisSlice->SetMarkerStyle(21);
      iasPredHistThisSlice->Draw("same");
      // draw 3rd hist?
      string thirdHistName = histNameBeg;
      thirdHistName+="iasPredictionFixedSlopeMinusOneSigmaLimitsHist";
      TH1D* iasPredThirdHist = (TH1D*)_file0->Get(thirdHistName.c_str());
      iasPredThirdHist->SetLineColor(kRed);
      iasPredThirdHist->SetMarkerColor(kRed);
      iasPredThirdHist->SetMarkerStyle(23);
      iasPredThirdHist->Draw("same");
      // draw 4th hist?
      string fourthHistName = histNameBeg;
      fourthHistName+="iasPredictionFixedSlopeAddOneSigmaLimitsHist";
      TH1D* iasPredFourthHist = (TH1D*)_file0->Get(fourthHistName.c_str());
      iasPredFourthHist->SetLineColor(kGreen+2);
      iasPredFourthHist->SetMarkerColor(kGreen+2);
      iasPredFourthHist->SetMarkerStyle(22);
      iasPredFourthHist->Draw("same");

      if(iasDataHistThisSlice->GetMaximum() <= 0)
        iasDataHistThisSlice->GetYaxis()->SetRangeUser(1e-10,iasPredHistThisSlice->GetMaximum()*1.2);
      else
        iasDataHistThisSlice->GetYaxis()->SetRangeUser(1e-10,iasDataHistThisSlice->GetMaximum()*1.2);

      string name="nom";
      name+=intToString(nom);
      name+="to";
      if(nom==5)
        name+=intToString(nom+3);
      else
        name+=intToString(nom+1);
      name+="eta";
      name+=intToString(etaSlice);
      name+="to";
      name+=intToString(etaSlice+2);
      name+="iasDataPred";
      t.Print((name+".png").c_str());
    }
    if(nom==5) nom+=2;
  }
}
