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


makeIasPredictionPngsFixedBins()
{
  TCanvas t("myCanvas","myCanvas",100,100,800,600);
  t.cd();
  t.SetLogy();

  //for(int nom=5; nom < 23; nom+=2)
  for(int nom=5; nom < 18; nom+=2)
  {
    for(int etaSlice=0; etaSlice < 16; etaSlice+=2)
    {
      string histName="nom";
      histName+=intToString(nom);
      histName+="to";
      if(nom==5)
        histName+=intToString(nom+2);
      else
        histName+=intToString(nom+1);
      histName+="eta";
      histName+=intToString(etaSlice);
      histName+="to";
      histName+=intToString(etaSlice+2);
      histName+="iasPointsAndFitHist";

      string fullPath = "iasPredictionsFixedBins/";
      fullPath+=histName;
      TH1D* iasPredHistThisSlice = (TH1D*)_file0->Get(fullPath.c_str());

      iasPredHistThisSlice->SetLineColor(kBlue);
      iasPredHistThisSlice->SetMarkerColor(kBlue);
      iasPredHistThisSlice->SetMarkerStyle(22);
      iasPredHistThisSlice->Draw();
      string name="nom";
      name+=intToString(nom);
      name+="to";
      if(nom==5)
        name+=intToString(nom+2);
      else
        name+=intToString(nom+1);
      name+="eta";
      name+=intToString(etaSlice);
      name+="to";
      name+=intToString(etaSlice+2);
      name+="iasPred";
      t.Print((name+".png").c_str());
    }
    if(nom==5) nom+=2;
  }
}
