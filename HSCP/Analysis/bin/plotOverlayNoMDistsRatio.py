#/usr/bin/env python

import array
import math
import datetime

import ROOT
from ROOT import TFile,kBlue,kBlack,TCanvas,TPad,TGraph,gROOT,TPaveText,TLegend

myMCFile = TFile('may7.mcDYToMuMu.makeHSCParticlePlots.root')
myDataFile = TFile('may7.data2011fewRuns.makeHSCParticlePlots.root')

def DrawPreliminary(Lumi,  X1=0.2,  Y1=0.995,  X2=0.5, Y2=0.93):
  text = TPaveText(X1,Y1,X2,Y2, "NDC")
  text.SetFillColor(0)
  text.SetTextSize(0.0275)
  text.SetBorderSize(0)
  if(Lumi<0):
    text.AddText("CMS Preliminary 2011,  #sqrt{s} 7 TeV")
  if(Lumi>0):
    text.AddText("CMS Preliminary 2011, "+str(Lumi)+" fb^{-1}    #sqrt{s} = 7 TeV   ")
  return text

def DrawLegend(dataHist,mcHist):
  leg = TLegend(0.75,0.7,0.95,0.9)
  leg.SetBorderSize(0)
  leg.AddEntry(dataHist,"data (control sample)","p")
  leg.AddEntry(mcHist,"MC (control sample)","f")
  return leg

# Load/apply TDR style
gROOT.ProcessLine('.L tdrStyle.C')
gROOT.ProcessLine('setTDRStyle()')

nomHistMC = myMCFile.Get("nomCRegionHist")
nomHistMC.Sumw2()
integral = nomHistMC.Integral()
nomHistMC.Scale(1/integral)
#nomHistMC.SetMarkerColor(kBlue)
#nomHistMC.SetMarkerStyle(23)
#nomHistMC.SetLineColor(kBlue)
nomHistMC.SetTitle("")
nomHistMC.GetXaxis().SetRangeUser(5,21)

nomHistData = myDataFile.Get("nomCRegionHist")
nomHistData.Sumw2()
nomHistData.SetTitle("")
integralData = nomHistData.Integral()
nomHistData.Scale(1/integralData)
nomHistData.SetMarkerStyle(22)

canvas2 = TCanvas()
canvas2.cd()
pad1 = TPad("pad1","pad1",0,0.3,1,0.945)
pad1.SetBottomMargin(0)
pad1.Draw()
pad1.cd()
nomHistMC.SetStats(0)
nomHistData.SetStats(0)
nomHistData.SetMarkerStyle(20)
nomHistData.SetMarkerSize(0.8)
nomHistData.SetMarkerColor(kBlack)
nomHistData.SetLineColor(kBlack)
nomHistMC.SetFillColor(26)
nomHistMC.SetLineColor(kBlack)
nomHistMC.GetYaxis().SetTitleSize(0.05)
nomHistMC.GetYaxis().SetTitleOffset(1.2)
nomHistMC.GetYaxis().SetTitle("A.U.")
nomHistMC.Draw("HIST")
histDataCopy = nomHistData.DrawCopy("elpx0SAME")
nomHistMC.SetMinimum(-0.005)
leg = DrawLegend(nomHistData,nomHistMC)
leg.Draw("same")
canvas2.cd()
pad2 = TPad("pad2","pad2",0,0,1,0.3)
pad2.SetTopMargin(0)
pad2.Draw()
pad2.cd()
nomHistData.Divide(nomHistMC)
nomHistData.GetXaxis().SetTitle("Number of measurements")
nomHistData.GetYaxis().SetTitle("Data/MC")
nomHistData.GetYaxis().SetTitleSize(0.1)
nomHistData.GetYaxis().SetTitleOffset(0.4)
nomHistData.Draw("ex0")
canvas2.cd()
text = DrawPreliminary(4976)
text.Draw()

# print
now = datetime.datetime.now()
Date = now.strftime("%b%d")
fileNamePng = Date+'.mcData.overlayNoM.png'
fileNameEps = Date+'.mcData.overlayNoM.eps'
fileNameC = Date+'.mcData.overlayNoM.C'
canvas2.Print(fileNamePng)
canvas2.Print(fileNameEps)
canvas2.Print(fileNameC)

### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]

myMCFile.Close()
myDataFile.Close()
