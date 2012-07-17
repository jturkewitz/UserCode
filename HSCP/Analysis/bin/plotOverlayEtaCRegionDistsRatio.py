#/usr/bin/env python

import array
import math
import datetime

import ROOT
from ROOT import TFile,kBlue,kBlack,TCanvas,TPad,TGraph,gROOT,TPaveText,TLegend

myMCFile = TFile('may7.2.mcDYToMuMu.makeHSCParticlePlots.root')
myDataFile = TFile('may7.2.data2011fewRuns.makeHSCParticlePlots.root')

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

etaHistMC = myMCFile.Get("etaCRegionHist")
etaHistMC.Sumw2()
integral = etaHistMC.Integral()
etaHistMC.Scale(1/integral)
etaHistMC.SetTitle("")
etaHistMC.GetXaxis().SetRangeUser(-1.6,1.6)
etaHistMC.GetYaxis().SetRangeUser(0.015,0.033)

etaHistData = myDataFile.Get("etaCRegionHist")
etaHistData.Sumw2()
etaHistData.SetTitle("")
integralData = etaHistData.Integral()
etaHistData.Scale(1/integralData)
etaHistData.SetMarkerStyle(22)
etaHistData.GetXaxis().SetRangeUser(-1.6,1.6)

canvas2 = TCanvas()
canvas2.cd()
pad1 = TPad("pad1","pad1",0,0.3,1,0.945)
pad1.SetBottomMargin(0)
pad1.GetFrame().SetBorderSize(4)
pad1.Draw()
pad1.cd()
etaHistMC.SetStats(0)
etaHistData.SetStats(0)
etaHistData.SetMarkerStyle(20)
etaHistData.SetMarkerSize(0.8)
etaHistData.SetMarkerColor(kBlack)
etaHistData.SetLineColor(kBlack)
etaHistMC.SetFillColor(26)
etaHistMC.SetLineColor(kBlack)
etaHistMC.GetYaxis().SetTitleSize(0.05)
etaHistMC.GetYaxis().SetTitleOffset(1.25)
etaHistMC.GetYaxis().SetTitle("A.U.")
etaHistMC.Draw("HIST")
histDataCopy = etaHistData.DrawCopy("elpx0SAME")
leg = DrawLegend(etaHistData,etaHistMC)
leg.Draw("same")
canvas2.cd()
pad2 = TPad("pad2","pad2",0,0,1,0.3)
pad2.SetTopMargin(0)
pad2.GetFrame().SetBorderSize(4)
pad2.Draw()
pad2.cd()
etaHistData.Divide(etaHistMC)
etaHistData.GetXaxis().SetTitle("#eta")
etaHistData.GetYaxis().SetTitle("Data/MC")
etaHistData.GetYaxis().SetTitleSize(0.1)
etaHistData.GetYaxis().SetTitleOffset(0.4)
etaHistData.SetMaximum(1.135)
etaHistData.Draw("ex0")
canvas2.cd()
text = DrawPreliminary(4976)
text.Draw()

# print
now = datetime.datetime.now()
Date = now.strftime("%b%d")
fileNamePng = Date+'.mcData.overlayEta.png'
fileNameEps = Date+'.mcData.overlayEta.eps'
fileNameC = Date+'.mcData.overlayEta.C'
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
