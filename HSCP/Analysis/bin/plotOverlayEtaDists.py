#/usr/bin/env python

import array
import math
import datetime

import ROOT
from ROOT import TFile, kBlue, TCanvas, TGraph, TLegend, gROOT, TPaveText

myZFile = TFile('may4.gluino600z.makeHSCParticlePlots.root')
myNormalFile = TFile('may4.gluino600.makeHSCParticlePlots.root')
modelName = "gluino600"

def DrawPreliminary(Lumi,  X1=0.15,  Y1=0.995,  X2=0.5, Y2=0.945):
  text = TPaveText(X1,Y1,X2,Y2, "NDC")
  text.SetFillColor(0)
  text.SetTextSize(0.0275)
  text.SetBorderSize(0)
  if(Lumi<0):
    text.AddText("CMS Preliminary 2011,  #sqrt{s} 7 TeV")
  if(Lumi>0):
    text.AddText("CMS Preliminary 2011, "+str(Lumi)+" fb^{-1}    #sqrt{s} = 7 TeV   ")
  return text

def DrawLegend(mczHist,mcHist):
  leg = TLegend(0.75,0.75,0.97,0.95)
  leg.SetBorderSize(0)
  leg.AddEntry(mczHist,"Gluino 600 (Z2 tune)","pl")
  leg.AddEntry(mcHist,"Gluino 600 (D6T tune)","pl")
  return leg

# Load/apply TDR style
gROOT.ProcessLine('.L tdrStyle.C')
gROOT.ProcessLine('setTDRStyle()')


canvas = TCanvas()
canvas.cd()

etaHistZ = myZFile.Get("hscpEtaHist")
etaHistZ.Sumw2()
etaHistZ.Rebin()
integral = etaHistZ.Integral()
etaHistZ.Scale(1/integral)
etaHistZ.SetMarkerColor(kBlue)
etaHistZ.SetMarkerStyle(23)
etaHistZ.SetMarkerSize(1)
etaHistZ.SetLineColor(kBlue)
etaHistZ.SetStats(0)
title = "HSCP, "+modelName+" and "+modelName+"Z (blue);#eta"
etaHistZ.SetTitle(title)
etaHistZ.SetTitle("")
etaHistZ.GetYaxis().SetTitle("A.U.")
etaHistZ.GetYaxis().SetTitleSize(0.045)
etaHistZ.GetYaxis().SetTitleOffset(1.75)
etaHistZ.Draw()

etaHistNormal = myNormalFile.Get("hscpEtaHist")
etaHistNormal.Sumw2()
etaHistNormal.Rebin()
integralNormal = etaHistNormal.Integral()
etaHistNormal.Scale(1/integralNormal)
etaHistNormal.SetMarkerStyle(22)
etaHistNormal.SetMarkerSize(1)
etaHistNormal.SetStats(0)
etaHistNormal.Draw("sames")
leg = DrawLegend(etaHistZ,etaHistNormal)
leg.Draw("same")
text = DrawPreliminary(4976)
text.Draw()

#canvas2 = TCanvas()
#canvas2.cd()
#
#etaVals = []
#diffOverErrs = []
#diffOverNorms = []
#
#for bin in range(1,etaHistNormal.GetNbinsX()+1):
#  etaVal = etaHistZ.GetBinCenter(bin)
#  if(etaVal < -1.6 or etaVal > 1.6):
#    continue
#  bincZ = etaHistZ.GetBinContent(bin)
#  if(bincZ <= 0):
#    continue
#  bincZerr = etaHistZ.GetBinError(bin)
#  bincNorm = etaHistNormal.GetBinContent(bin)
#  bincNormErr = etaHistNormal.GetBinError(bin)
#  diffOverErr = (bincZ-bincNorm)/bincZerr
#  diffOverErrs.append(diffOverErr)
#  diffOverNorm = math.fabs((bincZ-bincNorm)/bincNorm)
#  diffOverNorms.append(diffOverNorm)
#  etaVals.append(etaVal)
#
#graph = TGraph(len(etaVals),array.array("f",etaVals),array.array("f",diffOverNorms))
#graph.SetTitle("abs(z-normal)/normal;#eta")
#graph.Draw("ap")

### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]

# print
now = datetime.datetime.now()
Date = now.strftime("%b%d")
fileNamePng = Date+'.mcMCZ.overlayEta.png'
fileNameEps = Date+'.mcMCZ.overlayEta.eps'
fileNameC = Date+'.mcMCZ.overlayEta.C'
canvas.Print(fileNamePng)
canvas.Print(fileNameEps)
canvas.Print(fileNameC)


myZFile.Close()
myNormalFile.Close()
