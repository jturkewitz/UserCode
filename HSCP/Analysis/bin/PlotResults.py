#!/usr/bin/env python

import string
import os
import sys
import glob
import math
import array
from SignalDefinitions import *
# run in batch
sys.argv.append('-b')
import ROOT
#
from ROOT import TGraph,TMultiGraph,TCanvas,TPaveText,gROOT,TLegend,TCutG,kGreen,TFile,RooStats,Math,vector,gROOT,TLine,kRed,gPad

# at UMN, must use /usr/bin/python (2.4) and my compiled ROOT version (5.32.02)
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_allModels_cutPt50GeVcutIas0.1_allSlices_Jun18'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_50largestDeltaBStatErrors_halfIasHalfNoMAsNominal_highToyStatTestSomeModels_cutPt50GeVcutIas0.1_allSlices_Jun25'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_second50largestDeltaBStatErrors_halfIasHalfNoMAsNominal_highToyStatTestSomeModels_cutPt50GeVcutIas0.1_allSlices_Jun25'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_discBGMinusOneSigma_allModels_cutPt50GeVcutIas0.1_allSlices_Jul3'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_nominal_allModels_cutPt50GeVcutIas0.1_allSlices_Jul4'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_iasMinusOneSigmaShapeShift_allModels_cutPt50GeVcutIas0.1_allSlices_Jul4'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVary_allModels_cutPt50GeVcutIas0.1_allSlices_Jul5'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVaryFIX_allModels_cutPt50GeVcutIas0.1_allSlices_Jul5'
#BaseDir = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVaryFIX_allModelsMoreToys_cutPt50GeVcutIas0.1_allSlices_Jul9'
#BaseDir = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_moreModels_Jul10'
#BaseDir = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_noBGShapeSyst_moreModels_Jul10'
#BaseDir = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_moreModels_Jul10'
BaseDir = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_moreModels_moreToys_Jul11'

runCERN = False
PlotMinScale = 0.0005
PlotMaxScale = 3
GluinoXSecFile = 'gluino_XSec.txt'
StopXSecFile = 'stop_XSec.txt'
StauXSecFile = 'stau_XSec.txt'

gROOT.Reset()

# Load/apply TDR style
gROOT.ProcessLine('.L tdrStyle.C')
gROOT.ProcessLine('setTDRStyle()')

if(runCERN):
  fileEnd = '.log'
else:
  fileEnd = '.out'


# STD ANA INFO
StdAnaInfo = {}
# expLim, obsLim, 5-sigma xsec
# gluino f=0.1
StdAnaInfo['Gluino300'] = [3.8e-3,3.7e-3,8.5e-3]
StdAnaInfo['Gluino400'] = [3.3e-3,3.3e-3,5.5e-3]
StdAnaInfo['Gluino500'] = [3.1e-3,3.1e-3,5.2e-3]
StdAnaInfo['Gluino600'] = [3.0e-3,4.1e-3,4.8e-3]
StdAnaInfo['Gluino700'] = [3.0e-3,4.0e-3,4.8e-3]
StdAnaInfo['Gluino800'] = [3.0e-3,4.1e-3,4.9e-3]
StdAnaInfo['Gluino900'] = [3.2e-3,4.3e-3,5.2e-3]
StdAnaInfo['Gluino1000'] = [3.5e-3,3.4e-3,5.7e-3]
StdAnaInfo['Gluino1100'] = [4.0e-3,3.9e-3,6.6e-3]
StdAnaInfo['Gluino1200'] = [4.8e-3,6.4e-3,7.6e-3]
# Stop
StdAnaInfo['Stop130'] = [1.5e-2,1.9e-2,3.8e-2]
StdAnaInfo['Stop200'] = [5.8e-3,1.1e-2,1.4e-2]
StdAnaInfo['Stop300'] = [3.4e-3,3.4e-3,6.7e-3]
StdAnaInfo['Stop400'] = [2.8e-3,2.8e-3,4.8e-3]
StdAnaInfo['Stop500'] = [2.4e-3,2.3e-3,4.0e-3]
StdAnaInfo['Stop600'] = [2.4e-3,3.2e-3,3.8e-3]
StdAnaInfo['Stop700'] = [2.3e-3,3.0e-3,3.6e-3]
StdAnaInfo['Stop800'] = [2.1e-3,2.8e-3,3.4e-3]
# GMStau
StdAnaInfo['GMStau100'] = [1.2e-2,1.3e-2,3.1e-2]
StdAnaInfo['GMStau126'] = [6.5e-3,8.1e-3,1.7e-2]
StdAnaInfo['GMStau156'] = [4.1e-3,7.2e-3,1.0e-2]
StdAnaInfo['GMStau200'] = [2.3e-3,4.1e-3,5.7e-3]
StdAnaInfo['GMStau247'] = [1.4e-3,1.7e-3,3.7e-3]
StdAnaInfo['GMStau308'] = [1.2e-3,1.2e-3,2.3e-3]
StdAnaInfo['GMStau370'] = [1.1e-3,1.1e-3,1.8e-3]
StdAnaInfo['GMStau432'] = [1.0e-3,1.0e-3,1.7e-3]
StdAnaInfo['GMStau494'] = [9.3e-4,9.3e-4,1.5e-3]
# gluinoN f=0.1
StdAnaInfo['Gluino300N'] = [3.0e-2,3.0e-2,6.8e-2]
StdAnaInfo['Gluino400N'] = [1.8e-2,1.7e-2,3.0e-2]
StdAnaInfo['Gluino500N'] = [1.6e-2,1.6e-2,2.7e-2]
StdAnaInfo['Gluino600N'] = [1.3e-2,1.3e-2,2.2e-2]
StdAnaInfo['Gluino700N'] = [1.4e-2,1.8e-2,2.2e-2]
StdAnaInfo['Gluino800N'] = [1.4e-2,1.9e-2,2.2e-2]
StdAnaInfo['Gluino900N'] = [1.4e-2,1.9e-2,2.2e-2]
StdAnaInfo['Gluino1000N'] = [1.6e-2,2.2e-2,2.6e-2]
StdAnaInfo['Gluino1100N'] = [1.9e-2,2.6e-2,3.1e-2]
StdAnaInfo['Gluino1200N'] = [2.3e-2,3.0e-2,3.6e-2]
# StopN
StdAnaInfo['Stop130N'] = [1.5e-1,2.1e-1,4.2e-1]
StdAnaInfo['Stop200N'] = [4.0e-2,7.2e-2,1.0e-1]
StdAnaInfo['Stop300N'] = [1.5e-2,1.5e-2,3.4e-2]
StdAnaInfo['Stop400N'] = [1.2e-2,1.1e-2,1.9e-2]
StdAnaInfo['Stop500N'] = [9.1e-3,9.1e-3,1.5e-2]
StdAnaInfo['Stop600N'] = [8.3e-3,8.2e-3,1.4e-2]
StdAnaInfo['Stop700N'] = [8.3e-3,1.1e-2,1.3e-2]
StdAnaInfo['Stop800N'] = [8.1e-3,1.1e-2,1.3e-2]
# gluino f=0.5
StdAnaInfo['Gluino300f5'] = [6.9e-3,6.8e-3,1.6e-2]
StdAnaInfo['Gluino400f5'] = [6.0e-3,6.0e-3,1.0e-2]
StdAnaInfo['Gluino500f5'] = [5.5e-3,5.5e-3,9.3e-3]
StdAnaInfo['Gluino600f5'] = [5.4e-3,7.3e-3,8.6e-3]
StdAnaInfo['Gluino700f5'] = [5.3e-3,7.1e-3,8.5e-3]
StdAnaInfo['Gluino800f5'] = [5.4e-3,7.3e-3,8.7e-3]
StdAnaInfo['Gluino900f5'] = [5.8e-3,7.7e-3,9.2e-3]
StdAnaInfo['Gluino1000f5'] = [5.9e-3,5.9e-3,1.0e-2]
StdAnaInfo['Gluino1100f5'] = [7.0e-3,6.9e-3,1.2e-2]
StdAnaInfo['Gluino1200f5'] = [8.6e-3,1.2e-2,1.4e-2]

class LimitResult:
  def __init__(self,name,mass,thCrossSection,ptCut,iasCut,massCut,sigEff,expBg,expBgStatErr,obsEvts,expLimit,obsLimit,fiveSigmaXSec,obsSignalSignif):
    self.name = name
    self.mass = mass
    self.thCrossSection = thCrossSection
    self.ptCut = ptCut
    self.iasCut = iasCut
    self.massCut = massCut
    self.sigEff = sigEff
    self.expBg = expBg
    self.expBgStatErr = expBgStatErr
    self.obsEvts = obsEvts
    self.expLimit = expLimit
    self.obsLimit = obsLimit
    self.fiveSigmaXSec = fiveSigmaXSec
    self.obsSignalSignif = obsSignalSignif

  def Print(self):
    print 'LimitResults object information:'
    print '\tname=',self.name
    print '\tmass=',self.mass
    print '\tthCrossSection=',self.thCrossSection
    print '\tptCut=',self.ptCut
    print '\tiasCut=',self.iasCut
    print '\tmassCut=',self.massCut
    print '\tsigEff=',self.sigEff
    print '\texpBg=',self.expBg
    print '\texpBgStatErr=',self.expBgStatErr
    print '\tobsEvts=',self.obsEvts
    print '\texpLimit=',self.expLimit
    print '\tobsLimit=',self.obsLimit
    print '\tfiveSigmaXSec=',self.fiveSigmaXSec
    print '\tobsSignalSignif=',self.obsSignalSignif

  def LatexTableLine(self):
    modelSymbol = ''
    if('Gluino' in self.name):
      modelSymbol = '$\\tilde{g}$'
    elif('GMStau' in self.name):
      modelSymbol = '$\\tilde\\tau$'
    elif('Stop' in self.name):
      modelSymbol = '$\\tilde{t}$'
    latexLine = modelSymbol+' '+str(self.mass)
    if 'N' in self.name:
      latexLine+='N'
    latexLine+='&'+self.massCut+'&'+str(round(self.sigEff,2))+'&'+str(round(self.expBg,1))+'$\pm$'+str(round(self.expBgStatErr,1))+'&'
    latexLine+=str(int(self.obsEvts))+'&'
    latexLine+='%.1E'%float(self.thCrossSection)+'&'
    latexLine+='%.1E'%self.expLimit+'&'
    latexLine+='%.1E'%self.obsLimit+'&'
    latexLine+='%.1E'%self.fiveSigmaXSec
    #latexLine+='&'+str(round(self.obsSignalSignif,2))
    latexLine+='\\\\'
    return latexLine

  def StringTableLine(self):
    tableString = string.ljust(self.name,15)+string.ljust(self.ptCut,6)+string.ljust(self.iasCut,6)
    tableString+=string.center(str(self.massCut),8)+string.center(str(round(self.sigEff,4)),10)
    backExpString = str(round(self.expBg,4))+' +/- '+str(round(self.expBgStatErr,4))
    tableString+=string.center(backExpString,22)+string.center(str(self.obsEvts),10)
    tableString+=string.center(str(round(float(self.thCrossSection),6)),10)+string.center(str(round(self.expLimit,6)),10)
    tableString+=string.center(str(round(self.obsLimit,6)),10)+string.center(str(round(self.fiveSigmaXSec,6)),12)
    tableString+=string.center(str(round(self.obsSignalSignif,6)),10)
    return tableString



def GetModelName(fileName):
  massCutPos = fileName.find("_massCut")
  fileNameReduced = fileName[0:massCutPos]
  lastUndscPos = fileNameReduced.rfind("_")
  signalName = fileNameReduced[lastUndscPos+1:len(fileNameReduced)]
  return signalName

def GetMassCut(fileName):
  ptCutPos = fileName.find("_ptCut")
  substr = fileName[fileName.find("_massCut")+len("_massCut"):ptCutPos]
  return substr

def GetPtCut(fileName):
  if(runCERN):
    dotOutPos = fileName.rfind(".log")
  else:
    dotOutPos = fileName.rfind(".out")
  substr = fileName[fileName.find("_ptCut")+len("_ptCut"):dotOutPos]
  return substr

def GetIasCut(filePath):
  for line in open(filePath):
    if "Ias" in line and "Running" in line:
      break
  #print "GetIasCut ",line
  return line[line.find("(")+1:len(line)-2]

def GetBackExp(filePath):
  for line in open(filePath):
    if "found" in line and not "WARNING" in line:
       break
  lineSplit = line.split(' ')
  if "found" in line and not "WARNING" in line:
    return float(lineSplit[1])
  else:
    return 0

def GetBackExpError(filePath):
  for line in open(filePath):
    if "found" in line and not "WARNING" in line:
      break
  lineSplit = line.split(' ')
  if "found" in line and not "WARNING" in line:
    return float(lineSplit[3])
  else:
    return 0

def GetSigEff(filePath):
  for line in open(filePath):
    if "event level" in line:
      break
  lineSplit = line.split(' ')
  if "event level" in line:
    return float(lineSplit[0])
  else:
    return 0

def GetBkOverIas(filePath):
  for line in open(filePath):
    if "Bk" in line:
      break
  lineSplit = line.split(' ')
  if "Bk" in line:
    return float(lineSplit[0])
  else:
    return 0.0

def GetFile(fileList,modelName):
  for fileName in fileList:
    if fileEnd in fileName and modelName in fileName:
      return fileName
  return ''

def GetRootFile(fileList,modelName,limitFile):
  for fileName in fileList:
    if modelName in fileName:
      if limitFile and 'combined' in fileName:
        return fileName
      elif not limitFile:
        return fileName
  return ''

def GetExpDiscSignif(filePath):
  for line in open(filePath):
    if "significance at 0" in line:
      value = line[line.rfind("significance")+13:line.rfind("sigma")-2]
      if(value[0]=='-'):
        return 0.0
      else:
        return float(value)
  return 0.0

def GetSignalSignificanceFiveSigma(doSignificanceDir,fileList,modelName):
  if(len(fileList) == 0):
    return -666
  poiSigDict = {}
  for fileName in fileList:
    if modelName in fileName and fileEnd in fileName:
      if 'asymptotic' in fileName:
        if 'N' in fileName and not 'N' in modelName:
          continue
        if 'f5' in fileName and not 'f5' in modelName:
          continue
        medianSignificance = GetExpDiscSignif(doSignificanceDir+fileName)
        #print 'modelName=',modelName
        #print 'fileName=',fileName
        #print 'poi=',fileName[fileName.find(modelName)+len(modelName)+1:fileName.rfind("asymptotic")-1]
        poi = float(fileName[fileName.find(modelName)+len(modelName)+1:fileName.rfind("asymptotic")-1])
        if(medianSignificance > 0):
          poiSigDict[poi] = medianSignificance
  medianSignifs = vector("double")()
  poiValues = vector("double")()
  for poi,sig in sorted(poiSigDict.items()):
    poiValues.push_back(poi)
    medianSignifs.push_back(sig)
  if(len(medianSignifs) < 3):
    print 'ERROR: Not enough (nonzero) points to do interpolation!'
    return -666
  interp = Math.Interpolator(medianSignifs,poiValues)
  return interp.Eval(5)

def GetObsDiscSignif(filePath):
  for line in open(filePath):
    if "Significance =" in line:
      signif = line[line.rfind("Significance =")+15:line.rfind("sigma")-1]
      if "-" in signif:
        signif = signif.lstrip("-")
      return float(signif)
  return 0.0

def GetSignalSignificance(doSignificanceDir,fileList,modelName):
  significance = -999
  for fileName in fileList:
    if modelName in fileName and fileEnd in fileName:
      if 'asymptotic' in fileName:
        if 'N' in fileName and not 'N' in modelName:
          continue
        if 'f5' in fileName and not 'f5' in modelName:
          continue
        significance = GetObsDiscSignif(doSignificanceDir+fileName)
        break
  return significance

def FindIntersection(obsGraph, thGraph, min, max, step, thUncertainty, debug):
  intersection = -1
  thShift = 1.0-thUncertainty
  previousX = min
  previousV = obsGraph.Eval(previousX, 0, "") - (thShift * thGraph.Eval(previousX, 0, ""))
  if previousV > 0 :
     return -1
  for x in xrange(min,max,step):
      v = obsGraph.Eval(x, 0, "") - (thShift * thGraph.Eval(x, 0, "") )
      #if debug:
      #  print "%7.2f --> Obs=%6.2E  Th=%6.2E",x,obs->Eval(x, 0, ""),ThShift * th->Eval(x, 0, "")
      #  if(V>=0)printf("   X\n");
      #  else printf("\n");
      if v < 0:
        previousX = x
        previousV = v
      else:
        intersection = previousX
  return intersection

def DrawPreliminary(Lumi,  X=0.40,  Y=0.995,  W=0.82, H=0.945):
  text = TPaveText(X,Y,W,H, "NDC");
  text.SetFillColor(0);
  text.SetTextAlign(22);
  if(Lumi<0):
    text.AddText("CMS Preliminary   #sqrt{s} 7 TeV");
  if(Lumi>0):
    text.AddText("CMS Preliminary   #sqrt{s} = 7 TeV   "+str(Lumi)+" pb^{-1}");
  text.Draw("same");

def GetErrorBand(name, Mass, Low, High, MinLow=PlotMinScale,MaxHigh=PlotMaxScale):
   N = len(Mass)
   cutg = TCutG(name,2*N)
   cutg.SetFillColor(kGreen-7)
   for i in xrange(0,N):
     Min = max(Low[i],MinLow)
     cutg.SetPoint(i,Mass[i], Min)
   for i in xrange(0,N):
     Max = min(High[N-1-i],MaxHigh)
     cutg.SetPoint(N+i,Mass[N-1-i],Max)
   return cutg

def ReadXSection(InputFile):
  Mass = []
  XSec = []
  Low = []
  High = []
  for line in open(InputFile):
     line = line.strip()
     split = line.split(' ')
     splitString = []
     for i in split:
       if len(i) > 0:
         splitString.append(i)
     Mass.append(float(splitString[0]))
     XSec.append(float(splitString[1]))
     High.append(float(splitString[2]))
     Low.append(float(splitString[3]))
  return {'mass':Mass,'xsec':XSec,'high':High,'low':Low}

def MakePValuePlot(htip):
  gobs = htip.MakePlot()
  gobs.Draw("APL")   
  gplot = gobs
  gplot.GetHistogram().SetTitle( htip.GetTitle() )
  # exp
  gexp = htip.MakeExpectedPlot()
  gexp.Draw()
  # line for CL
  alpha = 0.05
  x1 = gplot.GetXaxis().GetXmin()
  x2 = gplot.GetXaxis().GetXmax()
  line = TLine(x1, alpha, x2,alpha)
  line.SetLineColor(kRed)
  line.Draw()
  # put axis labels 
  gplot.GetXaxis().SetTitle("#sigma (pb)")
  gplot.GetXaxis().SetNdivisions(508)
  gplot.GetYaxis().SetTitle("p value")
  # draw again observed values otherwise will be covered by the bands
  gobs.Draw("PL");
  # redraw the axis 
  gPad.RedrawAxis()
  #
  y0 = 0.6
  verticalSize = 0.3
  y1 = y0 + verticalSize
  l = TLegend(0.6,y0,0.8,y1)
  l.AddEntry(gobs,"","PEL")
  # loop in reverse order (opposite to drawing one)
  #ngraphs = gexp.GetListOfGraphs().GetSize()
  #for i in range(ngraphs-1,-1,-1):
  #  obj =  gexp.GetListOfGraphs().At(i)
  #  lopt = "F"
  #  if(i == ngraphs-1):
  #    lopt = "L"
  #  if(obj != 0):
  #    l.AddEntry(obj,"",lopt);
  l.Draw()
  gPad.Update()


#XXX Run
makeScaledPredictionsDir = BaseDir+'/logs/makeScaledPredictions/'
makeScaledPredictionsLogFileList = os.listdir(makeScaledPredictionsDir)
makeScaledPredictionsRootDir = BaseDir+'/outputs/makeScaledPredictions/'
makeScaledPredictionsRootFileList = os.listdir(makeScaledPredictionsRootDir)
doLimitsDir = BaseDir+'/logs/doLimits/'
try:
  doLimitsLogFileList = os.listdir(doLimitsDir)
except OSError:
  doLimitsLogFileList = []
doLimitsRootDir = BaseDir+'/outputs/doLimits/'
try:
  doLimitsRootFileList = os.listdir(doLimitsRootDir)
except OSError:
  doLimitsRootFileList = []
doSignificanceDir = BaseDir+'/logs/doSignificance/'
try:
  doSignificanceLogFileList = os.listdir(doSignificanceDir)
except OSError:
  doSignificanceLogFileList = []

LimitResultsGluinos = []
LimitResultsStops = []
LimitResultsStaus = []
LimitResultsGluinoNs = []
LimitResultsStopNs = []

# MKDIR
plotsDir = BaseDir+'/results/'
if os.path.isdir(plotsDir) == False:
  os.mkdir(plotsDir)
myCanvas = TCanvas("myCanvas", "myCanvas",1000,800)
##########
#
##########
for model in modelList:
    print 'Do model',model.name
    predLogFile = GetFile(makeScaledPredictionsLogFileList,model.name)
    limitLogFile = GetFile(doLimitsLogFileList,model.name)
    predRootFile = GetRootFile(makeScaledPredictionsRootFileList,model.name,False)
    limitRootFile = GetRootFile(doLimitsRootFileList,model.name,True)
    if predLogFile=='' or limitLogFile=='':
      print 'Warning: Found predLogFile ',predLogFile,' and limitLogFile ',limitLogFile,' one of which is empty'
    if predRootFile=='' or limitRootFile=='':
      print 'Warning: Found predRootFile ',predRootFile,' and limitRootFile ',limitRootFile,' one of which is empty'

    backExp = GetBackExp(makeScaledPredictionsDir+predLogFile)
    backExpError = GetBackExpError(makeScaledPredictionsDir+predLogFile)
    iasCut  = GetIasCut(makeScaledPredictionsDir+predLogFile)
    sigEff = GetSigEff(makeScaledPredictionsDir+predLogFile)
    bk = GetBkOverIas(makeScaledPredictionsDir+predLogFile)
    massCut = GetMassCut(predLogFile)
    ptCut = GetPtCut(predLogFile)
    signalName = GetModelName(predLogFile)
    if(len(limitRootFile) > 0):
      limitTFile = TFile(doLimitsRootDir+limitRootFile)
      htr = limitTFile.Get("result_SigXsec")
      expLimit = htr.GetExpectedUpperLimit(0)
      obsLimit = htr.UpperLimit()
      #pValuePlotTitle = "Cross section limit (pb) for "+model.name
      pValuePlotTitle = ""
      #pValuePlotName = "crossSectionLimit_pValue_"+model.name
      pValueFileName = "crossSectionLimit_pValue_"+model.name
      pValuePlotName = ""
      myCanvas.cd()
      pValuePlot = RooStats.HypoTestInverterPlot("htrResultPlot",pValuePlotName,htr)
      pValuePlot.Draw()
      #MakePValuePlot(pValuePlot)
      Lumi=4976
      X=0.40
      Y=0.995
      W=0.82
      H=0.945
      text = TPaveText(X,Y,W,H, "NDC")
      text.SetFillColor(0)
      text.SetTextAlign(22)
      text.AddText("CMS Preliminary   #sqrt{s} = 7 TeV   "+str(Lumi)+" pb^{-1}")
      text.SetBorderSize(0)
      text.Draw("same")
      myCanvas.Update()
      myCanvas.Print(plotsDir+pValueFileName+".C")
      myCanvas.Print(plotsDir+pValueFileName+".pdf")
      myCanvas.Print(plotsDir+pValueFileName+".png")
      ## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
      #if __name__ == '__main__':
      #   rep = ''
      #   while not rep in [ 'q', 'Q' ]:
      #      rep = raw_input( 'enter "q" to quit: ' )
      #      if 1 < len(rep):
      #         rep = rep[0]

      limitTFile.Close()
    else:
      expLimit = -99
      obsLimit = -99
    predTFile = TFile(makeScaledPredictionsRootDir+predRootFile)
    dataHist = predTFile.Get("dataAllNoMAllEtaUnrolledHist")
    obsEvts = dataHist.Integral()
    predTFile.Close()
    expSignalSignifRaw = GetSignalSignificanceFiveSigma(doSignificanceDir,doSignificanceLogFileList,model.name)
    # Multiply by 1.5, since we are underestimating 5-sigma cross section by using asymptotic method vs. toys
    # see spreadsheet compareToysAsymptoticDiscovery_GMStau100_Gluino1200.xls
    # 1.5 is the spread between the interpolations for the gluino1200, where the asymptotic approx should be the worst
    expSignalSignif = 1.5*expSignalSignifRaw
    obsSignalSignif = GetSignalSignificance(doSignificanceDir,doSignificanceLogFileList,model.name)

    if 'Gluino' in model.name:
      thisLR = LimitResult(model.name,model.mass,model.crossSection,ptCut,iasCut,massCut,sigEff,backExp,backExpError,obsEvts,expLimit,obsLimit,expSignalSignif,obsSignalSignif)
      if 'N' in model.name:
        LimitResultsGluinoNs.append(thisLR)
      else:
        LimitResultsGluinos.append(thisLR)
      #print thisLR.StringTableLine()
    if 'Stop' in model.name:
      thisLR = LimitResult(model.name,model.mass,model.crossSection,ptCut,iasCut,massCut,sigEff,backExp,backExpError,obsEvts,expLimit,obsLimit,expSignalSignif,obsSignalSignif)
      if 'N' in model.name:
        LimitResultsStopNs.append(thisLR)
      else:
        LimitResultsStops.append(thisLR)
      #print thisLR.StringTableLine()
    if 'Stau' in model.name and not 'PPStau' in model.name:
      thisLR = LimitResult(model.name,model.mass,model.crossSection,ptCut,iasCut,massCut,sigEff,backExp,backExpError,obsEvts,expLimit,obsLimit,expSignalSignif,obsSignalSignif)
      LimitResultsStaus.append(thisLR)
      #print thisLR.StringTableLine()

# TABLE OUTPUT
print
print
titleString = string.ljust('Model',15)+string.ljust('Pt',6)+string.ljust('Ias',6)+string.center('Mass',8)
titleString+=string.center('SigEff',10)+string.center('BackExpOverIas',22)+string.center('Obs.',10)
titleString+=string.center('ThCrossSec',10)+string.center('ExpLim',10)+string.center('ObsLim',10)
titleString+=string.center('Exp5SigDisc',12)
titleString+=string.center('ObsSignf.',10)
print titleString
# write text file
txt_file=open(plotsDir+'results.txt','w')
txt_file.write(titleString+'\n')
# print table lines
for limitRes in LimitResultsGluinos:
  print limitRes.StringTableLine()
  txt_file.write(limitRes.StringTableLine()+'\n')
for limitRes in LimitResultsGluinoNs:
  print limitRes.StringTableLine()
  txt_file.write(limitRes.StringTableLine()+'\n')
for limitRes in LimitResultsStops:
  print limitRes.StringTableLine()
  txt_file.write(limitRes.StringTableLine()+'\n')
for limitRes in LimitResultsStopNs:
  print limitRes.StringTableLine()
  txt_file.write(limitRes.StringTableLine()+'\n')
for limitRes in LimitResultsStaus:
  print limitRes.StringTableLine()
  txt_file.write(limitRes.StringTableLine()+'\n')
txt_file.close()

# write tbl fragment file
titleTableString = '\\begin{tabular}{|l|c|c|c|}\n'
titleTableString+= '\\hline\n'
    #latexLine = self.name+'&'+self.massCut+'&'+str(round(self.sigEff,4))+'&'+str(round(self.expBg,4))+'&'+str(round(self.expBgStatErr,4))+'&'
    #latexLine+=self.obsEvts+'&'+str(round(self.ThCrossSec,6))+'&'+str(round(self.expLimit,6))+'&'
    #latexLine+=str(round(self.obsLimit,6))+'&'+str(round(self.fiveSigmaXSec,6))+'&'+str(round(self.obsSignalSignif,6))
titleTableString+= '\\textbf{Model}&$\\mathbf{M_{reco}}$&\\textbf{Sig. Eff.}&\\textbf{Back. Exp.}&'
titleTableString+= '\\textbf{Obs.}&\\textbf{Th. $\sigma$}&\\textbf{Exp. $\sigma$}&\\textbf{Obs. $\sigma$}&'
titleTableString+= '\\textbf{Disc. $\sigma$}\\'
#\textbf{Model and mass}&\\textbf{Mass Cut}&\\textbf{Signal Eff.}&\\textbf{Background Exp.}&'
#\textbf{Obs.}&\\textbf{Th. Cross Sec.}&\\textbf{Exp. Lim.}&\\textbf{Obs. Lim.}&\\textbf{5$\sigma$ Disc. Cross Section}\n'#&\\textbf{Obs. Sig. Eff.}\n'
titleTableString+='\\hline\n'
txt_file=open(plotsDir+'results.tbl','w')
txt_file.write(titleTableString+'\n')
print titleTableString
for limitRes in LimitResultsGluinos+LimitResultsGluinoNs+LimitResultsStops+LimitResultsStopNs+LimitResultsStaus:
  txt_file.write(limitRes.LatexTableLine()+'\n')
  print limitRes.LatexTableLine()
txt_file.write('\\hline\n')
print '\\hline'
txt_file.write('\\end{tabular}\n')
print '\\end{tabular}'
txt_file.close()

# Ratios
print
print 'Tables of ratios with standard analysis'
print
for limitRes in LimitResultsGluinos+LimitResultsGluinoNs+LimitResultsStops+LimitResultsStopNs+LimitResultsStaus:
  expLimStd = StdAnaInfo[limitRes.name][0]
  obsLimStd = StdAnaInfo[limitRes.name][1]
  fiveSigmaStd = StdAnaInfo[limitRes.name][2]
  ratioLimExp = round(expLimStd/limitRes.expLimit,2)
  ratioLimObs = round(obsLimStd/limitRes.obsLimit,2)
  ratioFiveSig = round(fiveSigmaStd/limitRes.fiveSigmaXSec,2)
  modelSymbol = ''
  if('Gluino' in limitRes.name):
    modelSymbol = '$\\tilde{g}$'
  elif('GMStau' in limitRes.name):
    modelSymbol = '$\\tilde\\tau$'
  elif('Stop' in limitRes.name):
    modelSymbol = '$\\tilde{t}$'
  latexLine = modelSymbol+' '+str(limitRes.mass)+'&'
  latexLine+='%.1E'%expLimStd+'&'
  latexLine+='%.1E'%limitRes.expLimit+'&'
  #latexLine+='%.1E'%obsLimStd+'&'
  #latexLine+='%.1E'%limitRes.obsLimit+'&'
  latexLine+='%.1E'%fiveSigmaStd+'&'
  latexLine+='%.1E'%limitRes.fiveSigmaXSec+'&'
  latexLine+=str(ratioLimExp)+'&'
  #latexline+=str(ratioLimObs)+'&'
  latexLine+=str(ratioFiveSig)
  latexLine+='\\\\'
  print latexLine
print
print
titleString = string.ljust('Model',15)+string.center('ExpLimStd',12)+string.center('ObsLimStd',12)+string.center('5XsecStd',12)
titleString+=string.center('ExpLimShape',12)+string.center('ObsLimShape',12)+string.center('5XsecShape',12)
titleString+=string.center('ExpLimStd/sh',14)+string.center('ObsLimStd/sh',14)+string.center('5XsecStd/sh',14)
print titleString
for limitRes in LimitResultsGluinos+LimitResultsGluinoNs+LimitResultsStops+LimitResultsStopNs+LimitResultsStaus:
  expLimStd = StdAnaInfo[limitRes.name][0]
  obsLimStd = StdAnaInfo[limitRes.name][1]
  fiveSigmaStd = StdAnaInfo[limitRes.name][2]
  ratioLimExp = round(expLimStd/limitRes.expLimit,2)
  ratioLimObs = round(obsLimStd/limitRes.obsLimit,2)
  ratioFiveSig = round(fiveSigmaStd/limitRes.fiveSigmaXSec,2)
  printString = string.ljust(limitRes.name,15)+string.center(str(expLimStd),12)+string.center(str(obsLimStd),12)+string.center(str(fiveSigmaStd),12)
  printString+=string.center('%.1E'%limitRes.expLimit,12)+string.center('%.1E'%limitRes.obsLimit,12)+string.center('%.1E'%limitRes.fiveSigmaXSec,12)
  printString+=string.center(str(ratioLimExp),14)+string.center(str(ratioLimObs),14)+string.center(str(ratioFiveSig),14)
  print printString

# TESTING
#sys.exit()

# gluinos
massesGluinos = []
expCrossSectionsGluinos = []
obsCrossSectionsGluinos = []
stdShapeExpRatioGluinos = []
stdShapeDiscRatioGluinos = []
for lr in LimitResultsGluinos:
  massesGluinos.append(lr.mass)
  expCrossSectionsGluinos.append(lr.expLimit)
  obsCrossSectionsGluinos.append(lr.obsLimit)
  expLimStd = StdAnaInfo[lr.name][0]
  fiveSigmaStd = StdAnaInfo[lr.name][2]
  ratioLimExp = round(expLimStd/lr.expLimit,2)
  ratioFiveSig = round(fiveSigmaStd/lr.fiveSigmaXSec,2)
  stdShapeExpRatioGluinos.append(ratioLimExp)
  stdShapeDiscRatioGluinos.append(ratioFiveSig)
# stops
massesStops = []
expCrossSectionsStops = []
obsCrossSectionsStops = []
stdShapeExpRatioStops = []
stdShapeDiscRatioStops = []
for lr in LimitResultsStops:
  massesStops.append(lr.mass)
  expCrossSectionsStops.append(lr.expLimit)
  obsCrossSectionsStops.append(lr.obsLimit)
  expLimStd = StdAnaInfo[lr.name][0]
  fiveSigmaStd = StdAnaInfo[lr.name][2]
  ratioLimExp = round(expLimStd/lr.expLimit,2)
  ratioFiveSig = round(fiveSigmaStd/lr.fiveSigmaXSec,2)
  stdShapeExpRatioStops.append(ratioLimExp)
  stdShapeDiscRatioStops.append(ratioFiveSig)
# staus
massesStaus = []
expCrossSectionsStaus = []
obsCrossSectionsStaus = []
stdShapeExpRatioStaus = []
stdShapeDiscRatioStaus = []
for lr in LimitResultsStaus:
  massesStaus.append(lr.mass)
  expCrossSectionsStaus.append(lr.expLimit)
  obsCrossSectionsStaus.append(lr.obsLimit)
  expLimStd = StdAnaInfo[lr.name][0]
  fiveSigmaStd = StdAnaInfo[lr.name][2]
  ratioLimExp = round(expLimStd/lr.expLimit,2)
  ratioFiveSig = round(fiveSigmaStd/lr.fiveSigmaXSec,2)
  stdShapeExpRatioStaus.append(ratioLimExp)
  stdShapeDiscRatioStaus.append(ratioFiveSig)

# now for charge suppression models
# gluinoNs
massesGluinoNs = []
expCrossSectionsGluinoNs = []
obsCrossSectionsGluinoNs = []
stdShapeExpRatioGluinoNs = []
stdShapeDiscRatioGluinoNs = []
for lr in LimitResultsGluinoNs:
  massesGluinoNs.append(lr.mass)
  expCrossSectionsGluinoNs.append(lr.expLimit)
  obsCrossSectionsGluinoNs.append(lr.obsLimit)
  expLimStd = StdAnaInfo[lr.name][0]
  fiveSigmaStd = StdAnaInfo[lr.name][2]
  ratioLimExp = round(expLimStd/lr.expLimit,2)
  ratioFiveSig = round(fiveSigmaStd/lr.fiveSigmaXSec,2)
  stdShapeExpRatioGluinoNs.append(ratioLimExp)
  stdShapeDiscRatioGluinoNs.append(ratioFiveSig)
# stopNs
massesStopNs = []
expCrossSectionsStopNs = []
obsCrossSectionsStopNs = []
stdShapeExpRatioStopNs = []
stdShapeDiscRatioStopNs = []
for lr in LimitResultsStopNs:
  massesStopNs.append(lr.mass)
  expCrossSectionsStopNs.append(lr.expLimit)
  obsCrossSectionsStopNs.append(lr.obsLimit)
  expLimStd = StdAnaInfo[lr.name][0]
  fiveSigmaStd = StdAnaInfo[lr.name][2]
  ratioLimExp = round(expLimStd/lr.expLimit,2)
  ratioFiveSig = round(fiveSigmaStd/lr.fiveSigmaXSec,2)
  stdShapeExpRatioStopNs.append(ratioLimExp)
  stdShapeDiscRatioStopNs.append(ratioFiveSig)

    
# gluinos
gluinoThInfo = ReadXSection(GluinoXSecFile)
gluinoXSecErr = GetErrorBand("gluinoErr",gluinoThInfo['mass'],gluinoThInfo['low'],gluinoThInfo['high'])
theoryGraphGluinos = TGraph(len(gluinoThInfo['mass']),array.array("f",gluinoThInfo['mass']),array.array("f",gluinoThInfo['xsec']))
theoryGraphGluinos.SetTitle('')
theoryGraphGluinos.GetYaxis().SetTitle("CrossSection ( pb )")
theoryGraphGluinos.GetYaxis().SetTitleOffset(1.70)
theoryGraphGluinos.SetLineColor(4)
theoryGraphGluinos.SetMarkerColor(4)
theoryGraphGluinos.SetLineWidth(1)
theoryGraphGluinos.SetLineStyle(3)
theoryGraphGluinos.SetMarkerStyle(1)
expLimGraphGluinos = TGraph(len(massesGluinos),array.array("f",massesGluinos),array.array("f",expCrossSectionsGluinos))
expLimGraphGluinos.SetTitle('')
expLimGraphGluinos.GetYaxis().SetTitle("CrossSection ( pb )")
expLimGraphGluinos.GetYaxis().SetTitleOffset(1.70)
expLimGraphGluinos.SetLineColor(4)
expLimGraphGluinos.SetMarkerColor(4)
expLimGraphGluinos.SetLineWidth(2)
expLimGraphGluinos.SetLineStyle(1)
expLimGraphGluinos.SetMarkerStyle(22)
obsLimGraphGluinos = TGraph(len(massesGluinos),array.array("f",massesGluinos),array.array("f",obsCrossSectionsGluinos))
obsLimGraphGluinos.SetTitle('')
obsLimGraphGluinos.GetYaxis().SetTitle("CrossSection ( pb )")
obsLimGraphGluinos.GetYaxis().SetTitleOffset(1.70)
obsLimGraphGluinos.SetLineColor(4)
obsLimGraphGluinos.SetMarkerColor(4)
obsLimGraphGluinos.SetLineWidth(2)
obsLimGraphGluinos.SetLineStyle(1)
obsLimGraphGluinos.SetMarkerStyle(22)
stdShapeExpLimRatioGraphGluinos = TGraph(len(massesGluinos),array.array("f",massesGluinos),array.array("f",stdShapeExpRatioGluinos))
stdShapeExpLimRatioGraphGluinos.SetTitle('')
#stdShapeExpLimRatioGraphGluinos.GetYaxis().SetTitle("")
#stdShapeExpLimRatioGraphGluinos.GetYaxis().SetTitleOffset(1.70)
stdShapeExpLimRatioGraphGluinos.SetLineColor(4)
stdShapeExpLimRatioGraphGluinos.SetMarkerColor(4)
stdShapeExpLimRatioGraphGluinos.SetLineWidth(2)
stdShapeExpLimRatioGraphGluinos.SetLineStyle(1)
stdShapeExpLimRatioGraphGluinos.SetMarkerStyle(22)
stdShapeDiscRatioGraphGluinos = TGraph(len(massesGluinos),array.array("f",massesGluinos),array.array("f",stdShapeDiscRatioGluinos))
stdShapeDiscRatioGraphGluinos.SetTitle('')
#stdShapeDiscRatioGraphGluinos.GetYaxis().SetTitle("")
#stdShapeDiscRatioGraphGluinos.GetYaxis().SetTitleOffset(1.70)
stdShapeDiscRatioGraphGluinos.SetLineColor(4)
stdShapeDiscRatioGraphGluinos.SetMarkerColor(4)
stdShapeDiscRatioGraphGluinos.SetLineWidth(2)
stdShapeDiscRatioGraphGluinos.SetLineStyle(1)
stdShapeDiscRatioGraphGluinos.SetMarkerStyle(22)
# gluinoNs
expLimGraphGluinoNs = TGraph(len(massesGluinoNs),array.array("f",massesGluinoNs),array.array("f",expCrossSectionsGluinoNs))
expLimGraphGluinoNs.SetTitle('')
expLimGraphGluinoNs.GetYaxis().SetTitle("CrossSection ( pb )")
expLimGraphGluinoNs.GetYaxis().SetTitleOffset(1.70)
expLimGraphGluinoNs.SetLineColor(4)
expLimGraphGluinoNs.SetMarkerColor(4)
expLimGraphGluinoNs.SetLineWidth(2)
expLimGraphGluinoNs.SetLineStyle(1)
expLimGraphGluinoNs.SetMarkerStyle(26)
obsLimGraphGluinoNs = TGraph(len(massesGluinoNs),array.array("f",massesGluinoNs),array.array("f",obsCrossSectionsGluinoNs))
obsLimGraphGluinoNs.SetTitle('')
obsLimGraphGluinoNs.GetYaxis().SetTitle("CrossSection ( pb )")
obsLimGraphGluinoNs.GetYaxis().SetTitleOffset(1.70)
obsLimGraphGluinoNs.SetLineColor(4)
obsLimGraphGluinoNs.SetMarkerColor(4)
obsLimGraphGluinoNs.SetLineWidth(2)
obsLimGraphGluinoNs.SetLineStyle(1)
obsLimGraphGluinoNs.SetMarkerStyle(26)
stdShapeExpLimRatioGraphGluinoNs = TGraph(len(massesGluinoNs),array.array("f",massesGluinoNs),array.array("f",stdShapeExpRatioGluinoNs))
stdShapeExpLimRatioGraphGluinoNs.SetTitle('')
#stdShapeExpLimRatioGraphGluinoNs.GetYaxis().SetTitle("")
#stdShapeExpLimRatioGraphGluinoNs.GetYaxis().SetTitleOffset(1.70)
stdShapeExpLimRatioGraphGluinoNs.SetLineColor(4)
stdShapeExpLimRatioGraphGluinoNs.SetMarkerColor(4)
stdShapeExpLimRatioGraphGluinoNs.SetLineWidth(2)
stdShapeExpLimRatioGraphGluinoNs.SetLineStyle(1)
stdShapeExpLimRatioGraphGluinoNs.SetMarkerStyle(26)
stdShapeDiscRatioGraphGluinoNs = TGraph(len(massesGluinoNs),array.array("f",massesGluinoNs),array.array("f",stdShapeDiscRatioGluinoNs))
stdShapeDiscRatioGraphGluinoNs.SetTitle('')
#stdShapeDiscRatioGraphGluinoNs.GetYaxis().SetTitle("")
#stdShapeDiscRatioGraphGluinoNs.GetYaxis().SetTitleOffset(1.70)
stdShapeDiscRatioGraphGluinoNs.SetLineColor(4)
stdShapeDiscRatioGraphGluinoNs.SetMarkerColor(4)
stdShapeDiscRatioGraphGluinoNs.SetLineWidth(2)
stdShapeDiscRatioGraphGluinoNs.SetLineStyle(1)
stdShapeDiscRatioGraphGluinoNs.SetMarkerStyle(26)
# Stops
stopThInfo = ReadXSection(StopXSecFile)
stopXSecErr = GetErrorBand("stopErr",stopThInfo['mass'],stopThInfo['low'],stopThInfo['high'])
theoryGraphStops = TGraph(len(stopThInfo['mass']),array.array("f",stopThInfo['mass']),array.array("f",stopThInfo['xsec']))
theoryGraphStops.SetTitle('')
theoryGraphStops.GetYaxis().SetTitle("CrossSection ( pb )")
theoryGraphStops.GetYaxis().SetTitleOffset(1.70)
theoryGraphStops.SetLineColor(2)
theoryGraphStops.SetMarkerColor(2)
theoryGraphStops.SetLineWidth(1)
theoryGraphStops.SetLineStyle(2)
theoryGraphStops.SetMarkerStyle(1);
expLimGraphStops = TGraph(len(massesStops),array.array("f",massesStops),array.array("f",expCrossSectionsStops))
expLimGraphStops.SetTitle('')
expLimGraphStops.GetYaxis().SetTitle("CrossSection ( pb )")
expLimGraphStops.GetYaxis().SetTitleOffset(1.70)
expLimGraphStops.SetLineColor(2)
expLimGraphStops.SetMarkerColor(2)
expLimGraphStops.SetLineWidth(2)
expLimGraphStops.SetLineStyle(1)
expLimGraphStops.SetMarkerStyle(21)
obsLimGraphStops = TGraph(len(massesStops),array.array("f",massesStops),array.array("f",obsCrossSectionsStops))
obsLimGraphStops.SetTitle('')
obsLimGraphStops.GetYaxis().SetTitle("CrossSection ( pb )")
obsLimGraphStops.GetYaxis().SetTitleOffset(1.70)
obsLimGraphStops.SetLineColor(2)
obsLimGraphStops.SetMarkerColor(2)
obsLimGraphStops.SetLineWidth(2)
obsLimGraphStops.SetLineStyle(1)
obsLimGraphStops.SetMarkerStyle(21)
stdShapeExpLimRatioGraphStops = TGraph(len(massesStops),array.array("f",massesStops),array.array("f",stdShapeExpRatioStops))
stdShapeExpLimRatioGraphStops.SetTitle('')
#stdShapeExpLimRatioGraphStops.GetYaxis().SetTitle("")
#stdShapeExpLimRatioGraphStops.GetYaxis().SetTitleOffset(1.70)
stdShapeExpLimRatioGraphStops.SetLineColor(2)
stdShapeExpLimRatioGraphStops.SetMarkerColor(2)
stdShapeExpLimRatioGraphStops.SetLineWidth(2)
stdShapeExpLimRatioGraphStops.SetLineStyle(1)
stdShapeExpLimRatioGraphStops.SetMarkerStyle(21)
stdShapeDiscRatioGraphStops = TGraph(len(massesStops),array.array("f",massesStops),array.array("f",stdShapeDiscRatioStops))
stdShapeDiscRatioGraphStops.SetTitle('')
#stdShapeDiscRatioGraphStops.GetYaxis().SetTitle("")
#stdShapeDiscRatioGraphStops.GetYaxis().SetTitleOffset(1.70)
stdShapeDiscRatioGraphStops.SetLineColor(2)
stdShapeDiscRatioGraphStops.SetMarkerColor(2)
stdShapeDiscRatioGraphStops.SetLineWidth(2)
stdShapeDiscRatioGraphStops.SetLineStyle(1)
stdShapeDiscRatioGraphStops.SetMarkerStyle(21)
# stopNs
expLimGraphStopNs = TGraph(len(massesStopNs),array.array("f",massesStopNs),array.array("f",expCrossSectionsStopNs))
expLimGraphStopNs.SetTitle('')
expLimGraphStopNs.GetYaxis().SetTitle("CrossSection ( pb )")
expLimGraphStopNs.GetYaxis().SetTitleOffset(1.70)
expLimGraphStopNs.SetLineColor(2)
expLimGraphStopNs.SetMarkerColor(2)
expLimGraphStopNs.SetLineWidth(2)
expLimGraphStopNs.SetLineStyle(1)
expLimGraphStopNs.SetMarkerStyle(25)
obsLimGraphStopNs = TGraph(len(massesStopNs),array.array("f",massesStopNs),array.array("f",obsCrossSectionsStopNs))
obsLimGraphStopNs.SetTitle('')
obsLimGraphStopNs.GetYaxis().SetTitle("CrossSection ( pb )")
obsLimGraphStopNs.GetYaxis().SetTitleOffset(1.70)
obsLimGraphStopNs.SetLineColor(2)
obsLimGraphStopNs.SetMarkerColor(2)
obsLimGraphStopNs.SetLineWidth(2)
obsLimGraphStopNs.SetLineStyle(1)
obsLimGraphStopNs.SetMarkerStyle(25)
stdShapeExpLimRatioGraphStopNs = TGraph(len(massesStopNs),array.array("f",massesStopNs),array.array("f",stdShapeExpRatioStopNs))
stdShapeExpLimRatioGraphStopNs.SetTitle('')
#stdShapeExpLimRatioGraphStopNs.GetYaxis().SetTitle("")
#stdShapeExpLimRatioGraphStopNs.GetYaxis().SetTitleOffset(1.70)
stdShapeExpLimRatioGraphStopNs.SetLineColor(2)
stdShapeExpLimRatioGraphStopNs.SetMarkerColor(2)
stdShapeExpLimRatioGraphStopNs.SetLineWidth(2)
stdShapeExpLimRatioGraphStopNs.SetLineStyle(1)
stdShapeExpLimRatioGraphStopNs.SetMarkerStyle(25)
stdShapeDiscRatioGraphStopNs = TGraph(len(massesStopNs),array.array("f",massesStopNs),array.array("f",stdShapeDiscRatioStopNs))
stdShapeDiscRatioGraphStopNs.SetTitle('')
#stdShapeDiscRatioGraphStopNs.GetYaxis().SetTitle("")
#stdShapeDiscRatioGraphStopNs.GetYaxis().SetTitleOffset(1.70)
stdShapeDiscRatioGraphStopNs.SetLineColor(2)
stdShapeDiscRatioGraphStopNs.SetMarkerColor(2)
stdShapeDiscRatioGraphStopNs.SetLineWidth(2)
stdShapeDiscRatioGraphStopNs.SetLineStyle(1)
stdShapeDiscRatioGraphStopNs.SetMarkerStyle(25)
# staus
stauThInfo = ReadXSection(StauXSecFile)
stauXSecErr = GetErrorBand("stauErr",stauThInfo['mass'],stauThInfo['low'],stauThInfo['high'])
theoryGraphStaus = TGraph(len(stauThInfo['mass']),array.array("f",stauThInfo['mass']),array.array("f",stauThInfo['xsec']))
theoryGraphStaus.SetTitle('')
theoryGraphStaus.GetYaxis().SetTitle("CrossSection ( pb )")
theoryGraphStaus.GetYaxis().SetTitleOffset(1.70)
theoryGraphStaus.SetLineColor(1)
theoryGraphStaus.SetMarkerColor(1)
theoryGraphStaus.SetLineWidth(1)
theoryGraphStaus.SetLineStyle(1)
theoryGraphStaus.SetMarkerStyle(1);
expLimGraphStaus = TGraph(len(massesStaus),array.array("f",massesStaus),array.array("f",expCrossSectionsStaus))
expLimGraphStaus.SetTitle('')
expLimGraphStaus.GetYaxis().SetTitle("CrossSection ( pb )")
expLimGraphStaus.GetYaxis().SetTitleOffset(1.70)
expLimGraphStaus.SetLineColor(1)
expLimGraphStaus.SetMarkerColor(1)
expLimGraphStaus.SetLineWidth(2)
expLimGraphStaus.SetLineStyle(1)
expLimGraphStaus.SetMarkerStyle(20)
obsLimGraphStaus = TGraph(len(massesStaus),array.array("f",massesStaus),array.array("f",obsCrossSectionsStaus))
obsLimGraphStaus.SetTitle('')
obsLimGraphStaus.GetYaxis().SetTitle("CrossSection ( pb )")
obsLimGraphStaus.GetYaxis().SetTitleOffset(1.70)
obsLimGraphStaus.SetLineColor(1)
obsLimGraphStaus.SetMarkerColor(1)
obsLimGraphStaus.SetLineWidth(2)
obsLimGraphStaus.SetLineStyle(1)
obsLimGraphStaus.SetMarkerStyle(20)
stdShapeExpLimRatioGraphStaus = TGraph(len(massesStaus),array.array("f",massesStaus),array.array("f",stdShapeExpRatioStaus))
stdShapeExpLimRatioGraphStaus.SetTitle('')
#stdShapeExpLimRatioGraphStaus.GetYaxis().SetTitle("")
#stdShapeExpLimRatioGraphStaus.GetYaxis().SetTitleOffset(1.70)
stdShapeExpLimRatioGraphStaus.SetLineColor(1)
stdShapeExpLimRatioGraphStaus.SetMarkerColor(1)
stdShapeExpLimRatioGraphStaus.SetLineWidth(2)
stdShapeExpLimRatioGraphStaus.SetLineStyle(1)
stdShapeExpLimRatioGraphStaus.SetMarkerStyle(20)
stdShapeDiscRatioGraphStaus = TGraph(len(massesStaus),array.array("f",massesStaus),array.array("f",stdShapeDiscRatioStaus))
stdShapeDiscRatioGraphStaus.SetTitle('')
#stdShapeDiscRatioGraphStaus.GetYaxis().SetTitle("")
#stdShapeDiscRatioGraphStaus.GetYaxis().SetTitleOffset(1.70)
stdShapeDiscRatioGraphStaus.SetLineColor(1)
stdShapeDiscRatioGraphStaus.SetMarkerColor(1)
stdShapeDiscRatioGraphStaus.SetLineWidth(2)
stdShapeDiscRatioGraphStaus.SetLineStyle(1)
stdShapeDiscRatioGraphStaus.SetMarkerStyle(20)

c1 = TCanvas("c1", "c1",600,600)
c1.SetLogy()
gluinoXSecErr.Draw("f")
stopXSecErr.Draw("f")
stauXSecErr.Draw("f")
#PPStauXSecErr  ->Draw("f")
#MGTk.Draw("same")
MGTk = TMultiGraph()
MGTk.Add(theoryGraphGluinos      ,"L")
MGTk.Add(theoryGraphStops      ,"L")
MGTk.Add(theoryGraphStaus      ,"L")
#MGTk.Add(expLimGraphGluinos      ,"LP")
#MGTk.Add(expLimGraphStops      ,"LP")
#MGTk.Add(expLimGraphStaus      ,"LP")
MGTk.Add(obsLimGraphGluinos      ,"LP")
MGTk.Add(obsLimGraphGluinoNs      ,"LP")
MGTk.Add(obsLimGraphStops      ,"LP")
MGTk.Add(obsLimGraphStopNs      ,"LP")
MGTk.Add(obsLimGraphStaus      ,"LP")
MGTk.Draw("A")
MGTk.SetTitle("")
MGTk.GetXaxis().SetTitle("Mass (GeV/c^{2})")
MGTk.GetYaxis().SetTitle("#sigma (pb)")
MGTk.GetYaxis().SetTitleOffset(1.0)
MGTk.GetYaxis().SetRangeUser(PlotMinScale,PlotMaxScale)
# lumi text
Lumi=4976
X=0.40
Y=0.995
W=0.82
H=0.945
text = TPaveText(X,Y,W,H, "NDC")
text.SetFillColor(0)
text.SetTextAlign(22)
text.AddText("CMS Preliminary   #sqrt{s} = 7 TeV   "+str(Lumi)+" pb^{-1}")
text.SetBorderSize(0)
text.Draw("same")
# legend th
LEGTh = TLegend(0.2,0.7,0.42,0.9)
LEGTh.SetHeader("Theoretical Prediction")
LEGTh.SetFillColor(0)
LEGTh.SetBorderSize(0)
GlThLeg = theoryGraphGluinos.Clone("GluinoThLeg")
GlThLeg.SetFillColor(gluinoXSecErr.GetFillColor())
LEGTh.AddEntry(GlThLeg, "gluino (NLO+NLL)" ,"LF")
StThLeg = theoryGraphStops.Clone("StopThLeg")
StThLeg.SetFillColor(stopXSecErr.GetFillColor())
LEGTh.AddEntry(StThLeg   ,"stop   (NLO+NLL)" ,"LF")
StauThLeg = theoryGraphStaus.Clone("StauThLeg")
StauThLeg.SetFillColor(stauXSecErr.GetFillColor())
LEGTh.AddEntry(StauThLeg   ,"GMSB stau   (NLO)" ,"LF")
LEGTh.Draw()
# leg limits
#LEGTk = TLegend(0.7,0.7,0.9,0.9)
LEGTk = TLegend(0.6,0.6,0.9,0.9)
#LEGTk->SetHeader("95% C.L. Limits")
LEGTk.SetHeader("")
LEGTk.SetFillColor(0)
LEGTk.SetBorderSize(0)
#LEGTk.AddEntry(Tk_Obs_GluinoF5 , "gluino; 50% #tilde{g}g"    ,"LP")
LEGTk.AddEntry(obsLimGraphGluinos , "gluino; 10% #tilde{g}g"    ,"LP")
LEGTk.AddEntry(obsLimGraphGluinoNs , "gluino; 10% #tilde{g}g, ch. supp."    ,"LP")
#LEGTk.AddEntry(Tk_Obs_GluinoNF1, "gluino; 10% #tilde{g}g; ch. suppr.","LP")
LEGTk.AddEntry(obsLimGraphStops     , "stop"            ,"LP")
LEGTk.AddEntry(obsLimGraphStopNs     , "stop; ch. supp."            ,"LP")
#LEGTk.AddEntry(Tk_Obs_StopN    , "stop; ch. suppr.","LP")
#LEGTk.AddEntry(Tk_Obs_PPStau   , "Pair Prod. stau"       ,"LP")
LEGTk.AddEntry(obsLimGraphStaus   , "GMSB stau"       ,"LP")
LEGTk.Draw()

limit = FindIntersection(obsLimGraphGluinos,theoryGraphGluinos,300,1200,1,0.00,False)
print 'Found Gluino (f=0.1) mass limit = ',limit
limit = FindIntersection(obsLimGraphGluinoNs,theoryGraphGluinos,300,1200,1,0.00,False)
print 'Found GluinoN (f=0.1) mass limit = ',limit
limit = FindIntersection(obsLimGraphStops,theoryGraphStops,300,1200,1,0.00,False)
print 'Found Stop mass limit = ',limit
limit = FindIntersection(obsLimGraphStopNs,theoryGraphStops,300,1200,1,0.00,False)
print 'Found StopN mass limit = ',limit
limit = FindIntersection(obsLimGraphStaus,theoryGraphStaus,100,1000,1,0.00,False)
print 'Found GMStau mass limit = ',limit
c1.Update()
c1.Print(plotsDir+'gluinosStopsStaus.png')
c1.Print(plotsDir+'gluinosStopsStaus.pdf')

# Ratio of exp lim graph
c2 = TCanvas("c2", "c2",600,600)
#c2.SetLogy()
MGTk = TMultiGraph()
MGTk.Add(stdShapeExpLimRatioGraphGluinos      ,"LP")
MGTk.Add(stdShapeExpLimRatioGraphGluinoNs      ,"LP")
MGTk.Add(stdShapeExpLimRatioGraphStops      ,"LP")
MGTk.Add(stdShapeExpLimRatioGraphStopNs      ,"LP")
MGTk.Add(stdShapeExpLimRatioGraphStaus      ,"LP")
MGTk.Draw("A")
MGTk.SetTitle("")
MGTk.GetXaxis().SetTitle("Mass (GeV/c^{2})")
#MGTk.GetYaxis().SetTitle("#sigma (pb)")
MGTk.GetYaxis().SetTitleOffset(1.0)
MGTk.GetYaxis().SetRangeUser(1.0,1.65)
# lumi text
Lumi=4976
X=0.40
Y=0.995
W=0.82
H=0.945
text = TPaveText(X,Y,W,H, "NDC")
text.SetFillColor(0)
text.SetTextAlign(22)
text.AddText("CMS Preliminary   #sqrt{s} = 7 TeV   "+str(Lumi)+" pb^{-1}")
text.SetBorderSize(0)
text.Draw("same")
# leg limits
LEGTk = TLegend(0.6,0.6,0.9,0.9)
LEGTk.SetHeader("")
LEGTk.SetFillColor(0)
LEGTk.SetBorderSize(0)
LEGTk.AddEntry(stdShapeExpLimRatioGraphGluinos , "gluino; 10% #tilde{g}g"    ,"LP")
LEGTk.AddEntry(stdShapeExpLimRatioGraphGluinoNs , "gluino; 10% #tilde{g}g, ch. supp."    ,"LP")
LEGTk.AddEntry(stdShapeExpLimRatioGraphStops     , "stop"            ,"LP")
LEGTk.AddEntry(stdShapeExpLimRatioGraphStopNs     , "stop; ch. supp."            ,"LP")
LEGTk.AddEntry(stdShapeExpLimRatioGraphStaus   , "GMSB stau"       ,"LP")
LEGTk.Draw()
c2.Update()
c2.Print(plotsDir+'gluinosStopsStausRatioExpLim.png')
c2.Print(plotsDir+'gluinosStopsStausRatioExpLim.pdf')
c2.Print(plotsDir+'gluinosStopsStausRatioExpLim.C')


# Ratio of exp disc 5-sigma xsec graph
c3 = TCanvas("c3", "c3",600,600)
#c3.SetLogy()
MGTk = TMultiGraph()
MGTk.Add(stdShapeDiscRatioGraphGluinos      ,"LP")
MGTk.Add(stdShapeDiscRatioGraphGluinoNs      ,"LP")
MGTk.Add(stdShapeDiscRatioGraphStops      ,"LP")
MGTk.Add(stdShapeDiscRatioGraphStopNs      ,"LP")
MGTk.Add(stdShapeDiscRatioGraphStaus      ,"LP")
MGTk.Draw("A")
MGTk.SetTitle("")
MGTk.GetXaxis().SetTitle("Mass (GeV/c^{2})")
#MGTk.GetYaxis().SetTitle("#sigma (pb)")
MGTk.GetYaxis().SetTitleOffset(1.0)
MGTk.GetYaxis().SetRangeUser(1.5,3.2)
# lumi text
Lumi=4976
X=0.40
Y=0.995
W=0.82
H=0.945
text = TPaveText(X,Y,W,H, "NDC")
text.SetFillColor(0)
text.SetTextAlign(22)
text.AddText("CMS Preliminary   #sqrt{s} = 7 TeV   "+str(Lumi)+" pb^{-1}")
text.SetBorderSize(0)
text.Draw("same")
# leg limits
LEGTk = TLegend(0.25,0.65,0.55,0.95)
LEGTk.SetHeader("")
LEGTk.SetFillColor(0)
LEGTk.SetBorderSize(0)
LEGTk.AddEntry(stdShapeDiscRatioGraphGluinos , "gluino; 10% #tilde{g}g"    ,"LP")
LEGTk.AddEntry(stdShapeDiscRatioGraphGluinoNs , "gluino; 10% #tilde{g}g, ch. supp."    ,"LP")
LEGTk.AddEntry(stdShapeDiscRatioGraphStops     , "stop"            ,"LP")
LEGTk.AddEntry(stdShapeDiscRatioGraphStopNs     , "stop; ch. supp."            ,"LP")
LEGTk.AddEntry(stdShapeDiscRatioGraphStaus   , "GMSB stau"       ,"LP")
LEGTk.Draw()
c3.Update()
c3.Print(plotsDir+'gluinosStopsStausRatioDiscXSec.png')
c3.Print(plotsDir+'gluinosStopsStausRatioDiscXSec.pdf')
c3.Print(plotsDir+'gluinosStopsStausRatioDiscXSec.C')

