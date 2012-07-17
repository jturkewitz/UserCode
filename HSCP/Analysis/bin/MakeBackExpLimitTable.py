#!/usr/bin/env python

import string
import os
import sys
import glob
import math
from subprocess import call

BaseDirPt50Ias0p1 = 'FARM_dReg_NL_50IasBins_emptyBins1e-25_iasPredScalingFix_cutPt50GeVcutIas0.1_allSlices_Apr30'
BaseDirPt50IasStd = ''
BaseDirPtStdIasStd = ''
runCERN = False

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
    #print line
    return float(lineSplit[0])
  else:
    return 0.0

def GetExpLimit(filePath):
  line = ''
  for thisLine in open(filePath):
    if "median" in thisLine:
      line = thisLine
      break
  limitFind = line[line.find("median")+8:len(line)-1]
  if(len(limitFind) > 0):
    return float(limitFind)
  else:
    return 0.0

def GetObsLimit(filePath):
  line = ''
  for thisLine in open(filePath):
    if "computed" in thisLine:
      line = thisLine
      break
  limitFind = line[line.find("is:")+4:line.find("+")-1]
  if(len(limitFind) > 0):
    return float(limitFind)
  else:
    return 0.0
  
doPt50IasStd = True
doPt50Ias0p1 = True
doPtStdIasStd = True
makeScaledPredictionsDirPtStdIasStd = BaseDirPtStdIasStd+'/logs/makeScaledPredictions/'
try:
  fileListPtStdIasStd = os.listdir(makeScaledPredictionsDirPtStdIasStd)
  fileListPtStdIasStd.sort()
except OSError:
  doPtStdIasStd = False
makeScaledPredictionsDirPt50IasStd = BaseDirPt50IasStd+'/logs/makeScaledPredictions/'
try:
  fileListPt50IasStd = os.listdir(makeScaledPredictionsDirPt50IasStd)
  fileListPt50IasStd.sort()
except OSError:
  doPt50IasStd = False
makeScaledPredictionsDirPt50Ias0p1 = BaseDirPt50Ias0p1+'/logs/makeScaledPredictions/'
try:
  fileListPt50Ias0p1 = os.listdir(makeScaledPredictionsDirPt50Ias0p1)
  fileListPt50Ias0p1.sort()
except OSError:
  doPt50Ias0p1 = False
#
doLimitsDirPtStdIasStd = BaseDirPtStdIasStd+'/logs/doLimits/'
if(doPtStdIasStd):
  fileListLimitsPtStdIasStd = os.listdir(doLimitsDirPtStdIasStd)
  fileListLimitsPtStdIasStd.sort()
doLimitsDirPt50IasStd = BaseDirPt50IasStd+'/logs/doLimits/'
if(doPt50IasStd): 
  fileListLimitsPt50IasStd = os.listdir(doLimitsDirPt50IasStd)
  fileListLimitsPt50IasStd.sort()
doLimitsDirPt50Ias0p1 = BaseDirPt50Ias0p1+'/logs/doLimits/'
if(doPt50Ias0p1):
  try:
    fileListLimitsPt50Ias0p1 = os.listdir(doLimitsDirPt50Ias0p1)
    fileListLimitsPt50Ias0p1.sort()
  except OSError:
    fileListLimitsPt50Ias0p1 = []
    
print
#print 'FarmDirectory =',BaseDir
titleString = string.ljust('Model',15)+string.ljust('Pt',8)+string.ljust('Ias',10)+string.center('Mass',8)+string.center('BackExpOverIas',30)
titleString+=string.center('BkOverIas',10)+string.center('SigEff',10)+string.center('ExpLim',10)+string.center('ObsLim',10)
print
print titleString
if(runCERN):
  fileEnd = '.log'
else:
  fileEnd = '.out'
for index,file in enumerate(fileListPt50Ias0p1):
  if fileEnd in file:
    backExpPt50Ias0p1 = GetBackExp(makeScaledPredictionsDirPt50Ias0p1+file)
    backExpErrorPt50Ias0p1 = GetBackExpError(makeScaledPredictionsDirPt50Ias0p1+file)
    iasCutPt50Ias0p1  = GetIasCut(makeScaledPredictionsDirPt50Ias0p1+file)
    sigEffPt50Ias0p1 = GetSigEff(makeScaledPredictionsDirPt50Ias0p1+file)
    bkPt50Ias0p1 = GetBkOverIas(makeScaledPredictionsDirPt50Ias0p1+file)
    massCutPt50Ias0p1 = GetMassCut(file)
    ptCutPt50Ias0p1 = GetPtCut(file)
    limIndex = index
    signalNamePt50Ias0p1 = GetModelName(file)
    signalNameLimPt50Ias0p1 = GetModelName(fileListLimitsPt50Ias0p1[limIndex])
    if(signalNameLimPt50Ias0p1 != signalNamePt50Ias0p1):
      #print "Pt50Ias0p1: signal names don't match: ",signalNameLimPt50Ias0p1," for limits vs ",signalNamePt50Ias0p1
      while signalNameLimPt50Ias0p1 != signalNamePt50Ias0p1 and index < len(fileListLimitsPt50Ias0p1):
        limIndex+=1
        while not fileEnd in fileListLimitsPt50Ias0p1[limIndex]:
          limIndex+=1
        signalNameLimPt50Ias0p1 = GetModelName(fileListLimitsPt50Ias0p1[limIndex])
    if(signalNameLimPt50Ias0p1 != signalNamePt50Ias0p1):
      print "signal names don't match"
      exit(1)
    expLimitPt50Ias0p1 = GetExpLimit(doLimitsDirPt50Ias0p1+fileListLimitsPt50Ias0p1[limIndex])
    printStringLoose = string.ljust(signalNamePt50Ias0p1,15)+string.ljust(ptCutPt50Ias0p1,8)+string.ljust(iasCutPt50Ias0p1,10)
    printStringLoose+=string.center(massCutPt50Ias0p1,8)
    printStringLoose+=string.center(str(round(backExpPt50Ias0p1,30))+' +/- '+str(round(backExpErrorPt50Ias0p1,4)),30)
    printStringLoose+=string.center(str(bkPt50Ias0p1),10)
    printStringLoose+=string.center(str(round(sigEffPt50Ias0p1,4)),10)
    printStringLoose+=string.center(str(round(expLimitPt50Ias0p1,6)),10)
    if(doPtStdIasStd):
      file = fileListPtStdIasStd[index]
      backExpPtStdIasStd = GetBackExp(makeScaledPredictionsDirPtStdIasStd+file)
      backExpErrorPtStdIasStd = GetBackExpError(makeScaledPredictionsDirPtStdIasStd+file)
      iasCutPtStdIasStd  = GetIasCut(makeScaledPredictionsDirPtStdIasStd+file)
      sigEffPtStdIasStd = GetSigEff(makeScaledPredictionsDirPtStdIasStd+file)
      bkPtStdIasStd = GetBkOverIas(makeScaledPredictionsDirPtStdIasStd+file)
      massCutPtStdIasStd = GetMassCut(file)
      ptCutPtStdIasStd = GetPtCut(file)
      signalNamePtStdIasStd = GetModelName(file)
      signalNameLimPtStdIasStd = GetModelName(fileListLimitsPtStdIasStd[limIndex])
      if(signalNameLimPtStdIasStd != signalNamePtStdIasStd):
        print "PtStdIasStd: signal names don't match: ",signalNameLimPtStdIasStd," for limits vs ",signalNamePtStdIasStd
        exit(1)
      expLimitPtStdIasStd = GetExpLimit(doLimitsDirPtStdIasStd+fileListLimitsPtStdIasStd[limIndex])
      printStringStd = string.ljust(signalNamePtStdIasStd,15)+string.ljust(ptCutPtStdIasStd,8)+string.ljust(iasCutPtStdIasStd,10)
      printStringStd+=string.center(massCutPtStdIasStd,8)
      printStringStd+=string.center(str(round(backExpPtStdIasStd,30))+' +/- '+str(round(backExpErrorPtStdIasStd,4)),30)
      printStringStd+=string.center(str(bkPtStdIasStd),10)
      printStringStd+=string.center(str(round(sigEffPtStdIasStd,4)),10)
      printStringStd+=string.center(str(round(expLimitPtStdIasStd,6)),10)
    if(doPt50IasStd):
      file = fileListPt50IasStd[index]
      backExpPt50IasStd = GetBackExp(makeScaledPredictionsDirPt50IasStd+file)
      backExpErrorPt50IasStd = GetBackExpError(makeScaledPredictionsDirPt50IasStd+file)
      iasCutPt50IasStd  = GetIasCut(makeScaledPredictionsDirPt50IasStd+file)
      sigEffPt50IasStd = GetSigEff(makeScaledPredictionsDirPt50IasStd+file)
      bkPt50IasStd = GetBkOverIas(makeScaledPredictionsDirPt50IasStd+file)
      massCutPt50IasStd = GetMassCut(file)
      ptCutPt50IasStd = GetPtCut(file)
      signalNamePt50IasStd = GetModelName(file)
      signalNameLimPt50IasStd = GetModelName(fileListLimitsPt50IasStd[limIndex])
      if(signalNameLimPt50IasStd != signalNamePt50IasStd):
        print "signal names don't match"
        exit(1)
      expLimitPt50IasStd = GetExpLimit(doLimitsDirPt50IasStd+fileListLimitsPt50IasStd[limIndex])
      printString50 = string.ljust(signalNamePt50IasStd,15)+string.ljust(ptCutPt50IasStd,8)+string.ljust(iasCutPt50IasStd,10)
      printString50+=string.center(massCutPt50IasStd,8)
      printString50+=string.center(str(round(backExpPt50IasStd,30))+' +/- '+str(round(backExpErrorPt50IasStd,4)),30)
      printString50+=string.center(str(bkPt50IasStd),10)
      printString50+=string.center(str(round(sigEffPt50IasStd,4)),10)
      printString50+=string.center(str(round(expLimitPt50IasStd,6)),10)
    if(doPtStdIasStd):
      print printStringStd
    if(doPt50IasStd):
      print printString50
    print printStringLoose


