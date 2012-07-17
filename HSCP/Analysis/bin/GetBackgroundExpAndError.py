#!/usr/bin/env python

import string
import os
import sys
import glob
import math
from subprocess import call

#BaseDir2 = "FARM_moreModels_cutPt50GeV_Feb21/logs/makeScaledPredictions/"
#BaseDir1 = "FARM_moreModels_usePtOver50GeVOnly_Feb22/logs/makeScaledPredictions/"
#BaseDir2 = "FARM_moreModels_cutPt50GeVNewError_Feb28/logs/makeScaledPredictions/"
#BaseDir1 = "FARM_moreModels_usePtOver50GeVOnlyNewError_Feb28/logs/makeScaledPredictions/"
BaseDir2 = "FARM_moreModels_cutPt50GeVNewBGPred_Mar04/logs/makeScaledPredictions/"
BaseDir1 = "FARM_moreModels_cutPtOver50GeVOnlyNewBGPred_Mar04/logs/makeScaledPredictions/"

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
  dotOutPos = fileName.rfind(".out")
  substr = fileName[fileName.find("_ptCut")+len("_ptCut"):dotOutPos]
  return substr

def GetIasCut(filePath):
  for line in open(filePath):
    if "Ias" in line:
      break
  return line[line.find("(")+1:len(line)-2]

def GetBackExp(filePath):
  for line in open(filePath):
    if "found" in line:
      break
  lineSplit = line.split(' ')
  if "found" in line:
    return float(lineSplit[1])
  else:
    return 0

def GetBackExpError(filePath):
  for line in open(filePath):
    if "found" in line:
      break
  lineSplit = line.split(' ')
  if "found" in line:
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
    return float(lineSplit[3])
  else:
    return 0.0



fileList1 = os.listdir(BaseDir1)
print
print '1'
print 'Over 50 GeV Only'
print 'Model         PtCut   IasCut     MassCut    BackExpOverIas   BkOverIas   SignalEff'
print
for file in fileList1:
  if ".out" in file:
    backExp = GetBackExp(BaseDir1+file)
    backExpError = GetBackExpError(BaseDir1+file)
    iasCut  = GetIasCut(BaseDir1+file)
    sigEff = GetSigEff(BaseDir1+file)
    bk = GetBkOverIas(BaseDir1+file)
    massCut = GetMassCut(file)
    ptCut = GetPtCut(file)
    signalName = GetModelName(file)
    printString = string.ljust(signalName,15)+string.ljust(ptCut,8)+string.ljust(iasCut,10)
    printString+=string.center(massCut,8)
    printString+=string.center(str(round(backExp,4))+' +/- '+str(round(backExpError,4)),20)
    printString+=string.center(str(bk),5)
    printString+=string.center(str(round(sigEff,4)),10)
    print printString


fileList2 = os.listdir(BaseDir2)
print
print '2'
print 'Up to 50 GeV'
print 'Model         PtCut   IasCut     MassCut    BackExpOverIas   BkOverIas   SignalEff'
print
for file in fileList2:
  if ".out" in file:
    backExp = GetBackExp(BaseDir2+file)
    backExpError = GetBackExpError(BaseDir2+file)
    iasCut  = GetIasCut(BaseDir2+file)
    sigEff = GetSigEff(BaseDir2+file)
    bk = GetBkOverIas(BaseDir2+file)
    massCut = GetMassCut(file)
    ptCut = GetPtCut(file)
    signalName = GetModelName(file)
    printString = string.ljust(signalName,15)+string.ljust(ptCut,8)+string.ljust(iasCut,10)
    printString+=string.center(massCut,8)
    printString+=string.center(str(round(backExp,4))+' +/- '+str(round(backExpError,4)),20)
    printString+=string.center(str(bk),5)
    printString+=string.center(str(round(sigEff,4)),10)
    print printString

# Now do ratios/error
print
print 'Model           BackExp1/BackExp2    SigEff1/SigEff2'
print
for i,file in enumerate(fileList1):
  if ".out" in file:
    backExp1 = GetBackExp(BaseDir1+file)
    backExpError1 = GetBackExpError(BaseDir1+file)
    sigEff1 = GetSigEff(BaseDir1+file)
    signalName1 = GetModelName(file)
    for i,file2 in enumerate(fileList2):
      if ".out" in file2:
        signalName2 = GetModelName(BaseDir2+file2)
        if signalName2==signalName1:
          break
    if signalName1 != signalName2:
      print 'Error: signalNames do not match: '+signalName1+' '+signalName2
      continue
    backExp2 = GetBackExp(BaseDir2+file2)
    backExpError2 = GetBackExpError(BaseDir2+file2)
    sigEff2 = GetSigEff(BaseDir2+file2)
    if backExp1==0:
      if backExp2==0:
        backExpErrorComb = 0
      else:
        backExpErrorComb = backExpError2/backExp2
    else:
      if backExp2==0:
        backExpErrorComb = backExpError1/backExp1
      else:
        backExpErrorComb = math.sqrt(math.pow(backExpError1/backExp1,2)+math.pow(backExpError2/backExp2,2))
    if backExp2 != 0:
      backExpErrorComb*=(backExp1/backExp2)
    printString = string.ljust(signalName1,15)
    if backExp2 != 0:
      printString+=string.center(str(round(backExp1/backExp2,4))+' +/- '+str(round(backExpErrorComb,4)),20)
    else:
      printString+=string.center("0",2)+' +/- '+str(0)
    if sigEff2 != 0:
      printString+=string.center(str(round(sigEff1/sigEff2,4)),10)
    else:
      printString+=string.center("0",10)
    print printString

