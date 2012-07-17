#!/usr/bin/env python

import string
import os
import sys
import glob
import math
import array
from SignalDefinitions import *

def GetExpObsLim(filePath, modelName):
  expObsLim = []
  for line in open(filePath):
    if modelName in line:
      lineSplit = line.split()
      expObsLim.append(str(round(float(lineSplit[10]),4)))
      expObsLim.append(str(round(float(lineSplit[11]),4)))
  return expObsLim


#BaseDir1 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatErr_nominalAsNominal_cutPt50GeVcutIas0.1_allSlices_May28'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatErr_nomShiftAsNominal_cutPt50GeVcutIas0.1_allSlices_May26'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatErr_halfIasShiftAsNominal_cutPt50GeVcutIas0.1_allSlices_May28'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystWithStatErr400SortedBins_nominalAsNominal_cutPt50GeVcutIas0.1_allSlices_May28'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystWithStatErrIgnoreFirst400SortedBins_nominalAsNominal_cutPt50GeVcutIas0.1_allSlices_May28/'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_fullIasShiftAsNominal_cutPt50GeVcutIas0.1_allSlices_Jun05'
#BaseDir1 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_nominalAsNominal_highToyStatTestFewModels_cutPt50GeVcutIas0.1_allSlices_Jun05'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystWithStatErr400SortedBins_nominalAsNominal_highToyStatTestFewModels_cutPt50GeVcutIas0.1_allSlices_Jun01'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystWithStatErrIgnoreFirst400SortedBins_nominalAsNominal_highToyStatTestFewModels_cutPt50GeVcutIas0.1_allSlices_Jun04'
BaseDir1 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_nominalAsNominal_highToyStatAllModels_cutPt50GeVcutIas0.1_allSlices_Jun06'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_halfIasShiftAsNominalPlusShapeVariation_cutPt50GeVcutIas0.1_allSlices_Jun06'
#BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_nomShiftAsNominal_highToyStatsAllModels_cutPt50GeVcutIas0.1_allSlices_Jun06'
BaseDir2 = 'FARM_chk50IasB_emptyB1e-25_1gausSystNoStatError_halfIasShiftAsNominalNoShapeVariation_cutPt50GeVcutIas0.1_allSlices_Jun07'

BaseFile1 = BaseDir1+'/results/results.txt'
BaseFile2 = BaseDir2+'/results/results.txt'

print
print 'Using 1:',BaseFile1
print 'Using 2:',BaseFile2
print

titleString = string.ljust('Model',15)+string.center('ExpLim1',10)+string.center('ExpLim2',10)+string.center('ExpLim1/ExpLim2',20)
titleString+=string.center('ObsLim1',10)+string.center('ObsLim2',10)+string.center('ObsLim1/ObsLim2',10)
print titleString

avgExpLimRatio = 0.0
avgObsLimRatio = 0.0
modelCounter = 0
for model in modelList:
  expObsLim1 = GetExpObsLim(BaseFile1,model.name)
  expObsLim2 = GetExpObsLim(BaseFile2,model.name)
  if(len(expObsLim1) < 2 or len(expObsLim2) < 2):
    print 'ERROR with model',model.name
    continue
  printString = string.ljust(model.name,15)+string.center(expObsLim1[0],10)+string.center(expObsLim2[0],10)+string.center(str(round(float(expObsLim1[0])/float(expObsLim2[0]),4)),20)
  printString+=string.center(expObsLim1[1],10)+string.center(expObsLim2[1],10)+string.center(str(round(float(expObsLim1[1])/float(expObsLim2[1]),4)),10)
  print printString
  if(float(expObsLim1[0]) > 0 and float(expObsLim1[1]) > 0 and float(expObsLim2[0]) > 0 and float(expObsLim2[1]) > 0): 
    avgExpLimRatio+=float(expObsLim1[0])/float(expObsLim2[0])
    avgObsLimRatio+=float(expObsLim1[1])/float(expObsLim2[1])
    modelCounter+=1

print 'Average ExpLim1/ExpLim2 =',avgExpLimRatio/modelCounter
print 'Average ObsLim1/ObsLim2 =',avgObsLimRatio/modelCounter

