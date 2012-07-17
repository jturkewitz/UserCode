#!/usr/bin/env python

import urllib
import string
import os
import sys
import glob
from subprocess import call
from SignalDefinitions import *

OutputIasPredictionDir = os.getcwd()
OutputIasPredictionDir+="/FARM_moreModels_cutPt50GeV_Feb21/"
OutputIasPredictionDir+="outputs/makeIasPredictions/"

#for model in modelList:
#    thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*massCut'+str(model.massCut)+'_pt'+str(model.ptCut)+'_ias'+str(model.iasCut)+'_eta*')
#    thisMassCutFiles.insert(0,OutputIasPredictionDir+'makeIasPredictionsCombined_massCut'+str(model.massCut)+'_ptCut'+str(model.ptCut)+'_ias'+str(model.iasCut)+'.root')
#    thisMassCutFiles.insert(0,"hadd")
#    print 'Merging files for mass cut = ' + str(model.massCut) + ' pt cut = ' + str(model.ptCut) + ' ias cut = ' + str(model.iasCut)
#    call(thisMassCutFiles)

# all Pt = 50
for model in modelList:
    thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*massCut'+str(model.massCut)+'_pt50_ias'+str(model.iasCut)+'_eta*')
    thisMassCutFiles.insert(0,OutputIasPredictionDir+'makeIasPredictionsCombined_massCut'+str(model.massCut)+'_ptCut50_ias'+str(model.iasCut)+'.root')
    thisMassCutFiles.insert(0,"hadd")
    print 'Merging files for mass cut = ' + str(model.massCut) + ' pt cut = 50 ias cut = ' + str(model.iasCut)
    call(thisMassCutFiles)
