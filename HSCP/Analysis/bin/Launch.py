#!/usr/bin/env python

import urllib
import string
import os
import sys
import HSCPMakeIasPredictionsLaunch
import HSCPMakeScaledPredictionsLaunch
import HSCPDoLimitsLaunch
import HSCPDoSignificanceLaunch
import glob
from SignalDefinitions import *
from subprocess import call
import datetime

now = datetime.datetime.now()
Date = now.strftime("%b%d")

# Set things here
AllSlices = True
RunCondor = True
BayesianLimit = False
DoSignificanceToys = False
DoMass = False
ModelListToys = []
PoiList = [1e-1,8e-2,5e-2,2e-2,8e-3,6e-3,5e-3,4e-3,3e-3,2e-3,1e-3,7.5e-4,5e-4,2.5e-4]
QueueName = '8nh' # LSF
## FOR NOW, CutPt=Std means CutIas Std also!
CutPt = 50
CutIas = 0.1
#CutPt = 'Std'
#CutIas = 'Std'

FarmDirectory = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_stop200N_moreToys_Jul12'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_glu300N400N_moreToys_Jul12'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_moreModels_moreToys_Jul11'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_pureFreqAndInitialFit_overallSigSystOnly_moreModels_Jul10'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_noBGShapeSyst_moreModels_Jul10'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_stdSysts_pureFreqAndInitialFit_moreModels_Jul10'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVaryFIX_allModelsMoreToys_cutPt50GeVcutIas0.1_allSlices_Jul9'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVaryFIX_moreModels_cutPt50GeVcutIas0.1_allSlices_Jul5'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVaryFIX_allModels_cutPt50GeVcutIas0.1_allSlices_Jul5'
#FarmDirectory = 'FARM_50IasB_emptyB1e-25_1gausSyst_noStatErr_halfIasHalfNoMAsNominal_handleZeroBGSlices_bgShapeVary_allModels_cutPt50GeVcutIas0.1_allSlices_Jul5'

##
#FarmDirectory+='cutPt'
#FarmDirectory+=str(CutPt)
#FarmDirectory+='GeVcutIas'
#FarmDirectory+=str(CutIas)
#if(AllSlices):
#  FarmDirectory+='_allSlices_'
#else:
#  FarmDirectory+='_oneSlice_'
#FarmDirectory+=Date

InputDataRootFile = os.getcwd()
InputDataRootFile+='/FARM_MakeHSCParticlePlots_Data_Apr28/outputs/makeHSCParticlePlots_Data2011_all.root'
#InputDataRootFile+='/FARM_MakeHSCParticlePlots_dataWithTightRPCTrig_absEta_ptErrorPresel_feb6/outputs/makeHSCParticlePlots_feb1_Data2011_all.root'
#InputDataRootFile+='/FARM_MakeHSCParticlePlots_data_absEta_ptErrorPresel_Mar08/outputs/makeHSCParticlePlots_Data2011_all.root'
#InputDataRootFile+='/FARM_MakeHSCParticlePlots_absEta_ptErrorPresel_Mar21/outputs/makeHSCParticlePlots_Mar21_Data2011_all.root'
sigInput = os.getcwd()
#sigInput+='/FARM_MakeHSCParticlePlots_Signals_deDxSystematics_May25/outputs/makeHSCParticlePlots_May25_'
sigInput+='/FARM_MakeHSCParticlePlots_Signals_Jul06/outputs/makeHSCParticlePlots_Jul06_'
IntLumi = 4976 # 1/pb (2011, new pixel measurement)
#

OutputIasPredictionDir = os.getcwd()
OutputIasPredictionDir+="/"
OutputIasPredictionDir+=FarmDirectory
OutputIasPredictionDir+="/outputs/makeIasPredictions/"
bgInput = OutputIasPredictionDir+'makeIasPredictionsCombined_'
BaseCfgMakeScaled = "makeScaledPredictionHistograms_template_cfg.py"
BaseCfgMakePred = "makeIasPredictions_template_cfg.py"
if(DoMass):
  BaseChXML = "hscp_dataDriven_ch_mass_template.xml"
else:
  BaseChXML = "hscp_dataDriven_ch_template.xml"
BaseChXMLSig = "hscp_dataDriven_ch_sig_template.xml"
BaseCombXML = "hscp_dataDriven_template.xml"
BaseCombXMLSig = "hscp_dataDriven_sig_template.xml"
BaseMacro = os.getcwd() + "/TestStandardHypoTestInvDemo.C"
BaseMacroBayesian = os.getcwd() + "/StandardBayesianNumericalDemo.C"
#BaseMacroBayesian = os.getcwd() + "/StandardBayesianMCMCDemo.C"
BaseMacroSignificance = os.getcwd() + "/StandardHypoTestDemo.C"


def printConfiguration():
  if(CutPt=='Std'):
    print 'Using standard analysis optimized Ias and Pt cuts'
  elif(CutIas=='Std'):
    print 'Using standard analysis optimized Ias cut only (pt=50)'
  else:
    print 'Using Pt=50 Ias=',CutIas
  if(AllSlices):
    print 'Doing all eta/NoM slices'
  else:
    print 'Doing single eta/NoM slice only'


def checkForIasPredictions(modelName,massCut,ptCut,iasCut):
  # have to change this if slicing is changed in HSCPMakeIasPredictionsLaunch.py
  etaCuts = [0.0,0.4,0.8,1.2]
  #nomCuts = [5,11]
  nomCuts = [5,13]
  for minEta in etaCuts:
    for minNom in nomCuts:
      fileName = os.getcwd()+'/'+FarmDirectory+'/outputs/makeIasPredictions/makeIasPredictions_'+modelName+'_massCut'+str(massCut)+'_pt'+str(ptCut)+'_ias'+str(iasCut)+'_etaMin'+str(minEta)+'_nomMin'+str(minNom)+'.root'
      try:
        open(fileName)
      except IOError:
        print 'File: ',fileName,' does not exist.'
        return False
  return True


def checkForScaledPredictions(modelName, checkForLimitFile):
  if(checkForLimitFile):
    fileName = os.getcwd()+'/'+FarmDirectory+'/outputs/makeScaledPredictions/'+modelName+'_Limits/hscp_combined_hscp_model.root'
  else:
    fileName = os.getcwd()+'/'+FarmDirectory+'/outputs/makeScaledPredictions/'+modelName+'_Significance/hscp_combined_hscp_model.root'
  try:
    open(fileName)
  except IOError:
    print 'File: ',fileName,' does not exist.'
    return False
  return True

def checkForLimitResult(model):
  #fileName = os.getcwd()+'/'+FarmDirectory+'/outputs/doLimits/doLimits_'+model.name+'_massCut'+str(model.massCut)+'_ptCut50_index0.root'
  outputsDir = os.getcwd()+'/'+FarmDirectory+'/outputs/doLimits'
  doLimitsFiles = os.listdir(outputsDir)
  for fileName in doLimitsFiles:
    if model.name in fileName:
      try:
        open(outputsDir+'/'+fileName)
      except IOError:
        print 'File: ',fileName,' does not exist.'
        return False
      return True
  return False


def doIasPredictions():
  # in the real analysis we will use the same ih/pt cuts for sidebands
  # so make a set of the mass cuts (no duplicates)
  printConfiguration()
  JobName = "makeIasPredictions_"
  for model in modelList:
    if(CutPt=='Std'):
      HSCPMakeIasPredictionsLaunch.SendCluster_Create(FarmDirectory,JobName+model.name+"_",InputDataRootFile,
                                                BaseCfgMakePred, model.massCut, model.ptCut, model.iasCut, AllSlices, RunCondor, QueueName)
    elif(CutIas=='Std'):
      HSCPMakeIasPredictionsLaunch.SendCluster_Create(FarmDirectory,JobName+model.name+"_",InputDataRootFile,
                                                BaseCfgMakePred, model.massCut, 50, model.iasCut, AllSlices, RunCondor, QueueName)
    else:
      HSCPMakeIasPredictionsLaunch.SendCluster_Create(FarmDirectory,JobName+model.name+"_",InputDataRootFile,
                                                BaseCfgMakePred, model.massCut, 50, CutIas, AllSlices, RunCondor, QueueName)

  HSCPMakeIasPredictionsLaunch.SendCluster_Submit()


def doMergeIasPredictions():
  printConfiguration()
  fnull = open(os.devnull, 'w')
  if(CutPt=='Std'):
    for model in modelList:
      if(checkForIasPredictions(model.name,model.massCut,model.ptCut,model.iasCut)):
        thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt'+str(model.ptCut)+'_ias'+model.iasCut+'_eta*')
        thisMassCutFiles.insert(0,OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut'+str(model.ptCut)+'_ias'+model.iasCut+'.root')
        thisMassCutFiles.insert(0,"-f")
        thisMassCutFiles.insert(0,"hadd")
        print 'Merging files for ' + model.name + ' mass cut = ' + str(model.massCut) + ' pt cut = ' + str(model.ptCut) + ' ias cut = ' + model.iasCut
        call(thisMassCutFiles,stdout = fnull)
      else:
        print 'ERROR: Ias predictions for',model.name,'not found!'
        break
  elif(CutIas=='Std'):
    for model in modelList:
      if(checkForIasPredictions(model.name,model.massCut,CutPt,model.iasCut)):
        thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt50_ias'+model.iasCut+'_eta*')
        thisMassCutFiles.insert(0,OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut50_ias'+model.iasCut+'.root')
        thisMassCutFiles.insert(0,"-f")
        thisMassCutFiles.insert(0,"hadd")
        print 'Merging files for ' + model.name + ': mass cut = ' + str(model.massCut) + ' pt cut = 50 ias cut = ' + model.iasCut
        call(thisMassCutFiles,stdout = fnull)
      else:
        print 'ERROR: Ias predictions for',model.name,'not found!'
        break
  else:
    for model in modelList:
      if(checkForIasPredictions(model.name,model.massCut,CutPt,CutIas)):
        thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt50_ias'+str(CutIas)+'_eta*')
        thisMassCutFiles.insert(0,OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut50_ias'+str(CutIas)+'.root')
        thisMassCutFiles.insert(0,"-f")
        thisMassCutFiles.insert(0,"hadd")
        print 'Merging files for ' + model.name + ' mass cut = ' + str(model.massCut) + ' pt cut = 50 ias cut = ' + str(CutIas)
        call(thisMassCutFiles,stdout = fnull)
      else:
        print 'ERROR: Ias predictions for',model.name,'not found!'
        break
  fnull.close()

def doCopyIasPredictions():
  printConfiguration()
  if(CutPt=='Std'):
    for model in modelList:
      thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt'+str(model.ptCut)+'_ias'+model.iasCut+'_eta*')
      thisMassCutFiles.append(OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut'+str(model.ptCut)+'_ias'+model.iasCut+'.root')
      thisMassCutFiles.insert(0,"cp")
      print 'Copy file for ' + model.name + ' mass cut = ' + str(model.massCut) + ' pt cut = ' + str(model.ptCut) + ' ias cut = ' + model.iasCut
      call(thisMassCutFiles)
  elif(CutIas=='Std'):
    for model in modelList:
      thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt50_ias'+model.iasCut+'_eta*')
      thisMassCutFiles.append(OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut50_ias'+model.iasCut+'.root')
      thisMassCutFiles.insert(0,"cp")
      print 'Copy file for ' + model.name + ': mass cut = ' + str(model.massCut) + ' pt cut = 50 ias cut = ' + model.iasCut
      call(thisMassCutFiles)
  else:
    for model in modelList:
      thisMassCutFiles = glob.glob(OutputIasPredictionDir+'*'+model.name+'_massCut'+str(model.massCut)+'_pt50_ias'+str(CutIas)+'_eta*')
      thisMassCutFiles.append(OutputIasPredictionDir+'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut50_ias'+str(CutIas)+'.root')
      thisMassCutFiles.insert(0,"cp")
      print 'Copy file for ' + model.name + ' mass cut = ' + str(model.massCut) + ' pt cut = 50 ias cut = ' + str(CutIas)
      call(thisMassCutFiles)
    
def doScaledPredictions():
  printConfiguration()
  # BG
  JobName = "makeScaledPredictions_"
  HSCPMakeScaledPredictionsLaunch.SendCluster_Create(FarmDirectory,JobName,
                                              IntLumi, BaseCfgMakeScaled, BaseChXML, 
                                                           BaseCombXML, BaseChXMLSig, BaseCombXMLSig, AllSlices, RunCondor, QueueName)
                                                           
  for model in modelList:
    looseRPCFileName = sigInput+model.name+'BX1.root'
    tightRPCFileName = sigInput+model.name+'.root'
    if(CutPt=='Std'):
      HSCPMakeScaledPredictionsLaunch.SendCluster_Push(
          bgInput+model.name+'_massCut',looseRPCFileName,tightRPCFileName,model.name,model.massCut,model.iasCut,model.ptCut,model.crossSection)
    elif(CutIas=='Std'):
      HSCPMakeScaledPredictionsLaunch.SendCluster_Push(
          bgInput+model.name+'_massCut',looseRPCFileName,tightRPCFileName,model.name,model.massCut,model.iasCut,50,model.crossSection)
    else:
      HSCPMakeScaledPredictionsLaunch.SendCluster_Push(
          bgInput+model.name+'_massCut',looseRPCFileName,tightRPCFileName,model.name,model.massCut,CutIas,50,model.crossSection)

  HSCPMakeScaledPredictionsLaunch.SendCluster_Submit()


def doLimits():
  printConfiguration()
  JobName = "doLimits_"
  if(BayesianLimit):
    print 'Doing Bayesian limit computation'
    print 'Run expected limits'
    HSCPDoLimitsLaunch.SendCluster_Create(FarmDirectory, JobName, IntLumi, BaseMacroBayesian, True, True, RunCondor, QueueName)
    for model in modelList:
      fileExists = checkForScaledPredictions(model.name,True)
      if not fileExists:
        return
      else:
        HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+model.name+'.root',model.massCut,model.iasCut,model.ptCut)
    #HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+'Gluino1000'+'.root',Gluino1000.massCut,0.1,50)
    HSCPDoLimitsLaunch.SendCluster_Submit()
    print 'Run observed limits'
    HSCPDoLimitsLaunch.SendCluster_Create(FarmDirectory, JobName, IntLumi, BaseMacroBayesian, True, False, RunCondor, QueueName)
    for model in modelList:
      fileExists = checkForScaledPredictions(model.name,True)
      if not fileExists:
        return
      else:
        HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+model.name+'.root',model.massCut,model.iasCut,model.ptCut)
    HSCPDoLimitsLaunch.SendCluster_Submit()
  else:
    print 'Doing CLs limit computation'
    #print 'Run expected limits'
    #HSCPDoLimitsLaunch.SendCluster_Create(FarmDirectory, JobName, IntLumi, BaseMacro, False, True, RunCondor, QueueName)
    ##for model in modelList:
    ##  fileExists = checkForScaledPredictions(model.name,True)
    ##  if not fileExists:
    ##    return
    ##  else:
    ##    HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+model.name+'.root',model.massCut,model.iasCut,model.ptCut)
    #HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+'Gluino1000'+'.root',Gluino1000.massCut,0.1,50)
    #HSCPDoLimitsLaunch.SendCluster_Submit()
    print 'Run observed limits'
    HSCPDoLimitsLaunch.SendCluster_Create(FarmDirectory, JobName, IntLumi, BaseMacro, False, False, RunCondor, QueueName)
    for model in modelList:
      fileExists = checkForScaledPredictions(model.name,True)
      if not fileExists:
        return
      else:
        HSCPDoLimitsLaunch.SendCluster_Push(bgInput,sigInput+model.name+'.root',model.massCut,CutIas,CutPt)
    HSCPDoLimitsLaunch.SendCluster_Submit()


def combineLimitResults():
  # used to combine limit files from CLs
  # must execute makeCombineLimitResults to build the binary first
  doLimitsOutputDir = FarmDirectory+'/outputs/doLimits/'
  for model in modelList:
    atLeastOneFile = checkForLimitResult(model)
    if(not atLeastOneFile):
      print 'ERROR: Files for model:',model.name,'not found!'
      continue
    print 'Combining limit results from ',model.name
    logFileName = FarmDirectory+'/logs/doLimits/doLimitsCombined_'+model.name+'.out'
    shell_file=open('/tmp/combine.sh','w')
    shell_file.write('#!/bin/sh\n')
    # CERN
    #shell_file.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    #shell_file.write('source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh\n')
    #shell_file.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
    # UMN
    shell_file.write('source /local/cms/user/cooper/root/bin/thisroot.sh\n')
    shell_file.write(os.getcwd()+'/combineLimitResults '+doLimitsOutputDir+' '+model.name+' |tee '+logFileName+'\n')
    shell_file.close()
    os.system("chmod 777 /tmp/combine.sh")
    os.system("source /tmp/combine.sh")
   

def doSignificance():
  printConfiguration()
  for model in modelList:
    fileExists = checkForScaledPredictions(model.name,False)
    if not fileExists:
      return
  JobName = "doSignificance_"
  inputWorkspaceFile = 'hscp_combined_hscp_model.root'
  inputWorkspacePathBeg = os.getcwd()+'/'+FarmDirectory+'/outputs/makeScaledPredictions/'
  # use observed data as "black line" or test point
  # first do asymptotic for all
  print 'Running asymptotic significance test for all models'
  ThrowToysSignificance = False
  HSCPDoSignificanceLaunch.SendCluster_Create(FarmDirectory, JobName, BaseMacroSignificance, ThrowToysSignificance, False, RunCondor, QueueName)
  for model in modelList:
    for poi in PoiList:
      HSCPDoSignificanceLaunch.SendCluster_Push(model.name,inputWorkspacePathBeg+model.name+'_Significance/'+inputWorkspaceFile,poi)
  HSCPDoSignificanceLaunch.SendCluster_Submit()
  if(DoSignificanceToys):
    # then throw toys for a few test models
    print
    print 'Running toys for selected models'
    ThrowToysSignificance = True
    HSCPDoSignificanceLaunch.SendCluster_Create(FarmDirectory, JobName, BaseMacroSignificance, ThrowToysSignificance, False, RunCondor, QueueName)
    for model in ModelListToys:
      poi = 0.001
      HSCPDoSignificanceLaunch.SendCluster_Push(model.name,inputWorkspacePathBeg+model.name+'_Significance/'+inputWorkspaceFile,poi)
      #for poi in PoiList:
      #  HSCPDoSignificanceLaunch.SendCluster_Push(model.name,inputWorkspacePathBeg+model.name+'_Significance/'+inputWorkspaceFile,poi)
    HSCPDoSignificanceLaunch.SendCluster_Submit()

def combineHypoTestResults():
  # for combining significance test with toys -- CLs
  doSignificanceOutputDir = FarmDirectory+'/outputs/doSignificance'
  for model in ModelListToys:
    for poi in PoiList:
      doSignificanceMerge = []
      doSignificanceMerge.append(os.getcwd()+'/combineHypoTestResults')
      doSignificanceMerge.append(doSignificanceOutputDir)
      doSignificanceMerge.append(model.name)
      doSignificanceMerge.append(str(poi))
      print 'Merging files for ' + model.name + ' poi = ' + str(poi)
      print doSignificanceMerge
      call(doSignificanceMerge)
  
def plotIasPredictionHists():
  fnull = open(os.devnull, 'w')
  IasPredictionDir = FarmDirectory+'/outputs/makeIasPredictions/'
  if os.path.isdir(FarmDirectory+'/plots') == False:
    os.mkdir(FarmDirectory+'/plots')
  for model in modelList:
    fileName = 'makeIasPredictionsCombined_'+model.name+'_massCut'+str(model.massCut)+'_ptCut'+str(CutPt)+'_ias'+str(CutIas)+'.root'
    filePath = IasPredictionDir+fileName
    try:
      open(filePath)
    except IOError:
      print 'File: ',filePath,' does not exist.'
      break
    print 'Making plots for ',model.name
    dirName = 'makeIasPredictionPlots.'+model.name+'_'+Date
    command = ['./makeIasPredictionPlotsHtml.sh','-d',dirName,'-r',filePath]
    call(command,stdout = fnull)
    command = ['mv',dirName,FarmDirectory+'/plots/']
    call(command)
  fnull.close()


def Usage():
  print 'Usage: Launch.py [arg] where arg can be:'
  print '    0 --> makeHSCParticlePlots (process into RooDataSets)'
  print '    1 --> mergeHSCParticlePlots for data'
  print '    2 --> makeIasPredictions'
  print '    3 --> mergeIasPredictions'
  print '    4 --> makeScaledPredictionPlots (process input for limit-setting macro)'
  print '    5 --> do limits (expected and observed)'
  print '    6 --> combine limit results (when using CLs limits)'
  print '    7 --> significance test using CLs'
  print '    8 --> combine significance test results for toys'
  print '    9 --> make ias prediction plots and html'



# Running the functions

if len(sys.argv)==1:
  Usage()
  sys.exit()

elif sys.argv[1]=='0':
  print 'Step 0: MakeHSCParticlePlots'
  print 'TODO: Implement this'

elif sys.argv[1]=='1':
  print 'Step 1: MergeHSCParticlePlots for data'
  print 'TODO: Implement this'

elif sys.argv[1]=='2':
  print 'Step 2: Make Ias predictions for background'
  doIasPredictions()

elif sys.argv[1]=='3':
  print 'Step 3: Merge Ias predictions'
  if(AllSlices):
    doMergeIasPredictions()
  else:
    doCopyIasPredictions()

elif sys.argv[1]=='4':
  print 'Step 4: Make scaled prediction plots'
  doScaledPredictions()

elif sys.argv[1]=='5':
  print 'Step 5: Compute limits'
  doLimits()

elif sys.argv[1]=='6':
  print 'Step 6: Merge CLs limit results'
  combineLimitResults()

elif sys.argv[1]=='7':
  print 'Step 7: Calculate significance using CLs'
  doSignificance()

elif sys.argv[1]=='8':
  print 'Step 8: Merge significance test toys'
  combineHypoTestResults()

elif sys.argv[1]=='9':
  print 'Step 9: Make all ias prediction plots and html'
  plotIasPredictionHists()
  

else:
  print 'Did not understand input.'
  Usage()
  sys.exit()


