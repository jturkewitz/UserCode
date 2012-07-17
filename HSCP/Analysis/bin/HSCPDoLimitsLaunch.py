#!/usr/bin/env python

import urllib
import string
import os
import sys
import glob
from decimal import *

CopyRights  = '####################################################\n'
CopyRights += '#           HSCPDoLimitsLaunch script              #\n'
CopyRights += '#               seth.cooper@cern.ch                #\n'
CopyRights += '#                    Mar 2012                      #\n'
CopyRights += '#         modified from script by Loic             #\n'
CopyRights += '####################################################\n'

Path_Cfg          = ''
Jobs_Count        = 0
Jobs_Name         = ''
Jobs_Index        = ''
#
Farm_Directories  = []
Path_Cmd          = ''
Path_Shell        = ''
Path_Log          = ''
Base_Cfg = ''
Int_Lumi = 0
Base_xml_channel = ''
Base_xml_combined = ''
Do_Bayesian = False
Condor = True
Queue_Name = ''
All_Slices = True
Expected_Limits = True


def CreateTheShellFile(bgInputFile,sigInputFile,massCut,iasCut,ptCut,xSecMin,xSecMax,index):
    global Path_Shell
    global CopyRights
    global Jobs_Name
    global Base_macro
    global Do_Bayesian
    global Expected_Limits
    global Condor
    posLastUndsc = sigInputFile.rfind("_")
    signalName = sigInputFile[posLastUndsc+1:len(sigInputFile)-5]
    #endRemoved = sigInputFile[0:posLastUndsc]
    #posLastUndsc = endRemoved.rfind("_")
    #signalName = endRemoved[posLastUndsc+1:len(endRemoved)]
    #print signalName
    if(Expected_Limits):
      dataName = 'asimovData'
      prepend = 'expected_'
    else:
      dataName = 'obsData'
      prepend = ''
    if(Do_Bayesian):
      outputFile = 'doLimits_'+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_bayesian.root'
      Path_Shell = Farm_Directories[1]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_bayesian.sh'
      shell_file=open(Path_Shell,'w')
      shell_file.write('#!/bin/sh\n')
      shell_file.write(CopyRights + '\n')
      shell_file.write('# use later root\n')
      # UMN
      if(Condor):
        shell_file.write('source /local/cms/user/cooper/root/bin/thisroot.sh\n')
      else:
        # CERN
        shell_file.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
        shell_file.write('source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh\n')
        shell_file.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
      shell_file.write('cd ' + os.getcwd() + '/' + Farm_Directories[0] + 'outputs/makeScaledPredictions/' + signalName + '_Limits\n')
      shell_file.write('root -l -b -q "'+Base_macro
        +'(\\"hscp_combined_hscp_model.root\\",\\"combined\\",\\"ModelConfig\\",\\"'+dataName+'\\",\\"'+outputFile+'\\")"'+'\n')
    else:
      # CLs
      # variables for limit setting
      #nPoints = 8 # per job
      nPoints = 1 # per job
      nToys = 2000 # per xSec trial point
      #nToys = 5000 # for high stats test; since each gets run twice we end up with nToysx2 total toys per point
      # end variables for limit setting
      outputFile = 'doLimits_'+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)+'.root'
      Path_Shell = Farm_Directories[1]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)+'.sh'
      shell_file=open(Path_Shell,'w')
      shell_file.write('#!/bin/sh\n')
      shell_file.write(CopyRights + '\n')
      shell_file.write('# use later root\n')
      # UMN
      if(Condor):
        shell_file.write('source /local/cms/user/cooper/root/bin/thisroot.sh\n')
      else:
        # CERN
        shell_file.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
        shell_file.write('source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh\n')
        shell_file.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
      shell_file.write('cd ' + os.getcwd() + '/' + Farm_Directories[0] + 'outputs/makeScaledPredictions/' + signalName + '_Limits\n')
      shell_file.write('root -l -b -q "'+Base_macro
          # Frequentist calculator
          +'(\\"hscp_combined_hscp_model.root\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",0,3,true,'
          # Hybrid calculator
          #+'(\\"hscp_combined_hscp_model.root\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",1,3,true,'
          #asymptCLs#+'(\\"hscp_combined_hscp_model.root\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",2,3,true,'
          +str(nPoints)+','
          +str(xSecMin)+','
          +str(xSecMax)+','
          +str(nToys)+','
          +'false,\\"\\",\\"'
          +outputFile
          +'\\")"'
          + '\n')
    shell_file.write('mv ' + outputFile + ' ' + os.getcwd() + '/' + Farm_Directories[4] + '\n\n')
    shell_file.close()
    os.system("chmod 777 "+Path_Shell)


def CreateTheCmdFile():
    global Path_Cmd
    global CopyRights
    global Jobs_Name
    global Farm_Directories
    global Condor
    global Queue_Name
    global Expected_Limits
    if(Expected_Limits):
      Path_Cmd   = Farm_Directories[0]+Jobs_Name+'expected_submit.sh'
    else:
      Path_Cmd   = Farm_Directories[0]+Jobs_Name+'submit.sh'
    cmd_file=open(Path_Cmd,'w')
    if(Condor):
      cmd_file.write(CopyRights + '\n')
      #cmd_file.write('cd ' + os.getcwd() + Farm_Directories[1] + '\n')
      #cmd_file.write('# init root and prepare histfactory\n')
      #cmd_file.write('source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh\n')
      #cmd_file.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.32.00/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
      #cmd_file.write('prepareHistFactory\n')
      #cmd_file.write('cd -\n\n')
      cmd_file.write('Universe                = vanilla\n')
      cmd_file.write('Environment             = CONDORJOBID=$(Process)\n')
      cmd_file.write('notification            = Error\n')
      cmd_file.write('requirements            = (Memory > 1024)&&(Arch=?="X86_64")&&(Machine=!="zebra01.spa.umn.edu")&&(Machine=!="zebra02.spa.umn.edu")&&(Machine=!="zebra03.spa.umn.edu")&&(Machine=!="caffeine.spa.umn.edu")\n')
      #cmd_file.write('+CondorGroup            = "cmsfarm"\n')
      cmd_file.write('+CondorGroup            = "twins"\n')
      cmd_file.write('should_transfer_files   = NO\n')
      cmd_file.write('Notify_user = cooper@physics.umn.edu\n')
      #cmd_file.write('should_transfer_files   = YES\n')
      #cmd_file.write('when_to_transfer_output = ON_EXIT\n')
    else:
      cmd_file.write('#!/bin/bash\n')
      cmd_file.write(CopyRights + '\n')
      cmd_file.write('# list all bsub commands' + '\n')
    cmd_file.close()

def AddJobToCmdFile(massCut,ptCut,sigInputFile,index):
    global Path_Shell
    global Path_Cmd
    global Jobs_Name
    global Do_Bayesian
    global Condor
    global Queue_Name
    global Expected_Limits
    #posLastUndsc = sigInputFile.rfind("_")
    #endRemoved = sigInputFile[0:posLastUndsc]
    #posLastUndsc = endRemoved.rfind("_")
    #signalName = endRemoved[posLastUndsc+1:len(endRemoved)]
    posLastUndsc = sigInputFile.rfind("_")
    signalName = sigInputFile[posLastUndsc+1:len(sigInputFile)-5]
    if(Expected_Limits):
      prepend = 'expected_'
    else:
      prepend = ''
    if(Do_Bayesian):
      Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_bayesian'
      Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_bayesian'
    else:
      Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)
      Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+prepend+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)
    cmd_file=open(Path_Cmd,'a')
    cmd_file.write('\n')
    if(Condor):
      cmd_file.write('Executable              = %s\n'     % Path_Shell)
      cmd_file.write('output                  = %s.out\n' % Path_Log)
      cmd_file.write('error                   = %s.err\n' % Path_Error)
      cmd_file.write('log                     = %s.log\n' % Path_Log)
      cmd_file.write('Queue 1\n')
      #cmd_file.write(os.getcwd() + '/' + Path_Shell + '\n')
    else:
      Path_Log+='.log'
      Path_Error+='.err'
      cmd_file.write('bsub -q "%s" '     % Queue_Name)
      cmd_file.write(' -o %s' % Path_Log)
      cmd_file.write(' -e %s ' % Path_Error)
      cmd_file.write(os.getcwd() + '/' + Path_Shell + '\n')
    cmd_file.close()

def CreateDirectoryStructure(FarmDirectory):
    global Jobs_Name
    global Farm_Directories
    Farm_Directories = [FarmDirectory+'/', FarmDirectory+'/inputs/doLimits/',
                        FarmDirectory+'/logs/doLimits/',
                        FarmDirectory+'/errors/doLimits/',
                        FarmDirectory+'/outputs/doLimits/']
    for i in range(0,len(Farm_Directories)):
        if os.path.isdir(Farm_Directories[i]) == False:
            os.system('mkdir -p ' + Farm_Directories[i])

def SendCluster_Create(FarmDirectory, jobName, intLumi, baseMacro, doBayesian, doExpLimits,
                        doCondor, queueName):
    global Jobs_Name
    global Jobs_Count
    global Base_macro
    global Int_Lumi
    global Farm_Directories
    global Do_Bayesian
    global Expected_Limits
    global Condor
    global Queue_Name
    Jobs_Name  = jobName
    Jobs_Count = 0
    Base_macro = baseMacro
    Int_Lumi=intLumi
    Do_Bayesian=doBayesian
    Expected_Limits=doExpLimits
    Condor = doCondor
    Queue_Name = queueName
    CreateDirectoryStructure(FarmDirectory)
    CreateTheCmdFile()


def SendCluster_Push(bgInputFilesBase,sigInputFile,massCut,iasCut,ptCut):
    global Jobs_Count
    global Jobs_Index
    global Do_Bayesian
    Jobs_Index = "%04i" % Jobs_Count
    bgInputFile = bgInputFilesBase+`massCut`+'_ptCut'+`ptCut`+'.root'
    posLastUndsc = sigInputFile.rfind("_")
    signalName = sigInputFile[posLastUndsc+1:len(sigInputFile)-5]
    if(Do_Bayesian):
      CreateTheShellFile(bgInputFile,sigInputFile,massCut,iasCut,ptCut,0,0,0)
      AddJobToCmdFile(massCut,ptCut,sigInputFile,0)
      Jobs_Count = Jobs_Count+1
    else:
      xSecMin = 0.0001
      #xSecMax = 0.01
      xSecMax = 0.1
      # zoom in on limit area
      if(signalName=='Gluino300'):
        xSecMin = 0.002
        xSecMax = 0.004
      elif(signalName=='Gluino400'):
        xSecMin = 0.002
        xSecMax = 0.004
      elif(signalName=='Gluino500'):
        xSecMin = 0.0015
        xSecMax = 0.0035
      elif(signalName=='Gluino600'):
        xSecMin = 0.0015
        xSecMax = 0.0035
      elif(signalName=='Gluino700'):
        xSecMin = 0.0015
        xSecMax = 0.0035
      elif(signalName=='Gluino800'):
        xSecMin = 0.002
        xSecMax = 0.003
      elif(signalName=='Gluino900'):
        xSecMin = 0.002
        xSecMax = 0.004
      elif(signalName=='Gluino1000'):
        xSecMin = 0.002
        xSecMax = 0.005
      elif(signalName=='Gluino1100'):
        xSecMin = 0.002
        xSecMax = 0.005
      elif(signalName=='Gluino1200'):
        xSecMin = 0.003
        xSecMax = 0.006
      ## GluinoN
      if(signalName=='Gluino300N'):
        xSecMin = 0.01
        xSecMax = 0.05
      elif(signalName=='Gluino400N'):
        xSecMin = 0.01
        xSecMax = 0.05
      elif(signalName=='Gluino500N'):
        xSecMin = 0.009
        xSecMax = 0.02
      elif(signalName=='Gluino600N'):
        xSecMin = 0.009
        xSecMax = 0.02
      elif(signalName=='Gluino700N'):
        xSecMin = 0.008
        xSecMax = 0.02
      elif(signalName=='Gluino800N'):
        xSecMin = 0.008
        xSecMax = 0.02
      elif(signalName=='Gluino900N'):
        xSecMin = 0.008
        xSecMax = 0.02
      elif(signalName=='Gluino1000N'):
        xSecMin = 0.008
        xSecMax = 0.02
      elif(signalName=='Gluino1100N'):
        xSecMin = 0.009
        xSecMax = 0.02
      elif(signalName=='Gluino1200N'):
        xSecMin = 0.009
        xSecMax = 0.02
      ## Stop
      elif(signalName=='Stop200'):
        xSecMin = 0.003
        xSecMax = 0.005
      elif(signalName=='Stop300'):
        xSecMin = 0.002
        xSecMax = 0.004
      elif(signalName=='Stop400'):
        xSecMin = 0.0015
        xSecMax = 0.0035
      elif(signalName=='Stop500'):
        xSecMin = 0.0015
        xSecMax = 0.0035
      elif(signalName=='Stop600'):
        xSecMin = 0.001
        xSecMax = 0.003
      elif(signalName=='Stop700'):
        xSecMin = 0.001
        xSecMax = 0.003
      elif(signalName=='Stop800'):
        xSecMin = 0.001
        xSecMax = 0.003
      ## StopN
      elif(signalName=='Stop200N'):
        xSecMin = 0.008
        xSecMax = 0.025
      elif(signalName=='Stop300N'):
        xSecMin = 0.008
        xSecMax = 0.015
      elif(signalName=='Stop400N'):
        xSecMin = 0.008
        xSecMax = 0.012
      elif(signalName=='Stop500N'):
        xSecMin = 0.006
        xSecMax = 0.01
      elif(signalName=='Stop600N'):
        xSecMin = 0.005
        xSecMax = 0.009
      elif(signalName=='Stop700N'):
        xSecMin = 0.006
        xSecMax = 0.008
      elif(signalName=='Stop800N'):
        xSecMin = 0.005
        xSecMax = 0.008
      ## GMStau
      elif(signalName=='GMStau200'):
        xSecMin = 0.001
        xSecMax = 0.002
      elif(signalName=='GMStau247'):
        xSecMin = 0.001
        xSecMax = 0.002
      elif(signalName=='GMStau308'):
        xSecMin = 0.0007
        xSecMax = 0.0015
      elif(signalName=='GMStau370'):
        xSecMin = 0.0007
        xSecMax = 0.0015
      elif(signalName=='GMStau432'):
        xSecMin = 0.0005
        xSecMax = 0.0015
      elif(signalName=='GMStau494'):
        xSecMin = 0.0005
        xSecMax = 0.0015
      # too low mass
      elif(signalName=='GMStau100'):
        xSecMin = 0.003
        xSecMax = 0.005
      elif(signalName=='GMStau126'):
        xSecMin = 0.003
        xSecMax = 0.005
      elif(signalName=='GMStau156'):
        xSecMin = 0.0015
        xSecMax = 0.0025

      nSteps = 12 # gets multiplied by points/job later
      stepSize = (xSecMax-xSecMin)/nSteps

      for index in range(0,nSteps):
        xSecMinPt = xSecMin+index*stepSize
        xSecMaxPt = xSecMin+(index+1)*stepSize
        # repeat same toys a few times
        repeatTimes = 20
        for count in range(0,repeatTimes):
          CreateTheShellFile(bgInputFile,sigInputFile,massCut,iasCut,ptCut,xSecMinPt,xSecMaxPt,Jobs_Count+count)
          AddJobToCmdFile(massCut,ptCut,sigInputFile,Jobs_Count+count)
        Jobs_Count = Jobs_Count+repeatTimes

def SendCluster_Submit():
    global CopyRights
    global Jobs_Count
    global Path_Cmd
    if(Condor):
      os.system("condor_submit " + Path_Cmd)
      print '\n'+CopyRights
      print '%i job(s) has/have been submitted to condor' % Jobs_Count
    else:
      os.system("source " + Path_Cmd)
      print '\n'+CopyRights
      print '%i job(s) has/have been submitted to lxbatch' % Jobs_Count

