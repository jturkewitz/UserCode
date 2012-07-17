#!/usr/bin/env python

import urllib
import string
import os
import sys
import glob
from decimal import *

CopyRights  = '####################################################\n'
CopyRights += '#      HSCPMakeIasPredictionsLaunch script         #\n'
CopyRights += '#            seth.cooper@cern.ch                   #\n'
CopyRights += '#                 Mar 2012                         #\n'
CopyRights += '#         modified from script by Loic             #\n'
CopyRights += '####################################################\n'

Path_Cfg          = ''
Jobs_Count        = 0
Jobs_Name         = ''
#
Farm_Directories  = []
Path_Cmd          = ''
Path_Shell        = ''
Path_Log          = ''
Base_Cfg = ''
Input_File = ''
Condor = True
Queue_Name = ''
All_Slices = True
#
Init = True

def GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh):
    toReturn='massCut'+`massCut`+'_pt'+str(ptThresh)+'_ias'+str(iasThresh)+'_etaMin'+str(etaMin)+'_nomMin'+`nomMin`
    return toReturn

def CreateTheConfigFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh):
    global Jobs_Name
    global Jobs_Count
    global Path_Cfg
    global CopyRights
    global Base_cfg
    global Input_File
    Path_Cfg = Farm_Directories[1]+Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)+'_cfg.py'
    outputFile = Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)+'.root'
    config_file=open(Base_Cfg,'r')
    config_txt   = '\n\n' + CopyRights + '\n\n'
    config_txt  += config_file.read()
    config_file.close()

    #Default Replacements
    config_txt = config_txt.replace("XXX_INPUTFILE_XXX", Input_File)
    config_txt = config_txt.replace("XXX_OUTPUTFILE_XXX"    , outputFile)
    config_txt = config_txt.replace("XXX_MASSCUT_XXX"       , `massCut`)
    config_txt = config_txt.replace("XXX_ETAMIN_XXX"        , str(etaMin))
    config_txt = config_txt.replace("XXX_ETAMAX_XXX"        , str(etaMax))
    config_txt = config_txt.replace("XXX_NOMMIN_XXX"        , `nomMin`)
    config_txt = config_txt.replace("XXX_NOMMAX_XXX"        , `nomMax`)
    config_txt = config_txt.replace("XXX_PTTHRESH_XXX"      , str(ptThresh))
    config_txt = config_txt.replace("XXX_IASTHRESH_XXX"     , str(iasThresh))

    config_file=open(Path_Cfg,'w')
    config_file.write(config_txt)
    config_file.close()


def CreateTheShellFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh):
    global Path_Shell
    global CopyRights
    global Jobs_Name
    global All_Slices
    global Condor
    CreateTheConfigFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh)
    outputFile = Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)+'.root'
    Path_Shell = Farm_Directories[1]+Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)+'.sh'
    shell_file=open(Path_Shell,'w')
    shell_file.write('#!/bin/sh\n')
    shell_file.write(CopyRights + '\n')
    # if at UMN
    if(Condor):
      shell_file.write('source /local/cms/sw/cmsset_default.sh\n')
    else:
      # CERN
      shell_file.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
      shell_file.write('source /afs/cern.ch/cms/sw/cmsset_default.sh\n')
    #
    shell_file.write('cd ' + os.getcwd() + '\n')
    shell_file.write('eval `scramv1 runtime -sh`\n')
    #shell_file.write('cd -\n')
    shell_file.write('cd ' + Farm_Directories[0] + '\n')
    if(All_Slices):
      shell_file.write('makeIasPredictions ' + os.getcwd() + '/'+Path_Cfg + '\n')
    else:
      shell_file.write('makeIasPredictionsOneSlice ' + os.getcwd() + '/'+Path_Cfg + '\n')
    #shell_file.write('cd -\n')
    shell_file.write('mv ' + outputFile + ' ' + os.getcwd() + '/' + Farm_Directories[4] + '\n')
    shell_file.close()
    os.system('chmod 777 '+Path_Shell)

def CreateTheCmdFile():
    global Path_Cmd
    global CopyRights
    global Jobs_Name
    global Condor
    Path_Cmd = Farm_Directories[0]+Jobs_Name+'submit.sh'
    cmd_file=open(Path_Cmd,'w')
    if(Condor):
      cmd_file.write(CopyRights + '\n')
      cmd_file.write('Universe                = vanilla\n')
      cmd_file.write('Environment             = CONDORJOBID=$(Process)\n')
      cmd_file.write('notification            = Error\n')
      cmd_file.write('requirements            = (Memory > 1024)&&(Arch=?="X86_64")&&(Machine=!="caffeine.spa.umn.edu")\n')
      cmd_file.write('+CondorGroup            = "cmsfarm"\n')
      #cmd_file.write('+CondorGroup            = "twins"\n')
      cmd_file.write('should_transfer_files   = NO\n')
      cmd_file.write('Notify_user = cooper@physics.umn.edu\n')
      #cmd_file.write('should_transfer_files   = YES\n')
      #cmd_file.write('when_to_transfer_output = ON_EXIT\n')
    else:
      cmd_file.write('#!/bin/bash\n')
      cmd_file.write(CopyRights + '\n')
      cmd_file.write('# File listing all bsub commands' + '\n')
    cmd_file.close()


def AddJobToCmdFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh):
    global Path_Shell
    global Path_Cmd
    global Jobs_Name
    global Queue_Name
    global Condor
    Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)
    Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+GetSelectionString(massCut,etaMin,nomMin,ptThresh,iasThresh)
    cmd_file=open(Path_Cmd,'a')
    cmd_file.write('\n')
    if(Condor):
      cmd_file.write('Executable              = %s\n'     % Path_Shell)
      cmd_file.write('output                  = %s.out\n' % Path_Log)
      cmd_file.write('error                   = %s.err\n' % Path_Error)
      cmd_file.write('log                     = %s.log\n' % Path_Log)
      cmd_file.write('Queue 1\n')
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
    Farm_Directories = [FarmDirectory+'/', FarmDirectory+'/inputs/makeIasPredictions/', 
                        FarmDirectory+'/logs/makeIasPredictions/',
                        FarmDirectory+'/errors/makeIasPredictions/',
                        FarmDirectory+'/outputs/makeIasPredictions/']
    for i in range(0,len(Farm_Directories)):
        if os.path.isdir(Farm_Directories[i]) == False:
            os.system('mkdir -p ' + Farm_Directories[i])

def SendCluster_Create(farmDirectory, jobName, inputRootFile, baseCfg,
                       massCut, ptCut, iasCut, allSlices, doCondor, queueName):
    global Jobs_Name
    global Jobs_Count
    global Input_File
    global Base_Cfg
    global Init
    global Queue_Name
    global All_Slices
    global Condor
    Jobs_Name  = jobName
    if(Init):
      Jobs_Count = 0
      Input_File = inputRootFile
      Base_Cfg = baseCfg
      Condor = doCondor
      Queue_Name = queueName
      All_Slices = allSlices
      Init = False
      CreateDirectoryStructure(farmDirectory)
      CreateTheCmdFile()
    # if this splitting is changed, have to change checkForIasPredictions function in Launch.py
    if(allSlices):
      for etaIndex in range(0,16,4):
        #SendCluster_Push(massCut,etaIndex/10.0,(etaIndex+4)/10.0,5,10,ptCut,iasCut)
        #SendCluster_Push(massCut,etaIndex/10.0,(etaIndex+4)/10.0,11,18,ptCut,iasCut)
        SendCluster_Push(massCut,etaIndex/10.0,(etaIndex+4)/10.0,5,12,ptCut,iasCut)
        SendCluster_Push(massCut,etaIndex/10.0,(etaIndex+4)/10.0,13,18,ptCut,iasCut)
    else:
      SendCluster_Push(massCut,0.0,1.5,5,21,ptCut,iasCut)


def SendCluster_Push(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh):
    global Jobs_Count
    CreateTheShellFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh)
    AddJobToCmdFile(massCut,etaMin,etaMax,nomMin,nomMax,ptThresh,iasThresh)
    Jobs_Count = Jobs_Count+1

def SendCluster_Submit():
    global CopyRights
    global Jobs_Count
    global Path_Cmd
    global Condor
    if(Condor):
      os.system("condor_submit " + Path_Cmd)
      print '\n'+CopyRights
      print '%i job(s) has/have been submitted to condor' % Jobs_Count
    else:
      os.system("source " + Path_Cmd)
      print '\n'+CopyRights
      print '%i job(s) has/have been submitted to lxbatch' % Jobs_Count

