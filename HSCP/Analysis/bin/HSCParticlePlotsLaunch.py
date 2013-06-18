#!/usr/bin/env python

import urllib
import string
import os
import sys
import glob
from decimal import *

CopyRights  = '########################################\n'
CopyRights += '#     HSCParticlePlotsLaunch script    #\n'
CopyRights += '#        seth.cooper@cern.ch           #\n'
CopyRights += '#             Feb 2012                 #\n'
CopyRights += '#     modified from script by Loic     #\n'
CopyRights += '########################################\n'

Path_Cfg          = ''
Jobs_Index        = ''
Jobs_Inputs       = []
Jobs_Name         = ''
Jobs_Count        = 0
Farm_Directories  = []
Path_Log          = ''
Path_Shell        = ''
RunCondor         = False


def CreateTheConfigFile(jobName, baseCfg, inputFile, crossSection, massCut, intLumi, intLumiBefTrigChange,
                        is8TeV, isMC, isHSCP):
    global CopyRights
    global Farm_Directories
    global Jobs_Name
    path_Cfg   = Farm_Directories[1]+Jobs_Name+jobName+'_cfg.py'
    outputFile = Jobs_Name+jobName+'.root'
    config_file=open(baseCfg,'r')
    config_txt   = '\n\n' + CopyRights + '\n\n'
    config_txt  += config_file.read()
    config_file.close()

    #print 'inputFiles= ',inputFiles
    #Default Replacements
    #config_txt = config_txt.replace("XXX_INPUTFILES_XXX", "    '"+inputFile+"'")
    config_txt = config_txt.replace("XXX_INPUTFILES_XXX", inputFile)
    config_txt = config_txt.replace("XXX_OUTPUTFILE_XXX", outputFile)
    config_txt = config_txt.replace("XXX_CROSSSECTION_XXX" , str(crossSection))
    config_txt = config_txt.replace("XXX_MASSCUT_XXX", str(massCut))
    config_txt = config_txt.replace("XXX_IS8TEV_XXX", is8TeV)
    config_txt = config_txt.replace("XXX_ISMC_XXX", isMC)
    config_txt = config_txt.replace("XXX_ISHSCP_XXX", isHSCP)
    config_txt = config_txt.replace("XXX_INTEGRATEDLUMI_XXX" , str(intLumi))
    config_txt = config_txt.replace("XXX_INTEGRATEDLUMIBEFTRIGCHANGE_XXX" , str(intLumiBefTrigChange))

    config_file=open(path_Cfg,'w')
    config_file.write(config_txt)
    config_file.close()

def CreateTheShellFile(jobName,baseCfg,inputFile,crossSection,massCut,intLumi,intLumiBefTrigChange,is8TeV,isMC,isHSCP):
    global Path_Shell
    global CopyRights
    global Jobs_Name
    CreateTheConfigFile(jobName,baseCfg,inputFile,crossSection,massCut,intLumi,intLumiBefTrigChange,is8TeV,isMC,isHSCP)
    outputFile = Jobs_Name+jobName+'.root'
    path_Cfg   = Farm_Directories[1]+Jobs_Name+jobName+'_cfg.py'
    Path_Shell = Farm_Directories[1]+Jobs_Name+jobName+'.sh'
    file = inputFile.lstrip(" 'file:")
    file = file.strip("'")
    shell_file=open(Path_Shell,'w')
    shell_file.write('#!/bin/sh\n')
    shell_file.write(CopyRights + '\n')
    #shell_file.write('export SCRAM_ARCH=slc5_amd64_gcc434\n')
    #shell_file.write('source /local/cms/sw/cmsset_default.sh\n')
    shell_file.write('source /local/cms/sw/cmsset_CMSSW5X.sh\n')
    shell_file.write('cd ' + os.getcwd() + '\n')
    shell_file.write('eval `scramv1 runtime -sh`\n')
    shell_file.write('cd -\n')
    shell_file.write('head %s >& /dev/null\n' % file)
    shell_file.write('sleep 1\n')
    shell_file.write('head %s >& /dev/null\n' % file)
    shell_file.write('sleep 1\n')
    shell_file.write('head %s >& /dev/null\n' % file)
    shell_file.write('sleep 1\n')
    #shell_file.write('echo Making the plots for %s\n' % jobName)
    shell_file.write('makeHSCParticlePlots %s\n' % path_Cfg)
    shell_file.write('mv ' + outputFile + ' ' + os.getcwd() + '/' + Farm_Directories[4] )
    shell_file.close()
    os.system('chmod 777 '+Path_Shell)

def CreateTheCmdFile():
    global Path_Cmd
    global CopyRights
    global Jobs_Name
    global RunCondor
    Path_Cmd   = Farm_Directories[0]+Jobs_Name+'submit.sh'
    cmd_file=open(Path_Cmd,'w')
    cmd_file.write(CopyRights + '\n')
    if(RunCondor):
      cmd_file.write('Universe                = vanilla\n')
      cmd_file.write('Environment             = CONDORJOBID=$(Process)\n')
      cmd_file.write('notification            = Error\n')
      cmd_file.write('requirements            = (Memory > 1024)&&(Arch=?="X86_64")&&(Machine=!="zebra01.spa.umn.edu")&&(Machine=!="zebra02.spa.umn.edu")&&(Machine=!="zebra03.spa.umn.edu")&&(Machine=!="caffeine.spa.umn.edu")&&(Machine=!="gc1-ce-sl5.spa.umn.edu.spa.umn.edu")&&(Machine=!="hadoop-test.spa.umn.edu")\n')
      cmd_file.write('+CondorGroup            = "cmsfarm"\n')
      cmd_file.write('should_transfer_files   = NO\n')
      cmd_file.write('Notify_user = turkewitz@physics.umn.edu\n')
    else:
      cmd_file.write('# List of bash commands\n')
    cmd_file.close()

def AddJobToCmdFile(jobName):
    global Path_Shell
    global Path_Cmd
    global Jobs_Name
    global RunCondor
    Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+jobName
    Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+jobName
    cmd_file=open(Path_Cmd,'a')
    cmd_file.write('\n')
    if(RunCondor):
      cmd_file.write('Executable              = %s\n'     % Path_Shell)
      cmd_file.write('output                  = %s.out\n' % Path_Log)
      cmd_file.write('error                   = %s.err\n' % Path_Error)
      cmd_file.write('log                     = %s.log\n' % Path_Log)
      cmd_file.write('Queue 1\n')
    else:
      cmd_file.write('echo "Executing: '+Path_Shell+' >& '+Path_Log+'"\n')
      cmd_file.write('source '+Path_Shell+' >& '+Path_Log+'\n')
    cmd_file.close()

def CreateDirectoryStructure(FarmDirectory, jobsName):
    global Farm_Directories
    global Jobs_Name
    Jobs_Name = jobsName
    Farm_Directories = [FarmDirectory+'/', FarmDirectory+'/inputs/', FarmDirectory+'/logs/',
                        FarmDirectory+'/errors/', FarmDirectory+'/outputs/']
    for i in range(0,len(Farm_Directories)):
        if os.path.isdir(Farm_Directories[i]) == False:
            os.system('mkdir ' + Farm_Directories[i])

def SendCluster_Create(FarmDirectory,jobsName,runCondor):
    global RunCondor
    RunCondor = runCondor
    CreateDirectoryStructure(FarmDirectory,jobsName)
    CreateTheCmdFile()

def SendCluster_Push(jobName,baseCfg,inputFile,crossSection,massCut,intLumi,intLumiBefTrigChange,is8TeV,isMC,isHSCP):
    global Jobs_Count
    global Jobs_Index
    Jobs_Index = "%04i" % Jobs_Count
    CreateTheShellFile(jobName,baseCfg,inputFile,crossSection,massCut,intLumi,intLumiBefTrigChange,is8TeV,isMC,isHSCP)
    AddJobToCmdFile(jobName)
    Jobs_Count = Jobs_Count+1

def SendCluster_Submit():
    global CopyRights
    global Jobs_Count
    global Path_Cmd
    global RunCondor
    if(RunCondor):
      os.system("condor_submit " + Path_Cmd)
      print '\n'+CopyRights
      print '%i job(s) has/have been submitted to condor' % Jobs_Count
    else:
     print 'Run everything locally'
     os.system("source " + Path_Cmd)

