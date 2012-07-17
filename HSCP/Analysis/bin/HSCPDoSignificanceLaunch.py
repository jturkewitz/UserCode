#!/usr/bin/env python

import urllib
import string
import os
import sys
import glob
from decimal import *

CopyRights  = '####################################################\n'
CopyRights += '#        HSCPDoSignificanceLaunch script           #\n'
CopyRights += '#               seth.cooper@cern.ch                #\n'
CopyRights += '#                    Apr 2012                      #\n'
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
Throw_Toys = False
Expected_Significance = True
Condor = True
Queue_Name = ''
All_Slices = True


def CreateTheShellFile(modelName,inputFile,poiValue,index):
    global Path_Shell
    global CopyRights
    global Jobs_Name
    global Base_macro
    global Throw_Toys
    global Expected_Significance
    global Condor
    if(Expected_Significance):
      dataName = 'asimovData'
      prepend = 'expected_'
    else:
      dataName = 'obsData'
      prepend = ''
    if(Throw_Toys):
      outputFile = 'doSignificance_'+prepend+modelName+'_poi'+str(poiValue)+'_'+str(index)+'.root'
      Path_Shell = Farm_Directories[1]+Jobs_Name+prepend+modelName+'_poi'+str(poiValue)+'_'+str(index)+'.sh'
    else:
      outputFile = 'doSignificance_'+prepend+modelName+'_poi'+str(poiValue)+'_asymptotic.root'
      Path_Shell = Farm_Directories[1]+Jobs_Name+prepend+modelName+'_poi'+str(poiValue)+'_asymptotic.sh'
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
      shell_file.write('source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.33.02/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh\n')
    #numToys = 25000
    numToys = 5000
    shell_file.write('cd ' + os.getcwd() + '/' + Farm_Directories[4] + '\n')
    if(Throw_Toys):
      shell_file.write('root -l -b -q "'+Base_macro
        +'(\\"'
        +inputFile
        #+'\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",0,3,'
        +'\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",1,3,'
        +str(numToys)
        +',false,\\"\\",\\"'
        +outputFile
        +'\\",'
        +str(poiValue)
        +')"'
        +'\n')
    else:
      shell_file.write('root -l -b -q "'+Base_macro
        +'(\\"'
        +inputFile
        +'\\",\\"combined\\",\\"ModelConfig\\",\\"\\",\\"'+dataName+'\\",2,3,'
        +str(numToys)
        +',false,\\"\\",\\"'
        +outputFile
        +'\\",'
        +str(poiValue)
        +')"'
        +'\n')
    #shell_file.write(Base_macro+' '+inputFile+' '+str(numToys)+' '+outputFile+' '+str(poiValue)+'\n')
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
    global Expected_Significance
    if(Expected_Significance):
      Path_Cmd   = Farm_Directories[0]+Jobs_Name+'expected_submit.sh'
    else:
      Path_Cmd   = Farm_Directories[0]+Jobs_Name+'submit.sh'
    cmd_file=open(Path_Cmd,'w')
    if(Condor):
      cmd_file.write(CopyRights + '\n')
      cmd_file.write('Universe                = vanilla\n')
      cmd_file.write('Environment             = CONDORJOBID=$(Process)\n')
      cmd_file.write('notification            = Error\n')
      cmd_file.write('requirements            = (Memory > 512)&&(Arch=?="X86_64")&&(Machine=!="zebra01.spa.umn.edu")&&(Machine=!="zebra02.spa.umn.edu")&&(Machine=!="zebra03.spa.umn.edu")&&(Machine=!="caffeine.spa.umn.edu")\n')
      cmd_file.write('+CondorGroup            = "cmsfarm"\n')
      #cmd_file.write('+CondorGroup            = "twins"\n')
      cmd_file.write('should_transfer_files   = NO\n')
      cmd_file.write('Notify_user = cooper@physics.umn.edu\n')
      #cmd_file.write('should_transfer_files   = YES\n')
      #cmd_file.write('when_to_transfer_output = ON_EXIT\n')
    else:
      cmd_file.write('#!/bin/bash\n')
      cmd_file.write(CopyRights + '\n')
      cmd_file.write('# list all bsub commands' + '\n')
    cmd_file.close()

def AddJobToCmdFile(modelName,inputFile,poiValue,index):
    global Path_Shell
    global Path_Cmd
    global Jobs_Name
    global Condor
    global Queue_Name
    global Throw_Toys
    global Expected_Significance
    #posLastUndsc = sigInputFile.rfind("_")
    #signalName = sigInputFile[posLastUndsc+1:len(sigInputFile)-5]
    #Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)
    #Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+signalName+'_massCut'+`massCut`+'_ptCut'+`ptCut`+'_index'+str(index)
    if(Expected_Significance):
      prepend = 'expected_'
    else:
      prepend = ''
    if(Throw_Toys):
      Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+prepend+modelName+'_'+str(poiValue)+'_index'+str(index)
      Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+prepend+modelName+'_'+str(poiValue)+'_index'+str(index)
    else:
      Path_Log   = os.getcwd()+'/'+Farm_Directories[2]+Jobs_Name+prepend+modelName+'_'+str(poiValue)+'_asymptotic'
      Path_Error   = os.getcwd()+'/'+Farm_Directories[3]+Jobs_Name+prepend+modelName+'_'+str(poiValue)+'_asymptotic'
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
      #cmd_file.write('sleep 15' + '\n')
    cmd_file.close()

def CreateDirectoryStructure(FarmDirectory):
    global Jobs_Name
    global Farm_Directories
    Farm_Directories = [FarmDirectory+'/', FarmDirectory+'/inputs/doSignificance/',
                        FarmDirectory+'/logs/doSignificance/',
                        FarmDirectory+'/errors/doSignificance/',
                        FarmDirectory+'/outputs/doSignificance/']
    for i in range(0,len(Farm_Directories)):
        if os.path.isdir(Farm_Directories[i]) == False:
            os.system('mkdir -p ' + Farm_Directories[i])

def SendCluster_Create(FarmDirectory, jobName, baseMacro, throwToys, expSignificance, doCondor, queueName):
    global Jobs_Name
    global Jobs_Count
    global Base_macro
    global Throw_Toys
    global Expected_Significance
    global Farm_Directories
    global Condor
    global Queue_Name
    Jobs_Name  = jobName
    Jobs_Count = 0
    Base_macro = baseMacro
    Throw_Toys = throwToys
    Expected_Significance = expSignificance
    Condor = doCondor
    Queue_Name = queueName
    CreateDirectoryStructure(FarmDirectory)
    CreateTheCmdFile()


def SendCluster_Push(modelName,inputFile,poiValue):
    global Jobs_Count
    global Jobs_Index
    global Throw_Toys
    Jobs_Index = "%04i" % Jobs_Count
    if(Throw_Toys):
      #numToyJobs = 400
      numToyJobs = 1
      for index in range(0,numToyJobs):
        CreateTheShellFile(modelName,inputFile,poiValue,index)
        AddJobToCmdFile(modelName,inputFile,poiValue,index)
        Jobs_Count = Jobs_Count+1
    else:
        CreateTheShellFile(modelName,inputFile,poiValue,0)
        AddJobToCmdFile(modelName,inputFile,poiValue,0)
        Jobs_Count = Jobs_Count+1

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

