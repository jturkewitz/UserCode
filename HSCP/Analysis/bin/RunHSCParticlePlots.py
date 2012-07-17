#!/usr/bin/env python

import urllib
import string
import os
import sys
import HSCParticlePlotsLaunch
import glob
from subprocess import call
import datetime
from SignalDefinitions import *

RunCondor = True
FarmDirBeg = "FARM_MakeHSCParticlePlots_data_2012_"
now = datetime.datetime.now()
Date = now.strftime("%b%d")
BaseCfg = "makeHSCParticlePlots_template_cfg.py"
JobsName = "makeHSCParticlePlots_"
JobsName+=Date
JobsName+="_"
FarmDirBeg+=Date

IntLumiBefTrigChange = 352.067 #Total luminosity taken before RPC L1 trigger change (went into effect on run 165970)
##IntLumi = 4976 # 1/pb (2011, new pixel measurement)
IntLumi = 77 # 1/pb (2012, few runs)

# initialize dirs and shell file
HSCParticlePlotsLaunch.SendCluster_Create(FarmDirBeg,JobsName,RunCondor)

### 2012 data
DataCrossSection = 0.0
DataMassCut = 0.0
DataIsMC = "False"
DataIsHSCP = "False"
DataInputFiles=[]
DataInputFiles.append("    'file:/local/cms/user/turkewitz/HSCP/data/Data_RunA_190645-191100.root'")
for index,item in enumerate(DataInputFiles):
  HSCParticlePlotsLaunch.SendCluster_Push("Data2012_"+str(index), BaseCfg, item, DataCrossSection,
                      DataMassCut, IntLumi,IntLumiBefTrigChange,DataIsMC, DataIsHSCP)


# SIGNAL MC
HSCPisMC = "True"
HSCPisHSCP = "True"
#
for model in modelList:
  fileName = SignalBasePath+model.name+"BX1.root'"
  HSCParticlePlotsLaunch.SendCluster_Push(model.name+'BX1',BaseCfg,fileName,
                                          model.crossSection,model.massCut,
                                          IntLumi,IntLumiBefTrigChange,HSCPisMC,HSCPisHSCP)
  fileName = SignalBasePath+model.name+".root'"
  HSCParticlePlotsLaunch.SendCluster_Push(model.name,BaseCfg,fileName,
                                          model.crossSection,model.massCut,
                                          IntLumi,IntLumiBefTrigChange,HSCPisMC,HSCPisHSCP)



HSCParticlePlotsLaunch.SendCluster_Submit()
