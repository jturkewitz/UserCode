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
##Is8TeV = "True"
##FarmDirBeg = "FARM_MakeHSCParticlePlots_data_2012_"
FarmDirBeg = "FARM_MakeHSCParticlePlots_data_2012_testingMC_Shift_"
now = datetime.datetime.now()
Date = now.strftime("%b%d")
BaseCfg = "makeHSCParticlePlots_template_cfg.py"
JobsName = "makeHSCParticlePlots_"
JobsName+=Date
JobsName+="_"
FarmDirBeg+=Date

IntLumiBefTrigChange = 352.067 #Total luminosity taken before RPC L1 trigger change (went into effect on run 165970)
##IntLumiBefTrigChange = 0.0 #Total luminosity taken before RPC L1 trigger change (went into effect on run 165970)
##IntLumi = 4976 # 1/pb (2011, new pixel measurement, i think this should be 0.0 for 2012, but it may not matter, should test)
IntLumi = 18823 
##IntLumi = 4801 # 1/pb (2012, few runs estimate - needs to be revised)

# initialize dirs and shell file
HSCParticlePlotsLaunch.SendCluster_Create(FarmDirBeg,JobsName,RunCondor)

### 2012 data
DataCrossSection = 0.0
DataMassCut = 0.0
DataIs8TeV = "True"
DataIsMC = "False"
DataIsHSCP = "False"
DataInputFiles=[]
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012A_190645_190999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012A_190645_190999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012A_191000_191999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012A_193000_193621.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012B_193622_193999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012B_194000_194999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012B_195000_195999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012B_196000_196531.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_198000_198345.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_198488_198919.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_198920_198999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_199000_199999.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_200000_200532.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_200533_202016.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012C_202017_203002.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_203601_203900.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_203901_204200.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_204201_204500.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_204501_204800.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_204801_205100.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_205101_205400.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_205401_205700.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_205701_206000.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_206001_206300.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_206301_206600.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_206601_206900.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_206901_207200.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_207201_207500.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_207501_207800.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_207801_208100.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_208101_208357.root'")
##DataInputFiles.append("    'file:/hdfs/cms/user/turkewitz/hscp/data/Run2012D_208358_208686.root'")
for index,item in enumerate(DataInputFiles):
  HSCParticlePlotsLaunch.SendCluster_Push("Data2012_"+str(index), BaseCfg, item, DataCrossSection,
                      DataMassCut, IntLumi, IntLumiBefTrigChange, DataIs8TeV, DataIsMC, DataIsHSCP)

# SIGNAL MC
HSCPis8TeV = "True"
HSCPisMC = "True"
HSCPisHSCP = "True"
#
for model in modelList:
##  fileName = SignalBasePath+model.name+"BX1.root'"
##  HSCParticlePlotsLaunch.SendCluster_Push(model.name+'BX1',BaseCfg,fileName,
##                                          model.crossSection,model.massCut,
##                                          IntLumi,IntLumiBefTrigChange,HSCPis8TeV,HSCPisMC,HSCPisHSCP)
  fileName = SignalBasePath+model.name+".root'"
  HSCParticlePlotsLaunch.SendCluster_Push(model.name,BaseCfg,fileName,
                                          model.crossSection,model.massCut,
                                          IntLumi,IntLumiBefTrigChange,HSCPis8TeV,HSCPisMC,HSCPisHSCP)



HSCParticlePlotsLaunch.SendCluster_Submit()
