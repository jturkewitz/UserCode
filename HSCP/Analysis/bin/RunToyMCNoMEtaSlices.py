#!/usr/bin/env python

import urllib
import string
import os
import sys
import HSCPToyMCLaunchOnLxBatch
import glob
from subprocess import call


if len(sys.argv[1])==0:
        print "Usage:"
        print "  RunIasNoMEtaSlicedSBDists.py nTrials"
        print "          nTrials is the number of trials to use"
        sys.exit()

else:
        nTrials = sys.argv[1]
        QueueName = "8nh"
        BaseCfg = "doToyMC_cfg.py"
        bgRootFile = "1.nov9.allSlices.EPSrunPostEPSprod.makeIasPredictions.root"
        sigRootFile = "nov1.withRooDataSet.makeHSCParticlePlots.GMstau200.root"
        nSigTracksD = 11
        nBackgroundTracksD = 71994
        Date = "nov18"

        # For both
        JobName = "toyMC_"
        JobName+=Date
        JobName+="_"
        CastorPathBeg = "/castor/cern.ch/user/s/scooper/hscp/428/ToyMC/"
        CastorPathBeg+=Date
        CastorPathBeg+="_"
        FarmDirBeg = "FARM_ToyMC_"
        FarmDirBeg+=Date
        FarmDirBeg+="_"

        # S+B
        FarmDirectory = FarmDirBeg+`nSigTracksD`
        FarmDirectory+="signalTracks_"
        FarmDirectory+=`nBackgroundTracksD`
        FarmDirectory+="backgroundTracks"
        backOnly = "False"
        sbCastorPath = CastorPathBeg+`nSigTracksD`
        sbCastorPath+="signalTracks_"
        sbCastorPath+=`nBackgroundTracksD`
        sbCastorPath+="backgroundTracks"
        HSCPToyMCLaunchOnLxBatch.SendCluster_Create(FarmDirectory,JobName,bgRootFile,sigRootFile,nTrials,
                                                    nSigTracksD,nBackgroundTracksD,backOnly,QueueName,BaseCfg,sbCastorPath)
        HSCPToyMCLaunchOnLxBatch.SendCluster_Submit()

        # B only
        FarmDirectory = FarmDirBeg+"0signalTracks_"
        FarmDirectory+=`nBackgroundTracksD`
        FarmDirectory+="backgroundTracks"
        backOnly = "True"
        bCastorPath = CastorPathBeg+"0signalTracks_"
        bCastorPath+=`nBackgroundTracksD`
        bCastorPath+="backgroundTracks"
        HSCPToyMCLaunchOnLxBatch.SendCluster_Create(FarmDirectory,JobName,bgRootFile,sigRootFile,nTrials,
                                                    nSigTracksD,nBackgroundTracksD,backOnly,QueueName,BaseCfg,bCastorPath)
        HSCPToyMCLaunchOnLxBatch.SendCluster_Submit()


