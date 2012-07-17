import FWCore.ParameterSet.Config as cms

process = cms.Process("GenerateHSCPToyMC")

# shaping our Message logger to suit our needs
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.generateHSCPToyMC = cms.EDAnalyzer("GenerateHSCPToyMC",
    GenerateSignalAndBackground = cms.bool(True),
    SignalFraction = cms.double(0.001828), #51./27893
    #SignalFraction = cms.double(0.000501918), #14./27893
    #SignalFraction = cms.double(0.000430215), #12./27893
    #SignalFraction = cms.double(0.000358513), #10./27893
    #SignalFraction = cms.double(0.00028681), #8./27893
    #SignalFraction = cms.double(0.000215108), #6./27893
    #SignalFraction = cms.double(0.000143405), #4./27893
    NumberOfTrials = cms.int32(10),
    EventsPerTrial = cms.int32(27944) # 27893+51
    #EventsPerTrial = cms.int32(27907) # 27893+14
    #EventsPerTrial = cms.int32(27905) # 27893+12
    #EventsPerTrial = cms.int32(27903) # 27893+10
    #EventsPerTrial = cms.int32(27901) # 27893+8
    #EventsPerTrial = cms.int32(27899) # 27893+6
    #EventsPerTrial = cms.int32(27897) # 27893+4
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource",
       numberEventsInRun = cms.untracked.uint32(1),
       firstRun = cms.untracked.uint32(1)
)


process.p = cms.Path(process.generateHSCPToyMC)

