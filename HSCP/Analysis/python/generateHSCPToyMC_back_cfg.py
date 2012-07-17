import FWCore.ParameterSet.Config as cms

process = cms.Process("GenerateHSCPToyMC")

# shaping our Message logger to suit our needs
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
    categories = cms.untracked.vstring('*'),
    destinations = cms.untracked.vstring('cout')
)

process.generateHSCPToyMC = cms.EDAnalyzer("GenerateHSCPToyMC",
    GenerateSignalAndBackground = cms.bool(False),
    SignalFraction = cms.double(0.0),
    NumberOfTrials = cms.int32(50),
    #EventsPerTrial = cms.int32(27907) # 27893+14
    EventsPerTrial = cms.int32(27944) # 27893+51
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )
process.source = cms.Source("EmptySource",
       numberEventsInRun = cms.untracked.uint32(1),
       firstRun = cms.untracked.uint32(1)
)


process.p = cms.Path(process.generateHSCPToyMC)

