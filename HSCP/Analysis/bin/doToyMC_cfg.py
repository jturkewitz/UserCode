import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

process.doToyMC = cms.PSet(
    analysisCommonParameters,
    BackgroundPredictionInputRootFile = cms.string("XXX_BACKGROUNDINPUTFILE_XXX"),
    SignalInputRootFile = cms.string("XXX_SIGNALINPUTFILE_XXX"),
    OutputRootFile = cms.string("XXX_OUTPUT_XXX"),
    NumberOfBackgroundTracksInDRegion = cms.int32(XXX_NBACKTRACKSD_XXX), # after all skim/presels
    NumberOfSignalTracks = cms.int32(XXX_NSIGTRACKS_XXX), # for this int lumi, after all skim/presels
    LowerNoMOfSlice = cms.int32(XXX_LOWERNOM_XXX),
    LowerEtaOfSlice = cms.double(XXX_LOWERETA_XXX),
    NumberOfTrials = cms.int32(XXX_NUMTRIALS_XXX),
    GenerateBackgroundOnlySamples = cms.bool(XXX_BACKGROUNDONLY_XXX),
    Verbose = cms.bool(False),
    ReportEvery = cms.int32(1000)
)

print 'Running with p threshold: ',
print analysisCommonParameters.PSidebandThreshold
print 'Running with Ih threshold:',
print analysisCommonParameters.IhSidebandThreshold
print 'Running on lowerNoM',
print process.doToyMC.LowerNoMOfSlice
print 'Running on lowerEta',
print process.doToyMC.LowerEtaOfSlice
#print 'Doing:',
#print process.doToyMC.NumberOfTrials,
#print ' trials for',
#print process.doToyMC.NumberOfSignalTracks,
#print ' signal tracks and ',
#print process.doToyMC.NumberOfBackgroundTracksInDRegion,
#print ' background tracks in D region (after skim/preselections)'

