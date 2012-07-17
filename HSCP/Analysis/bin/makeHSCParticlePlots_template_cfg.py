import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),
    outputEvery = cms.uint32(100000),
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('XXX_OUTPUTFILE_XXX'),

)

process.makeHSCParticlePlots = cms.PSet(
    analysisCommonParameters,
    ScaleFactor = cms.double(1),
    MatchToHSCP = cms.bool(XXX_ISHSCP_XXX),
    IsMC = cms.bool(XXX_ISMC_XXX),
    MassCut = cms.double(XXX_MASSCUT_XXX),
    SignalEventCrossSection = cms.double(XXX_CROSSSECTION_XXX), #NLO, pb
    IntegratedLumi = cms.double(XXX_INTEGRATEDLUMI_XXX),
    IntegratedLumiBeforeTriggerChange = cms.double(XXX_INTEGRATEDLUMIBEFTRIGCHANGE_XXX),
)

#print 'Running with p threshold: ',
#print analysisCommonParameters.PSidebandThreshold
#print 'Running with Ih threshold:',
#print analysisCommonParameters.IhSidebandThreshold

process.fwliteInput.fileNames.extend([
    XXX_INPUTFILES_XXX
]);

