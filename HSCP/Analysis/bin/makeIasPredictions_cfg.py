import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),
)

process.fwliteOutput = cms.PSet(
    #fileName  = cms.string('1.dec1.allSlices.pt40.ih3p5.EPSrunPostEPSprod.makeIasPredictions.root'),
    #fileName  = cms.string('2.nov9.testOneSlice.EPSrunPostEPSprod.makeIasPredictions.root'),
    #fileName  = cms.string('dec19.oneSliceSmall.mass410.EPSrunPostEPSprod.makeIasPredictions.root'),
    fileName  = cms.string('test.dec21.EPSrunPostEPSprod.makeIasPredictions'),
)

analysisCommonParameters.MassCut = 180

process.makeIasPredictions = cms.PSet(
    analysisCommonParameters,
    UsePtForSideband = cms.bool(False),
    PtSidebandThreshold = cms.double(120.0),
)

if process.makeIasPredictions.UsePtForSideband==cms.bool(True):
    print 'Running with pt threshold: ',
    print process.makeIasPredictions.PtSidebandThreshold
else:
    print 'Running with p threshold: ',
    print analysisCommonParameters.PSidebandThreshold

print 'Running with Ih threshold:',
print analysisCommonParameters.IhSidebandThreshold
print 'Running with mass cut:',
print analysisCommonParameters.MassCut

process.fwliteInput.fileNames.extend([
    #'file:nov2.withRooDataSet.makeHSCParticlePlots.EPSrunPostEPSprod.root'
    #'file:nov30.withRooDataSet.pt40.makeHSCParticlePlots.EPSrunPostEPSprod.root'
    'file:dec20.p100.ih3p5.makeHSCParticlePlots.EPSrunPostEPSprod.root'
]);

