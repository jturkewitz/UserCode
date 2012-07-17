import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('XXX_OUTPUTFILE_XXX'),
)


process.makeIasPredictions = cms.PSet(
    analysisCommonParameters,
    MassCut = cms.double(XXX_MASSCUT_XXX),
    EtaMax = cms.double(XXX_ETAMAX_XXX),
    EtaMin = cms.double(XXX_ETAMIN_XXX),
    NoMMax = cms.int32(XXX_NOMMAX_XXX),
    NoMMin = cms.int32(XXX_NOMMIN_XXX),
    UsePtForSideband = cms.bool(True),
    PtSidebandThreshold = cms.double(XXX_PTTHRESH_XXX),
    UseIasForSideband = cms.bool(True),
    IasSidebandThreshold = cms.double(XXX_IASTHRESH_XXX),
)

if process.makeIasPredictions.UsePtForSideband==cms.bool(True):
    print 'Running with pt threshold: ',
    print process.makeIasPredictions.PtSidebandThreshold
else:
    print 'Running with p threshold: ',
    print analysisCommonParameters.PSidebandThreshold

if process.makeIasPredictions.UseIasForSideband==cms.bool(True):
  print 'Running with Ias threshold:',
  print process.makeIasPredictions.IasSidebandThreshold
else:
  print 'Running with Ih threshold:',
  print analysisCommonParameters.IhSidebandThreshold

print 'Running with mass cut:',
print process.makeIasPredictions.MassCut
print 'nom range: ',
print process.makeIasPredictions.NoMMin,
print '-',
print process.makeIasPredictions.NoMMax
print 'eta range: ',
print process.makeIasPredictions.EtaMin,
print '-',
print process.makeIasPredictions.EtaMax

process.fwliteInput.fileNames.extend([
    'file:XXX_INPUTFILE_XXX'
]);

