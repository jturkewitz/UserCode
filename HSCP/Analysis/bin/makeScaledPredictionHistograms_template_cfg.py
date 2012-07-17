import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

## MASS CUT HERE

process.makeScaledPredictionHistograms = cms.PSet(
    analysisCommonParameters,
    MassCut = cms.double(XXX_MASSCUT_XXX),
    IntegratedLuminosity= cms.double(XXX_INTLUMI_XXX),
    BackgroundPredictionInputRootFile = cms.string("XXX_BGINPUTFILE_XXX"),
    SignalLooseRPCInputRootFile = cms.string("XXX_SIGLOOSERPCINPUTFILE_XXX"),
    SignalTightRPCInputRootFile = cms.string("XXX_SIGTIGHTRPCINPUTFILE_XXX"),
    OutputRootFile = cms.string("XXX_OUTPUTFILE_XXX"),
    SignalCrossSection = cms.double(XXX_CROSSSECTION_XXX),

    IasCutForEfficiency = cms.double(XXX_IASCUT_XXX),
    UsePtForSideband = cms.bool(True),
    PtSidebandThreshold = cms.double(XXX_PTTHRESH_XXX),
    UseIasForSideband = cms.bool(True),
)

if process.makeScaledPredictionHistograms.UsePtForSideband==cms.bool(True):
    print 'Running with pt threshold: ',
    print process.makeScaledPredictionHistograms.PtSidebandThreshold
else:
    print 'Running with p threshold: ',
    print analysisCommonParameters.PSidebandThreshold

if process.makeScaledPredictionHistograms.UseIasForSideband==cms.bool(True):
  print 'Running with Ias threshold:',
  print process.makeScaledPredictionHistograms.IasCutForEfficiency
else:
  print 'Running with Ih threshold:',
  print analysisCommonParameters.IhSidebandThreshold
print 'Running with mass cut:',
print process.makeScaledPredictionHistograms.MassCut
