import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

## MASS CUT HERE
analysisCommonParameters.MassCut = 180.0

process.makeScaledPredictionHistograms = cms.PSet(
    analysisCommonParameters,
    IntegratedLuminosity= cms.double(1091),
    #BackgroundPredictionInputRootFile = cms.string("1.nov9.allSlices.EPSrunPostEPSprod.makeIasPredictions.root"),
    #BackgroundPredictionInputRootFile = cms.string("dec19.mass120.EPSrunPostEPSprod.makeIasPredictions.root"),
    #SignalInputRootFile = cms.string("nov1.withRooDataSet.makeHSCParticlePlots.GMstau200.root"),
    #OutputRootFile = cms.string("dec20.makeScaledPredictionHistograms.gmstau200.mass120.unrolledHists.root"),

    BackgroundPredictionInputRootFile = cms.string("dec19.mass180.EPSrunPostEPSprod.makeIasPredictions.root"),
    #BackgroundPredictionInputRootFile = cms.string("dec21.mass180.pt120.EPSrunPostEPSprod.makeIasPredictions.root"),
    SignalInputRootFile = cms.string("dec15.makeHSCParticlePlots.GMStau308.root"),
    OutputRootFile = cms.string("test.dec21.makeScaledPredictionHistograms.gmstau308"),

    #BackgroundPredictionInputRootFile = cms.string("dec15.mass420.EPSrunPostEPSprod.makeIasPredictions.root"),
    #BackgroundPredictionInputRootFile = cms.string("dec16.setOneSliceSmall.mass410.EPSrunPostEPSprod.makeIasPredictions.root"),
    #BackgroundPredictionInputRootFile = cms.string("dec16.normal.mass410.EPSrunPostEPSprod.makeIasPredictions.root"),
    #BackgroundPredictionInputRootFile = cms.string("dec19.oneSliceSmall.mass410.EPSrunPostEPSprod.makeIasPredictions.root"),
    #SignalInputRootFile = cms.string("dec9.makeHSCParticlePlots.gluino1000.tof.root"),
    #OutputRootFile = cms.string("dec20.makeScaledPredictionHistograms.gluino1000.mass410.unrolledHists.normal.root"),
    #OutputRootFile = cms.string("dec20.makeScaledPredictionHistograms.gluino1000.mass410.unrolledHists.oneSliceSmall.root"),

    IasCutForEfficiency = cms.double(0.4),
    UsePtForSideband = cms.bool(False),
    PtSidebandThreshold = cms.double(120.0),
)

if process.makeScaledPredictionHistograms.UsePtForSideband==cms.bool(True):
    print 'Running with pt threshold: ',
    print process.makeScaledPredictionHistograms.PtSidebandThreshold
else:
    print 'Running with p threshold: ',
    print analysisCommonParameters.PSidebandThreshold

print 'Running with Ih threshold:',
print analysisCommonParameters.IhSidebandThreshold
print 'Running with mass cut:',
print analysisCommonParameters.MassCut
