import FWCore.ParameterSet.Config as cms

analysisCommonParameters = cms.PSet(
  #MassCut = cms.double(100),
  # dE/dx calibration
  dEdx_k = cms.double(2.529),
  dEdx_c = cms.double(2.772),
  # search region
  PSidebandThreshold = cms.double(100),
  IhSidebandThreshold = cms.double(3.5),
  # preselections
  #minTrackNoH = 11;
  #minTrackValidFraction = 0.8;
  #minValidPixelHits = 2;
  #minIasNoM = 5;
  #minTofNdof = 8;
  #minTofNdofDt = 6;
  #minTofNdofCsc = 6;
  #minTrackQualityMask = 2;
  #maxTrackChi2OverNdf = 5.0;
  #minTrackPt = 35.0;
  #minIas = 0.0;
  #minIh = 3.0;
  #minTofInvBeta = 1.0;
  #maxTofIntBetaErr = 0.07;
  #minNumVertices = 1;
  #maxV3D = 2.0;
  #maxTrackEtIso = 50;
  #maxCalEOverPIso = 0.3;
  #maxTrackPtErr = 0.25;
  #maxTrackEta = 2.5;
)
