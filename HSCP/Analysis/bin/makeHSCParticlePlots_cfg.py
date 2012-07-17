import FWCore.ParameterSet.Config as cms

process = cms.PSet()

from HSCP2011.ToyMC.makeHSCParticlePlotsAndDoLikelihoodFitsCommonParameters_cfi import *

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(),
    maxEvents   = cms.int32(-1),
    outputEvery = cms.uint32(100000),
)

process.fwliteOutput = cms.PSet(
    #fileName  = cms.string('dec20.p100.ih3p5.makeHSCParticlePlots.EPSrunPostEPSprod.root'),
    #fileName  = cms.string('dec5.p100.ih3p5.makeHSCParticlePlots.GMstau200.root'),
    #fileName  = cms.string('dec9.makeHSCParticlePlots.GMStau200.tof.root'),
    #fileName  = cms.string('dec9.makeHSCParticlePlots.gluino1000.tof.root'),
    fileName  = cms.string('dec22.makeHSCParticlePlots.GMStau308.root'),

)

process.makeHSCParticlePlots = cms.PSet(
    analysisCommonParameters,
    MatchToHSCP = cms.bool(True),
    IsMC = cms.bool(True),
    #MatchToHSCP = cms.bool(False),
    #IsMC = cms.bool(False),
    ScaleFactor = cms.double(1),
    MassCut = cms.double(180.0),
    SignalEventCrossSection = cms.double(0.00098447), #NLO, pb
)

#print 'Running with p threshold: ',
#print analysisCommonParameters.PSidebandThreshold
#print 'Running with Ih threshold:',
#print analysisCommonParameters.IhSidebandThreshold

process.fwliteInput.fileNames.extend([
    # EPS runs -- from post EPS production (trigger thresh different?)
    # ==> 1.09/fb
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_160404_163869.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_165001_166033.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_166034_166500.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_166894_167151.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_166501_166893.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Data_RunA_167153_167913.root',

    # Jet-triggered data; should be same run range as EPS runs
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_May10ReReco_Jet_160404_163869/filtered_merged.root',
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_PromptRecoV4_Jet_165001_166033/filtered_merged.root',
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_PromptRecoV4_Jet_166034_166500/filtered_merged.root',
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_PromptRecoV4_Jet_166894_167151/filtered_merged.root',
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_PromptRecoV4_Jet_166501_166893/filtered_merged.root',
    #'file:/data3/scooper/data/hscp/428/Data/MyProcessing_SkimAndHSCParticles/RunA_PromptRecoV4_Jet_167153_167913/filtered_merged.root',

    # Z+Jet --> MuMu MonteCarlo Summer11
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-0to15.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-15to20.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-20to30.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-30to50.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-50to80.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-80to120.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-120to170.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-170to230.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-230to300.root',
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_ZJetToMuMu_Pt-300.root',
    # DY (m20) to MuMu
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/MC_DYToMuMu.root',
    # GMSB STAU 100
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/GMStau100.root'
    # GMSB STAU 200
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/GMStau200.root'
    # GMSB STAU 308
    'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/GMStau308.root'
    # Stop 200
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Stop200.root'
    # Gluino 1000
    #'file:/data3/scooper/data/hscp/428/Data/HSCParticles_HSCP2011_PostEPS_Sep30/Gluino1000.root'

]);

