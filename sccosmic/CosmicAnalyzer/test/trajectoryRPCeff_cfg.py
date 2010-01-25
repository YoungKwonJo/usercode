# The following comments couldn't be translated into the new config version:
#keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("TrackingTools.TrackRefitter.cosmicMuonTrajectories_cff")

# Conditions (Global Tag is used here):
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
process.GlobalTag.globaltag = "GR09_P_V8::All"#"GR09_P_V8_34X::All"#"GR09_P_V8_34X::All"#"GR09_P_V8::All"#"CRUZET4_V4P::All"
#CRAFT_V4P::All , CRUZET4_V4P::All
#process.GlobalTag.globaltag = 'COSMMC_22X_TK::All'
#process.GlobalTag.globaltag = "CRZT210_V3P::All"
process.prefer("GlobalTag")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

#process.load("Configuration.StandardSequences.MagneticField_0T_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")

process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("TrackingTools.TrackRefitter.cosmicMuonTrajectories_cff")

process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('traRPCeff.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(


        'file:///afs/cern.ch/user/y/youngjo/public/reco.root'


     )
)

process.muonAnalyzer = cms.EDAnalyzer("TrajectoryRPCEff",
    process.MuonServiceProxy,
    InputLabel = cms.InputTag("cosmicMuons"),#"standAloneMuons"), #"standAloneMuons"), #"globalMuons"),#"standAloneSETMuons"),
#    trajectoryInput = cms.untracked.string('cosmicMuons'),

    CSCRecSegmentLabel = cms.untracked.InputTag("csc2DRecHits"),
    RPCRecSegmentLabel = cms.untracked.InputTag("rpcRecHits"),
    DTRecSegmentLabel = cms.untracked.InputTag("dt1DRecHits"),
    cscSegments = cms.untracked.InputTag("cscSegments"),
    dt4DSegments = cms.untracked.InputTag("dt4DSegments"),
    TrackTransformer = cms.PSet(
          DoPredictionsOnly = cms.bool(False),
          Fitter = cms.string('KFFitterForRefitInsideOut'),
          TrackerRecHitBuilder = cms.string('WithTrackAngle'),
          Smoother = cms.string('KFSmootherForRefitInsideOut'),
          MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
          RefitDirection = cms.string('alongMomentum'),
          RefitRPCHits = cms.bool(False),
          Propagator = cms.string('SmartPropagatorAnyRKOpposite')
    )

)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root')
)

process.MessageLogger = cms.Service("MessageLogger")

#process.glbMuons = cms.Sequence(process.cosmicMuons*process.muonAnalyzer)
process.p = cms.Path(process.muonAnalyzer)
#process.outpath = cms.EndPath(process.out)


