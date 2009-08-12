import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCEff")

# Conditions (Global Tag is used here):
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_31X_GLOBALTAG"
process.GlobalTag.globaltag = "GR09_31X_V3P::All"
process.prefer("GlobalTag")

#Geometry
process.load("Configuration.StandardSequences.Geometry_cff")

# reconstruction sequence for Cosmics
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")

##------------ Cosmics ----------------------------------------------------
process.load("RecoLocalMuon.Configuration.RecoLocalMuonCosmics_cff")
process.load("RecoMuon.Configuration.RecoMuonCosmics_cff")
process.load("RecoMuon.MuonSeedGenerator.CosmicMuonSeedProducer_cfi")
process.load("RecoMuon.CosmicMuonProducer.cosmicMuons_cfi")
process.cosmicMuons.TrajectoryBuilderParameters.EnableRPCMeasurement = False
##--------------------------------------------------------------------------

process.Timing = cms.Service("Timing")

process.load("DQMServices.Core.DQM_cfg")

process.source = cms.Source ("PoolSource",
                             fileNames = cms.untracked.vstring(
                                       'file:/afs/cern.ch/user/y/youngjo/public/reco.root'
                                                               )
                             
                             )
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)


process.ModuleWebRegistry = cms.Service("ModuleWebRegistry")

process.LockService = cms.Service("LockService",
    labels = cms.untracked.vstring('source')
)

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.rpcefftrack = cms.EDFilter("RPCEffTrackExtrapolationNew",
                                   MuonServiceProxy,
                                   dt4DSegments = cms.InputTag("dt4DSegments"),
                                   cscSegments =  cms.InputTag("cscSegments"),
                                   RPCRecHits = cms.InputTag("rpcRecHits"),
                                   NavigationType = cms.string("Direct"),
                                   trajectoryInput = cms.untracked.string('cosmicMuons'),
                                   ##trajectoryInput = cms.untracked.string('standAloneMuons'),
                                   EffRootFileName = cms.untracked.string('RPCEffTrackExtrapolationNew.root'),
                                   PropagatorName = cms.string('SteppingHelixPropagatorAny'),
                                   TrackTransformer = cms.PSet(TrackerRecHitBuilder = cms.string('WithTrackAngle'), # added for trajectory of cosmicmuon
                                                               MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
                                                               RefitRPCHits = cms.bool(False)
                                                               )
                                   
                                   )

##process.DQM.collectorHost = 'myhost'
##process.DQM.collectorPort = 9090
##process.DQM.debug = False

##process.p1 = cms.Path(process.rpcefftrack)
process.p1 = cms.Path(process.cosmicMuons*process.rpcefftrack)



