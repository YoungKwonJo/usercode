import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCEff")

# Conditions (Global Tag is used here):
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_31X_GLOBALTAG"
process.GlobalTag.globaltag = "GR09_P_V6::All"
process.prefer("GlobalTag")

#Geometry
process.load("Configuration.StandardSequences.Geometry_cff")

# reconstruction sequence for Cosmics
process.load("Configuration.StandardSequences.ReconstructionCosmics_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")

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
                                       '/store/data/BeamCommissioning09/MinimumBias/RECO/v2/000/123/977/34CFFC67-0DE6-DE11-B867-003048D37456.root'
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
                                   trajectoryInput = cms.untracked.string('standAloneMuons'),
                                   ##trajectoryInput = cms.untracked.string('standAloneMuons'),
                                   EffRootFileName = cms.untracked.string('RPCEffTrackExtrapolationNew.root'),
                                   PropagatorName = cms.string('SteppingHelixPropagatorAny'),
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

##process.DQM.collectorHost = 'myhost'
##process.DQM.collectorPort = 9090
##process.DQM.debug = False

##process.p1 = cms.Path(process.rpcefftrack)
process.p1 = cms.Path(process.cosmicMuons*process.rpcefftrack)
