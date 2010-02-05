import FWCore.ParameterSet.Config as cms

##process = cms.Process("RPCPointProducer")
process = cms.Process("OWNPARTICLES")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")


process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR09_P_V8::All"
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       "/store/data/BeamCommissioning09/PhysicsMuonBkg/RECO/v2/000/124/056/FC1CCF05-B3E7-DE11-B062-001617E30CC8.root"
    )
)

process.rpcPointProducer = cms.EDProducer('RPCPointProducer',
  incldt = cms.untracked.bool(True),
  inclcsc = cms.untracked.bool(True),
  incltrack =  cms.untracked.bool(True),

  debug = cms.untracked.bool(True),

  rangestrips = cms.untracked.double(4.),
  rangestripsRB4 = cms.untracked.double(4.),
  MinCosAng = cms.untracked.double(0.85),
  MaxD = cms.untracked.double(80.0),
  MaxDrb4 = cms.untracked.double(150.0),
  ExtrapolatedRegion = cms.untracked.double(0.6), #in stripl/2 in Y and stripw*nstrips/2 in X

  cscSegments = cms.untracked.string('cscSegments'),
  dt4DSegments = cms.untracked.string('dt4DSegments'),
  tracks = cms.untracked.string('standAloneMuons'),
  TrackTransformer = cms.PSet(
      DoPredictionsOnly = cms.bool(False),
      Fitter = cms.string('KFFitterForRefitInsideOut'),
      TrackerRecHitBuilder = cms.string('WithTrackAngle'),
      Smoother = cms.string('KFSmootherForRefitInsideOut'),
      MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
      RefitDirection = cms.string('alongMomentum'),
      RefitRPCHits = cms.bool(False),
      Propagator = cms.string('SmartPropagatorAnyRKOpposite')
  ),
# cscSegments = cms.untracked.string('hltCscSegments'),
# dt4DSegments = cms.untracked.string('hltDt4DSegments'),
)

process.out = cms.OutputModule("PoolOutputModule",
  outputCommands = cms.untracked.vstring('drop *',
        'keep *_dt4DSegments_*_*',
        'keep *_cscSegments_*_*',
        'keep *_rpcPointProducer_*_*',
        'keep *_rpcRecHits_*_*'),
  fileName = cms.untracked.string('/tmp/youngjo/output.root')
)
  
process.p = cms.Path(process.rpcPointProducer)

process.e = cms.EndPath(process.out)
