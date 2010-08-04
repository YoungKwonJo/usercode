import FWCore.ParameterSet.Config as cms

##process = cms.Process("RPCPointProducer")
process = cms.Process("OWNPARTICLES")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_36X_V7::All"

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'/store/data/Commissioning10/Cosmics/RECO/v3/000/127/155/005F9301-2E16-DF11-B60B-0030487CD6B4.root'
        "rfio:/castor/cern.ch/user/s/stupputi/Beams/Data/SecondEra/PromptReco140127140330/rpcSkim_1_1_2By.root"

        )                           
)


process.load("RecoLocalMuon.RPCRecHit.rpcPointProducer_cfi")
#process.rpcPointProducer.tracks = cms.InputTag("cosmicMuons") # for cosmicMuons
process.rpcPointProducer.tracks = cms.InputTag("standAloneMuonsNoRPC") # for collision

process.out = cms.OutputModule("PoolOutputModule",
  outputCommands = cms.untracked.vstring('drop *',
        'keep *_dt4DSegments_*_*',
        'keep *_cscSegments_*_*',
        'keep *_rpcPointProducer_*_*',
        'keep *_rpcRecHits_*_*',
        'keep *_standAloneMuons_*_*',
        'keep *_cosmicMuons_*_*',
        'keep *_globalMuons_*_*',
        'keep *_*NoRPC_*_*',
        'keep L1MuRegionalCand*_*_*_*',
        'keep *_simMuonRPCDigis_*_*',
        'keep L1GlobalTriggerObjectMapRecord*_*_*_*',
        'keep *_g4SimHits_*_*',
        'keep *_muonRPCDigis_*_*'),
 fileName = cms.untracked.string('/tmp/carrillo/outs/output.root')
)
  
process.p = cms.Path(process.rpcPointProducer)

process.e = cms.EndPath(process.out)
