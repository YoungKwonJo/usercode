import FWCore.ParameterSet.Config as cms

##process = cms.Process("RPCPointProducer")
process = cms.Process("OWNPARTICLES")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_R_36X_V7::All"#"GR10_P_V7::All"#"GR_R_35X_V6::All" 

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(



        '/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/802/F6FB8B26-8563-DF11-8488-000423D98950.root',



        )                           
)

#process.load("UserCode.Fabozzi.rpcDataSkimAndProduce_cff")
process.load("RecoLocalMuon.RPCRecHit.rpcPointProducer_cfi")
#process.load("UserCode.Fabozzi.rpcSkimComplDetTrigger_cff")
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('132440:157-132440:999','132569:169-132569:999','132596:317-132596:999','132606:92-132606:999','132599:317-132599:999')
#process.source.eventsToSkip = cms.untracked.VEventRange('132440:1-132440:MAX','132569:1-132569:MAX','132596:1-132596:MAX','132606:1-132606:MAX','132599:1-132599:MAX')
process.out = cms.OutputModule("PoolOutputModule",
  outputCommands = cms.untracked.vstring('drop *',
        'keep *_gtDigis_*_*',
        'keep *_dt4DSegments_*_*',
        'keep *_cscSegments_*_*',
        'keep *_rpcPointProducer_*_*',
        'keep *_rpcRecHits_*_*',
        'keep *_standAloneMuonsNoRPC_*_*',
        'keep *_standAloneMuons_*_*',
        'keep *_cosmicMuons_*_*',
        'keep *_globalMuons_*_*',
        'keep *_*NoRPC_*_*',
        'keep L1MuRegionalCand*_*_*_*',
        'keep *_simMuonRPCDigis_*_*',
        'keep L1GlobalTriggerObjectMapRecord*_*_*_*',
        'keep *_g4SimHits_*_*',
        'keep *_muonRPCDigis_*_*'),
#  fileName = cms.untracked.string('/tmp/carrillo/outs/output.root')
 fileName = cms.untracked.string('output.root')
)

process.options = cms.untracked.PSet(
  SkipEvent = cms.untracked.vstring('ProductNotFound')
)
  
process.p = cms.Path(process.rpcPointProducer)

process.e = cms.EndPath(process.out)
