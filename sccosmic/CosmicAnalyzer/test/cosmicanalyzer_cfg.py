# The following comments couldn't be translated into the new config version:
#keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
# Conditions (Global Tag is used here):
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.connect = "frontier://PromptProd/CMS_COND_21X_GLOBALTAG"
#process.GlobalTag.globaltag = "CRAFT_V4P::All"

process.GlobalTag.globaltag = 'COSMMC_22X_TK::All'
#process.GlobalTag.globaltag = "CRZT210_V3P::All"
process.prefer("GlobalTag")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.MagneticField_0T_cff")
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.MagneticField_40T_cff")

process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")

process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")

process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('sccosmicana.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(


       "/store/data/Commissioning08/Cosmics/RECO/v1/000/068/021/CE84A0AC-2AA6-DD11-91BC-000423D98BE8.root"


     )
)

process.muonAnalyzer = cms.EDAnalyzer("CosmicAnalyzer",
    process.MuonServiceProxy,
    TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dREcal = cms.double(9999.0),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
        crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        accountForTrajectoryChangeCalo = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(True),
        useCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
    ),
    Tracks3 = cms.untracked.InputTag("cosmicMuonsNoDriftBarrelOnly"), #"cosmicMuonsNoRPC"),"cosmicMuonsNoDriftBarrelOnly"
    #Tracks3 = cms.untracked.InputTag("cosmicMuonsNoRPC"), #"cosmicMuonsNoDriftBarrelOnly"
    Tracks2 = cms.untracked.InputTag("cosmicMuonsEndCapsOnly"),
    TrackAssociatorParameterBlock = cms.PSet(
        TrackAssociatorParameters = cms.PSet(
            muonMaxDistanceSigmaX = cms.double(0.0),
            muonMaxDistanceSigmaY = cms.double(0.0),
            CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
            dRHcal = cms.double(9999.0),
            dREcal = cms.double(9999.0),
            CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
            useEcal = cms.bool(True),
            dREcalPreselection = cms.double(0.05),
            HORecHitCollectionLabel = cms.InputTag("horeco"),
            dRMuon = cms.double(9999.0),
            crossedEnergyType = cms.string('SinglePointAlongTrajectory'),
            propagateAllDirections = cms.bool(True),
            muonMaxDistanceX = cms.double(5.0),
            muonMaxDistanceY = cms.double(5.0),
            useHO = cms.bool(True),
            accountForTrajectoryChangeCalo = cms.bool(False),
            DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
            EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
            dRHcalPreselection = cms.double(0.2),
            useMuon = cms.bool(True),
            useCalo = cms.bool(False),
            EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
            dRMuonPreselection = cms.double(0.2),
            truthMatch = cms.bool(False),
            HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
            useHcal = cms.bool(True)
        )
    ),
    CSCRecSegmentLabel = cms.untracked.InputTag("csc2DRecHits"),
    Tracks = cms.untracked.InputTag("cosmicMuonsBarrelOnly"),
    Tracks4 = cms.untracked.InputTag("muonDTDigis"),
    RPCRecSegmentLabel = cms.untracked.InputTag("rpcRecHits"),
    debug = cms.untracked.bool(True),
    Muons = cms.untracked.InputTag("GLBMuons"), #"muons"), "GLBMuons"
    #Muons = cms.untracked.InputTag("muons"), #"muons"), "GLBMuons"
    STAMuons = cms.untracked.InputTag("STAMuons"), # "muonsNoRPC", "STAMuons"
    #STAMuons = cms.untracked.InputTag("muonsNoRPC"), # "muonsNoRPC", "STAMuons"
    DTRecSegmentLabel = cms.untracked.InputTag("dt1DRecHits")
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('output.root')
)

process.MessageLogger = cms.Service("MessageLogger")

process.p = cms.Path(process.muonAnalyzer)
#process.outpath = cms.EndPath(process.out)
process.muonAnalyzer.TrackAssociatorParameters.dRMuonPreselection = 0.6
process.muonAnalyzer.TrackAssociatorParameters.muonMaxDistanceX = 3.0
process.muonAnalyzer.TrackAssociatorParameters.muonMaxDistanceY = 3.0


