// -*- C++ -*-
//
// Package:    CosmicAnalyzer
// Class:      CosmicAnalyzer
//
/**\class CosmicAnalyzer CosmicAnalyzer.cc sccosmic/CosmicAnalyzer/src/CosmicAnalyzer.cc

 Description: <one line class summary>

// This code is run with : Tracking Tools(http://higgs.skku.ac.kr/CMS/TrackingTools_CMSSW_2_1_10.tgz) 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Su Yong Choi
//         Created:  Wed Jul  9 17:03:35 CEST 2008
// $Id$
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"

#include <DataFormats/MuonDetId/interface/DTChamberId.h>
#include <DataFormats/DTDigi/interface/DTLocalTrigger.h>
#include <DataFormats/MuonData/interface/MuonDigiCollection.h>

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TAMuonChamberMatch.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfoCollection.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
//
// class decleration
//

using namespace edm;
using namespace std;
typedef MuonDigiCollection<DTChamberId, DTLocalTrigger> DTLocalTriggerCollection;

bool findmatch(const  edm::Handle<DTRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual);
bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum);
bool findmatch(const edm::Handle<CSCRecHit2DCollection> &hitcoll, int detid, float locx, float locy, float &residual);

double maxdist = 7777777.0;

edm::ESHandle<GlobalTrackingGeometry> theG;

struct dethitstruct {
    int  count;
    int  subdet;
    double locx;
    double locy;
    double locz;
    double glbx;
    double glby;
    double glbz;
};

class CosmicAnalyzer : public edm::EDAnalyzer
{
    public:
        explicit CosmicAnalyzer(const edm::ParameterSet&);
        ~CosmicAnalyzer();

    private:
        virtual void beginJob(const edm::EventSetup&) ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        edm::InputTag MuonTags_;
        edm::InputTag STAMuonTags_;
        edm::InputTag TrackTags_;
        edm::InputTag TrackTags2_;
        edm::InputTag TrackTags3_;
        edm::InputTag TrackTags4_;

// ----------member data ---------------------------
        int icn; 
        TH1F *hnstamuon;
        TH1F *hstamuon_eta;
        TH2F *hstamuon_etaphi;
        TH1F *hstamuon_pt;

        TH1F *hnglbmuon;
        TH1F *hglbmuon_eta;
        TH2F *hglbmuon_etaphi;
        TH1F *hglbmuon_pt;

        TH1F *hntracks;
        TH2F *htracks_etaphi;
        TH2F *hnvalidhits;
        TH2F *hchi2prob;
        TH2F *hchi2ndof;
        TH2F *hqoverppull;
        TH2F *htrackinnerxy;
        TH2F *htrackouterxy;
        TH1F *htrackinnerz;
        TH1F *htrackouterz;
        TH2F *htrackq;

        TH2F *htrackrechitsxy;
        TH1F *htrackrechitsz;

        TH2F *htrackrechitsDTxy;
        TH1F *htrackrechitsDTz;

        TH2F *htrackrechitsRPCxy;
        TH1F *htrackrechitsRPCz;

        TH2F *htrackrechitsCSCxy;
        TH1F *htrackrechitsCSCz;

        TH1F *hntracks_2;
        TH2F *htracks_etaphi_2;
        TH2F *hnvalidhits_2;
        TH2F *hchi2prob_2;
        TH2F *hchi2ndof_2;
        TH2F *hqoverppull_2;
        TH2F *htrackinnerxy_2;
        TH2F *htrackouterxy_2;
        TH1F *htrackinnerz_2;
        TH1F *htrackouterz_2;
        TH2F *htrackq_2;

        TH2F *htrackrechitsxy_2;
        TH1F *htrackrechitsz_2;

        TH2F *htrackrechitsDTxy_2;
        TH1F *htrackrechitsDTz_2;

        TH2F *htrackrechitsRPCxy_2;
        TH1F *htrackrechitsRPCz_2;

        TH2F *htrackrechitsCSCxy_2;
        TH1F *htrackrechitsCSCz_2;

        TH1F *hntracks_3;
        TH2F *htracks_etaphi_3;
        TH2F *hnvalidhits_3;
        TH2F *hchi2prob_3;
        TH2F *hchi2ndof_3;
        TH2F *hqoverppull_3;
        TH2F *htrackinnerxy_3;
        TH2F *htrackouterxy_3;
        TH1F *htrackinnerz_3;
        TH1F *htrackouterz_3;
        TH2F *htrackq_3;

        TH2F *hntrackcorr;

        TH2F *hdttriggerhitxy[5];
        TH1F *hdttriggerhitz;

        TH2F *hdttrajectoryxy[5];
        TH1F *hdttrajectoryz;
        TH2F *hcsctrajectoryxy_e1[4];
        TH2F *hcsctrajectoryxy_e2[4];
        TH1F *hcsctrajectoryz;
        TH2F *hrpctrajectoryxy_r0[5];
        TH1F *hrpctrajectoryz;
        TH2F *hrpctrajectoryxy_rm1[3];
        TH2F *hrpctrajectoryxy_rp1[3];

        TH2F *hdtmatchrechitxy[5];
        TH1F *hdtmatchrechitz;
        TH2F *hcscmatchrechitxy_e1[4];
        TH2F *hcscmatchrechitxy_e2[4];
        TH1F *hcscmatchrechitz;
        TH2F *hrpcmatchrechitxy_r0[5];
        TH1F *hrpcmatchrechitz;
        TH2F *hrpcmatchrechitxy_rm1[3];
        TH2F *hrpcmatchrechitxy_rp1[3];

        TH2F *hdtmatchrechitontrackxy[5];
        TH1F *hdtmatchrechitontrackz;
        TH2F *hcscmatchrechitontrackxy_e1[4];
        TH2F *hcscmatchrechitontrackxy_e2[4];
        TH1F *hcscmatchrechitontrackz;
        TH2F *hrpcmatchrechitontrackxy_r0[5];
        TH1F *hrpcmatchrechitontrackz;
        TH2F *hrpcmatchrechitontrackxy_rm1[3];
        TH2F *hrpcmatchrechitontrackxy_rp1[3];

        MuonServiceProxy *theMuonService_;
        edm::ESHandle<Propagator> thePropagator;

/// the name of the DT rec hits collection
        edm::InputTag theDTRecSegmentLabel;
/// the name of the CSC rec hits collection
        edm::InputTag theCSCRecSegmentLabel;
/// the name of the RPC rec hits collection
        edm::InputTag theRPCRecSegmentLabel;

        TH1F *hntracks_4;

        TrackDetectorAssociator trackAssociator_;
        TrackAssociatorParameters parameters_;
        Chi2MeasurementEstimator trchi2est;
        MuonDetLayerMeasurements muonMeasurements;
        std::map<int, int> allmuonhits;
        std::map<int, int> matchedmuonhitsontrack;
        std::map<int, int> matchedmuonhitsondet;
        bool _debug;

        std::map<int, TH2F*> detLptMap_;
        std::map<int, TH2F*> detLpmMap_;
        std::map<int, TH2F*> detLprMap_;

        std::map<int, TH1F*> detMapResidureX_;
        std::map<int, TH1F*> detMapResidureY_;
        std::map<int, TH1F*> detMapResidureZ_;
        std::map<int, TH1F*> detMapResidureXadd;

        std::map<int, TH1F*> detMapResidureGX_;
        std::map<int, TH1F*> detMapResidureGY_;
        std::map<int, TH1F*> detMapResidureGZ_;

       int dtlsize;
       int csclsize;
       int rpclsize;

    edm::Service<TFileService> fs;

    //std::map<int, dethitstruct> rpcrechits_;

};


//
// constructors and destructor
//
CosmicAnalyzer::CosmicAnalyzer(const edm::ParameterSet& iConfig)
:MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
STAMuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("STAMuons")),
TrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks")),
TrackTags2_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks2")),
TrackTags3_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks3")),
TrackTags4_(iConfig.getUntrackedParameter<edm::InputTag>("Tracks4")),
_debug(iConfig.getUntrackedParameter<bool>("debug")),
theDTRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("DTRecSegmentLabel")),
theCSCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("CSCRecSegmentLabel")),
theRPCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("RPCRecSegmentLabel")),
trchi2est(100,3.0),
muonMeasurements(theDTRecSegmentLabel,theCSCRecSegmentLabel, theRPCRecSegmentLabel, true, true, true)
{
//now do what ever initialization is needed
    icn = 1; 
    dtlsize=140;
    csclsize=110;
    rpclsize=200;

    edm::ParameterSet serviceParameters
        = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
    theMuonService_ = new MuonServiceProxy(serviceParameters);

// TrackAssociator parameters
    edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
    parameters_.loadParameters( parameters );
    trackAssociator_.useDefaultPropagator();

    edm::Service<TFileService> fs;

    hnstamuon = fs->make<TH1F>("hnsta", "Number of standalone muons", 20, 0.0, 20.0);
    hstamuon_eta = fs->make<TH1F>("hsta_eta", "|#eta| standalone muons", 50, -2.5, 2.5);
    hstamuon_etaphi = fs->make<TH2F>("hsta_etaphi", "|#eta| #phi standalone muons", 50, -2.5, 2.5, 50, -3.1415, 3.14158);
    hstamuon_pt = fs->make<TH1F>("hsta_pt", "Pt standalone muons", 600, 0, 300);

    hnglbmuon = fs->make<TH1F>("hnglb", "Number of Global muons", 20, 0.0, 20.0);
    hglbmuon_eta = fs->make<TH1F>("hglb_eta", "|#eta| global muons", 50, -2.5, 2.5);
    hglbmuon_etaphi = fs->make<TH2F>("hglb_etaphi", "|#eta| #phi global muons", 50, -2.5, 2.5, 50, -3.1415, 3.14158);
    hglbmuon_pt = fs->make<TH1F>("hglb_pt", "Pt global muons", 600, 0, 300);

    hntracks = fs->make<TH1F>("hntracks", "Number of tracks", 20, 0.0, 20.0);
    htracks_etaphi = fs->make<TH2F>("htracks_etaphi", "|#eta| #phi tracks", 50, -2.5, 2.5, 50, -3.1415, 3.14158);
    htrackinnerxy = fs->make<TH2F>("htrackinnerxy", "y vs x of inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxy = fs->make<TH2F>("htrackouterxy", "y vs x of outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackinnerz = fs->make<TH1F>("htrackinnerz", "z of inner hit ", 100, -1500.0, 1500.0);
    htrackouterz = fs->make<TH1F>("htrackouterz", "z of outer hit ", 100, -1500.0, 1500.0);
    hnvalidhits = fs->make<TH2F>("hnvalidhits", "Number of hits on track", 50, -3.1415, 3.1415, 50, 0.0, 50.0);
    hchi2prob = fs->make<TH2F>("hchi2prob", "#chi^2 probability of tracks", 50, -3.1415, 3.1415, 50, 0.0, 1.0);
    hchi2ndof = fs->make<TH2F>("hchi2ndof", "Normalized #chi^2 of tracks", 50, -3.1415, 3.1415, 50, 0.0, 10.0);
    hqoverppull = fs->make<TH2F>("hqoverppull", "qoverp pull of tracks", 50, -3.1415, 3.1415, 50, -1.0, 1.0);
    htrackq = fs->make<TH2F>("htrackq", "q of tracks", 50, -3.1415, 3.1415, 50, -1.5, 1.5);

    htrackrechitsxy = fs->make<TH2F>("htrackrechitsxy", "y vs x of rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsz = fs->make<TH1F>("htrackrechitsz", "z of inner hit ", 500, -1500.0, 1500.0);
    htrackrechitsDTxy = fs->make<TH2F>("htrackrechitsDTxy", "y vs x of DT rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsDTz = fs->make<TH1F>("htrackrechitsDTz", "z of DT rechits ", 500, -1500.0, 1500.0);
    htrackrechitsRPCxy = fs->make<TH2F>("htrackrechitsRPCxy", "y vs x of RPC rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsRPCz = fs->make<TH1F>("htrackrechitsRPCz", "z of RPC rechits ", 500, -1500.0, 1500.0);
    htrackrechitsCSCxy = fs->make<TH2F>("htrackrechitsCSCxy", "y vs x of CSC rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsCSCz = fs->make<TH1F>("htrackrechitsCSCz", "z of CSC rechits ", 500, -1500.0, 1500.0);

    hntracks_2 = fs->make<TH1F>("hntracks_2", "Number of tracks", 20, 0.0, 20.0);
    htracks_etaphi_2 = fs->make<TH2F>("htracks_etaphi_2", "|#eta| #phi tracks", 50, -2.5, 2.5, 50, -3.1415, 3.14158);
    htrackinnerxy_2 = fs->make<TH2F>("htrackinnerxy_2", "y vs x of inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxy_2 = fs->make<TH2F>("htrackouterxy_2", "y vs x of outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackinnerz_2 = fs->make<TH1F>("htrackinnerz_2", "z of inner hit ", 100, -1500.0, 1500.0);
    htrackouterz_2 = fs->make<TH1F>("htrackouterz_2", "z of outer hit ", 100, -1500.0, 1500.0);
    hnvalidhits_2 = fs->make<TH2F>("hnvalidhits_2", "Number of hits on track", 50, -3.1415, 3.1415, 50, 0.0, 50.0);
    hchi2prob_2 = fs->make<TH2F>("hchi2prob_2", "#chi^2 probability of tracks", 50, -3.1415, 3.1415, 50, 0.0, 1.0);
    hchi2ndof_2 = fs->make<TH2F>("hchi2ndof_2", "Normalized #chi^2 of tracks", 50, -3.1415, 3.1415, 50, 0.0, 10.0);
    hqoverppull_2 = fs->make<TH2F>("hqoverppull_2", "qoverp pull of tracks", 50, -3.1415, 3.1415, 50, -1.0, 1.0);
    htrackq_2 = fs->make<TH2F>("htrackq_2", "q of tracks", 50, -3.1415, 3.1415, 50, -1.5, 1.5);

    htrackrechitsxy_2 = fs->make<TH2F>("htrackrechitsxy_2", "y vs x of rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsz_2 = fs->make<TH1F>("htrackrechitsz_2", "z of inner hit ", 500, -1500.0, 1500.0);
    htrackrechitsDTxy_2 = fs->make<TH2F>("htrackrechitsDTxy_2", "y vs x of DT rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsDTz_2 = fs->make<TH1F>("htrackrechitsDTz_2", "z of DT rechits ", 500, -1500.0, 1500.0);
    htrackrechitsRPCxy_2 = fs->make<TH2F>("htrackrechitsRPCxy_2", "y vs x of RPC rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsRPCz_2 = fs->make<TH1F>("htrackrechitsRPCz_2", "z of RPC rechits ", 500, -1500.0, 1500.0);
    htrackrechitsCSCxy_2 = fs->make<TH2F>("htrackrechitsCSCxy_2", "y vs x of CSC rechits", 250, -1000.0, 1000.0, 250, -1000.0, 1000.0);
    htrackrechitsCSCz_2 = fs->make<TH1F>("htrackrechitsCSCz_2", "z of CSC rechits ", 500, -1500.0, 1500.0);

    hntracks_3 = fs->make<TH1F>("hntracks_3", "Number of tracks", 20, 0.0, 20.0);
    htracks_etaphi_3 = fs->make<TH2F>("htracks_etaphi_3", "|#eta| #phi tracks", 50, -2.5, 2.5, 50, -3.1415, 3.14158);
    htrackinnerxy_3 = fs->make<TH2F>("htrackinnerxy_3", "y vs x of inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxy_3 = fs->make<TH2F>("htrackouterxy_3", "y vs x of outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackinnerz_3 = fs->make<TH1F>("htrackinnerz_3", "z of inner hit ", 100, -1500.0, 1500.0);
    htrackouterz_3 = fs->make<TH1F>("htrackouterz_3", "z of outer hit ", 100, -1500.0, 1500.0);
    hnvalidhits_3 = fs->make<TH2F>("hnvalidhits_3", "Number of hits on track", 50, -3.1415, 3.1415, 50, 0.0, 100.0);
    hchi2prob_3 = fs->make<TH2F>("hchi2prob_3", "#chi^2 probability of tracks", 50, -3.1415, 3.1415, 50, 0.0, 1.0);
    hchi2ndof_3 = fs->make<TH2F>("hchi2ndof_3", "Normalized #chi^2 of tracks", 50, -3.1415, 3.1415, 50, 0.0, 10.0);
    hqoverppull_3 = fs->make<TH2F>("hqoverppull_3", "qoverp pull of tracks", 50, -3.1415, 3.1415, 50, -1.0, 1.0);
    htrackq_3 = fs->make<TH2F>("htrackq_3", "q of tracks", 50, -3.1415, 3.1415, 50, -1.5, 1.5);

    hntrackcorr = fs->make<TH2F>("hntrackcorr", "track correlations", 10, 0.0, 10.0, 10, 0.0, 10.0);
   
    int nbins = 1000;
// ///////////////////////////
    hdttriggerhitxy[0] = fs->make<TH2F>("hdttriggerhitxy_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttriggerhitxy[1] = fs->make<TH2F>("hdttriggerhitxy_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttriggerhitxy[2] = fs->make<TH2F>("hdttriggerhitxy_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttriggerhitxy[3] = fs->make<TH2F>("hdttriggerhitxy_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttriggerhitxy[4] = fs->make<TH2F>("hdttriggerhitxy_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttriggerhitz = fs->make<TH1F>("hdttriggerhitz", "z demesion", nbins, -1000.0, 1000.0);

    // DT chambers

    hdttrajectoryxy[0] = fs->make<TH2F>("hdttrajectory_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttrajectoryxy[1] = fs->make<TH2F>("hdttrajectory_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttrajectoryxy[2] = fs->make<TH2F>("hdttrajectory_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttrajectoryxy[3] = fs->make<TH2F>("hdttrajectory_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttrajectoryxy[4] = fs->make<TH2F>("hdttrajectory_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdttrajectoryz = fs->make<TH1F>("hdttrajectoryz", "z demesion", nbins, -1000.0, 1000.0);

    // CSC chambers
    hcsctrajectoryxy_e1[0] = fs->make<TH2F>("hcsctrajectory_e1_s1", "y vs x endcap 1 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e1[1] = fs->make<TH2F>("hcsctrajectory_e1_s2", "y vs x endcap 1 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e1[2] = fs->make<TH2F>("hcsctrajectory_e1_s3", "y vs x endcap 1 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e1[3] = fs->make<TH2F>("hcsctrajectory_e1_s4", "y vs x endcap 1 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e2[0] = fs->make<TH2F>("hcsctrajectory_e2_s1", "y vs x endcap 2 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e2[1] = fs->make<TH2F>("hcsctrajectory_e2_s2", "y vs x endcap 2 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e2[2] = fs->make<TH2F>("hcsctrajectory_e2_s3", "y vs x endcap 2 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcsctrajectoryxy_e2[3] = fs->make<TH2F>("hcsctrajectory_e2_s4", "y vs x endcap 2 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    hcsctrajectoryz = fs->make<TH1F>("hcsctrajectoryz", "z demesion", nbins, -1500.0, 1500.0);

    // RPC chambers
    hrpctrajectoryxy_r0[0] = fs->make<TH2F>("hrpctrajectory_r0_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_r0[1] = fs->make<TH2F>("hrpctrajectory_r0_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_r0[2] = fs->make<TH2F>("hrpctrajectory_r0_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_r0[3] = fs->make<TH2F>("hrpctrajectory_r0_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_r0[4] = fs->make<TH2F>("hrpctrajectory_r0_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryz = fs->make<TH1F>("hrpctrajectoryz", "z demesion", nbins, -1500.0, 1500.0);

    hrpctrajectoryxy_rm1[0] = fs->make<TH2F>("hrpctrajectory_rm1_1", "y vs x endcap -1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_rm1[1] = fs->make<TH2F>("hrpctrajectory_rm1_2", "y vs x endcap -1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_rm1[2] = fs->make<TH2F>("hrpctrajectory_rm1_3", "y vs x endcap -1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    hrpctrajectoryxy_rp1[0] = fs->make<TH2F>("hrpctrajectory_rp1_1", "y vs x endcap 1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_rp1[1] = fs->make<TH2F>("hrpctrajectory_rp1_2", "y vs x endcap 1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpctrajectoryxy_rp1[2] = fs->make<TH2F>("hrpctrajectory_rp1_3", "y vs x endcap 1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    // DT chambers
    hdtmatchrechitxy[0] = fs->make<TH2F>("hdtmatchrechit_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitxy[1] = fs->make<TH2F>("hdtmatchrechit_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitxy[2] = fs->make<TH2F>("hdtmatchrechit_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitxy[3] = fs->make<TH2F>("hdtmatchrechit_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitxy[4] = fs->make<TH2F>("hdtmatchrechit_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitz = fs->make<TH1F>("hdtmatchrechitz", "z demesion", nbins, -1000.0, 1000.0);

    // CSC chambers
    hcscmatchrechitxy_e1[0] = fs->make<TH2F>("hcscmatchrechit_e1_s1", "y vs x endcap 1 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e1[1] = fs->make<TH2F>("hcscmatchrechit_e1_s2", "y vs x endcap 1 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e1[2] = fs->make<TH2F>("hcscmatchrechit_e1_s3", "y vs x endcap 1 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e1[3] = fs->make<TH2F>("hcscmatchrechit_e1_s4", "y vs x endcap 1 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e2[0] = fs->make<TH2F>("hcscmatchrechit_e2_s1", "y vs x endcap 2 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e2[1] = fs->make<TH2F>("hcscmatchrechit_e2_s2", "y vs x endcap 2 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e2[2] = fs->make<TH2F>("hcscmatchrechit_e2_s3", "y vs x endcap 2 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitxy_e2[3] = fs->make<TH2F>("hcscmatchrechit_e2_s4", "y vs x endcap 2 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    hcscmatchrechitz = fs->make<TH1F>("hcscmatchrechitz", "z demesion", nbins, -1500.0, 1500.0);

    // RPC chambers
    hrpcmatchrechitxy_r0[0] = fs->make<TH2F>("hrpcmatchrechit_r0_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_r0[1] = fs->make<TH2F>("hrpcmatchrechit_r0_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_r0[2] = fs->make<TH2F>("hrpcmatchrechit_r0_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_r0[3] = fs->make<TH2F>("hrpcmatchrechit_r0_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_r0[4] = fs->make<TH2F>("hrpcmatchrechit_r0_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitz = fs->make<TH1F>("hrpcmatchrechitz", "z demesion", nbins, -1500.0, 1500.0);

    hrpcmatchrechitxy_rm1[0] = fs->make<TH2F>("hrpcmatchrechit_rm1_1", "y vs x endcap -1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_rm1[1] = fs->make<TH2F>("hrpcmatchrechit_rm1_2", "y vs x endcap -1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_rm1[2] = fs->make<TH2F>("hrpcmatchrechit_rm1_3", "y vs x endcap -1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    hrpcmatchrechitxy_rp1[0] = fs->make<TH2F>("hrpcmatchrechit_rp1_1", "y vs x endcap 1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_rp1[1] = fs->make<TH2F>("hrpcmatchrechit_rp1_2", "y vs x endcap 1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitxy_rp1[2] = fs->make<TH2F>("hrpcmatchrechit_rp1_3", "y vs x endcap 1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

 // DT chambers
    hdtmatchrechitontrackxy[0] = fs->make<TH2F>("hdtmatchrechitontrack_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitontrackxy[1] = fs->make<TH2F>("hdtmatchrechitontrack_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitontrackxy[2] = fs->make<TH2F>("hdtmatchrechitontrack_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitontrackxy[3] = fs->make<TH2F>("hdtmatchrechitontrack_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitontrackxy[4] = fs->make<TH2F>("hdtmatchrechitontrack_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hdtmatchrechitontrackz = fs->make<TH1F>("hdtmatchrechitontrackz", "z demesion", nbins, -1000.0, 1000.0);

    // CSC chambers
    hcscmatchrechitontrackxy_e1[0] = fs->make<TH2F>("hcscmatchrechitontrack_e1_s1", "y vs x endcap 1 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e1[1] = fs->make<TH2F>("hcscmatchrechitontrack_e1_s2", "y vs x endcap 1 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e1[2] = fs->make<TH2F>("hcscmatchrechitontrack_e1_s3", "y vs x endcap 1 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e1[3] = fs->make<TH2F>("hcscmatchrechitontrack_e1_s4", "y vs x endcap 1 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e2[0] = fs->make<TH2F>("hcscmatchrechitontrack_e2_s1", "y vs x endcap 2 station 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e2[1] = fs->make<TH2F>("hcscmatchrechitontrack_e2_s2", "y vs x endcap 2 station 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e2[2] = fs->make<TH2F>("hcscmatchrechitontrack_e2_s3", "y vs x endcap 2 station 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackxy_e2[3] = fs->make<TH2F>("hcscmatchrechitontrack_e2_s4", "y vs x endcap 2 station 4", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hcscmatchrechitontrackz = fs->make<TH1F>("hcscmatchrechitontrackz", "z demesion", nbins, -1500.0, 1500.0);

    // RPC chambers
    hrpcmatchrechitontrackxy_r0[0] = fs->make<TH2F>("hrpcmatchrechitontrack_r0_m2", "y vs x wheel -2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_r0[1] = fs->make<TH2F>("hrpcmatchrechitontrack_r0_m1", "y vs x wheel -1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_r0[2] = fs->make<TH2F>("hrpcmatchrechitontrack_r0_0", "y vs x wheel 0", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_r0[3] = fs->make<TH2F>("hrpcmatchrechitontrack_r0_p1", "y vs x wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_r0[4] = fs->make<TH2F>("hrpcmatchrechitontrack_r0_p2", "y vs x wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackz = fs->make<TH1F>("hrpcmatchrechitontrackz", "z demesion", nbins, -1500.0, 1500.0);

    hrpcmatchrechitontrackxy_rm1[0] = fs->make<TH2F>("hrpcmatchrechitontrack_rm1_1", "y vs x endcap -1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_rm1[1] = fs->make<TH2F>("hrpcmatchrechitontrack_rm1_2", "y vs x endcap -1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_rm1[2] = fs->make<TH2F>("hrpcmatchrechitontrack_rm1_3", "y vs x endcap -1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

    hrpcmatchrechitontrackxy_rp1[0] = fs->make<TH2F>("hrpcmatchrechitontrack_rp1_1", "y vs x endcap 1 wheel 1", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_rp1[1] = fs->make<TH2F>("hrpcmatchrechitontrack_rp1_2", "y vs x endcap 1 wheel 2", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);
    hrpcmatchrechitontrackxy_rp1[2] = fs->make<TH2F>("hrpcmatchrechitontrack_rp1_3", "y vs x endcap 1 wheel 3", nbins, -1000.0, 1000.0, nbins, -1000.0, 1000.0);

// ////////////////////////////////
    hntracks_4 = fs->make<TH1F>("hntracks_4", "Number of tracks", 20, 0.0, 20.0);
}


CosmicAnalyzer::~CosmicAnalyzer()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
CosmicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using reco::MuonCollection;
    using reco::TrackCollection;
    using namespace GeomDetEnumerators;

    edm::Handle<MuonCollection> MuCollection;
    edm::Handle<MuonCollection> STACollection;
    edm::Handle<TrackCollection> TrCollection;
    edm::Handle<TrackCollection> TrCollection2;
    edm::Handle<TrackCollection> TrCollection3;
    edm::Handle<DTLocalTriggerCollection> TrCollection4;
    edm::Handle<reco::MuonTrackLinksCollection> muHandle;
    edm::ESHandle<TransientTrackBuilder> ttrackBuilder;

    edm::Handle<RPCRecHitCollection> allRPChits;
    edm::Handle<DTRecHitCollection> allDThits;
    edm::Handle<CSCRecHit2DCollection> allCSChits;

// Get the magnetic field
    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

// Get the global tracking geometry
    iSetup.get<GlobalTrackingGeometryRecord>().get(theG);

    edm::ESHandle<MuonDetLayerGeometry> theMuonLayers;
    iSetup.get<MuonRecoGeometryRecord>().get(theMuonLayers);
// get the Muon layers
    vector<DetLayer*> dtLayers = theMuonLayers->allDTLayers();
    vector<DetLayer*> rpcLayers = theMuonLayers->allRPCLayers();

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",thePropagator);

    iEvent.getByLabel(MuonTags_, MuCollection);
    iEvent.getByLabel(STAMuonTags_, STACollection);
    iEvent.getByLabel(TrackTags_, TrCollection);
    iEvent.getByLabel(TrackTags2_, TrCollection2);
    iEvent.getByLabel(TrackTags3_, TrCollection3);
    iEvent.getByLabel(TrackTags4_, TrCollection4);
    iEvent.getByLabel(theRPCRecSegmentLabel, allRPChits);
    iEvent.getByLabel(theDTRecSegmentLabel, allDThits);
    iEvent.getByLabel(theCSCRecSegmentLabel, allCSChits);
//   iEvent.getByLabel(MuonTags_,muHandle);
//   const reco::MuonTrackLinksCollection muColl = *(muHandle.product());

    if (_debug) 
    {
        cout << "number of global muons " << MuCollection->size() << endl;
        cout << "number of STA muons " << STACollection->size() << endl;
        cout << "number of tracks 1 " << TrCollection->size() << endl;
        cout << "number of tracks 2 " << TrCollection2->size() << endl;
        cout << "number of tracks 3 " << TrCollection3->size() << endl;
    //    cout << "number of tracks 4 " << TrCollection4->size() << endl << endl;
    
   }

    hnstamuon->Fill(STACollection->size());
    hnglbmuon->Fill(MuCollection->size());
    hntracks->Fill(TrCollection->size());
    hntracks_2->Fill(TrCollection2->size());
    hntracks_3->Fill(TrCollection3->size());
  //  hntracks_4->Fill(TrCollection4->get());

    hntrackcorr->Fill(TrCollection->size(), TrCollection3->size());

    bool enabledt=true;
    bool enablecsc=false;
    bool enablerpc=true;

    muonMeasurements.setEvent(iEvent);
    MuonTransientTrackingRecHit::MuonRecHitContainer allHits;


    std::vector<const TrackingRecHit *> alltrackrechits;
    //TrackDetMatchInfo trAcInfo;
    TrackDetMatchInfoCollection trAcInfo;
    std::map<int, dethitstruct> detectorscrossed;
    std::map<int, bool> detectorscrossed_histfilled;
    std::map<int, int> matchedrechits;
    std::map<int, int> rechitsontrack;
    std::map<int, int> allrechits;
    //std::map<int, dethitstruct> rpcrechits_;

    for (MuonCollection::const_iterator sta = STACollection->begin(); sta!=STACollection->end(); sta++)
    {
        hstamuon_etaphi->Fill(sta->eta(), sta->phi());
        hstamuon_pt->Fill(sta->pt());
        cout << "standalon muon pt :" << sta->pt()<< endl;
    }
    for (MuonCollection::const_iterator glbmu = MuCollection->begin(); glbmu !=MuCollection->end(); glbmu ++)
    {
        hglbmuon_etaphi->Fill(glbmu->eta(), glbmu->phi());
        hglbmuon_pt->Fill(glbmu->pt());
    }
 
    if (_debug)
    {
        cout << "Number of RPC Rechits " << allRPChits->size() << endl;
    }
    for (RPCRecHitCollection::const_iterator irpchit =  allRPChits->begin(); irpchit != allRPChits->end(); irpchit++)
    {
        const GeomDet *whichdet = theG->idToDet(irpchit->geographicalId());
        const GlobalPoint &p = whichdet->surface().toGlobal(irpchit->localPosition());
        if (_debug) printf("collection RPChit detid=%d global position: %6.2f %6.2f %6.2f\n", irpchit->geographicalId().rawId(), p.x(), p.y(), p.z());
        allrechits[irpchit->geographicalId().rawId()]++;

        const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
        if (aroll)
        {
            // cout << "-----PRC-Roll-- " << irpchit->geographicalId().rawId() << " " << endl;
        }
        int rawid = irpchit->geographicalId().rawId();

    //    rpcrechits_[rawid].locx = irpchit->localPosition().x();
     //   rpcrechits_[rawid].locy = irpchit->localPosition().y();
     //   rpcrechits_[rawid].locz = irpchit->localPosition().z();
        // global position
    //    rpcrechits_[rawid].glbx = p.x();
     //   rpcrechits_[rawid].glby = p.y();
     //   rpcrechits_[rawid].glbz = p.z();

    }

    if (_debug)
    {
        cout << "Number of DT Rechits " << allDThits->size() << endl;
    }
    for (DTRecHitCollection::const_iterator idthit =  allDThits->begin(); idthit != allDThits->end(); idthit++)
    {
        const GeomDet *whichdet = theG->idToDet(idthit->geographicalId());
        const GlobalPoint &p = whichdet->surface().toGlobal(idthit->localPosition());
        if (_debug) printf("collection DThit detid=%d global position: %6.2f %6.2f %6.2f\n", idthit->geographicalId().rawId(), p.x(), p.y(), p.z());
        allrechits[idthit->geographicalId().rawId()]++;
    }

    if (_debug)
    {
        cout << "Number of CSC Rechits " << allCSChits->size() << endl;
    }
    for (CSCRecHit2DCollection::const_iterator icschit =  allCSChits->begin(); icschit != allCSChits->end(); icschit++)
    {
        const GeomDet *whichdet = theG->idToDet(icschit->geographicalId());
        const GlobalPoint &p = whichdet->surface().toGlobal(icschit->localPosition());
        if (_debug) printf("collection CSChit detid=%d global position: %6.2f %6.2f %6.2f\n", icschit->geographicalId().rawId(), p.x(), p.y(), p.z());
        allrechits[icschit->geographicalId().rawId()]++;
    }

    bool muonwentthrough = false;

    for (TrackCollection::const_iterator track = TrCollection3->begin(); track !=TrCollection3->end(); track++)
    {
        htracks_etaphi_3->Fill(track->eta(), track->phi());
        htrackinnerxy_3->Fill(track->innerPosition().X(), track->innerPosition().Y());
        htrackouterxy_3->Fill(track->outerPosition().X(), track->outerPosition().Y());
        htrackinnerz_3->Fill(track->innerPosition().Z());
        htrackouterz_3->Fill(track->outerPosition().Z());
        const reco::HitPattern& p = track->hitPattern();

        // a muon with more than 50 hits went through the detector
        if (p.numberOfHits()>50) muonwentthrough = true;

        hnvalidhits_3->Fill(track->innerPosition().phi(), p.numberOfHits());
        hchi2prob_3->Fill(track->innerPosition().phi(), TMath::Prob(track->chi2(), track->ndof()));
        hchi2ndof_3->Fill(track->innerPosition().phi(), track->normalizedChi2());
        hqoverppull_3->Fill(track->innerPosition().phi(), track->qoverp()/track->qoverpError());
        htrackq_3->Fill(track->innerPosition().phi(), track->charge());
    }
////////////////////////////////////////////////////////////////////
    double endingypoint = 2000.0;
    //if (TrCollection->size()==2) muonwentthrough = true;
    for (TrackCollection::const_iterator track = TrCollection->begin(); track !=TrCollection->end(); track++)
    {

        double  endingypoinit_ = track->outerPosition().Y();
        if (endingypoint > endingypoinit_) endingypoint = track->outerPosition().Y();
 
        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                alltrackrechits.push_back(&(*(*hit)));
            }
        }
    }
    for (TrackCollection::const_iterator track = TrCollection2->begin(); track !=TrCollection2->end(); track++)
    {

        double  endingypoinit_ = track->outerPosition().Y();
        if (endingypoint > endingypoinit_) endingypoint = track->outerPosition().Y();

        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                alltrackrechits.push_back(&(*(*hit)));
            }
        }
    }    
///////////////////////////////////////////////////////////////////////
// Track = CosmicMuonBarrelOnly
    for (TrackCollection::const_iterator track = TrCollection->begin(); track !=TrCollection->end(); track++)
    {
        const reco::TransientTrack ttrack = ttrackBuilder->build(*track);
        const FreeTrajectoryState fts = ttrack.initialFreeState();
        const GlobalTrajectoryParameters trajpar = fts.parameters();

// the following all work OK
// method 1
        TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, fts, parameters_);

        trAcInfo.push_back(info);

        for(std::vector<TAMuonChamberMatch>::const_iterator chamber = info.chambers.begin();
            chamber!= info.chambers.end(); chamber++ )
        {
               int detRawId = chamber->id.rawId();
               if(_debug) cout << "subdetid : "<< chamber->id.subdetId() << endl;
                if ( !detLptMap_[detRawId] ) {

                        if ( chamber->id.subdetId() == 1 ) { // if DT hits
                                DTChamberId segId(detRawId);
                                int wheel = segId.wheel();
                                int station = segId.station(); 
                                int sector = segId.sector();
                                int lsize= dtlsize;

                                TFileDirectory subDir_DT = fs->mkdir( "DT" );
                                TFileDirectory dirWheel = subDir_DT.mkdir(Form("Wheel_%d", wheel));
                                TFileDirectory dirStation = dirWheel.mkdir(Form("Station_%d",station));
                                TFileDirectory subSub3DirSector=dirStation.mkdir( Form("sector_%d",sector));

                                detLptMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTt%d", detRawId),
                                                                        Form("DT Detid %d trajectory", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLpmMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTm%d", detRawId),
                                                                        Form("DT Detid %d match track->rec", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLprMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTr%d", detRawId),
                                                                        Form("DT Detid %d match rec ", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);

                                detMapResidureX_[detRawId] = subSub3DirSector.make<TH1F>(Form("DT_%d_Residure_X", detRawId),
                                                                        Form("DT Detid %d Residure X local", detRawId),
                                                                        400, -100, 100);

                        }
                        else if ( chamber->id.subdetId() == 2 ) { // if DT hits
                                CSCDetId segId(detRawId);
                                int endcap =  segId.endcap();
                                int station = segId.station();
                                int ring = segId.ring();
                                int chamber = segId.chamber();
                                int layer = segId.layer();
                                int lsize= csclsize;

                                TFileDirectory subDir_CSC = fs->mkdir( "CSC" );
                                TFileDirectory dirEndcap = subDir_CSC.mkdir(Form("Endcap_%d", endcap));
                                TFileDirectory dirStation = dirEndcap.mkdir(Form("Station_%d", station));
                                TFileDirectory dirRing = dirStation.mkdir(Form("Ring_%d", ring));
                                TFileDirectory dirChamber = dirRing.mkdir(Form("Chamber_%d", chamber));
                                TFileDirectory dirLayer = dirChamber.mkdir(Form("Layer_%d", layer));

                                detLptMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCt%d", detRawId),
                                                                        Form("CSC Detid %d trajectory", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLpmMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCm%d", detRawId),
                                                                        Form("CSC Detid %d match track->rec", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLprMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCr%d", detRawId),
                                                                        Form("CSC Detid %d match rec ", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detMapResidureX_[detRawId] = dirLayer.make<TH1F>(Form("CSC_%d_Residure_X", detRawId),
                                                                        Form("CSC Detid %d Residure X local", detRawId),
                                                                        400, -100, 100);

                        }
                        else if(chamber->id.subdetId() == 3) {
                                RPCDetId segId(detRawId);
                                RPCGeomServ servId(detRawId);

                                int region = segId.region(); 
	                        int ring = segId.ring();
	                        int station = segId. station();
	                        int sector = segId.sector();
	                        int layer = segId.layer();
	                        int subsector = segId.subsector();
	                        int roll  = segId.roll();
                                int lsize= rpclsize;
/*
                                string label;
                                string subs, rollN;

                                if(region == 0 )
                                {

                                    if((station ==1 || station ==2) && subsector==1) subs="in";
                                    else if((station ==1 || station ==2) && subsector==2) subs="out";

                                    if(station == 3 && subsector ==1) subs="-";
                                    else if(station == 3 && subsector ==2) subs="+";

                                    if(station == 4 && (sector==1 || sector==2 || sector==3
                                                   || sector==5 || sector==6
                                                   || sector==7 || sector==8
                                                   || sector==10             || sector==12) && subsector==2) subs = "+";
                                    else if(station==4 && sector==4 && subsector==1) subs="--";
                                    else if(station==4 && sector==4 && subsector==2) subs="-";
                                    else if(station==4 && sector==4 && subsector==3) subs="+";
                                    else if(station==4 && sector==4 && subsector==4) subs="++";
                                    else if(station==4 ) subs="-";

                                    if(roll==1) rollN="Backward";
                                    if(roll==2) rollN="Central";
                                    if(roll==3) rollN="Forward";
                                    if(roll==4) rollN="D";

                                    stringstream ss;
                                    if ( ring<1 && sector<10 ) ss << "W"  << ring << "_RB" << station << subs << "_S0" << sector << "_" << rollN;
                                    if ( ring<1 && sector>9  ) ss << "W"  << ring << "_RB" << station << subs << "_S"  << sector << "_" << rollN;
                                    if ( ring>0 && sector<10 ) ss << "W+" << ring << "_RB" << station << subs << "_S0" << sector << "_" << rollN;
                                    if ( ring>0 && sector>9  ) ss << "W+"  << ring << "_RB" << station << subs << "_S"  << sector << "_" << rollN;
                                    ss >> label;

                                    cout << label << "\n" << endl;

                                }
                                else 
                                {
                                    stringstream ss;
                                    ss << "detid" << detRawId;
                                    ss >> label;
                                }
*/                               // RPCGeomServ servId(detRawId); 
                                cout << "RPCGeomServ : " << servId.name() << " :: " << endl;

                                TFileDirectory subDir_RPC = fs->mkdir( "RPC" );
                                TFileDirectory dirRegion = subDir_RPC.mkdir(Form("Region_%d", region));
                                TFileDirectory dirRing = dirRegion.mkdir(Form("Ring_%d", ring));
                                TFileDirectory dirStation = dirRing.mkdir(Form("Station_%d", station));
                                TFileDirectory dirSector = dirStation.mkdir( Form("Sector_%d",sector));
                               // TFileDirectory dirLayer = dirSector.mkdir(Form("Layer_%d", layer));
                               // TFileDirectory dirSubsector = dirLayer.mkdir(Form("Subsector_%d", subsector));
                               // TFileDirectory dirRoll = dirSubsector.mkdir(Form("Roll_%d", roll));

                                detLptMap_[detRawId] = dirSector.make<TH2F>(Form("RPCt%d", detRawId),
                                                                        Form("%s trajectory", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLpmMap_[detRawId] = dirSector.make<TH2F>(Form("RPCm%d", detRawId),
                                                                        Form("%s match track->rec", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLprMap_[detRawId] = dirSector.make<TH2F>(Form("RPCr%d", detRawId),
                                                                        Form("%s match rec ", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);


                              //  detMapResidureXY_[detRawId] = dirRoll.make<TH2F>(Form("RPC_%d_Residure_X_Y", detRawId),
                              //                                          Form("RPC Detid %d Residure X local", detRawId),
                              //                                          100, -100, 100,100, -100, 100);

                                detMapResidureX_[detRawId] = dirSector.make<TH1F>(Form("RPC_%d_Residure_X", detRawId),
                                                                        Form("%s Residure X local", servId.name().c_str()),
                                                                        400, -100, 100);
                                detMapResidureXadd[detRawId] = dirSector.make<TH1F>(Form("RPC_%d_Residure_addX", detRawId),
                                                                        Form("%s Residure addX local", servId.name().c_str()),
                                                                        400, -100, 100);

                  /*              detMapResidureY_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_Y", detRawId),
                                                                        Form("RPC Detid %d Residure Y local", detRawId),
                                                                        400, -100, 100);
                                detMapResidureZ_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_Z", detRawId),
                                                                        Form("RPC Detid %d Residure Z local", detRawId),
                                                                        400, -100, 100);

                              //  detMapResidureGXY_[detRawId] = dirRoll.make<TH2F>(Form("RPC_%d_Residure_GX_GY", detRawId),
                              //                                          Form("RPC Detid %d Residure GX", detRawId),
                              //                                          100, -100, 100,100, -100, 100);
                                detMapResidureGX_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GX", detRawId),
                                                                        Form("RPC Detid %d Residure GX", detRawId),
                                                                        100, -200, 200);
                                detMapResidureGY_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GY", detRawId),
                                                                        Form("RPC Detid %d Residure GY", detRawId),
                                                                        100, -200, 200);
                                detMapResidureGZ_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GZ", detRawId),
                                                                        Form("RPC Detid %d Residure GZ", detRawId),
                                                                        100, -200, 200);
                   */
                        }
                }

// //////////////////////////
         if (chamber->tState.globalPosition().y() > endingypoint)
                detectorscrossed[chamber->id.rawId()].count++;
        }
    }
//////////////////////////////////////////////////////////////////////////
//Track = CosmicMuonEndCapOnly
    for (TrackCollection::const_iterator track = TrCollection2->begin(); track !=TrCollection2->end(); track++)
    {
        const reco::TransientTrack ttrack = ttrackBuilder->build(*track);
        const FreeTrajectoryState fts = ttrack.initialFreeState();
        const GlobalTrajectoryParameters trajpar = fts.parameters();

// the following all work OK
// method 1
        TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, fts, parameters_);

        trAcInfo.push_back(info);

        for(std::vector<TAMuonChamberMatch>::const_iterator chamber = info.chambers.begin();
            chamber!= info.chambers.end(); chamber++ )
        {
               int detRawId = chamber->id.rawId();
               int lsize = 200, csclsize = 110;
               if(_debug) cout << "subdetid : "<< chamber->id.subdetId() << endl;
                if ( !detLptMap_[detRawId] ) {

                        if ( chamber->id.subdetId() == 1 ) { // if DT hits
                                DTChamberId segId(detRawId);
                                int wheel = segId.wheel();
                                int station = segId.station(); 
                                int sector = segId.sector();
                                int lsize= dtlsize;

                                TFileDirectory subDir_DT = fs->mkdir( "DT" );
                                TFileDirectory dirWheel = subDir_DT.mkdir(Form("Wheel_%d", wheel));
                                TFileDirectory dirStation = dirWheel.mkdir(Form("Station_%d",station));
                                TFileDirectory subSub3DirSector=dirStation.mkdir( Form("sector_%d",sector));

                                detLptMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTt%d", detRawId),
                                                                        Form("DT Detid %d trajectory", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLpmMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTm%d", detRawId),
                                                                        Form("DT Detid %d match track->rec", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLprMap_[detRawId] = subSub3DirSector.make<TH2F>(Form("DTr%d", detRawId),
                                                                        Form("DT Detid %d match rec ", detRawId),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);

                                detMapResidureX_[detRawId] = subSub3DirSector.make<TH1F>(Form("DT_%d_Residure_X", detRawId),
                                                                        Form("DT Detid %d Residure X local", detRawId),
                                                                        400, -100, 100);

                        }
                        else if ( chamber->id.subdetId() == 2 ) { // if DT hits
                                CSCDetId segId(detRawId);
                                int endcap =  segId.endcap();
                                int station = segId.station();
                                int ring = segId.ring();
                                int chamber = segId.chamber();
                                int layer = segId.layer();
                                int lsize= csclsize;

                                TFileDirectory subDir_CSC = fs->mkdir( "CSC" );
                                TFileDirectory dirEndcap = subDir_CSC.mkdir(Form("Endcap_%d", endcap));
                                TFileDirectory dirStation = dirEndcap.mkdir(Form("Station_%d", station));
                                TFileDirectory dirRing = dirStation.mkdir(Form("Ring_%d", ring));
                                TFileDirectory dirChamber = dirRing.mkdir(Form("Chamber_%d", chamber));
                                TFileDirectory dirLayer = dirChamber.mkdir(Form("Layer_%d", layer));

                                detLptMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCt%d", detRawId),
                                                                        Form("CSC Detid %d trajectory", detRawId),
                                                                        25, -csclsize, csclsize, 25, -csclsize, csclsize);
                                detLpmMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCm%d", detRawId),
                                                                        Form("CSC Detid %d match track->rec", detRawId),
                                                                        25, -csclsize, csclsize, 25, -csclsize, csclsize);
                                detLprMap_[detRawId] = dirLayer.make<TH2F>(Form("CSCr%d", detRawId),
                                                                        Form("CSC Detid %d match rec ", detRawId),
                                                                        25, -csclsize, csclsize, 25, -csclsize, csclsize);

                                detMapResidureX_[detRawId] = dirLayer.make<TH1F>(Form("CSC_%d_Residure_X", detRawId),
                                                                        Form("CSC Detid %d Residure X local", detRawId),
                                                                        400, -100, 100);
                        }
                        else if(chamber->id.subdetId() == 3) {
                                RPCDetId segId(detRawId);
                                int region = segId.region(); 
	                        int ring = segId.ring();
	                        int station = segId. station();
	                        int sector = segId.sector();
	                        int layer = segId.layer();
	                        int subsector = segId.subsector();
	                        int roll  = segId.roll();
                                int lsize= rpclsize;

/*                                string label;
                                string subs, rollN;

                                if(region == 0 )
                                {

                                    if((station ==1 || station ==2) && subsector==1) subs="in";
                                    else if((station ==1 || station ==2) && subsector==2) subs="out";

                                    if(station == 3 && subsector ==1) subs="-";
                                    else if(station == 3 && subsector ==2) subs="+";

                                    if(station == 4 && (sector==1 || sector==2 || sector==3
                                                   || sector==5 || sector==6
                                                   || sector==7 || sector==8
                                                   || sector==10             || sector==12) && subsector==2) subs = "+";
                                    else if(station==4 && sector==4 && subsector==1) subs="--";
                                    else if(station==4 && sector==4 && subsector==2) subs="-";
                                    else if(station==4 && sector==4 && subsector==3) subs="+";
                                    else if(station==4 && sector==4 && subsector==4) subs="++";
                                    else if(station==4 ) subs="-";

                                    if(roll==1) rollN="Backward";
                                    if(roll==2) rollN="Central";
                                    if(roll==3) rollN="Forward";
                                    if(roll==4) rollN="D";

                                    stringstream ss;
                                    if ( ring<1 && sector<10 ) ss << "W"  << ring << "_RB" << station << subs << "_S0" << sector << "_" << rollN;
                                    if ( ring<1 && sector>9  ) ss << "W"  << ring << "_RB" << station << subs << "_S"  << sector << "_" << rollN;
                                    if ( ring>0 && sector<10 ) ss << "W+" << ring << "_RB" << station << subs << "_S0" << sector << "_" << rollN;
                                    if ( ring>0 && sector>9  ) ss << "W+"  << ring << "_RB" << station << subs << "_S"  << sector << "_" << rollN;
                                    ss >> label;

                                    cout << label << "\n" << endl;

                                }
                                else 
                                {
                                    stringstream ss;
                                    ss << "detid" << detRawId; 
                                    ss >> label;
                                }    
*/
                                RPCGeomServ servId(detRawId);
                                cout << "RPCGeomServ : " << servId.name() << " :: " << endl;


                                TFileDirectory subDir_RPC = fs->mkdir( "RPC" );
                                TFileDirectory dirRegion = subDir_RPC.mkdir(Form("Region_%d", region));
                                TFileDirectory dirRing = dirRegion.mkdir(Form("Ring_%d", ring));
                                TFileDirectory dirStation = dirRing.mkdir(Form("Station_%d", station));
                                TFileDirectory dirSector = dirStation.mkdir( Form("Sector_%d",sector));
                                //TFileDirectory dirLayer = dirSector.mkdir(Form("Layer_%d", layer));
                                //TFileDirectory dirSubsector = dirLayer.mkdir(Form("Subsector_%d", subsector));
                                //TFileDirectory dirRoll = dirSubsector.mkdir(Form("Roll_%d", roll));

                                detLptMap_[detRawId] = dirSector.make<TH2F>(Form("RPCt%d", detRawId),
                                                                        Form("%s trajectory", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLpmMap_[detRawId] = dirSector.make<TH2F>(Form("RPCm%d", detRawId),
                                                                        Form("%s match track->rec", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);
                                detLprMap_[detRawId] = dirSector.make<TH2F>(Form("RPCr%d", detRawId),
                                                                        Form("%s match rec ", servId.name().c_str()),
                                                                        25, -lsize, lsize, 25, -lsize, lsize);

                              //  detMapResidureXY_[detRawId] = dirRoll.make<TH2F>(Form("RPC_%d_Residure_X_Y", detRawId),
                              //                                          Form("RPC Detid %d Residure X local", detRawId),
                              //                                          100, -100, 100,100, -100, 100);

                                detMapResidureX_[detRawId] = dirSector.make<TH1F>(Form("RPC_%d_Residure_X", detRawId),
                                                                        Form("%s Residure X local", servId.name().c_str()),
                                                                        400, -100, 100);
                                detMapResidureXadd[detRawId] = dirSector.make<TH1F>(Form("RPC_%d_Residure_addX", detRawId),
                                                                        Form("%s Residure addX local", servId.name().c_str()),
                                                                        400, -100, 100);

                    /*          detMapResidureY_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_Y", detRawId),
                                                                        Form("RPC Detid %d Residure Y local", detRawId),
                                                                        400, -100, 100);
                                detMapResidureZ_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_Z", detRawId),
                                                                        Form("RPC Detid %d Residure Z local", detRawId),
                                                                        400, -100, 100);

                              //  detMapResidureGXY_[detRawId] = dirRoll.make<TH2F>(Form("RPC_%d_Residure_GX_GY", detRawId),
                              //                                          Form("RPC Detid %d Residure GX", detRawId),
                              //                                          100, -100, 100,100, -100, 100);
                                detMapResidureGX_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GX", detRawId),
                                                                        Form("RPC Detid %d Residure GX", detRawId),
                                                                        100, -200, 200);
                                detMapResidureGY_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GY", detRawId),
                                                                        Form("RPC Detid %d Residure GY", detRawId),
                                                                        100, -200, 200);
                                detMapResidureGZ_[detRawId] = dirRoll.make<TH1F>(Form("RPC_%d_Residure_GZ", detRawId),
                                                                        Form("RPC Detid %d Residure GZ", detRawId),
                                                                        100, -200, 200);
                     */
                        }
                }

// //////////////////////////
         if (chamber->tState.globalPosition().y() > endingypoint)
                detectorscrossed[chamber->id.rawId()].count++;
        }
    }


//////////////////////////////////////////////////////////////////////    
    for (TrackCollection::const_iterator track = TrCollection->begin(); track !=TrCollection->end(); track++)
    {
        htracks_etaphi->Fill(track->eta(), track->phi());
        htrackinnerxy->Fill(track->innerPosition().X(), track->innerPosition().Y());
        htrackouterxy->Fill(track->outerPosition().X(), track->outerPosition().Y());

        htrackinnerz->Fill(track->innerPosition().Z());
        htrackouterz->Fill(track->outerPosition().Z());
        const reco::HitPattern& p = track->hitPattern();
//hnvalidhits->Fill(p.numberOfValidTrackerHits());
        hnvalidhits->Fill(track->innerPosition().phi(), p.numberOfHits());
        hchi2prob->Fill(track->innerPosition().phi(), TMath::Prob(track->chi2(), track->ndof()));
        hchi2ndof->Fill(track->innerPosition().phi(), track->normalizedChi2());
        hqoverppull->Fill(track->innerPosition().phi(), track->qoverp()/track->qoverpError());
        htrackq->Fill(track->innerPosition().phi(), track->charge());

        int ihit=1;
        if (_debug) cout << "Number of Hits on track" << p.numberOfHits() << endl;
        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                htrackrechitsxy->Fill(p.x(), p.y());
                htrackrechitsz->Fill(p.z());
                if (_debug && (*hit)->geographicalId().subdetId()==1) printf("DT  Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());
                else if (_debug && (*hit)->geographicalId().subdetId()==2) printf("CSC Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());
                else if (_debug  && (*hit)->geographicalId().subdetId()==3) printf("RPC Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());

                if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==1)
                {
                    htrackrechitsDTxy->Fill(p.x(), p.y());
                    htrackrechitsDTz->Fill(p.z());
                }
                else if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==2)
                {
                    htrackrechitsCSCxy->Fill(p.x(), p.y());
                    htrackrechitsCSCz->Fill(p.z());
                }
                else if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==3)
                {
                    htrackrechitsRPCxy->Fill(p.x(), p.y());
                    htrackrechitsRPCz->Fill(p.z());
                }
            }
        }
    }
// ////////////////////////////////////////////////
    for (TrackCollection::const_iterator track = TrCollection2->begin(); track !=TrCollection2->end(); track++)
    {
        htracks_etaphi_2->Fill(track->eta(), track->phi());
        htrackinnerxy_2->Fill(track->innerPosition().X(), track->innerPosition().Y());
        htrackouterxy_2->Fill(track->outerPosition().X(), track->outerPosition().Y());
        htrackinnerz_2->Fill(track->innerPosition().Z());
        htrackouterz_2->Fill(track->outerPosition().Z());
        const reco::HitPattern& p = track->hitPattern();
//hnvalidhits->Fill(p.numberOfValidTrackerHits());
        hnvalidhits_2->Fill(track->innerPosition().phi(), p.numberOfHits());
        hchi2prob_2->Fill(track->innerPosition().phi(), TMath::Prob(track->chi2(), track->ndof()));
        hchi2ndof_2->Fill(track->innerPosition().phi(), track->normalizedChi2());
        hqoverppull_2->Fill(track->innerPosition().phi(), track->qoverp()/track->qoverpError());
        htrackq_2->Fill(track->innerPosition().phi(), track->charge());
// loop over hits associated with a track
        int ihit=1;
        if (_debug) cout << "Number of Hits on forward track" << p.numberOfHits() << endl;
        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                htrackrechitsxy_2->Fill(p.x(), p.y());
                htrackrechitsz_2->Fill(p.z());
             //   rechitsontrack[(*hit)->geographicalId().rawId()]++;
                if (_debug && (*hit)->geographicalId().subdetId()==1) printf("DT  Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());
                else if (_debug && (*hit)->geographicalId().subdetId()==2) printf("CSC Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());
                else if (_debug  && (*hit)->geographicalId().subdetId()==3) printf("RPC Rec Hit id=%d %6.2f %6.2f %6.2f\n", (*hit)->geographicalId().rawId(), p.x(), p.y(), p.z());

                if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==1)
                {
                    htrackrechitsDTxy_2->Fill(p.x(), p.y());
                    htrackrechitsDTz_2->Fill(p.z());
                }
                else if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==2)
                {
                    htrackrechitsCSCxy_2->Fill(p.x(), p.y());
                    htrackrechitsCSCz_2->Fill(p.z());
                }
                else if ((*hit)->geographicalId().det()==DetId::Muon
                    && (*hit)->geographicalId().subdetId()==3)
                {
                    htrackrechitsRPCxy_2->Fill(p.x(), p.y());
                    htrackrechitsRPCz_2->Fill(p.z());
                }
            }
        }
    }
// ////////////////////////////////////////////
    for(std::vector<TrackDetMatchInfo>::const_iterator tr = trAcInfo.begin();tr != trAcInfo.end(); tr++)
    {
        if (_debug) cout << "Muon detector matching details: "  << endl;
        for(std::vector<TAMuonChamberMatch>::const_iterator chamber = tr->chambers.begin();
            chamber!= tr->chambers.end(); chamber++ )
        {

            if (_debug) cout << chamber->info() << "\n\t(DetId, station, localx, localy): "
                << chamber->id.rawId() << ", "
                << chamber->station() << ", "
                << chamber->tState.localPosition().x() << " +- " << sqrt(chamber->tState.localError().positionError().xx()) << "  "
                << chamber->tState.localPosition().y() << " +- " << sqrt(chamber->tState.localError().positionError().yy())
                << endl;

            if (_debug) cout << "\t trajectory global point (x,y,z): "
                << chamber->tState.globalPosition().x() << ", "
                << chamber->tState.globalPosition().y() << ", "
                << chamber->tState.globalPosition().z() << endl;

            double hitz = chamber->tState.globalPosition().z();
            int zbin;

            int rawid = chamber->id.rawId();

            //detectorscrossed[rawid].count++;
            detectorscrossed_histfilled[rawid] = false;
            // local position
            detectorscrossed[rawid].locx = chamber->tState.localPosition().x();
            detectorscrossed[rawid].locy = chamber->tState.localPosition().y();
            detectorscrossed[rawid].locz = chamber->tState.localPosition().z();
            // global position
            detectorscrossed[rawid].glbx = chamber->tState.globalPosition().x();
            detectorscrossed[rawid].glby = chamber->tState.globalPosition().y();
            detectorscrossed[rawid].glbz = chamber->tState.globalPosition().z();


            const GeomDet *whichdet = theG->idToDet(chamber->id);
            const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
            const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(whichdet);
            const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(whichdet);

            if (aroll)
            {
                const float stripPredicted =aroll->strip(LocalPoint(chamber->tState.localPosition().x(),chamber->tState.localPosition().y(),0.));
                if (_debug) cout << "Expected strip # " <<  stripPredicted << " out of " << aroll->nstrips() << endl;

       /*         if (detectorscrossed[rawid].count>icn)
                {
                   for(map<int, dethitstruct>::iterator rpchit__= rpcrechits_.begin();
                             rpchit__ != rpcrechits_.end();++rpchit__ )
                   {
                       detMapResidureGX_[rawid]->Fill(detectorscrossed[rawid].glbx-(* rpchit__).second.glbx);
                       detMapResidureGY_[rawid]->Fill(detectorscrossed[rawid].glby-(* rpchit__).second.glby);
                       detMapResidureGZ_[rawid]->Fill(detectorscrossed[rawid].glbz-(* rpchit__).second.glbz);
                       if(_debug) cout << "RPC residual : " << rawid <<" " << (* rpchit__).first << endl;
                  }
         
                }
         */    
                detectorscrossed[rawid].subdet = 3;
                if (detectorscrossed[rawid].count>icn 
                        && !detectorscrossed_histfilled[rawid])
                {
                    float residual, sum;
                    if (findmatch(allRPChits, rawid, detectorscrossed[rawid].locx, detectorscrossed[rawid].locy, residual, sum))
                    {
                        RPCDetId segId(rawid);
                        int region = segId.region();
                        int ring = segId.ring()-1;
                        zbin = ring+3;
                       if(region==0) hrpctrajectoryxy_r0[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                       else if(region==-1) hrpctrajectoryxy_rm1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                       else hrpctrajectoryxy_rp1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                        hrpctrajectoryz->Fill(detectorscrossed[rawid].glbz);
                       if(region==0) hrpcmatchrechitxy_r0[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                       else if(region==-1) hrpcmatchrechitxy_rm1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                       else hrpcmatchrechitxy_rp1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                        hrpcmatchrechitz->Fill(detectorscrossed[rawid].glbz);
                        detectorscrossed_histfilled[rawid] = true;

                        matchedrechits[rawid]++;

                        detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                        detLprMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);


                        detMapResidureX_[rawid]->Fill(residual);
                        detMapResidureXadd[rawid]->Fill(sum);

                       // detMapResidureY_[rawid]->Fill(detectorscrossed[rawid].locy-rpcrechits_[rawid].locy);
                      //  detMapResidureZ_[rawid]->Fill(detectorscrossed[rawid].locz-rpcrechits_[rawid].locz );

                        cout << "RPC residual :" << residual << endl; //" localposition z: " << rpcrechits_[rawid].locz << endl;


/*                        for(map<int, dethitstruct>::iterator rpchit__= rpcrechits_.begin();
                                  rpchit__ != rpcrechits_.end();++rpchit__ )
                        {    
                             detMapResidureGX_[rawid]->Fill(detectorscrossed[rawid].glbx-(* rpchit__).second.glbx);
                             detMapResidureGY_[rawid]->Fill(detectorscrossed[rawid].glby-(* rpchit__).second.glby);
                             detMapResidureGZ_[rawid]->Fill(detectorscrossed[rawid].glbz-(* rpchit__).second.glbz );
                             if(_debug) cout << "RPC residual : " << rawid <<" " << (* rpchit__).first << endl;  
                        }
*/
                        bool matchfoundontrack = false;
                        for (std::vector<const TrackingRecHit *>::const_iterator hit=alltrackrechits.begin(); 
                                hit != alltrackrechits.end() && !matchfoundontrack; hit++)
                        {
                            //const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                            //const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                            if ((*hit)->geographicalId().rawId() == rawid
                                    && fabs( detectorscrossed[rawid].locx - (*hit)->localPosition().x()) < maxdist )
                            {

                                if(region==0) hrpcmatchrechitontrackxy_r0[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                                else if(region==-1) hrpcmatchrechitontrackxy_rm1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby); 
                                else hrpcmatchrechitontrackxy_rp1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                                hrpcmatchrechitontrackz->Fill(detectorscrossed[rawid].glbz);
                                rechitsontrack[rawid]++;
                                matchfoundontrack = true;
                                detLpmMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                            }
                        }
                    }
                }
            }
            else if (dtlayer)
            {
//                typedef std::vector<DTLocalTrigger>::const_iterator DigiConstIterator;
//                std::pair<DigiConstIterator,DigiConstIterator> range = TrCollection4->get(chamber->id.rawId());
                //  cout << "++++-----" << (double)(range) << " ++++" << endl;
//             if ( range.first == range.second )
//             {
                 //cout << "------ END " << endl;
                 //else  cout << " found -------------  "  <<  endl;

                const float stripPredicted =dtlayer->specificTopology().channel(LocalPoint(chamber->tState.localPosition().x(),chamber->tState.localPosition().y(),0.));
                if (_debug) cout << "Expected wire # " <<  stripPredicted 
                    << " out of " << dtlayer->specificTopology().channels() << endl;
                detectorscrossed[rawid].subdet = 1;
                if (detectorscrossed[rawid].count>icn 
                        && !detectorscrossed_histfilled[rawid])
                {
                    float residual;
                    if (findmatch(allDThits, rawid, detectorscrossed[rawid].locx, detectorscrossed[rawid].locy, residual))
                    {
                        DTChamberId segId(rawid);
                        int wheel = segId.wheel();
                        zbin = wheel+2;

                        hdttrajectoryxy[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        hdttrajectoryz->Fill(detectorscrossed[rawid].glbz);
                        hdtmatchrechitxy[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        hdtmatchrechitz->Fill(detectorscrossed[rawid].glbz);
                        detectorscrossed_histfilled[rawid] = true;

                        matchedrechits[rawid]++;
//////
                        detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                        detLprMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
//////
                        detMapResidureX_[rawid]->Fill(residual);

                        bool matchfoundontrack = false;
                        for (std::vector<const TrackingRecHit *>::const_iterator hit=alltrackrechits.begin(); 
                                hit != alltrackrechits.end() && !matchfoundontrack; hit++)
                        {
                            //const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                            //const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                            if ((*hit)->geographicalId().rawId() == rawid
                                    && fabs( detectorscrossed[rawid].locx - (*hit)->localPosition().x()) < maxdist )
                            {
                                hdtmatchrechitontrackxy[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                                hdtmatchrechitontrackz->Fill(detectorscrossed[rawid].glbz);
                                rechitsontrack[rawid]++;
                                matchfoundontrack = true;
                                detLpmMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                            }
                        }
                    }
//                }
              }
              else {
                  const float stripPredicted =dtlayer->specificTopology().channel(LocalPoint(chamber->tState.localPosition().x(),chamber->tState.localPosition().y(),0.));
                  if (_debug) cout << "trigger hit Expected wire # " <<  stripPredicted
                       << " out of " << dtlayer->specificTopology().channels() << endl;
                  detectorscrossed[rawid].subdet = 1;
                  if (detectorscrossed[rawid].count>icn
                        && !detectorscrossed_histfilled[rawid])
                  {
                     float residual;
                     if (findmatch(allDThits, rawid, detectorscrossed[rawid].locx, detectorscrossed[rawid].locy, residual))
                     {
                        DTChamberId segId(rawid);
                        int wheel = segId.wheel();
                        zbin = wheel+2;

                        hdttriggerhitxy[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        hdttriggerhitz->Fill(detectorscrossed[rawid].glbz);
                        cout << " DT trigger hit found. id:"<< chamber->id.rawId() << "  " <<  endl;
                        detectorscrossed_histfilled[rawid] = true;
                     }
                  }
              }
            }
            else if (csclayer)
            {
                        //std::cout << "1CSC-Global-position Z: " << detectorscrossed[rawid].glbz << "test" << std::endl;
                detectorscrossed[rawid].subdet = 2;
                if (detectorscrossed[rawid].count>icn 
                        && !detectorscrossed_histfilled[rawid])
                {
                    float residual;
                    if (findmatch(allCSChits, rawid, detectorscrossed[rawid].locx, detectorscrossed[rawid].locy, residual))
                    {
                        CSCDetId segId(rawid);
                        int endcap =  segId.endcap();
                        int station = segId.station();
                        zbin = station-1;

                        std::cout << "2CSC-Global-position Zbin: " << zbin << std::endl;
                      if(zbin<4 && zbin > -1)
                      {  
                        if(endcap == 1) hcsctrajectoryxy_e1[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        else hcsctrajectoryxy_e2[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                        hcsctrajectoryz->Fill(detectorscrossed[rawid].glbz);

                        if(endcap == 1) hcscmatchrechitxy_e1[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        else hcscmatchrechitxy_e2[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                        hcscmatchrechitz->Fill(detectorscrossed[rawid].glbz);
                        detectorscrossed_histfilled[rawid] = true;

                        matchedrechits[rawid]++;
                        detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                        detLprMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);

                        detMapResidureX_[rawid]->Fill(residual);

                        bool matchfoundontrack = false;
                        for (std::vector<const TrackingRecHit *>::const_iterator hit=alltrackrechits.begin(); 
                                hit != alltrackrechits.end() && !matchfoundontrack; hit++)
                        {
                            //const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                            //const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                            if ((*hit)->geographicalId().rawId() == rawid
                                    && fabs( detectorscrossed[rawid].locx - (*hit)->localPosition().x()) < maxdist )
                            {
                                if(endcap == 1) hcscmatchrechitontrackxy_e1[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                                else hcscmatchrechitontrackxy_e2[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                                hcscmatchrechitontrackz->Fill(detectorscrossed[rawid].glbz);
                                rechitsontrack[rawid]++;
                                matchfoundontrack = true;
                                detLpmMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                            }
                        }
                      }
                    }
                }
            }
/*
  LogVerbatim("CosmicAnalyzer") << "\t trajectory global point (z,perp,eta,phi): "
<< chamber->tState.globalPosition().z() << ", "
<< chamber->tState.globalPosition().perp() << ", "
<< chamber->tState.globalPosition().eta() << ", "
<< chamber->tState.globalPosition().phi() ;
*/
            if (_debug) 
            {
                cout << "\t trajectory local point (x,y): "
                << chamber->tState.localPosition().x() << ", "
                << chamber->tState.localPosition().y() <<endl;

                for(std::vector<TAMuonSegmentMatch>::const_iterator segment=chamber->segments.begin();
                    segment!=chamber->segments.end(); segment++)
                {
                    cout << "\t segment global point (x,y,z): "
                        << segment->segmentGlobalPosition.x() << ", "
                        << segment->segmentGlobalPosition.y() << ", "
                        << segment->segmentGlobalPosition.z() << endl;
    /*
     LogVerbatim("TrackAssociator") << "\t segment position (z,Rho,eta,phi,DetId): "
       << segment->segmentGlobalPosition.z() << ", "
       << segment->segmentGlobalPosition.Rho() << ", "
       << segment->segmentGlobalPosition.eta() << ", "
       << segment->segmentGlobalPosition.phi() << ", "
       << chamber->id.rawId();
       */
    //LogVerbatim("TrackAssociator") << "\t segment local position (x,y): "
    //  << segment->segmentLocalPosition.x() << ", "
    //  << segment->segmentLocalPosition.y();
                }
            } // debug
        }
// build a transient track with reverse trajectory
/*
        const CartesianTrajectoryError trajerr = fts.cartesianError();
        const CurvilinearTrajectoryError trajerr2 = fts.curvilinearError();
// build a global trajectory paramter with momentum inverted to travel the other way
        GlobalVector invertedmomentum=trajpar.momentum();
        invertedmomentum *= -1.0;                 // invert
                                                  // you need the - sign
        GlobalTrajectoryParameters trajparinv(trajpar.position(), invertedmomentum, -trajpar.charge(), &(*theMGField));
        FreeTrajectoryState fts2(trajparinv, trajerr, trajerr2);

        TrackDetMatchInfo info2 = trackAssociator_.associate(iEvent, iSetup, fts2, parameters_);
        if (_debug) cout << "Reverse direction muon detector matching details: "  << endl;
        for(std::vector<TAMuonChamberMatch>::const_iterator chamber = info2.chambers.begin();
            chamber!=info2.chambers.end(); chamber++)
        {
            if (_debug) cout << chamber->info() << "\n\t(DetId, station, localX, localY): "
                << chamber->id.rawId() << ", "
                << chamber->station() << ", "
                << chamber->tState.localPosition().x() << " +- " << sqrt(chamber->tState.localError().positionError().xx()) << "  "
                << chamber->tState.localPosition().y() << " +- " << sqrt(chamber->tState.localError().positionError().yy()) << endl;
            if (_debug) cout << "\t trajectory global point (x,y,z): "
                << chamber->tState.globalPosition().x() << ", "
                << chamber->tState.globalPosition().y() << ", "
                << chamber->tState.globalPosition().z() << endl;

            detectorscrossed[chamber->id.rawId()]++;
            const GeomDet *whichdet = theG->idToDet(chamber->id);
            const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
            const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(whichdet);
            if (aroll)
            {
                const float stripPredicted =aroll->strip(LocalPoint(chamber->tState.localPosition().x(),chamber->tState.localPosition().y(),0.));
                if (_debug) cout << "Expected strip # " <<  stripPredicted
                    << " out of " << aroll->nstrips() << endl;
            }
            else if (dtlayer)
            {
                const float stripPredicted =dtlayer->specificTopology().channel(LocalPoint(chamber->tState.localPosition().x(),chamber->tState.localPosition().y(),0.));
                if (_debug) cout << "Expected wire # " <<  stripPredicted 
                    << " out of " << dtlayer->specificTopology().channels() << endl;
            }
            if (_debug) cout << "\t trajectory local point (x,y): "
                << chamber->tState.localPosition().x() << ", "
                << chamber->tState.localPosition().y() << endl;

            for(std::vector<TAMuonSegmentMatch>::const_iterator segment=chamber->segments.begin();
                segment!=chamber->segments.end(); segment++)
            {
                if (_debug) cout << "\t segment global point (x,y,z): "
                    << segment->segmentGlobalPosition.x() << ", "
                    << segment->segmentGlobalPosition.y() << ", "
                    << segment->segmentGlobalPosition.z() << endl;
//LogVerbatim("TrackAssociator") << "\t segment local position (x,y): "
//  << segment->segmentLocalPosition.x() << ", "
//  << segment->segmentLocalPosition.y();
            }
        }
        */

// the following doesn't work  7/22/08
// BoundSurface is virtual :(
/*

for (vector<DetLayer*>::reverse_iterator irpclayer = rpcLayers.rbegin();
               irpclayer != rpcLayers.rend(); ++irpclayer) {
    MuonTransientTrackingRecHit::MuonRecHitContainer RHMB = muonMeasurements.recHits(*irpclayer);
     allHits.insert(allHits.end(),RHMB.begin(),RHMB.end());
     const BoundSurface rpcsurface = (*irpclayer)->surface();
     TrajectoryStateOnSurface tsosAtRPC = thePropagator->propagate(fts, (*irpclayer)->surface());

     if(tsosAtRPC.isValid()
&& irpclayer->bounds().inside(tsosAtRPC.localPosition(),tsosAtRPC.localError())
&& fabs(tsosAtRPC.localPosition().z()) < 10.0
&& fabs(tsosAtRPC.localPosition().x()) < 10.0
&& fabs(tsosAtRPC.localPosition().y()) < 10.0)
{
cout << "found RPC layer close to track" << endl;
cout << rpcsurface->position().x() << endl;
cout << rpcsurface->position().y() << endl;
cout << rpcsurface->position().z() << endl;
}
}
*/

    }                                             // end of loop over track collection 1


    //if (muonwentthrough)
    {
        for (map<int, dethitstruct>::const_iterator imh=detectorscrossed.begin(); imh != detectorscrossed.end(); imh++)
        {
            int rawid = imh->first;
            if (imh->second.count>icn || muonwentthrough) // forward and backward trajectories
            {
                allmuonhits[imh->first]++;
                if (_debug) cout << rawid << '\t';
                if (rechitsontrack.find(rawid)!=rechitsontrack.end())
                {
                    if (_debug) cout << 1 << '\t';
                    matchedmuonhitsontrack[rawid]++;
                }
                else
                {
                    if (_debug) cout << 0 << '\t';
                }

                if (allrechits.find(imh->first)!=allrechits.end())
                {
                    if (_debug) cout << 1 << endl;
                    matchedmuonhitsondet[rawid]++;
                }
                else
                {
                    if (_debug) cout << 0 << endl;
                }
            }

            // denominator histogram was not filled
            // need to fill it
            if (!detectorscrossed_histfilled[rawid] && imh->second.count>icn)
            {
                int det = detectorscrossed[rawid].subdet;
                
                int zbin;

                if (det==1) //DT
                {
                        DTChamberId segId(rawid);
                        int wheel = segId.wheel();
                        zbin = wheel+2;

                        hdttrajectoryxy[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                        hdttrajectoryz->Fill(detectorscrossed[rawid].glbz);
                        detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);

                }
                else if (det==2) // CSC
                {
                        CSCDetId segId(rawid);
                        int endcap =  segId.endcap();
                        int station = segId.station();
                        zbin = station-1;
                  if(zbin<4 && zbin > -1)
                  {   
                    if(endcap == 1) hcsctrajectoryxy_e1[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                    else hcsctrajectoryxy_e2[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);

                    hcsctrajectoryz->Fill(detectorscrossed[rawid].glbz);
                    detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                  }     
                }
                else if (det==3) // RPC
                {
                        RPCDetId segId(rawid);
                        int region = segId.region();
                        int ring = segId.ring()-1;
                        zbin = ring+3;

                    if(region==0) hrpctrajectoryxy_r0[zbin]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                    else if(region==-1) hrpctrajectoryxy_rm1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                    else hrpctrajectoryxy_rp1[ring]->Fill(detectorscrossed[rawid].glbx, detectorscrossed[rawid].glby);
                    hrpctrajectoryz->Fill(detectorscrossed[rawid].glbz);
                    detLptMap_[rawid]->Fill(detectorscrossed[rawid].locx,detectorscrossed[rawid].locy);
                }

            }
        }
    }

}


// ------------ method called once each job just before starting event loop  ------------
void
CosmicAnalyzer::beginJob(const edm::EventSetup&)
{
}


// ------------ method called once each job just after ending the event loop  ------------
void
CosmicAnalyzer::endJob()
{
    for (std::map<int,int>::const_iterator imh=allmuonhits.begin(); imh != allmuonhits.end(); imh++)
    {
        if (_debug) cout << imh->first << "\t" << imh->second << endl;
        if (matchedmuonhitsontrack.find(imh->first)!=0)
        {
            if (_debug) cout << matchedmuonhitsontrack[imh->first] << "\t efficiency = " << double(matchedmuonhitsontrack[imh->first])/double(imh->second) << endl;
        }
        else
        {
            if (_debug) cout <<  0 << "\t efficiency = " << 0 <<endl;
        }

        if (matchedmuonhitsondet.find(imh->first)!=0)
        {
            if (_debug) cout << matchedmuonhitsondet[imh->first] << "\t efficiency = " << double(matchedmuonhitsondet[imh->first])/double(imh->second) << endl;
        }
        else
        {
            if (_debug) cout <<  0 << "\t efficiency = " << 0 <<endl;
        }
    }
}

bool findmatch(const  edm::Handle<DTRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual)
{
    bool matchfound = false;
    double maxres = maxdist;

    for (DTRecHitCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = theG->idToDet(hit->geographicalId());
        //const GlobalPoint &p = whichdet->surface().toGlobal(hit->localPosition());

        if (hit->geographicalId().rawId() == detid && maxres>fabs(locx-hit->localPosition().x())) 
        {
            matchfound = true;
            residual = locx-hit->localPosition().x();
            maxres = fabs(residual);
        }
    }

    return matchfound;
}

bool findmatch(const edm::Handle<CSCRecHit2DCollection> &hitcoll, int detid, float locx, float locy, float &residual)
{
    bool matchfound = false;
    double maxres = maxdist;

    for (CSCRecHit2DCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = theG->idToDet(hit->geographicalId());
        //const GlobalPoint &p = whichdet->surface().toGlobal(hit->localPosition());

        if (hit->geographicalId().rawId() == detid && maxres>fabs(locx-hit->localPosition().x())) 
        {
            matchfound = true;
            residual = locx-hit->localPosition().x();
            maxres = fabs(residual);
        }
    }

    return matchfound;
}

bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum)
{
    bool matchfound = false;
    double maxres = maxdist;

    for (RPCRecHitCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = theG->idToDet(hit->geographicalId());
        //const GlobalPoint &p = whichdet->surface().toGlobal(hit->localPosition());

        if (hit->geographicalId().rawId() == detid && maxres>fabs(locx-hit->localPosition().x())) 
        {
            matchfound = true;
            residual = locx-hit->localPosition().x();
            maxres = fabs(residual);
            sum = locx+hit->localPosition().x();
        }
    }

    return matchfound;
}
//define this as a plug-in
DEFINE_FWK_MODULE(CosmicAnalyzer);
