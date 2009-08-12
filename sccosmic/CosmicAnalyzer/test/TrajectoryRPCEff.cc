// -*- C++ -*-
//
// Package:    TrajectoryRPCEff
// Class:      TrajectoryRPCEff
//
/**\class TrajectoryRPCEff TrajectoryRPCEff.cc sccosmic/CosmicAnalyzer/test/TrajectoryRPCEff.cc

 Description: <one line class summary>


 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Su Yong Choi
//         Created:  Wed Aug 5 19:02:58 2009 UTC (38 hours, 53 minutes ago) by youngjo 
// $Id: TrajectoryRPCEff.cc,v 1.1 2009/08/05 19:02:58 youngjo Exp $
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

#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TrackingTools/TrackAssociator/interface/DetIdAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TAMuonChamberMatch.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfoCollection.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"

#include "RecoMuon/DetLayers/interface/MuonDetLayerGeometry.h"
#include "RecoMuon/Records/interface/MuonRecoGeometryRecord.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerBase.h"

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
typedef std::vector<Trajectory> Trajectories;

bool findmatch(const  edm::Handle<DTRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual);
bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum, float &residualyy, float &recX, float &recY,  int &mstrip);
bool findmatch(const edm::Handle<CSCRecHit2DCollection> &hitcoll, int detid, float locx, float locy, float &residual);

double maxdist = 7777777.0;

edm::ESHandle<GlobalTrackingGeometry> theG;

class TrajectoryRPCEff : public edm::EDAnalyzer
{
    public:
        explicit TrajectoryRPCEff(const edm::ParameterSet&);
        ~TrajectoryRPCEff();

    protected:
        bool TrajectoryclosestMeasurement(const Trajectory &,   GlobalPoint , float &, float &, int &, int ) ;
        bool SetFolderMuonDir(int );
        std::set<DetId> nearRPCChamberFinder();

    private:
        virtual void beginJob(const edm::EventSetup&) ;
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        edm::InputTag theInputLabel;
        TrackTransformerBase *theTrackTransformer;

        MuonServiceProxy *theMuonService_;
        edm::ESHandle<Propagator> thePropagator;

/// the name of the DT rec hits collection
        edm::InputTag theDTRecSegmentLabel;
/// the name of the CSC rec hits collection
        edm::InputTag theCSCRecSegmentLabel;
/// the name of the RPC rec hits collection
        edm::InputTag theRPCRecSegmentLabel;

        bool _debug;

        std::map<int, TH2F*> detrtMap_;
        std::map<int, TH2F*> detrmMap_;
        std::map<int, TH1F*> detrrMap_;

        std::map<int, TH1F*> detsrtMap_;
        std::map<int, TH1F*> detsrmMap_;

        TH1F *htrtjdiff; 
        TH2F *htrjtrachits;
        TH1F *hsize; 

        TH1F *htracks_pt;
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

    edm::Service<TFileService> fs;
    edm::ESHandle<DetIdAssociator> detidAsso;

};


//
// constructors and destructor
//
TrajectoryRPCEff::TrajectoryRPCEff(const edm::ParameterSet& iConfig)
:theDTRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("DTRecSegmentLabel")),
theCSCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("CSCRecSegmentLabel")),
theRPCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("RPCRecSegmentLabel"))
{
//now do what ever initialization is needed

    MuonServiceProxy *theMuonService_;
    edm::ParameterSet serviceParameters
        = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
    theMuonService_ = new MuonServiceProxy(serviceParameters);

    ParameterSet trackTransformerParam = iConfig.getParameter<ParameterSet>("TrackTransformer");
    theTrackTransformer = new TrackTransformerForCosmicMuons(trackTransformerParam);

    edm::Service<TFileService> fs;
    theInputLabel = iConfig.getParameter<InputTag>("InputLabel");

    htrtjdiff = fs->make<TH1F>("htrtjdiff", "diff counts  of track and trajectory", 50, 0.0, 50.0);
    htrjtrachits = fs->make<TH2F>("htrjtrachits", "hits  of trajectory and track ", 80, 0.0, 80.0, 80, 0.0, 80.0);
    hsize = fs->make<TH1F>("hsize","event count", 4, 0.0, 4.0);

    htracks_pt = fs->make<TH1F>("htracks_pt", "Pt of track", 600, 0, 300);
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

}


TrajectoryRPCEff::~TrajectoryRPCEff()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

std::set<DetId> TrajectoryRPCEff::nearRPCChamberFinder()
{
   std::set<DetId> setOfValidIds;
   // RPC
   if (! theG->slaveGeometry(RPCDetId()) ) throw cms::Exception("FatalError") << "Cannnot RPCGeometry\n";
   std::vector<GeomDet*> geomDetsRPC = theG->slaveGeometry(RPCDetId())->dets();
   for(std::vector<GeomDet*>::const_iterator it = geomDetsRPC.begin(); it != geomDetsRPC.end(); ++it)
     if (RPCChamber* rpc = dynamic_cast< RPCChamber*>(*it)) setOfValidIds.insert(rpc->id());

   return setOfValidIds;
}

bool TrajectoryRPCEff::TrajectoryclosestMeasurement(const Trajectory &trajectory,
                                            GlobalPoint GlPt , float &locx, float &locy, int &detid, int dtcscid) 
{
     TrajectoryMeasurement  tMt = trajectory.closestMeasurement(GlPt);
     TrajectoryStateOnSurface upd2 = (tMt).updatedState();
     if(upd2.isValid())
     { 
          double tx=upd2.globalPosition().x();
          double ty=upd2.globalPosition().y();
          double tz=upd2.globalPosition().z();

          double gDxyz=999., gDxyz_;
          int mrpcid=0;
     
      //    std::set<DetId> detidS = detidAsso->getDetIdsCloseToAPoint( upd2.globalPosition(), 0.5);
          std::set<DetId> detidS = nearRPCChamberFinder();
          for(std::set<DetId>::const_iterator  detid2 = detidS.begin(); detid2 != detidS.end(); ++detid2)
          {
              DetId id(*detid2);
              RPCDetId rpcid(*detid2);     
              const GeomDet *whichdet = theG->idToDet(id.rawId());
              const RPCChamber *rpcChamber = dynamic_cast<const RPCChamber *>(whichdet);

              bool SameState = false;
              const GeomDet *geomDet = theG->idToDet(dtcscid);
              if(DetId(geomDet->geographicalId().rawId()).det() == DetId::Muon && DetId(geomDet->geographicalId().rawId()).subdetId() == MuonSubdetId::DT)
              {
                   DTChamberId chdtid(geomDet->geographicalId().rawId());
                   if (chdtid.wheel() == rpcid.ring() && chdtid.sector() == rpcid.sector() && chdtid.station() == rpcid.station())
                   {
                     SameState = true;
                   }
              }
              else if(DetId(geomDet->geographicalId().rawId()).det() == DetId::Muon && DetId(geomDet->geographicalId().rawId()).subdetId() == MuonSubdetId::CSC)
              {
                    CSCDetId chdtid(geomDet->geographicalId().rawId());
                    if (chdtid.ring() == rpcid.ring() && chdtid.station() == rpcid.station())
                    {
                      SameState = true;
                    }   
              }   
     
              if(rpcChamber && SameState == true)
              {
                     RPCGeomServ servId(id.rawId());
     
                     std::vector< const RPCRoll*> roles = (rpcChamber->rolls());
                     for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r)
                     {
     
                          RPCDetId rpcId = (*r)->id();
     
                          RPCGeomServ rpcsrv(rpcId);
                          std::string nameRoll = rpcsrv.name();
     
                          const GeomDet *whichdet1 = theG->idToDet(rpcId.rawId());
                          TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, theG->idToDet(rpcId)->surface());
                          if(ptss.isValid())
                          {
                                const GlobalPoint &p1 = whichdet1->surface().toGlobal(ptss.localPosition()); // LocalPoint(0,0,0));
                                double rx=p1.x(); double ry=p1.y(); double rz=p1.z();
                              
                 //               cout << "rpc roll : "<< nameRoll << " (x,y,z) "<< rx <<", "<<ry<<", " << rz << endl;
                                gDxyz_ = sqrt((tz-rz)*(tz-rz)+(tx-rx)*(tx-rx)+(ty-ry)*(ty-ry));
                 //               cout << "dist : " << gDxyz_ << endl;
                              
                                if(gDxyz>gDxyz_)
                                {
                                     mrpcid = (int)(*r)->id();
                                     gDxyz=gDxyz_;
                                }
                          } 
                     }
              
              }
          }  
      //      cout << "---" << endl;
          if(mrpcid!=0)
          {
              TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, theG->idToDet(mrpcid)->surface());
              if(ptss.isValid())
              {
                   double px=ptss.globalPosition().x();
                   double py=ptss.globalPosition().y();
                   double pz=ptss.globalPosition().z();
                  
                  // cout << "TrajStateOnSur: ptss glo posi " << " (x,y,z) : "<< px <<", "<<py<<", " << pz << endl;
                  // cout << "local position :(x,y,z) "<< ptss.localPosition().x() << ", "<< ptss.localPosition().y() << ", "<< ptss.localPosition().z() << endl;
                  
                   if(fabs(ptss.localPosition().x())<100 && fabs(ptss.localPosition().y())<50)
                   {
                         RPCDetId rpcid = mrpcid;
                         RPCGeomServ rpcsrv(rpcid);
                         std::string nameRoll = rpcsrv.name();
                     
                         const GeomDet *whichdet = theG->idToDet(rpcid.rawId());
                         const GlobalPoint &p = whichdet->surface().toGlobal(ptss.localPosition());
                         double rx=p.x(); double ry=p.y(); double rz=p.z();
                     
                        // cout << "RPC : "<< nameRoll << " (x,y,z) " << rx << ", " <<ry<< ", " << rz << endl;
                     
                         locx = ptss.localPosition().x();
                         locy = ptss.localPosition().y();
                         detid = rpcid.rawId();
                       //  const float stripPredicted =aroll->strip(LocalPoint(locx,locy,0.));
                     
                         return true;
                    } else return false;
               } else  return false;
          } else  return false;
     }else  return false;

}

bool TrajectoryRPCEff::SetFolderMuonDir(int detid)
{
     const GeomDet *whichdet = theG->idToDet(detid);
     const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
     const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(whichdet);
     const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(whichdet);

     if ( !detrtMap_[detid] ) {
      
      if(aroll)
      { // if RPC hits
          RPCDetId segId(detid);
          int region = segId.region();
          int ring = segId.ring();
          int station = segId. station();
          int sector = segId.sector();
          int layer = segId.layer();
          int subsector = segId.subsector();
          int roll  = segId.roll();
          int lsize= 200, bsize= lsize/4;

          RPCGeomServ servId(detid);
          if (_debug) cout << "RPCGeomServ : " << servId.name() << " :: " << endl;


          TFileDirectory subDir_RPC = fs->mkdir( "RPC" );
          TFileDirectory dirRegion = subDir_RPC.mkdir(Form("Region_%d", region));
          TFileDirectory dirRing = dirRegion.mkdir(Form("Ring_%d", ring));
          TFileDirectory dirStation = dirRing.mkdir(Form("Station_%d", station));
          TFileDirectory dirSector = dirStation.mkdir( Form("Sector_%d",sector));

           detrtMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrt%d", detid),
                                        Form("%s Trajectory hit", servId.name().c_str()),
                                   bsize, -lsize, lsize, bsize, -lsize, lsize);

           detrmMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrm%d", detid),
                                        Form("%s Trajectory matched hit", servId.name().c_str()),
                                   bsize, -lsize, lsize, bsize, -lsize, lsize);
    
           detrrMap_[detid] =  dirSector.make<TH1F>(Form("RPCTtrr%d_Residure", detid),
                                        Form("%s Residure local X", servId.name().c_str()),
                                   400, -20, 20);

           detsrtMap_[detid] = dirSector.make<TH1F>(Form("RPCSrt%d", detid),
                                        Form("%s strips trajectory", servId.name().c_str()),
                                   95, 0.5, 95.5);
           detsrmMap_[detid] = dirSector.make<TH1F>(Form("RPCSrm%d", detid),
                                        Form("%s strips matched trajectory", servId.name().c_str()),
                                   95, 0.5, 95.5);

    }
     return true;
   } else return false;

}

// ------------ method called to for each event  ------------
void TrajectoryRPCEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using reco::MuonCollection;
    using reco::TrackCollection;
    using namespace GeomDetEnumerators;

   theTrackTransformer->setServices(iSetup);

// Get the magnetic field
    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

// Get the global tracking geometry
    iSetup.get<GlobalTrackingGeometryRecord>().get(theG);

    // Get the RecTrack collection from the event
    Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(theInputLabel.label(),tracks);

    hsize->Fill(1);
    if(tracks->empty()) return;
     hntracks->Fill(tracks->size());

//    Get the Trajectory collection from the event
//    Handle<Trajectories> trajectories;
//    iEvent.getByLabel(theInputLabel,trajectories);

//    Handle<TrajTrackAssociationCollection> assoMap;
//    iEvent.getByLabel(theInputLabel,assoMap);

    edm::Handle<reco::MuonTrackLinksCollection> muHandle;
    edm::ESHandle<TransientTrackBuilder> ttrackBuilder;

    edm::Handle<RPCRecHitCollection> allRPChits;
    edm::Handle<DTRecHitCollection> allDThits;
    edm::Handle<CSCRecHit2DCollection> allCSChits;

    edm::ESHandle<MuonDetLayerGeometry> theMuonLayers;
    iSetup.get<MuonRecoGeometryRecord>().get(theMuonLayers);
// get the Muon layers
    vector<DetLayer*> dtLayers = theMuonLayers->allDTLayers();
    vector<DetLayer*> rpcLayers = theMuonLayers->allRPCLayers();

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder);
    iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",thePropagator);
    iSetup.get<DetIdAssociatorRecord>().get("MuonDetIdAssociator", detidAsso);

    iEvent.getByLabel(theRPCRecSegmentLabel, allRPChits);
    iEvent.getByLabel(theDTRecSegmentLabel, allDThits);
    iEvent.getByLabel(theCSCRecSegmentLabel, allCSChits);

    
    //muonMeasurements.setEvent(iEvent);
    MuonTransientTrackingRecHit::MuonRecHitContainer allHits;

    std::vector<const TrackingRecHit *> alltrackrechits;

    int detid;
    float locx, locy;
    std::map<int, int> rpccheck;

    for (TrackCollection::const_iterator track = tracks->begin(); track !=tracks->end(); track++)
    {
        htracks_pt->Fill(track->pt());          //  track transverse momentum

        htracks_etaphi->Fill(track->eta(), track->phi()); //  pseudorapidity and azimuthal angle of momentum vector

        htrackinnerxy->Fill(track->innerPosition().X(), track->innerPosition().Y()); //  position of the innermost hit
        htrackouterxy->Fill(track->outerPosition().X(), track->outerPosition().Y()); //  position of the outermost hit
        htrackinnerz->Fill(track->innerPosition().Z()); // position of the innermost hit
        htrackouterz->Fill(track->outerPosition().Z()); // position of the outermost hit

        const reco::HitPattern& p = track->hitPattern();
        hnvalidhits->Fill(track->innerPosition().phi(), p.numberOfHits()); // azimuthal angle of Innermost hit, number Of Hits

        hchi2prob->Fill(track->innerPosition().phi(), TMath::Prob(track->chi2(), track->ndof())); // 
        hchi2ndof->Fill(track->innerPosition().phi(), track->normalizedChi2()); // chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero)
        hqoverppull->Fill(track->innerPosition().phi(), track->qoverp()/track->qoverpError()); // error on signed transverse curvature
        htrackq->Fill(track->innerPosition().phi(), track->charge()); // track electric charge


        Trajectories trajectories = theTrackTransformer->transform(*track);

        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                 const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                 const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(whichdet);
                 const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(whichdet);

                 const GlobalPoint &p = theG->idToDet((*hit)->geographicalId())->surface().toGlobal((*hit)->localPosition());
                 htrackrechitsxy->Fill(p.x(), p.y());
                 htrackrechitsz->Fill(p.z());

                 if ((*hit)->geographicalId().det()==DetId::Muon && (*hit)->geographicalId().subdetId()==1)
                 {
                     htrackrechitsDTxy->Fill(p.x(), p.y());
                     htrackrechitsDTz->Fill(p.z());
                 }
                 else if ((*hit)->geographicalId().det()==DetId::Muon && (*hit)->geographicalId().subdetId()==2)
                 {
                     htrackrechitsCSCxy->Fill(p.x(), p.y());
                     htrackrechitsCSCz->Fill(p.z());
                 }
                 else if ((*hit)->geographicalId().det()==DetId::Muon && (*hit)->geographicalId().subdetId()==3)
                 {
                     htrackrechitsRPCxy->Fill(p.x(), p.y());
                     htrackrechitsRPCz->Fill(p.z());
                 }
                 if (dtlayer || csclayer )
                 {
             //         alltrackrechits.push_back(&(*(*hit)));
                      for(Trajectories::const_iterator trajectory = trajectories.begin();
                         trajectory != trajectories.end(); ++trajectory)
                      {
                           if(TrajectoryclosestMeasurement(*trajectory, whichdet->toGlobal((*hit)->localPosition()), locx, locy, detid, (*hit)->geographicalId()))
                           {
                            //    if(rpccheck[detid]<1)
                            //    {
                                     RPCGeomServ servId(detid);
                                     cout << "---tra"<<  servId.name()<<"-- " << detid <<", "<< locx << ", " << locy << endl;
                                     SetFolderMuonDir(detid);
                                     rpccheck[detid]++;
                               
                                     detrtMap_[detid]->Fill(locx, locy);

                                     const GeomDet *whichdet1 = theG->idToDet(detid);
                                     const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet1);
                                     const float stripPredicted =aroll->strip(LocalPoint(locx,locy,0.));
                                     detsrtMap_[detid]->Fill(stripPredicted);
                                     int mstrip;
                                     float residual, sum, residualy, recX, recY;
                                     if (findmatch(allRPChits, detid, locx, locy, residual, sum, residualy, recX, recY, mstrip))
                                     {
                                        cout << "residual : "<< residual << endl;
                                        detrmMap_[detid]->Fill(locx, locy);
                                        detrrMap_[detid]->Fill(residual);

                                      if(mstrip == stripPredicted)  detsrmMap_[detid]->Fill(mstrip);

                                     }

                            //    }
                           }
                      }
                
                 }
            }
        }
    }


}


// ------------ method called once each job just before starting event loop  ------------
void
TrajectoryRPCEff::beginJob(const edm::EventSetup&)
{
}


// ------------ method called once each job just after ending the event loop  ------------
void
TrajectoryRPCEff::endJob()
{

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

bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum,float &residualy, float &recX, float &recY, int &mstrip)
{
    bool matchfound = false;
    double maxres = maxdist;

    for (RPCRecHitCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = theG->idToDet(hit->geographicalId());
        //const GlobalPoint &p = whichdet->surface().toGlobal(hit->localPosition());
        const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
        RPCGeomServ servId(detid);

  //     if (hit->geographicalId().rawId()==detid)  cout << "!! " << maxres <<", "<< fabs(locx-hit->localPosition().x()) << " !!matched "<< detid << servId.name();

        if (hit->geographicalId().rawId() == detid && maxres>fabs(locx-hit->localPosition().x())) 
        {
            matchfound = true;
            residual = locx-hit->localPosition().x();
            residualy = locy-hit->localPosition().y();
            maxres = fabs(residual);
            sum = locx+hit->localPosition().x();
            recX = hit->localPosition().x();
            recY = hit->localPosition().y();

            mstrip = aroll->strip(LocalPoint( hit->localPosition().x(), hit->localPosition().y(),0.));

        }
    }

    return matchfound;
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrajectoryRPCEff);
