// -*- C++ -*-
//
// Package:    TrajectoryRPCEff
// Class:      TrajectoryRPCEff
//
/**\class TrajectoryRPCEff TrajectoryRPCEff.cc sccosmic/TrajectoryRPCEff/src/TrajectoryRPCEff.cc

 Description: <one line class summary>

// This code is run with : Tracking Tools(http://higgs.skku.ac.kr/CMS/TrackingTools_CMSSW_2_1_10.tgz) 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Su Yong Choi
//         Created:  Wed Jul  9 17:03:35 CEST 2008
// $Id: TrajectoryRPCEff.cc,v 1.10 2009/07/17 11:59:39 youngjo Exp $
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
bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum, float &residualyy, float &recX, float &recY);
bool findmatch(const edm::Handle<CSCRecHit2DCollection> &hitcoll, int detid, float locx, float locy, float &residual);

double maxdist = 7777777.0;

edm::ESHandle<GlobalTrackingGeometry> theG;

class TrajectoryRPCEff : public edm::EDAnalyzer
{
    public:
        explicit TrajectoryRPCEff(const edm::ParameterSet&);
        ~TrajectoryRPCEff();

    protected:
        void printTrajectoryRecHits(const Trajectory &,
                             edm::ESHandle<GlobalTrackingGeometry>, edm::Handle<RPCRecHitCollection> ) ;
        bool TrajectoryclosestMeasurement(const Trajectory &,   GlobalPoint , float &, float &, int & ) ;
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
        TH1F *htrtjdiff; 
        TH2F *htrjtrachits;
        TH1F *hsize; 

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

}


TrajectoryRPCEff::~TrajectoryRPCEff()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void TrajectoryRPCEff::printTrajectoryRecHits(const Trajectory &trajectory, 
                                              ESHandle<GlobalTrackingGeometry> theG, edm::Handle<RPCRecHitCollection> allRPChits) {

  TransientTrackingRecHit::ConstRecHitContainer rechits = trajectory.recHits();
//  cout << "Size of the RecHit container: " << rechits.size();

//  for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
//    if ((*hit)->isValid())
//  {
  std::map<int, int> rpchitcheck;

  for(TransientTrackingRecHit::ConstRecHitContainer::const_iterator recHit = rechits.begin();
      recHit != rechits.end(); ++recHit)
    if((*recHit)->isValid())
  {
      const GeomDet* geomDet = theG->idToDet((*recHit)->geographicalId());
//      double r = geomDet->surface().position().perp();
//      double z = geomDet->toGlobal((*recHit)->localPosition()).z();

  //    cout << " trajlogs - r:" << r << " Global position:" << geomDet->toGlobal((*recHit)->localPosition()) << endl;

      const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(geomDet);
      const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(geomDet);
      const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(geomDet);

      if(dtlayer || csclayer)
      {
            int detid; 
	    float locx, locy;
          //  std::cout << " trajlogs DT trj:  r:" << r << " Global position:" << geomDet->toGlobal((*recHit)->localPosition()) << std::endl;
/*            if(TrajectoryclosestMeasurement(trajectory,geomDet->toGlobal((*recHit)->localPosition()),locx,locy, detid))
            {
               if(rpchitcheck[detid]<1)
               {
                   RPCGeomServ servId(detid);
                   cout << "---tra"<<  servId.name()<<"-- " << detid <<", "<< locx << ", " << locy << endl;
               
                   SetFolderMuonDir(detid);
                   rpchitcheck[detid]++;
    
                   detrtMap_[detid]->Fill(locx, locy);
                   float residual, sum, residualy, recX, recY;
                   if (findmatch(allRPChits, detid, locx, locy, residual, sum, residualy, recX, recY))
                   {
                      cout << "residual : "<< residual << endl;
                      detrmMap_[detid]->Fill(locx, locy);
                      detrrMap_[detid]->Fill(residual);
                   }  
               }
            }
*/
      }
      if(aroll)
      {
          RPCGeomServ servId((*recHit)->geographicalId());
          std::cout << " trajlogs RPC trj: " << servId.name() << " Global position:" << geomDet->toGlobal((*recHit)->localPosition()) << std::endl;

      }
    
  }
}
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
                                            GlobalPoint GlPt , float &locx, float &locy, int &detid) 
{
     TrajectoryMeasurement  tMt = trajectory.closestMeasurement(GlPt);
     TrajectoryStateOnSurface upd2 = (tMt).updatedState();
     double tx=upd2.globalPosition().x();
     double ty=upd2.globalPosition().y();
     double tz=upd2.globalPosition().z();
////////////////
     PropagationDirection pdirec = trajectory.direction();
//   SteppingHelixPropagator a = SteppingHelixPropagator(const MagneticField* field, dir):

////////////////

     const TransientTrackingRecHit *hit = &(*tMt.recHit());

  //  cout << "TrajStateOnSur: updi2 glo posi " << " (x,y,z) : "<< tx <<", "<<ty<<", " << tz << endl;
 //   cout << "local position :(x,y,z) "<< upd2.localPosition().x() << ", "<< upd2.localPosition().y() << ", "<< upd2.localPosition().z() << endl;
//////////////////////////
     double gDxyz=50, gDxyz_;
     int mrpcid=0;

     std::set<DetId> detidS = detidAsso->getDetIdsCloseToAPoint( upd2.globalPosition(), 0.5);
//     std::set<DetId> detidS = nearRPCChamberFinder();
     for(std::set<DetId>::const_iterator  detid2 = detidS.begin(); detid2 != detidS.end(); ++detid2)
     {
         DetId id(*detid2);

         const GeomDet *whichdet = theG->idToDet(id.rawId());
         const RPCChamber *rpcChamber = dynamic_cast<const RPCChamber *>(whichdet);

         if(rpcChamber)
         {
          //      cout << "closet hit is RPC Chamber "<<  id.rawId()  << endl;
                RPCGeomServ servId(id.rawId());

                std::vector< const RPCRoll*> roles = (rpcChamber->rolls());

                const GlobalPoint &p = whichdet->surface().toGlobal(upd2.localPosition());

          //      cout << "RPC Chamber : " << servId.name() << "(x,y,z)" << p.x() << ", " << p.y() << ", " << p.z() << endl;

                for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r)
                {

                     RPCDetId rpcId = (*r)->id();

                     RPCGeomServ rpcsrv(rpcId);
                     std::string nameRoll = rpcsrv.name();

                     const GeomDet *whichdet1 = theG->idToDet(rpcId.rawId());
                     const GlobalPoint &p1 = whichdet1->surface().toGlobal(LocalPoint(0,0,0));
                     double rx=p1.x(); double ry=p1.y(); double rz=p1.z();

            //         cout << "rpc roll : "<< nameRoll << " (x,y,z) "<< rx <<", "<<ry<<", " << rz << endl;
                     gDxyz_ = sqrt((tz-rz)*(tz-rz)+(tx-rx)*(tx-rx)+(ty-ry)*(ty-ry));
            //         cout << "dist : " << gDxyz_ << endl;

                     if(gDxyz>gDxyz_)
                     {
                          mrpcid = (int)(*r)->id();
                          gDxyz=gDxyz_;
                     }

                }
         
         }
     }  
//      cout << "---" << endl;
     if(mrpcid!=0)
     {
         TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, theG->idToDet(mrpcid)->surface());
         double px=ptss.globalPosition().x();
         double py=ptss.globalPosition().y();
         double pz=ptss.globalPosition().z();

        // cout << "TrajStateOnSur: ptss glo posi " << " (x,y,z) : "<< px <<", "<<py<<", " << pz << endl;
        // cout << "local position :(x,y,z) "<< ptss.localPosition().x() << ", "<< ptss.localPosition().y() << ", "<< ptss.localPosition().z() << endl;

      if(fabs(ptss.localPosition().x())<100 && fabs(ptss.localPosition().y())<100)
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

         return true;
       }
       else
       {
       // cout << " I can not match a detetor!! ";
        detid = 0;
        locx= 0;
        locy= 0;
        return false;
       }
     }
     else
     {
       // cout << " I can not match a detetor!! ";
        detid = 0;
        locx= 0;
        locy= 0;
        return false;
     }

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
          int lsize= 100;

          RPCGeomServ servId(detid);
          if (_debug) cout << "RPCGeomServ : " << servId.name() << " :: " << endl;


          TFileDirectory subDir_RPC = fs->mkdir( "RPC" );
          TFileDirectory dirRegion = subDir_RPC.mkdir(Form("Region_%d", region));
          TFileDirectory dirRing = dirRegion.mkdir(Form("Ring_%d", ring));
          TFileDirectory dirStation = dirRing.mkdir(Form("Station_%d", station));
          TFileDirectory dirSector = dirStation.mkdir( Form("Sector_%d",sector));

           detrtMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrt%d", detid),
                                        Form("%s Trajectory hit", servId.name().c_str()),
                                   25, -lsize, lsize, 25, -lsize, lsize);

           detrmMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrm%d", detid),
                                        Form("%s Trajectory matched hit", servId.name().c_str()),
                                   25, -lsize, lsize, 25, -lsize, lsize);
    
           detrrMap_[detid] =  dirSector.make<TH1F>(Form("RPCTtrr%d_Residure", detid),
                                        Form("%s Residure local X", servId.name().c_str()),
                                   400, -20, 20);

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

    if(tracks->empty()) return;
  //  if(tracks->size()<20) return;
    
    // Get the Trajectory collection from the event
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
/*
    for(Trajectories::const_iterator trajectory = trajectories->begin();
       trajectory != trajectories->end(); ++trajectory)
    {
//         printTrajectoryRecHits(*trajectory, theG, allRPChits);
    }


    for(TrajTrackAssociationCollection::const_iterator it = assoMap->begin();
       it != assoMap->end(); ++it)
    {

         vector<Trajectory> traj = it->key;
         const reco::TrackRef tk = it->val;
   
         // Check the difference in Pt
         reco::TransientTrack track(tk,&*theMGField,theG);
   
         int diff = track.recHitsSize()- traj->recHits().size();
   
        // if (_debug)
         cout << "track and trajectory Difference: " << diff << endl;
         htrtjdiff->Fill(diff);
         htrjtrachits->Fill(traj->recHits().size(), track.recHitsSize());

    }
*/
    hsize->Fill(1);

    int detid;
    float locx, locy;
    std::map<int, int> rpccheck;

    for (TrackCollection::const_iterator track = tracks->begin(); track !=tracks->end(); track++)
    {

         Trajectories trajectories = theTrackTransformer->transform(*track);

        for (trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
        {
            if ((*hit)->isValid())
            {
                 const GeomDet *whichdet = theG->idToDet((*hit)->geographicalId());
                 const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(whichdet);
                 const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(whichdet);
                 if ((*hit)->isValid() &&(dtlayer || csclayer) )
                 {
        //             alltrackrechits.push_back(&(*(*hit)));
                       for(Trajectories::const_iterator trajectory = trajectories.begin();
                          trajectory != trajectories.end(); ++trajectory)
                       {
                            if(TrajectoryclosestMeasurement(*trajectory, whichdet->toGlobal((*hit)->localPosition()), locx, locy, detid))
                            {
                               if(rpccheck[detid]<1)
                               {
                                    RPCGeomServ servId(detid);
                                    cout << "---tra"<<  servId.name()<<"-- " << detid <<", "<< locx << ", " << locy << endl;
                                    SetFolderMuonDir(detid);
                                    rpccheck[detid]++;
                       
                                    detrtMap_[detid]->Fill(locx, locy);
                                    float residual, sum, residualy, recX, recY;
                                    if (findmatch(allRPChits, detid, locx, locy, residual, sum, residualy, recX, recY))
                                    {
                                       cout << "residual : "<< residual << endl;
                                       detrmMap_[detid]->Fill(locx, locy);
                                       detrrMap_[detid]->Fill(residual);
                                    }
                               }
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

bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum,float &residualy, float &recX, float &recY)
{
    bool matchfound = false;
    double maxres = maxdist;

    for (RPCRecHitCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = theG->idToDet(hit->geographicalId());
        //const GlobalPoint &p = whichdet->surface().toGlobal(hit->localPosition());
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
        }
    }

    return matchfound;
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrajectoryRPCEff);
