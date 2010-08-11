#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/CommonTopologies/interface/RectangularStripTopology.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include "FWCore/Framework/interface/ESHandle.h"
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include <DataFormats/RPCRecHit/interface/RPCRecHitCollection.h>
//#include "RecoLocalMuon/RPCRecHit/interface/DTSegtoRPC.h"
//#include "RecoLocalMuon/RPCRecHit/interface/CSCSegtoRPC.h"
#include <DataFormats/DetId/interface/DetId.h>
#include <RecoLocalMuon/RPCRecHit/interface/TracktoRPC.h>
#include <ctime>
#include <TMath.h>

ObjectMap2* ObjectMap2::mapInstance = NULL;

ObjectMap2* ObjectMap2::GetInstance(const edm::EventSetup& iSetup){
  if (mapInstance == NULL){
    mapInstance = new ObjectMap2(iSetup);
  }
  return mapInstance;
}

ObjectMap2::ObjectMap2(const edm::EventSetup& iSetup){
  edm::ESHandle<RPCGeometry> rpcGeo;
  edm::ESHandle<DTGeometry> dtGeo;
  
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  iSetup.get<MuonGeometryRecord>().get(dtGeo);
  
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
    if(dynamic_cast< RPCChamber* >( *it ) != 0 ){
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
      std::vector< const RPCRoll*> roles = (ch->rolls());
      for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
	RPCDetId rpcId = (*r)->id();
	int region=rpcId.region();
	if(region==0){
	  int wheel=rpcId.ring();
	  int sector=rpcId.sector();
	  int station=rpcId.station();
	  DTStationIndex2 ind(region,wheel,sector,station);
	  std::set<RPCDetId> myrolls;
	  if (rollstoreDT.find(ind)!=rollstoreDT.end()) myrolls=rollstoreDT[ind];
	  myrolls.insert(rpcId);
	  rollstoreDT[ind]=myrolls;
	}
      }
    }
  }
}

int distsector2(int sector1,int sector2){
  if(sector1==13) sector1=4;
  if(sector1==14) sector1=10;
  
  if(sector2==13) sector2=4;
  if(sector2==14) sector2=10;
  
  int distance = abs(sector1 - sector2);
  if(distance>6) distance = 12-distance;
  return distance;
}

int distwheel2(int wheel1,int wheel2){
  int distance = abs(wheel1 - wheel2);
  return distance;
}
ObjectMap2CSC* ObjectMap2CSC::mapInstance = NULL;

ObjectMap2CSC* ObjectMap2CSC::GetInstance(const edm::EventSetup& iSetup){
  if (mapInstance == NULL){
    mapInstance = new ObjectMap2CSC(iSetup);
  }
  return mapInstance;
}

ObjectMap2CSC::ObjectMap2CSC(const edm::EventSetup& iSetup){
  edm::ESHandle<RPCGeometry> rpcGeo;
  edm::ESHandle<CSCGeometry> cscGeo;
  
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  iSetup.get<MuonGeometryRecord>().get(cscGeo);
  
  for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
    if(dynamic_cast< RPCChamber* >( *it ) != 0 ){
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
      std::vector< const RPCRoll*> roles = (ch->rolls());
      for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
	RPCDetId rpcId = (*r)->id();
	int region=rpcId.region();
	if(region!=0){
	  int station=rpcId.station();
          int ring=rpcId.ring();
          int cscring=ring;
          int cscstation=station;
	  RPCGeomServ rpcsrv(rpcId);
	  int rpcsegment = rpcsrv.segment();
	  int cscchamber = rpcsegment; //FIX THIS ACCORDING TO RPCGeomServ::segment()Definition
          if((station==2||station==3)&&ring==3){//Adding Ring 3 of RPC to the CSC Ring 2
            cscring = 2;
          }
	  CSCStationIndex2 ind(region,cscstation,cscring,cscchamber);
          std::set<RPCDetId> myrolls;
	  if (rollstoreCSC.find(ind)!=rollstoreCSC.end()) myrolls=rollstoreCSC[ind];
	  myrolls.insert(rpcId);
          rollstoreCSC[ind]=myrolls;
	}
      }
    }
  }
}

bool TracktoRPC::ValidRPCSurface(RPCDetId rpcid, LocalPoint LocalP, const edm::EventSetup& iSetup)
{
  edm::ESHandle<RPCGeometry> rpcGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);

  const GeomDet *whichdet3 = rpcGeo->idToDet(rpcid.rawId());
  const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet3);
  float locx=LocalP.x(), locy=LocalP.y();//, locz=LocalP.z();

 cout << "test validregion---- ;";

  if(aroll->isBarrel())
  {
     const Bounds &rollbound = rpcGeo->idToDet((rpcid))->surface().bounds();
     float boundlength = rollbound.length()-10.;
     float boundwidth = rollbound.width()-10.;

     if(fabs(locx) < boundwidth/2 && fabs(locy) < boundlength/2 && locy > -boundlength/2) return true;
     else return false;

   }
   else if(aroll->isForward())
   {
     const Bounds &rollbound = rpcGeo->idToDet((rpcid))->surface().bounds();
     float boundlength = rollbound.length()-.10;
     float boundwidth = rollbound.width()-10.;

     float nminx = TMath::Pi()*(18*boundwidth/ TMath::Pi() - boundlength)/18;
     float ylimit = ((boundlength)/(boundwidth/2 - nminx/2))*fabs(locx) + boundlength/2 - ((boundlength)/(boundwidth/2 - nminx/2))*(boundwidth/2);
     if(ylimit < -boundlength/2 ) ylimit = -boundlength/2;

//     if(fabs(locx) < boundwidth/2 && fabs(locy) < boundlength/2 && locy > ylimit) return true;
     if(fabs(locx) < boundwidth/2 && fabs(locy) < boundlength/2 ) return true;
     else return false;
   } else return false;
}

TracktoRPC::TracktoRPC(edm::Handle<reco::TrackCollection> alltracks, const edm::EventSetup& iSetup,const edm::Event& iEvent,bool debug,const edm::ParameterSet& iConfig,edm::InputTag& tracklabel){ 

 _ThePoints = new RPCRecHitCollection();
// if(alltracks->empty()) return;

 if(tracklabel.label().find("cosmic")==0) theTrackTransformer = new TrackTransformerForCosmicMuons(iConfig);
 else if(tracklabel.label().find("globalCosmic")==0) theTrackTransformer = new TrackTransformerForCosmicMuons(iConfig);
 else theTrackTransformer = new TrackTransformer(iConfig);
 theTrackTransformer->setServices(iSetup);  

 edm::ESHandle<RPCGeometry> rpcGeo;
 edm::ESHandle<DTGeometry> dtGeo;
 edm::ESHandle<CSCGeometry> cscGeo;
 
 iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",thePropagator); 
 iSetup.get<MuonGeometryRecord>().get(rpcGeo);
 iSetup.get<MuonGeometryRecord>().get(dtGeo);
 iSetup.get<MuonGeometryRecord>().get(cscGeo);

std::vector<uint32_t> rpcput;
double MaxD=9999.;

for (TrackCollection::const_iterator track = alltracks->begin(); track !=alltracks->end(); track++)
{
 Trajectories trajectories = theTrackTransformer->transform(*track);
 if(debug) cout << "Building Trajectory from Track. " << endl;
 if(track->numberOfValidHits()<5) continue;
 if(track->normalizedChi2()>10) continue;

 std::vector<uint32_t> rpcrolls;
 std::vector<uint32_t> rpcrolls2; 
 std::map<uint32_t, int> rpcNdtcsc;
 std::map<uint32_t, int> rpcrollCounter;

if(debug) cout << "1. Search expeted RPC roll detid !!" << endl;
for(trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
 {
    if((*hit)->isValid())
    {
      DetId id = (*hit)->geographicalId();

       if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::DT)
       {
          const GeomDet *geomDet =  dtGeo->idToDet((*hit)->geographicalId());
          const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(geomDet);
         if(dtlayer) for(Trajectories::const_iterator trajectory = trajectories.begin(); trajectory != trajectories.end(); ++trajectory)
          {
             const BoundPlane & DTSurface = dtlayer->surface();
             const GlobalPoint dcPoint = DTSurface.toGlobal(LocalPoint(0.,0.,0.));

             TrajectoryMeasurement tMt = trajectory->closestMeasurement(dcPoint);
             TrajectoryStateOnSurface upd2 = (tMt).updatedState();
             if(upd2.isValid())
             {
                LocalPoint trajLP = upd2.localPosition();
                LocalPoint trackLP = (*hit)->localPosition();
                float dx = trajLP.x()-trackLP.x(), dy=trajLP.y()-trackLP.y();//, dz=trajLP.z()-trackLP.z();
               // if( dx>10. || dy>10. ) continue;

                DTChamberId dtid(geomDet->geographicalId().rawId());
                int dtW=dtid.wheel(), dtS=dtid.sector(), dtT=dtid.station();
                if(dtS==13) dtS=4; if(dtS==14) dtS=10;
                ObjectMap2* TheObject = ObjectMap2::GetInstance(iSetup);
                DTStationIndex2 theindex(0,dtW,dtS,dtT);
                std::set<RPCDetId> rollsForThisDT = TheObject->GetInstance(iSetup)->GetRolls(theindex);
                for(std::set<RPCDetId>::iterator iteraRoll = rollsForThisDT.begin();iteraRoll != rollsForThisDT.end(); iteraRoll++)
                {
	            const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
	            RPCDetId rpcId = rollasociated->id();

                    bool check = true;
                    vector<uint32_t>::iterator rpcput_;
                    for( rpcput_ = rpcput.begin() ; rpcput_ < rpcput.end(); rpcput_++ )
                    if(rollasociated->id().rawId()==*rpcput_) check = false;
                    if(!check) continue;
                    
                    TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(rollasociated->id())->surface());
                    if(ptss.isValid()) if(ValidRPCSurface(rpcId, ptss.localPosition(), iSetup))
                    {
                        float rpcGPX = ptss.globalPosition().x();
                        float rpcGPY = ptss.globalPosition().y();
                        float rpcGPZ = ptss.globalPosition().z();
                     
                        const GeomDet *geomDet2 = rpcGeo->idToDet(rollasociated->id().rawId());
                        const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(geomDet2);
                        const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(aroll->topology()));
                        LocalPoint xmin = top_->localPosition(0.);
                        LocalPoint xmax = top_->localPosition((float)aroll->nstrips());
                        float rsize = fabs( xmax.x()-xmin.x() );
                        float stripl = top_->stripLength();
                        //float stripw = top_->pitch();
                        float eyr=1;
                     
                        float locx = ptss.localPosition().x(), locy = ptss.localPosition().y(), locz= ptss.localPosition().z();
                        if( locx < rsize*eyr && locy < stripl*eyr && locz < 1. )
                        {
                           RPCRecHit RPCPoint(rollasociated->id().rawId(),0,LocalPoint(locx,locy,locz));
                     
                           RPCGeomServ servId(rollasociated->id().rawId());
                           if(debug) cout << "3\t Barrel Expected RPC " << servId.name().c_str() <<
                             " \tLocalposition X: " << locx << ", Y: "<< locy << " GlobalPosition(x,y,z) (" << rpcGPX <<", "<< rpcGPY <<", " << rpcGPZ << ")"<< endl;
                           RPCPointVector.clear();
                           RPCPointVector.push_back(RPCPoint);
                           _ThePoints->put(rollasociated->id().rawId(),RPCPointVector.begin(),RPCPointVector.end());
                           rpcput.push_back(rollasociated->id().rawId());
                        }
                    }
	      	}                                 
	     }
	  }
       }
       else if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::CSC) 
       {
          const GeomDet *geomDet =  cscGeo->idToDet((*hit)->geographicalId());
          const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(geomDet);

          CSCDetId cscid(geomDet->geographicalId().rawId());
          if(csclayer) for(Trajectories::const_iterator trajectory = trajectories.begin(); trajectory != trajectories.end(); ++trajectory)
          {
             const BoundPlane & CSCSurface = csclayer->surface();
             const GlobalPoint dcPoint = CSCSurface.toGlobal(LocalPoint(0.,0.,0.));

             TrajectoryMeasurement tMt = trajectory->closestMeasurement(dcPoint);
             TrajectoryStateOnSurface upd2 = (tMt).updatedState();

             if(upd2.isValid() && cscid.station()!=4 && cscid.ring()!=1 )
             {
                LocalPoint trajLP = upd2.localPosition();
                LocalPoint trackLP = (*hit)->localPosition();
                float dx = trajLP.x()-trackLP.x(), dy=trajLP.y()-trackLP.y();//, dz=trajLP.z()-trackLP.z();
               // if( dx>10. || dy>10.) continue;

                ObjectMap2CSC* TheObjectCSC = ObjectMap2CSC::GetInstance(iSetup);
	        int En = cscid.endcap(), St = cscid.station();// Ri = cscid.ring();
	        int rpcSegment = cscid.chamber(), rpcSegment3;
                if(En==2) En= -1; //if(Ri==4) Ri =1; 

                for(int Ri= 2;Ri<4;Ri++)
                {
                  for(int rpcSegment2=rpcSegment-1;rpcSegment2<rpcSegment+2;rpcSegment2++)
                  {
                     rpcSegment3=rpcSegment2;
                     if(rpcSegment2<1) rpcSegment3=36;
                     if(rpcSegment2>36) rpcSegment3=1;

                     CSCStationIndex2 theindex(En,St,Ri,rpcSegment3);
                     std::set<RPCDetId> rollsForThisCSC = TheObjectCSC->GetInstance(iSetup)->GetRolls(theindex);
                     for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++)
                     {
	                 const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
	                 RPCDetId rpcId = rollasociated->id();

                         bool check = true;
                         vector<uint32_t>::iterator rpcput_;
                         for( rpcput_ = rpcput.begin() ; rpcput_ < rpcput.end(); rpcput_++ )
                         if(rollasociated->id().rawId()==*rpcput_) check = false;
                         if(!check) continue;
                     
                         TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(rollasociated->id())->surface());
                         if(ptss.isValid()) if(ValidRPCSurface(rpcId, ptss.localPosition(), iSetup))
                         {
                              float rpcGPX = ptss.globalPosition().x();
                              float rpcGPY = ptss.globalPosition().y();
                              float rpcGPZ = ptss.globalPosition().z();
                              
                              RPCDetId rpcid(rollasociated->id().rawId());
                              const GeomDet *geomDet3 = rpcGeo->idToDet(rollasociated->id().rawId());
                              const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(geomDet3);
                              const TrapezoidalStripTopology* top_=dynamic_cast<const TrapezoidalStripTopology*>(&(aroll->topology()));
                              LocalPoint xmin = top_->localPosition(0.);
                              LocalPoint xmax = top_->localPosition((float)aroll->nstrips());
                              float rsize = fabs( xmax.x()-xmin.x() );
                              float stripl = top_->stripLength();
                              //float stripw = top_->pitch();
                              
                              float eyr=1;
                              float locx = ptss.localPosition().x(), locy = ptss.localPosition().y(), locz= ptss.localPosition().z();
                              if( locx < rsize*eyr && locy < stripl*eyr && locz < 1. )
                              {
                                 RPCRecHit RPCPoint(rollasociated->id().rawId(),0,LocalPoint(locx,locy,locz));
                                 RPCGeomServ servId(rollasociated->id().rawId());
                                 if(debug) cout << "3\t Forward Expected RPC " << servId.name().c_str() <<
                                   " \tLocalposition X: " << locx << ", Y: "<< locy << " GlobalPosition(x,y,z) (" << rpcGPX <<", "<< rpcGPY <<", " << rpcGPZ << ")"<< endl;
                                 RPCPointVector.clear();
                                 RPCPointVector.push_back(RPCPoint);
                                 _ThePoints->put(rollasociated->id().rawId(),RPCPointVector.begin(),RPCPointVector.end());
                                 rpcput.push_back(rollasociated->id().rawId());
                              
                              }
///////////////////////////////
                         }
	      	     }
                  } 
                }
             }
          }
       } else { if(debug) cout << "1\t The hit is not DT/CSC's.   " << endl;} 
    }
 }
}
}

TracktoRPC::~TracktoRPC(){
}
