// -*- C++ -*-
//
// Package:    RPCPointProducer
// Class:      RPCPointProducer
// 
/**\class RPCPointProducer RPCPointProducer.cc Analysis/RPCPointProducer/src/RPCPointProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Camilo Andres Carrillo Montoya
//         Created:  Wed Sep 16 14:56:18 CEST 2009
// $Id: RPCPointProducer.cc,v 1.1 2010/01/25 21:03:56 carrillo Exp $
//
//

#include "RecoLocalMuon/RPCRecHit/interface/RPCPointProducer.h"

// system include files

#include <memory>
#include <ctime>

// user include files

RPCPointProducer::RPCPointProducer(const edm::ParameterSet& iConfig)
{
  dt4DSegments=iConfig.getUntrackedParameter<std::string>("dt4DSegments","dt4DSegments");
  cscSegments=iConfig.getUntrackedParameter<std::string>("cscSegments","cscSegments");
  tracks=iConfig.getUntrackedParameter<std::string>("tracks","standAloneMuons");

  debug=iConfig.getUntrackedParameter<bool>("debug",false);
  incldt=iConfig.getUntrackedParameter<bool>("incldt",true);
  inclcsc=iConfig.getUntrackedParameter<bool>("inclcsc",true);
  incltrack=iConfig.getUntrackedParameter<bool>("incltrack",true);
//  checksegment==iConfig.getUntrackedParameter<bool>("checksegment",true);

  MinCosAng=iConfig.getUntrackedParameter<double>("MinCosAng",0.95);
  MaxD=iConfig.getUntrackedParameter<double>("MaxD",80.);
  MaxDrb4=iConfig.getUntrackedParameter<double>("MaxDrb4",150.);
  MaxDistanceBetweenSegments=iConfig.getUntrackedParameter<double>("",150.);
  ExtrapolatedRegion=iConfig.getUntrackedParameter<double>("ExtrapolatedRegion",0.5);

  produces<RPCRecHitCollection>("RPCDTExtrapolatedPoints");
  produces<RPCRecHitCollection>("RPCCSCExtrapolatedPoints");
  produces<RPCRecHitCollection>("RPCTrackExtrapolatedPoints"); 
  trackTransformerParam = iConfig.getParameter<ParameterSet>("TrackTransformer");
}


RPCPointProducer::~RPCPointProducer(){

}
//void RPCPointProducer::setServices(edm::Event& iEvent, const edm::EventSetup& iSetup){
//      produce( iEvent, iSetup);        
//}    
void RPCPointProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  /*
  struct timespec start_time, stop_time;
  time_t fs;
  time_t fn;
  time_t ls;
  time_t ln;
  clock_gettime(CLOCK_REALTIME, &start_time);  
  */

  if(incldt){
    edm::Handle<DTRecSegment4DCollection> all4DSegments;
    iEvent.getByLabel(dt4DSegments, all4DSegments);
    if(all4DSegments.isValid()){
      DTSegtoRPC DTClass(all4DSegments,iSetup,iEvent,debug,ExtrapolatedRegion);
      std::auto_ptr<RPCRecHitCollection> TheDTPoints(DTClass.thePoints());     
      iEvent.put(TheDTPoints,"RPCDTExtrapolatedPoints"); 
    }else{
      std::cout<<"RPCHLT Invalid DTSegments collection"<<std::endl;
    }
  }

  if(inclcsc){
    edm::Handle<CSCSegmentCollection> allCSCSegments;
    iEvent.getByLabel(cscSegments, allCSCSegments);
    if(allCSCSegments.isValid()){
      CSCSegtoRPC CSCClass(allCSCSegments,iSetup,iEvent,debug,ExtrapolatedRegion);
      std::auto_ptr<RPCRecHitCollection> TheCSCPoints(CSCClass.thePoints());  
      iEvent.put(TheCSCPoints,"RPCCSCExtrapolatedPoints"); 
    }else{
      std::cout<<"RPCHLT Invalid CSCSegments collection"<<std::endl;
    }
  }
  if(incltrack){
    edm::Handle<reco::TrackCollection> alltracks;
    iEvent.getByLabel(tracks,alltracks);

    if(!(alltracks->empty())){
      TracktoRPC TrackClass(alltracks,iSetup,iEvent,debug,trackTransformerParam,tracks);
      std::auto_ptr<RPCRecHitCollection> TheTrackPoints(TrackClass.thePoints());
      iEvent.put(TheTrackPoints,"RPCTrackExtrapolatedPoints");
    }else{
      std::cout<<"RPCHLT Invalid Tracks ("<<tracks.c_str()<<") collection"<<std::endl;
    }

  /*  edm::Handle<CSCSegmentCollection> allCSCSegments;
    iEvent.getByLabel(cscSegments, allCSCSegments);
    edm::Handle<DTRecSegment4DCollection> all4DSegments;
    iEvent.getByLabel(dt4DSegments, all4DSegments);

    if(!(alltracks->empty()) && (all4DSegments.isValid()) && (allCSCSegments.isValid())){
      TracktoRPC TrackClass(alltracks,iSetup,iEvent,debug,trackTransformerParam,checksegment,all4DSegments,allCSCSegments );
      std::auto_ptr<RPCRecHitCollection> TheTrackPoints(TrackClass.thePoints());
      iEvent.put(TheTrackPoints,"RPCTrackExtrapolatedPoints");
    }else{
      std::cout<<"RPCHLT Invalid Tracks ("<<tracks.c_str()<<") collection"<<std::endl; 
    } 
  */
  }
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
RPCPointProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RPCPointProducer::endJob() {
}


