#ifndef TrackAssociator_TrackDetectorAssociator_h
#define TrackAssociator_TrackDetectorAssociator_h 1

// -*- C++ -*-
//
// Package:    TrackAssociator
// Class:      TrackDetectorAssociator
// 
/*

 Description: main class of tools to associate a track to calorimeter and muon detectors

*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Fri Apr 21 10:59:41 PDT 2006
// $Id: TrackDetectorAssociator.h,v 1.15 2009/03/15 04:22:17 dmytro Exp $
//
//

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "ForSccosmic/TrackAssociator/interface/CaloDetIdAssociator.h"
#include "ForSccosmic/TrackAssociator/interface/EcalDetIdAssociator.h"
#include "ForSccosmic/TrackAssociator/interface/MuonDetIdAssociator.h"
#include "ForSccosmic/TrackAssociator/interface/HcalDetIdAssociator.h"
#include "ForSccosmic/TrackAssociator/interface/HODetIdAssociator.h"

#include "ForSccosmic/TrackAssociator/interface/TrackDetMatchInfo.h"
#include "ForSccosmic/TrackAssociator/interface/CachedTrajectory.h"

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/RecSegment.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "ForSccosmic/TrackAssociator/interface/TrackAssociatorParameters.h"

class TrackDetectorAssociator {
 public:
   TrackDetectorAssociator();
   ~TrackDetectorAssociator();
   
   typedef TrackAssociatorParameters AssociatorParameters;
   enum Direction { Any, InsideOut, OutsideIn };
   
   /// propagate a track across the whole detector and
   /// find associated objects. Association is done in
   /// two modes 
   ///  1) an object is associated to a track only if 
   ///     crossed by track
   ///  2) an object is associated to a track if it is
   ///     withing an eta-phi cone of some radius with 
   ///     respect to a track.
   ///     (the cone origin is at (0,0,0))
   /// Trajectory bending in eta-phi is taking into account
   /// when matching is performed
   ///
   /// associate using FreeTrajectoryState
   TrackDetMatchInfo            associate( const edm::Event&,
					   const edm::EventSetup&,
					   const FreeTrajectoryState&,
					   const AssociatorParameters& );
   /// associate using inner and outer most states of a track
   /// in the silicon tracker. 
   TrackDetMatchInfo            associate( const edm::Event& iEvent,
					   const edm::EventSetup& iSetup,
					   const AssociatorParameters& parameters,
					   const FreeTrajectoryState* innerState,
					   const FreeTrajectoryState* outerState=0);
   /// associate using reco::Track
   TrackDetMatchInfo            associate( const edm::Event&,
					   const edm::EventSetup&,
					   const reco::Track&,
					   const AssociatorParameters&,
					   Direction direction = Any );
   /// associate using a simulated track
   TrackDetMatchInfo            associate( const edm::Event&,
					   const edm::EventSetup&,
					   const SimTrack&,
					   const SimVertex&,
					   const AssociatorParameters& );
   /// associate using 3-momentum, vertex and charge
   TrackDetMatchInfo            associate( const edm::Event&,
					   const edm::EventSetup&,
					   const GlobalVector&,
					   const GlobalPoint&,
					   const int,
					   const AssociatorParameters& );
   
   /// use a user configured propagator
   void setPropagator( const Propagator* );
   
   /// use the default propagator
   void useDefaultPropagator();
   
   /// get FreeTrajectoryState from different track representations
   static FreeTrajectoryState getFreeTrajectoryState( const edm::EventSetup&, 
						      const reco::Track& );
   static FreeTrajectoryState getFreeTrajectoryState( const edm::EventSetup&, 
						      const SimTrack&, 
						      const SimVertex& );
   static FreeTrajectoryState getFreeTrajectoryState( const edm::EventSetup&,
						      const GlobalVector&,
						      const GlobalPoint&,
						      const int);
        
   static bool                crossedIP(const reco::Track& track);

 private:
   DetIdAssociator::MapRange getMapRange( const std::pair<float,float>& delta,
					  const float dR );

   void fillEcal(       const edm::Event&,
			TrackDetMatchInfo&, 
			const AssociatorParameters&);
   
   void fillCaloTowers( const edm::Event&,
			TrackDetMatchInfo&,
			const AssociatorParameters&);
   
   void fillHcal(       const edm::Event&,
			TrackDetMatchInfo&,
			const AssociatorParameters&);
   
   void fillHO(         const edm::Event&,
			TrackDetMatchInfo&,
			const AssociatorParameters&);
  
   void fillMuon(       const edm::Event&,
			TrackDetMatchInfo&,
			const AssociatorParameters&);
   
   void fillCaloTruth(  const edm::Event&,
			TrackDetMatchInfo&,
			const AssociatorParameters&);
   
   void addTAMuonSegmentMatch(TAMuonChamberMatch&,
			    const RecSegment*,
			    const AssociatorParameters&);
   
   void getTAMuonChamberMatches(std::vector<TAMuonChamberMatch>& matches,
			      const AssociatorParameters& parameters);
  
   void           init( const edm::EventSetup&);
   bool isclose(const DetId &detId, TAMuonChamberMatch &match
                            , const AssociatorParameters& parameters);
                            
   math::XYZPoint getPoint( const GlobalPoint& point)
     {
	return math::XYZPoint(point.x(),point.y(),point.z());
     }
   
   math::XYZPoint getPoint( const LocalPoint& point)
     {
	return math::XYZPoint(point.x(),point.y(),point.z());
     }
   
   math::XYZVector getVector( const GlobalVector& vec)
     {
	return math::XYZVector(vec.x(),vec.y(),vec.z());
     }
   
   math::XYZVector getVector( const LocalVector& vec)
     {
	return math::XYZVector(vec.x(),vec.y(),vec.z());
     }
   
   const Propagator* ivProp_;
   Propagator* defProp_;
   CachedTrajectory cachedTrajectory_;
   bool useDefaultPropagator_;
   
   edm::ESHandle<DetIdAssociator> ecalDetIdAssociator_;
   edm::ESHandle<DetIdAssociator> hcalDetIdAssociator_;
   edm::ESHandle<DetIdAssociator>   hoDetIdAssociator_;
   edm::ESHandle<DetIdAssociator> caloDetIdAssociator_;
   edm::ESHandle<DetIdAssociator> muonDetIdAssociator_;
   
   edm::ESHandle<CaloGeometry> theCaloGeometry_;
   edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry_;
   
   edm::ESWatcher<IdealMagneticFieldRecord>     theMagneticFeildWatcher_;
};
#endif