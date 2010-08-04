#ifndef TrackAssociator_MuonDetIdAssociator_h
#define TrackAssociator_MuonDetIdAssociator_h 1
// -*- C++ -*-
//
// Package:    TrackAssociator
// Class:      MuonDetIdAssociator
// 
/*

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dmytro Kovalskyi
//         Created:  Fri Apr 21 10:59:41 PDT 2006
// $Id: MuonDetIdAssociator.h,v 1.7 2009/05/12 13:27:09 jribnik Exp $
//
//

#include "ForSccosmic/TrackAssociator/interface/DetIdAssociator.h"
#include "ForSccosmic/TrackAssociator/interface/TAMuonChamberMatch.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "CondFormats/CSCObjects/interface/CSCBadChambers.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class MuonDetIdAssociator: public DetIdAssociator{
 public:
   MuonDetIdAssociator():DetIdAssociator(48, 48 , 0.125),geometry_(0),cscbadchambers_(0),includeBadChambers_(0){};
   MuonDetIdAssociator(const int nPhi, const int nEta, const double etaBinSize)
     :DetIdAssociator(nPhi, nEta, etaBinSize),geometry_(0),cscbadchambers_(0),includeBadChambers_(0){};

   MuonDetIdAssociator(const edm::ParameterSet& pSet)
     :DetIdAssociator(pSet.getParameter<int>("nPhi"),pSet.getParameter<int>("nEta"),pSet.getParameter<double>("etaBinSize")),geometry_(0),cscbadchambers_(0),includeBadChambers_(pSet.getParameter<bool>("includeBadChambers")){};
   
   virtual void setGeometry(const GlobalTrackingGeometry* ptr){ geometry_ = ptr; }

   virtual void setGeometry(const DetIdAssociatorRecord& iRecord){
      edm::ESHandle<GlobalTrackingGeometry> geometryH;
      iRecord.getRecord<GlobalTrackingGeometryRecord>().get(geometryH);
      setGeometry(geometryH.product());
   };

   virtual void setCSCBadChambers(const CSCBadChambers* ptr){ cscbadchambers_ = ptr; }

   virtual void setConditions(const DetIdAssociatorRecord& iRecord){
      edm::ESHandle<CSCBadChambers> cscbadchambersH;
      iRecord.getRecord<CSCBadChambersRcd>().get(cscbadchambersH);
      setCSCBadChambers(cscbadchambersH.product());
   };

   virtual const GeomDet* getGeomDet( const DetId& id ) const;

 protected:
   
   virtual void check_setup() const;
   
   virtual GlobalPoint getPosition(const DetId& id) const;
   
   virtual std::set<DetId> getASetOfValidDetIds() const;
   
   virtual std::vector<GlobalPoint> getDetIdPoints(const DetId& id) const;

   virtual bool insideElement(const GlobalPoint& point, const DetId& id) const;

   const GlobalTrackingGeometry* geometry_;

   const CSCBadChambers* cscbadchambers_;
   bool includeBadChambers_;
};
#endif