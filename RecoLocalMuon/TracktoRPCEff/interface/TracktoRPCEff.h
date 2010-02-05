
// system include files
#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/Common/interface/Ref.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerBase.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
//#include "RecoLocalMuon/RPCRecHit/interface/TracktoRPC.h"
//#include "RecoLocalMuon/RPCRecHit/interface/CSCSegtoRPC.h"
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace edm;
using namespace std;
using reco::MuonCollection;
using reco::TrackCollection;
typedef std::vector<Trajectory> Trajectories;
//edm::ESHandle<RPCGeometry> rpcGeo;
double maxdist = 7777777.0;

class TracktoRPCEff : public edm::EDAnalyzer {
public:
  explicit TracktoRPCEff(const edm::ParameterSet&);
  ~TracktoRPCEff();
//  RPCRecHitCollection* thePoints(){return _ThePoints;}
    protected:
        bool SetFolderMuonDir(int, const edm::EventSetup&);
        bool findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, 
                       float &residual,float &sum, float &residualyy, float &recX, float &recY,  float &mstrip,
                       int &BunchX, int &clusterSize, int &firstClusterStrip, const edm::EventSetup& iSetup);

        bool ValidRPCSurface(RPCDetId rpcid,LocalPoint LocalP, const edm::EventSetup& iSetup);
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//  RPCRecHitCollection* _ThePoints;
//  edm::OwnVector<RPCRecHit> RPCPointVector;
//  double MaxD;
 edm::InputTag theInputLabel;

 TrackTransformerBase *theTrackTransformer;
 edm::ESHandle<Propagator> thePropagator;
 ParameterSet trackTransformerParam;
/// the name of the DT rec hits collection
        //edm::InputTag theDTRecSegmentLabel;
/// the name of the CSC rec hits collection
        //edm::InputTag theCSCRecSegmentLabel;
/// the name of the RPC rec hits collection
        edm::InputTag theRPCRecSegmentLabel;

        edm::InputTag cscSegments;
        edm::InputTag dt4DSegments;

        bool _debug;

        std::map<int, TH2F*> detrtMap_;
        std::map<int, TH2F*> detrmMap_;
        std::map<int, TH1F*> detrrMap_;

        std::map<int, TH1F*> detsrtMap_;
        std::map<int, TH1F*> detsrmMap_;
        std::map<int, TH2F*> det2srmMap_;

        std::map<int, TH1F*> det1BunchX_;
        std::map<int, TH2F*> detmsBunchX_;
        std::map<int, TH2F*> detmsClustX_;
        std::map<int, TH2F*> detResiClustX_;

        // reco hits
        std::map<int, TH1F*> detRec_;
        std::map<int, TH2F*> detRecLoc_;
        std::map<int, TH2F*> detR2XLoc_;
        std::map<int, TH2F*> detR2YLoc_;

        std::map<int, TH1F*> detBunchX_;
        std::map<int, TH1F*> detClustX_;

        //std::map<int, TH1F*> detTangle_;
        //std::map<int, TH1F*> detAngle_;
        //std::map<int, TH1F*> detmsAngle_;

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

        TH2F *htrackinnerxyp1;
        TH2F *htrackouterxyp1;

        TH2F *htrackinnerxyp2;
        TH2F *htrackouterxyp2;
         
        TH2F *htrackinnerxyp3;
        TH2F *htrackouterxyp3;

        TH2F *htrackinnerxym1;
        TH2F *htrackouterxym1;

        TH2F *htrackinnerxym2;
        TH2F *htrackouterxym2;

        TH2F *htrackinnerxym3;
        TH2F *htrackouterxym3;

        TH2F *htrackrechitsxy;
        TH1F *htrackrechitsz;

        TH2F *htrackrechitsDTxy;
        TH1F *htrackrechitsDTz;

        TH2F *htrackrechitsRPCxy;
        TH1F *htrackrechitsRPCz;

        TH2F *htrackrechitsCSCxy;
        TH1F *htrackrechitsCSCz;

    edm::Service<TFileService> fs;
 
};

/////////////////////////
class DTStationIndex2{
public: 
  DTStationIndex2():_region(0),_wheel(0),_sector(0),_station(0){}
  DTStationIndex2(int region, int wheel, int sector, int station) : 
    _region(region),
    _wheel(wheel),
    _sector(sector),
    _station(station){}
  ~DTStationIndex2(){}
  int region() const {return _region;}
  int wheel() const {return _wheel;}
  int sector() const {return _sector;}
  int station() const {return _station;}
  bool operator<(const DTStationIndex2& dtind) const{
    if(dtind.region()!=this->region())
      return dtind.region()<this->region();
    else if(dtind.wheel()!=this->wheel())
      return dtind.wheel()<this->wheel();
    else if(dtind.sector()!=this->sector())
      return dtind.sector()<this->sector();
    else if(dtind.station()!=this->station())
      return dtind.station()<this->station();
    return false;
  }

private:
  int _region;
  int _wheel;
  int _sector;
  int _station; 
};

class ObjectMapB2{
public:
  static ObjectMapB2* GetInstance(const edm::EventSetup& iSetup);
  std::set<RPCDetId> GetRolls(DTStationIndex2 dtstationindex){return mapInstance->rollstoreDT[dtstationindex];}
//protected:
  std::map<DTStationIndex2,std::set<RPCDetId> > rollstoreDT;
  ObjectMapB2(const edm::EventSetup& iSetup);
private:
  static ObjectMapB2* mapInstance;
}; 
///////////////////////
class CSCStationIndex2{
public:
  CSCStationIndex2():_region(0),_station(0),_ring(0),_chamber(0){}
  CSCStationIndex2(int region, int station, int ring, int chamber):
    _region(region),
    _station(station),
    _ring(ring),
    _chamber(chamber){}
  ~CSCStationIndex2(){}
  int region() const {return _region;}
  int station() const {return _station;}
  int ring() const {return _ring;}
  int chamber() const {return _chamber;}
  bool operator<(const CSCStationIndex2& cscind) const{
    if(cscind.region()!=this->region())
      return cscind.region()<this->region();
    else if(cscind.station()!=this->station())
      return cscind.station()<this->station();
    else if(cscind.ring()!=this->ring())
      return cscind.ring()<this->ring();
    else if(cscind.chamber()!=this->chamber())
      return cscind.chamber()<this->chamber();
    return false;
  }

private:
  int _region;
  int _station;
  int _ring;  
  int _chamber;
};

class ObjectMapB2CSC{
public:
  static ObjectMapB2CSC* GetInstance(const edm::EventSetup& iSetup);
  std::set<RPCDetId> GetRolls(CSCStationIndex2 cscstationindex){return mapInstance->rollstoreCSC[cscstationindex];}
//protected:
  std::map<CSCStationIndex2,std::set<RPCDetId> > rollstoreCSC;
  ObjectMapB2CSC(const edm::EventSetup& iSetup);
private:
  static ObjectMapB2CSC* mapInstance;
}; 
/////////////////////////

