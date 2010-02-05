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
#include <DataFormats/DTRecHit/interface/DTRecHitCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

//#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
//#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"

#include <DataFormats/DetId/interface/DetId.h>
#include <RecoLocalMuon/TracktoRPCEff/interface/TracktoRPCEff.h>
#include <ctime>
#include <TMath.h>

ObjectMapB2* ObjectMapB2::mapInstance = NULL;

ObjectMapB2* ObjectMapB2::GetInstance(const edm::EventSetup& iSetup){
  if (mapInstance == NULL){
    mapInstance = new ObjectMapB2(iSetup);
  }
  return mapInstance;
}

ObjectMapB2::ObjectMapB2(const edm::EventSetup& iSetup){
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
ObjectMapB2CSC* ObjectMapB2CSC::mapInstance = NULL;

ObjectMapB2CSC* ObjectMapB2CSC::GetInstance(const edm::EventSetup& iSetup){
  if (mapInstance == NULL){
    mapInstance = new ObjectMapB2CSC(iSetup);
  }
  return mapInstance;
}

ObjectMapB2CSC::ObjectMapB2CSC(const edm::EventSetup& iSetup){
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

TracktoRPCEff::TracktoRPCEff(const edm::ParameterSet& iConfig)
:theRPCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("RPCRecSegmentLabel")),
//theDTRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("DTRecSegmentLabel")),
//theCSCRecSegmentLabel(iConfig.getUntrackedParameter<edm::InputTag>("CSCRecSegmentLabel")),
cscSegments(iConfig.getUntrackedParameter<edm::InputTag>("cscSegments")),
dt4DSegments(iConfig.getUntrackedParameter<edm::InputTag>("dt4DSegments"))
{ 

    edm::Service<TFileService> fs;

    theInputLabel = iConfig.getParameter<InputTag>("InputLabel");
    trackTransformerParam = iConfig.getParameter<ParameterSet>("TrackTransformer");

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

    htrackinnerxyp1 = fs->make<TH2F>("htrackinnerxyp1", "y vs x of +1 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxyp1 = fs->make<TH2F>("htrackouterxyp1", "y vs x of +1 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

    htrackinnerxyp2 = fs->make<TH2F>("htrackinnerxyp2", "y vs x of +2 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxyp2 = fs->make<TH2F>("htrackouterxyp2", "y vs x of +2 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

    htrackinnerxyp3 = fs->make<TH2F>("htrackinnerxyp3", "y vs x of +3 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxyp3 = fs->make<TH2F>("htrackouterxyp3", "y vs x of +3 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

    htrackinnerxym1 = fs->make<TH2F>("htrackinnerxym1", "y vs x of -1 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxym1 = fs->make<TH2F>("htrackouterxym1", "y vs x of -1 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

    htrackinnerxym2 = fs->make<TH2F>("htrackinnerxym2", "y vs x of -2 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxym2 = fs->make<TH2F>("htrackouterxym2", "y vs x of -2 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

    htrackinnerxym3 = fs->make<TH2F>("htrackinnerxym3", "y vs x of -3 inner hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);
    htrackouterxym3 = fs->make<TH2F>("htrackouterxym3", "y vs x of -3 outer hit ", 100, -1000.0, 1000.0, 100, -1000.0, 1000.0);

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

bool TracktoRPCEff::SetFolderMuonDir(int detid, const edm::EventSetup& iSetup)
{
     edm::ESHandle<RPCGeometry> rpcGeo;
     iSetup.get<MuonGeometryRecord>().get(rpcGeo);
     const GeomDet *whichdet = rpcGeo->idToDet(detid);
     const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);

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

          int spsize = aroll->nstrips()+5;
          int hsize1, wsize;

          if(aroll->isBarrel())
          {
               LocalPoint xmin; LocalPoint xmax;
               const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(aroll->topology()));

               xmin = top_->localPosition(0.);
               xmax = top_->localPosition((float)aroll->nstrips());
               wsize = (int)(fabs( xmax.x()-xmin.x())/2+10);
               hsize1 = (int)(top_->stripLength()/2+6);
          } 
          else
          {
               LocalPoint xmin; LocalPoint xmax;
               const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(aroll->topology()));

               xmin = top_->localPosition(0.);
               xmax = top_->localPosition((float)aroll->nstrips());
               wsize = (int)(fabs( xmax.x()-xmin.x())/2+10);
               hsize1 = (int)(top_->stripLength()/2+6);
          }  

 
          RPCGeomServ servId(detid);
//          if (_debug) cout << "RPCGeomServ : " << servId.name() << " :: " << endl;


          TFileDirectory subDir_RPC = fs->mkdir( "RPC" );
          TFileDirectory dirRegion = subDir_RPC.mkdir(Form("Region_%d", region));
          TFileDirectory dirRing = dirRegion.mkdir(Form("Ring_%d", ring));
          TFileDirectory dirStation = dirRing.mkdir(Form("Station_%d", station));
          TFileDirectory dirSector = dirStation.mkdir( Form("Sector_%d",sector));


           detrtMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrt%s", servId.name().c_str()),
                                        Form("%s Trajectory hit", servId.name().c_str()),
                                   wsize/2, -wsize, wsize, hsize1/2, -hsize1, hsize1);

           detrmMap_[detid] = dirSector.make<TH2F>(Form("RPCTtrm%s", servId.name().c_str()),
                                        Form("%s Trajectory matched hit", servId.name().c_str()),
                                   wsize/2, -wsize, wsize, hsize1/2, -hsize1, hsize1);
    
           detrrMap_[detid] =  dirSector.make<TH1F>(Form("RPCTtrr%s_Residure", servId.name().c_str()),
                                        Form("%s Residure local X", servId.name().c_str()),
                                   400, -20, 20);


           detsrtMap_[detid] = dirSector.make<TH1F>(Form("RPCSrt%s", servId.name().c_str()),
                                        Form("%s strips trajectory", servId.name().c_str()),
                                   spsize, 0.5, spsize+0.5);
           detsrmMap_[detid] = dirSector.make<TH1F>(Form("RPCSrm%s", servId.name().c_str()),
                                        Form("%s strips matched trajectory", servId.name().c_str()),
                                   spsize, 0.5, spsize+0.5);

           det2srmMap_[detid] = dirSector.make<TH2F>(Form("RPC2Srm%s", servId.name().c_str()),
                                        Form("%s strips matched trajectory vs strip", servId.name().c_str()),
                                   51, -10.5, 10.5, spsize, 0.5, spsize+0.5);


           det1BunchX_[detid] =  dirSector.make<TH1F>(Form("RPC1BunchX%s", servId.name().c_str()),
                                        Form("%s BunchX ", servId.name().c_str()),
                                   9, -4.5, 4.5);
           detmsBunchX_[detid] =  dirSector.make<TH2F>(Form("RPCmsBunchXstrip%s", servId.name().c_str()),
                                        Form("%s BunchX vs matched strip", servId.name().c_str()),
                                   9, -4.5, 4.5, spsize, -0.5, spsize+0.5);

           detmsClustX_[detid] = dirSector.make<TH2F>(Form("RPCmsClustX%s", servId.name().c_str()),
                                        Form("%s ClusterSize vs matchedStrip", servId.name().c_str()),
                                   30, 0, 30, spsize, -0.5, spsize+0.5);
           detResiClustX_[detid] =  dirSector.make<TH2F>(Form("RPCResiClustX%s", servId.name().c_str()),
                                        Form("%s Residure and Cluster Size", servId.name().c_str()),
                                   400, -20, 20, 20, 0, 20);


           detRec_[detid] = dirSector.make<TH1F>(Form("RPCRec%s", servId.name().c_str()),
                                        Form("%s rec hit strips", servId.name().c_str()),
                                   spsize, 0.5, spsize+0.5);
           detRecLoc_[detid] = dirSector.make<TH2F>(Form("RPCRLoc%s", servId.name().c_str()),
                                        Form("%s Reco hit localposition", servId.name().c_str()),
                                   wsize/2, -wsize, wsize, hsize1/2, -hsize1, hsize1);

           detBunchX_[detid] = dirSector.make<TH1F>(Form("RPCBunchX%s", servId.name().c_str()),
                                        Form("%s BunchX of rechits", servId.name().c_str()),
                                   9, -4.5, 4.5);
           detClustX_[detid] = dirSector.make<TH1F>(Form("RPCClustX%s", servId.name().c_str()),
                                        Form("%s ClusterSize of rechits", servId.name().c_str()),
                                   30, 0, 30);
           detR2XLoc_[detid] =  dirSector.make<TH2F>(Form("RPCR2XLoc%s", servId.name().c_str()),
                                        Form("%s Reco hit strip : loc x", servId.name().c_str()),
                                   spsize, 0.5, spsize+0.5, wsize/2, -wsize, wsize);

           detR2YLoc_[detid] =  dirSector.make<TH2F>(Form("RPCR2YLoc%s", servId.name().c_str()),
                                        Form("%s Reco hit strip : loc y", servId.name().c_str()),
                                   spsize, 0.5, spsize+0.5, hsize1/2, -hsize1, hsize1);

/*
     if(region != 0)
     {
           detTangle_[detid] = dirSector.make<TH1F>(Form("RPCTang%s", servId.name().c_str()),
                                        Form("%s anlge trajectory strip", servId.name().c_str()),
                                   300, -3, -3);
           detAngle_[detid] =  dirSector.make<TH1F>(Form("RPCAng%s", servId.name().c_str()),
                                        Form("%s angle strip ", servId.name().c_str()),
                                   300, -3, -3);
           detmsAngle_[detid] =  dirSector.make<TH1F>(Form("RPCmsAng%s", servId.name().c_str()),
                                        Form("%s matched angle strip", servId.name().c_str()),
                                   300, -3, -3);
      }
*/
    }
     return true;
   } else return false;

}

bool TracktoRPCEff::ValidRPCSurface(RPCDetId rpcid, LocalPoint LocalP, const edm::EventSetup& iSetup)
{
  edm::ESHandle<RPCGeometry> rpcGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);

  const GeomDet *whichdet3 = rpcGeo->idToDet(rpcid.rawId());
  const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet3);
  float locx=LocalP.x(), locy=LocalP.y(), locz=LocalP.z();
  if(aroll->isBarrel())
  {
     const Bounds &rollbound = rpcGeo->idToDet((rpcid))->surface().bounds();
     float boundlength = rollbound.length();
     float boundwidth = rollbound.width();

     if(fabs(locx)< boundwidth/2 && fabs(locy)<boundlength/2 && locy>-boundlength/2) return true;
     else return false;

   }
   else if(aroll->isForward())
   {
     const Bounds &rollbound = rpcGeo->idToDet((rpcid))->surface().bounds();
     float boundlength = rollbound.length();
     float boundwidth = rollbound.width();

     float nminx = TMath::Pi()*(18*boundwidth/ TMath::Pi() - boundlength)/18;
     float ylimit = ((boundlength)/(boundwidth/2 - nminx/2))*fabs(locx) + boundlength/2 - ((boundlength)/(boundwidth/2 - nminx/2))*(boundwidth/2);
     if(ylimit < -boundlength/2 ) ylimit = -boundlength/2;

     if(fabs(locx)< boundwidth/2 && fabs(locy)<boundlength/2 && locy> ylimit) return true;
     else return false;
   } else return false;
}

void TracktoRPCEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 using reco::TrackCollection;
 Handle<reco::TrackCollection> alltracks;
 iEvent.getByLabel(theInputLabel.label(),alltracks);

/////////////////////////////////////////////////////////
 hsize->Fill(1);
 if(alltracks->empty()) return;
 hntracks->Fill(alltracks->size());
///////////////////////////////////////////////////////////

 TrackTransformerBase *theTrackTransformer;
 if(theInputLabel.label()=="cosmicMuons") theTrackTransformer = new TrackTransformerForCosmicMuons(trackTransformerParam);
 else  theTrackTransformer = new TrackTransformer(trackTransformerParam);
//////////////////////////////////////////////////////////////////////////////////////
 //theTrackTransformer = new TrackTransformer(iConfig);
 theTrackTransformer->setServices(iSetup);  

 edm::ESHandle<RPCGeometry> rpcGeo;
 edm::ESHandle<DTGeometry> dtGeo;
 edm::ESHandle<CSCGeometry> cscGeo;
 
 iSetup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny",thePropagator); 
 iSetup.get<MuonGeometryRecord>().get(rpcGeo);
 iSetup.get<MuonGeometryRecord>().get(dtGeo);
 iSetup.get<MuonGeometryRecord>().get(cscGeo);

//////////////////////////////////////////////////////////
    edm::Handle<RPCRecHitCollection> allRPChits;
    //edm::Handle<DTRecHitCollection> allDThits;
    //edm::Handle<CSCRecHit2DCollection> allCSChits;

    iEvent.getByLabel(theRPCRecSegmentLabel, allRPChits);
    //iEvent.getByLabel(theDTRecSegmentLabel, allDThits);
    //iEvent.getByLabel(theCSCRecSegmentLabel, allCSChits);

    edm::Handle<DTRecSegment4DCollection> all4DSegments;
    iEvent.getByLabel(dt4DSegments, all4DSegments);

    edm::Handle<CSCSegmentCollection> allCSCSegments;
    iEvent.getByLabel(cscSegments, allCSCSegments);
///////////////////////////////////////////////////////////

std::vector<uint32_t> rpcput;
bool checksegment= false;//true;

cout << "Startting TrackExtrapolation method." << endl;
for (TrackCollection::const_iterator track = alltracks->begin(); track !=alltracks->end(); track++)
{
///////////////////////////////////////////////////////////
       htracks_pt->Fill(track->pt());          //  track transverse momentum

       htracks_etaphi->Fill(track->eta(), track->phi()); //  pseudorapidity and azimuthal angle of momentum vector

       htrackinnerxy->Fill(track->innerPosition().X(), track->innerPosition().Y()); //  position of the innermost hit
       htrackouterxy->Fill(track->outerPosition().X(), track->outerPosition().Y()); //  position of the outermost hit
       htrackinnerz->Fill(track->innerPosition().Z()); // position of the innermost hit
       htrackouterz->Fill(track->outerPosition().Z()); // position of the outermost hit

       float tInX = track->innerPosition().X(), tInY = track->innerPosition().Y(), tInZ = track->innerPosition().Z();
       float tOuX = track->outerPosition().X(), tOuY = track->outerPosition().Y(), tOuZ = track->outerPosition().Z();
       if(tInX > tOuX) { float temp=tOuX; tOuX=tInX; tInX=temp; }
       if(tInY > tOuY) { float temp=tOuY; tOuY=tInY; tInY=temp; }
       if(tInZ > tOuZ) { float temp=tOuZ; tOuZ=tInZ; tInZ=temp; }

       cout << "in (x,y,z): ("<< tInX <<", "<< tInY <<", "<< tInZ << ")" << endl;
       cout << "out (x,y,z): ("<< tOuX <<", "<< tOuY <<", "<< tOuZ << ")" << endl;

       if(track->innerPosition().Z()<750 && track->innerPosition().Z()>670) htrackinnerxyp1->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->innerPosition().Z()>-750  && track->innerPosition().Z()<-670) htrackinnerxym1->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->outerPosition().Z()<750 && track->outerPosition().Z()>670) htrackouterxyp1->Fill(track->outerPosition().X(), track->outerPosition().Y());
       if(track->outerPosition().Z()>-750 && track->outerPosition().Z()<-670) htrackouterxym1->Fill(track->outerPosition().X(), track->outerPosition().Y());

       if(track->innerPosition().Z()<900 && track->outerPosition().Z()>750) htrackinnerxyp2->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->innerPosition().Z()>-900 && track->outerPosition().Z()<-750) htrackinnerxym2->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->outerPosition().Z()<900 && track->outerPosition().Z()>750) htrackouterxyp2->Fill(track->outerPosition().X(), track->outerPosition().Y());
       if(track->outerPosition().Z()>-900 && track->outerPosition().Z()<-750) htrackouterxym2->Fill(track->outerPosition().X(), track->outerPosition().Y());

       if(track->innerPosition().Z()>900) htrackinnerxyp3->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->innerPosition().Z()<-900) htrackinnerxym3->Fill(track->innerPosition().X(), track->innerPosition().Y());
       if(track->outerPosition().Z()>900) htrackouterxyp3->Fill(track->outerPosition().X(), track->outerPosition().Y());
       if(track->outerPosition().Z()<-900) htrackouterxym3->Fill(track->outerPosition().X(), track->outerPosition().Y());

       const reco::HitPattern& p = track->hitPattern();
       hnvalidhits->Fill(track->innerPosition().phi(), p.numberOfHits()); // azimuthal angle of Innermost hit, number Of Hits

       hchi2prob->Fill(track->innerPosition().phi(), TMath::Prob(track->chi2(), track->ndof())); // 
       hchi2ndof->Fill(track->innerPosition().phi(), track->normalizedChi2()); // chi-squared divided by n.d.o.f. (or chi-squared * 1e6 if n.d.o.f. is zero)
       hqoverppull->Fill(track->innerPosition().phi(), track->qoverp()/track->qoverpError()); // error on signed transverse curvature
       htrackq->Fill(track->innerPosition().phi(), track->charge()); // track electric charge
///////////////////////////////////////////////////////////

 Trajectories trajectories = theTrackTransformer->transform(*track);
 cout << "Building Trajectory from Track. " << endl;
 std::vector<uint32_t> rpcrolls;
 std::vector<uint32_t> rpcrolls2; 
 std::map<uint32_t, int> rpcNdtcsc;
 std::map<uint32_t, int> rpcrollCounter;

cout << "1. Search expeted RPC roll detid !!" << endl;
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
                float dx = trajLP.x()-trackLP.x(), dy=trajLP.y()-trackLP.y(), dz=trajLP.z()-trackLP.z();
                if( dx>10. && dy>10.) continue;

                DTChamberId dtid(geomDet->geographicalId().rawId());
                int dtW=dtid.wheel(), dtS=dtid.sector(), dtT=dtid.station();
                if(dtS==13) dtS=4; if(dtS==14) dtS=10; 
                ObjectMapB2* TheObject = ObjectMapB2::GetInstance(iSetup);
                DTStationIndex2 theindex(0,dtW,dtS,dtT);
                std::set<RPCDetId> rollsForThisDT = TheObject->GetInstance(iSetup)->GetRolls(theindex);
                for(std::set<RPCDetId>::iterator iteraRoll = rollsForThisDT.begin();iteraRoll != rollsForThisDT.end(); iteraRoll++)
                {                                 
	            const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
	            RPCDetId rpcid = rollasociated->id();

                    TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(rollasociated->id())->surface());
                    if(ptss.isValid()) if(ValidRPCSurface(rpcid, ptss.localPosition(), iSetup))
                    {
                      rpcrollCounter[rollasociated->id().rawId()]++;  
                      bool check = true;
                      vector<uint32_t>::iterator rpcroll;
                      for( rpcroll=rpcrolls.begin() ; rpcroll < rpcrolls.end(); rpcroll++ )
                      if(rollasociated->id().rawId()== *rpcroll) check=false; 

                      if(check==true)
                      {
                        rpcrolls.push_back(rollasociated->id().rawId());
                        RPCGeomServ servId(rollasociated->id().rawId());
                        cout << "1\t Barrel RPC roll" << rollasociated->id().rawId() << " "<< servId.name().c_str() <<endl; 
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
                float dx = trajLP.x()-trackLP.x(), dy=trajLP.y()-trackLP.y(), dz=trajLP.z()-trackLP.z();
                if( dx>10. && dy>10.) continue;

                ObjectMapB2CSC* TheObjectCSC = ObjectMapB2CSC::GetInstance(iSetup);
	        int En = cscid.endcap(), St = cscid.station(), Ri = cscid.ring();
	        int rpcSegment = cscid.chamber();
                if(En==2) En= -1; if(Ri==4) Ri =1; 

                CSCStationIndex2 theindex(En,St,Ri,rpcSegment);
                std::set<RPCDetId> rollsForThisCSC = TheObjectCSC->GetInstance(iSetup)->GetRolls(theindex);
                for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++)
                {
	            const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
	            RPCDetId rpcid = rollasociated->id();

                    TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(rollasociated->id())->surface());
                    if(ptss.isValid()) if(ValidRPCSurface(rpcid, ptss.localPosition(), iSetup))
                    {
                      rpcrollCounter[rollasociated->id().rawId()]++;
                      bool check = true;
                      vector<uint32_t>::iterator rpcroll;
                      for( rpcroll=rpcrolls.begin() ; rpcroll < rpcrolls.end(); rpcroll++ )
                      if(rollasociated->id().rawId()==*rpcroll) check=false; 
                      if(check==true)
                      {
                        rpcrolls.push_back(rollasociated->id().rawId());
                        RPCGeomServ servId(rollasociated->id().rawId());
                        cout << "1\t Forward RPC roll" << rollasociated->id().rawId() << " "<< servId.name().c_str() <<endl;
                      }
                    }
	      	}
             }
          }
       } else { cout << "1\t The hit is not DT/CSC's.   " << endl;} 
    }
 }
 cout << "First step OK!!\n2. Search nearest DT/CSC sufrace!!" << endl;

//////////////////////// rpcrollCounter ///////////////////////////////


 vector<uint32_t>::iterator rpcroll;
 for( rpcroll=rpcrolls.begin() ; rpcroll < rpcrolls.end(); rpcroll++ )
 {
    RPCDetId rpcid(*rpcroll);
  //  if((rEn==0 && rSt ==2) && rpcrollCounter[*rpcroll]<2 ) continue;
  //  else if((rEn!=0 && rSt!=3) && rpcrollCounter[*rpcroll]<4) continue ;
   if( rpcrollCounter[*rpcroll]<4) continue ; 

    //const GeomDet *geomDet = rpcGeo->idToDet(*rpcroll);
    //const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(geomDet);
    const GlobalPoint &rGP = rpcGeo->idToDet(*rpcroll)->surface().toGlobal(LocalPoint(0,0,0));
    RPCGeomServ servId(rpcid);
    int rEn=rpcid.region(), rSe=rpcid.sector(), rWr=rpcid.ring(), rSt=rpcid.station(), rCh=servId.segment();
    uint32_t dtcscid=0; double distance=150, MaxD=150;

    //if(rSt ==2 ) MaxD=100; 
    //else if(rSt ==3 ) MaxD=100;
    //else if(rSt ==4 && rEn==0) MaxD =150;
    for(trackingRecHit_iterator hit=track->recHitsBegin(); hit != track->recHitsEnd(); hit++)
    {
        if((*hit)->isValid())
        {
            DetId id = (*hit)->geographicalId(); 
            if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::DT)
            {
               const GeomDet *geomDet =  dtGeo->idToDet((*hit)->geographicalId());
               //const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(geomDet);
               const GlobalPoint &dtGP = dtGeo->idToDet((*hit)->geographicalId())->surface().toGlobal(LocalPoint(0,0,0));
               double dx = rGP.x()-dtGP.x(), dy = rGP.y()-dtGP.y(), dz = rGP.z()-dtGP.z();
               double distanceN = sqrt(dx*dx+dy*dy+dz*dz);

               DTChamberId dtid(geomDet->geographicalId().rawId());
               int Se = dtid.sector(), Wh = dtid.wheel(), St = dtid.station();
               if(Se == 13) Se=4; if(Se ==14) Se=10;

               //cout << "check a distance between expected RPC and DT position." << endl;
               if( rEn==0&& (rSe-Se)==0 && (rWr-Wh) ==0 && (rSt-St)==0 && distanceN < distance)
               {
                   dtcscid=geomDet->geographicalId().rawId();
                   distance = distanceN;
                   cout << "2\t DT "<< dtcscid << " Wheel : " << Wh << " station : " << St << " sector : " << Se << endl;
               }
            }
            else if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::CSC)
            {
               const GeomDet *geomDet =  cscGeo->idToDet((*hit)->geographicalId());
               //const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(geomDet);
               const GlobalPoint &cscGP = cscGeo->idToDet((*hit)->geographicalId())->surface().toGlobal(LocalPoint(0,0,0));
               double dx = rGP.x()-cscGP.x(), dy = rGP.y()-cscGP.y(), dz = rGP.z()-cscGP.z();
               double distanceN = sqrt(dx*dx+dy*dy+dz*dz);

               CSCDetId cscid(geomDet->geographicalId().rawId());
               int En =cscid.endcap(), Ri=cscid.ring(), St=cscid.station(), Ch=cscid.chamber();
               if(En==2) En=-1; if(Ri==4) Ri=1; //??

               //cout << "check a distance between expected RPC and CSC position." << endl;
               if((rEn-En)==0 && (rSt-St)==0 && (Ch-rCh) ==0 && rWr!=1 && rSt!=4 && distanceN < distance)
               {
                  dtcscid=geomDet->geographicalId().rawId();
                  distance = distanceN;
                  cout << "2\t CSC " <<dtcscid <<" region : " << En << " station : " << St << " Ring : " << Ri << " chamber : " << Ch <<endl;
               }
            } 
         }
    }
    if(dtcscid != 0 && distance < MaxD)
    {
       rpcrolls2.push_back(*rpcroll);
       rpcNdtcsc[*rpcroll] = dtcscid;
    }
 }
 cout << "Second step OK!! \n3. Propagate to RPC from DT/CSC!!" << endl;
 //std::map<uint32_t, int> rpcput;
 vector<uint32_t>::iterator rpcroll2;
 for( rpcroll2=rpcrolls2.begin() ; rpcroll2 < rpcrolls2.end(); rpcroll2++ )
 {
    bool check = true;
    vector<uint32_t>::iterator rpcput_;
    for( rpcput_=rpcput.begin() ; rpcput_ < rpcput.end(); rpcput_++ )
    if(*rpcroll2==*rpcput_) check = false;

    if(check == true)
    {
        RPCDetId rpcid(*rpcroll2);
        uint32_t dtcscid = rpcNdtcsc[*rpcroll2];

        DetId id(dtcscid);
        if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::DT)
        {
           const GeomDet *geomDet =  dtGeo->idToDet(dtcscid);
           const DTLayer *dtlayer = dynamic_cast<const DTLayer *>(geomDet);
        
           if(dtlayer) for(Trajectories::const_iterator trajectory = trajectories.begin(); trajectory != trajectories.end(); ++trajectory)
           {
              const BoundPlane & DTSurface = dtlayer->surface();
              const GlobalPoint dcPoint = DTSurface.toGlobal(LocalPoint(0.,0.,0.));

              TrajectoryMeasurement tMt = trajectory->closestMeasurement(dcPoint);
              TrajectoryStateOnSurface upd2 = (tMt).updatedState();
              if(upd2.isValid())
              {
                 TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(*rpcroll2)->surface());
                 if(ptss.isValid()) if(ValidRPCSurface(rpcid, ptss.localPosition(), iSetup))
                 {
                     float rpcGPX = ptss.globalPosition().x();
                     float rpcGPY = ptss.globalPosition().y();
                     float rpcGPZ = ptss.globalPosition().z();

                     if(tInX > rpcGPX || tOuX < rpcGPX ) continue;
                     if(tInY > rpcGPY || tOuY < rpcGPY ) continue;
                     if(tInZ > rpcGPZ || tOuZ < rpcGPZ ) continue;

                     rpcput.push_back(*rpcroll2);
                     float locx=ptss.localPosition().x();
                     float locy=ptss.localPosition().y(); 
 
                     RPCGeomServ servId(*rpcroll2);
                     cout << "3\t Barrel Expected RPC " << servId.name().c_str() <<
                           " \tLocalposition X: " << locx << ", Y: "<< locy << " GlobalPosition(x,y,z) (" << rpcGPX <<", "<< rpcGPY <<", " << rpcGPZ << ")"<< endl;

                     SetFolderMuonDir(*rpcroll2, iSetup);
                     detrtMap_[*rpcroll2]->Fill(locx, locy);
                     
                     const GeomDet *whichdet1 = rpcGeo->idToDet(*rpcroll2);
                     const RPCRoll *aroll1 = dynamic_cast<const RPCRoll *>(whichdet1);

                     const float stripPredicted =aroll1->strip(LocalPoint(locx,locy,0.));
                     detsrtMap_[*rpcroll2]->Fill(stripPredicted);
                     float mstrip;
                     float residual, sum, residualy, recX, recY;
                     int BunchX,  clusterSize, firstClusterStrip;
                     if (findmatch(allRPChits, *rpcroll2, locx, locy, residual, sum, residualy, recX, recY, mstrip, BunchX, clusterSize, firstClusterStrip, iSetup))
                     {
                        cout << "3\t residual : "<< residual << endl;
                        if(fabs(residual) < 10.0)
                        {
                             detrmMap_[*rpcroll2]->Fill(locx, locy);
                             detrrMap_[*rpcroll2]->Fill(residual);
                             det2srmMap_[*rpcroll2]->Fill(stripPredicted-mstrip, stripPredicted);
   
                             detsrmMap_[*rpcroll2]->Fill(stripPredicted);
                             detmsClustX_[*rpcroll2]->Fill(clusterSize, stripPredicted);
                             detmsBunchX_[*rpcroll2]->Fill(BunchX,stripPredicted);
                             det1BunchX_[*rpcroll2]->Fill(BunchX);
                        }
                        detResiClustX_[*rpcroll2]->Fill(residual, clusterSize);
                     }
/////////////////////////////
                     // RPCRecHit RPCPoint(*rpcroll2,0,LocalPoint(locx,locy,locz));
                     // cout << "Track \tDT\t RPCDetId: "<<*rpcroll2<<" Local Position x: "<<ptss.localPosition().x()<<", y: "<< ptss.localPosition().y()<< endl;
        
                      //RPCPointVector.clear();
                      //RPCPointVector.push_back(RPCPoint);
                      //_ThePoints->put(*rpcroll2,RPCPointVector.begin(),RPCPointVector.end());
                      //rpcput.push_back(*rpcroll2);
                   
                 }
              }
           }
        }
        else if (id.det() == DetId::Muon  &&  id.subdetId() == MuonSubdetId::CSC)
        {
           const GeomDet *geomDet4 =  cscGeo->idToDet(dtcscid);
           const CSCLayer *csclayer = dynamic_cast<const CSCLayer *>(geomDet4);
        
           if(csclayer) for(Trajectories::const_iterator trajectory = trajectories.begin(); trajectory != trajectories.end(); ++trajectory)
           {
              const BoundPlane & CSCSurface = csclayer->surface();
              const GlobalPoint dcPoint = CSCSurface.toGlobal(LocalPoint(0.,0.,0.));

              TrajectoryMeasurement tMt = trajectory->closestMeasurement(dcPoint);
              TrajectoryStateOnSurface upd2 = (tMt).updatedState();
              if(upd2.isValid())
              {
                 TrajectoryStateOnSurface ptss =  thePropagator->propagate(upd2, rpcGeo->idToDet(*rpcroll2)->surface());
                 if(ptss.isValid()) if(ValidRPCSurface(rpcid, ptss.localPosition(), iSetup))
                 {
                     float rpcGPX = ptss.globalPosition().x();
                     float rpcGPY = ptss.globalPosition().y();
                     float rpcGPZ = ptss.globalPosition().z();

                     if(tInX > rpcGPX || tOuX < rpcGPX ) continue;
                     if(tInY > rpcGPY || tOuY < rpcGPY ) continue;
                     if(tInZ > rpcGPZ || tOuZ < rpcGPZ ) continue;

                     rpcput.push_back(*rpcroll2);
                     float locx=ptss.localPosition().x();
                     float locy=ptss.localPosition().y();
//////////////////////////////////////
                     RPCGeomServ servId(*rpcroll2);
                     cout << "3\t Forward Expected RPC " << servId.name().c_str() <<
                           " localposition X: " << locx << ", Y: "<< locy << endl;


                     SetFolderMuonDir(*rpcroll2,iSetup);
                     detrtMap_[*rpcroll2]->Fill(locx, locy);
      
                     const GeomDet *whichdet1 = rpcGeo->idToDet(*rpcroll2);
                     const RPCRoll *aroll1 = dynamic_cast<const RPCRoll *>(whichdet1);
                     const float stripPredicted =aroll1->strip(LocalPoint(locx,locy,0.));
                     detsrtMap_[*rpcroll2]->Fill(stripPredicted);
                     float mstrip, residual, sum, residualy, recX, recY;
                     int BunchX,  clusterSize, firstClusterStrip;
                     if (findmatch(allRPChits, *rpcroll2, locx, locy, residual, sum, residualy, recX, recY, mstrip, BunchX, clusterSize, firstClusterStrip, iSetup))
                     {
                        cout << "3\t residual : "<< residual << endl;
           	        if(fabs(residual) < 10.0)
                        {
                             detrmMap_[*rpcroll2]->Fill(locx, locy);
                             detrrMap_[*rpcroll2]->Fill(residual);
                             det2srmMap_[*rpcroll2]->Fill(stripPredicted-mstrip, stripPredicted);
   
                             detsrmMap_[*rpcroll2]->Fill(stripPredicted);
                             detmsClustX_[*rpcroll2]->Fill(clusterSize, stripPredicted);
                             detmsBunchX_[*rpcroll2]->Fill(BunchX,stripPredicted);
                             det1BunchX_[*rpcroll2]->Fill(BunchX);
                             detResiClustX_[*rpcroll2]->Fill(residual, clusterSize);
   
                        }
                     }
////////////////////////////////////
                      //RPCRecHit RPCPoint(*rpcroll2,0,LocalPoint(locx,locy,locz));
                      //cout << "Track \t-CSC\t RPCDetId: "<<*rpcroll2<<" Local Position x: "<<ptss.localPosition().x()<<", y: "<< ptss.localPosition().y()<< endl;
                      //RPCPointVector.clear();
                      //RPCPointVector.push_back(RPCPoint);
                      //_ThePoints->put(*rpcroll2,RPCPointVector.begin(),RPCPointVector.end());

        
                   
                 }
              }
           }       
        }
    }
 }
 cout << "last steps OK!! " << endl; 
}

/////////////////////////////////////////////////////////////////////////
    for (RPCRecHitCollection::const_iterator irpchit =  allRPChits->begin(); irpchit != allRPChits->end(); irpchit++)
    {
        const GeomDet *whichdet = rpcGeo->idToDet(irpchit->geographicalId());
        const RPCRoll *aroll = dynamic_cast<const RPCRoll *>(whichdet);
        if (aroll)
        {
             int detid = aroll->geographicalId().rawId(); 
             SetFolderMuonDir(detid, iSetup);
             detBunchX_[detid]->Fill(irpchit->BunchX());
             for(int cs=0;cs< irpchit->clusterSize()+1;cs++)  detRec_[detid]->Fill(irpchit->firstClusterStrip()+cs);
             detClustX_[detid]->Fill(irpchit->firstClusterStrip()); 
             detRecLoc_[detid]->Fill(irpchit->localPosition().x(),irpchit->localPosition().y());
             detR2XLoc_[detid]->Fill(irpchit->firstClusterStrip(),irpchit->localPosition().x());
             detR2YLoc_[detid]->Fill(irpchit->firstClusterStrip(),irpchit->localPosition().y());
///
             if(aroll->isBarrel());
             else     
             {
//                 const GlobalPoint &p = whichdet->surface().toGlobal(irpchit->localPosition());
                 const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(aroll->topology()));
                 float stripw = top_->localPitch(irpchit->localPosition());
                 float strips = aroll->nstrips();
                 cout << "stripsw : " << stripw*strips << " ---  : " << detid << endl;
             }     
//////////

        }
    }
//////////////////////////////////////////////////////////////////////////////////////
}


void TracktoRPCEff::beginJob(const edm::EventSetup&)
{
}

void TracktoRPCEff::endJob()
{

}

TracktoRPCEff::~TracktoRPCEff(){

}

bool TracktoRPCEff::findmatch(const edm::Handle<RPCRecHitCollection> &hitcoll, int detid, float locx, float locy, float &residual,float &sum,float &residualy, float &recX, float &recY, float &mstrip, int &BunchX, int &clusterSize, int &firstClusterStrip, const edm::EventSetup& iSetup)
{
    bool matchfound = false;
    double maxres = maxdist;
    edm::ESHandle<RPCGeometry> rpcGeo;
    iSetup.get<MuonGeometryRecord>().get(rpcGeo);

    for (RPCRecHitCollection::const_iterator hit =  hitcoll->begin(); hit != hitcoll->end(); hit++)
    {
        const GeomDet *whichdet = rpcGeo->idToDet(hit->geographicalId());
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
            BunchX = hit->BunchX();
            clusterSize = hit->clusterSize();
            firstClusterStrip = hit->firstClusterStrip(); 
        }
    }

    return matchfound;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TracktoRPCEff);
