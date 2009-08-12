/**********************************************
 *                                            *
 *          Raffaello Trentadue               *
 *         INFN, Sezione di Bari              *
 *      Via Amendola 173, 70126 Bari          *
 *         Phone: +390805442441               *
 *      raffaello.trentadue@ba.infn.it        *
 *                                            *
 *                                            *
 **********************************************/


#include "./RPCEffTrackExtrapolationNew.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include <DataFormats/MuonDetId/interface/DTChamberId.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include <DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h>

#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/CSCGeometry/interface/CSCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuRegionalCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTExtendedCand.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"

#include <FWCore/Utilities/interface/Exception.h>
#include "DataFormats/Provenance/interface/Timestamp.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"



#include <sys/time.h>
#include <algorithm>
#include <memory>
#include <cmath>
#include "math.h"
#include "TFile.h"
#include "TTree.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdlib>
#include <stdio.h>
#include "stdlib.h"
#include <string>
#include <memory>

using namespace edm;
using namespace reco;
using namespace std;

RPCEffTrackExtrapolationNew::RPCEffTrackExtrapolationNew(const edm::ParameterSet& iConfig){
  LogTrace("RPCEffTrackExtrapolation") <<"Dentro Costruttore"<<std::endl;

  RootFileName  = iConfig.getUntrackedParameter<std::string>("EffRootFileName", "RPCEffTrackExtrapolationNew.root"); 
  TjInput  = iConfig.getUntrackedParameter<std::string>("trajectoryInput");
  thePropagatorName = iConfig.getParameter<std::string>("PropagatorName");
  
  dt4DSegments = iConfig.getParameter< edm::InputTag >("dt4DSegments");
  cscSegments =  iConfig.getParameter< edm::InputTag >("cscSegments");
  RPCRecHits = iConfig.getParameter< edm::InputTag >("RPCRecHits");
  theNavigationType = iConfig.getParameter<std::string>("NavigationType");

  theService = new MuonServiceProxy(iConfig.getParameter<edm::ParameterSet>("ServiceParameters"));
  theMeasurementExtractor = new MuonDetLayerMeasurements(dt4DSegments,
							 cscSegments,
							 RPCRecHits,
							 false,
							 false,
							 true);
  
  theEstimator = new Chi2MeasurementEstimator(30.,3.0);
  themyHistoClassDbe = new MyHistoClassDbeNew();
    
  ParameterSet trackTransformerParam = iConfig.getParameter<ParameterSet>("TrackTransformer");
  theTrackTransformer = new TrackTransformerForCosmicMuons(trackTransformerParam);

  count_goodevent = 0;
  count_effevent = 0;

}

//void RPCEffTrackExtrapolationNew::beginRun( edm::Run& r, const edm::EventSetup& iSetup)
void RPCEffTrackExtrapolationNew::beginJob( const edm::EventSetup& iSetup)
{

  LogTrace("RPCEffTrackExtrapolation") <<"Dentro Begin Job"<<std::endl;
  dbe = Service<DQMStore>().operator->();

  edm::ESHandle<RPCGeometry> rpcGeo;
  edm::ESHandle<DTGeometry> dtGeo;
  edm::ESHandle<CSCGeometry> cscGeo;
  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  iSetup.get<MuonGeometryRecord>().get(dtGeo);
  iSetup.get<MuonGeometryRecord>().get(cscGeo);

  const RPCGeometry *pGeom = &*rpcGeo;

  themyHistoClassDbe->init(dbe,pGeom,RootFileName);
}

vector<const DetLayer*> RPCEffTrackExtrapolationNew::compatibleLayers(const DetLayer *initialLayer,
								   FreeTrajectoryState& fts,
								   PropagationDirection propDir){
  vector<const DetLayer*> detLayers;
  
  if(theNavigationType == "Standard"){
    detLayers = initialLayer->compatibleLayers(fts,propDir); 
    detLayers.insert(detLayers.begin(),initialLayer);
  }
  else if (theNavigationType == "Direct"){
    DirectMuonNavigation navigation(&*theService->detLayerGeometry());
    detLayers = navigation.compatibleLayers(fts,propDir);
  }
  else
    edm::LogError("Muon|RecoMuon|StandAloneMuonFilter") << "No Properly Navigation Selected!!"<<endl;
  
  return detLayers;
}


void RPCEffTrackExtrapolationNew::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  try{
    //    std::cout<<"Understanding ---> ----------------------- New event: "<<iEvent.id().event()<<"-----------------------"<<std::endl;

    LogTrace("RPCEffTrackExtrapolation") <<"Dentro Analyzer"<<std::endl;

    theTrackTransformer->setServices(iSetup); //add for trajectory of cosmicmuon
    theService->update(iSetup);
    theMeasurementExtractor->setEvent(iEvent);
    
    int bxE = iEvent.bunchCrossing();
    int thisOrbit = iEvent.orbitNumber();
    int orbitTime = 3564*thisOrbit + bxE;
    
    TimeValue_t time=iEvent.time().value(); 
    timeval *tmval=(timeval*)&time;
    int nevent = iEvent.id().event();
    int nrun = iEvent.id().run();
    int tempo = tmval->tv_usec;
    int lumisection=iEvent.luminosityBlock();
    
    edm::ESHandle<DTGeometry> dtGeo;
    iSetup.get<MuonGeometryRecord>().get(dtGeo);
    
    edm::ESHandle<CSCGeometry> cscGeo;
    iSetup.get<MuonGeometryRecord>().get(cscGeo);
    
    edm::ESHandle<RPCGeometry> rpcGeo;
    iSetup.get<MuonGeometryRecord>().get(rpcGeo);
    
    edm::Handle<RPCRecHitCollection> rpcHits;
    iEvent.getByLabel(RPCRecHits,rpcHits);
    
    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);
    
    ESHandle<GlobalTrackingGeometry> TrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(TrackingGeometry);

    Handle<reco::TrackCollection> staTracks;
    iEvent.getByLabel(TjInput, staTracks);
    
    reco::TrackCollection::const_iterator staTrack;
    
    LogTrace("RPCEffTrackExtrapolation") << "Before any TSOS"<<std::endl;  
    
    if(staTracks->size() > 0){
      int counttrack = 0;
      for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack){
	counttrack++;

	Trajectories trajectories = theTrackTransformer->transform(*staTrack);  //add for trajectory of cosmicmuon 
	reco::TransientTrack track(*staTrack,&*theMGField,TrackingGeometry);
	
	LogTrace("RPCEffTrackExtrapolation") << "Loops on the track!"<<std::endl;  

	GlobalVector mfield;
	const MagneticField* field = track.field();
	const math::XYZVector & innMom =track.track().innerMomentum();
	double innMomX = innMom.x();
	double innMomY = innMom.y();
	double innMomZ = innMom.z();
	const math::XYZPoint & outPos =track.track().outerPosition();
	double outPosX = outPos.x();
	double outPosY = outPos.y();
	double outPosZ = outPos.z();
	const math::XYZVector & outMom = track.track().outerMomentum();
	double outMomX = outMom.x();     
	double outMomY = outMom.y();
	double outMomZ = outMom.z();
	const math::XYZPoint & innPos = track.track().innerPosition();
	double innPosX = innPos.x();     
	double innPosY = innPos.y();
	double innPosZ = innPos.z();
	int charge = track.charge();
	float chi2 = track.chi2() ;
	int ndof = track.ndof();
	float normalizedChi2 = track.normalizedChi2();
	int numValidHits = track.numberOfValidHits();
	int numLostHits = track.numberOfLostHits();
	unsigned int outDetId = track.track().outerDetId();
	unsigned int innDetId = track.track().innerDetId();
	double outPx = track.track().outerPx();
	double outPy = track.track().outerPy();
	double outPz = track.track().outerPz();
	double outY =track.track().outerY();
	double outX = track.track().outerX();
	double outZ = track.track().outerZ();
	double outP = track.track().outerP();
	double outPt = track.track().outerPt();
	double outPhi = track.track().outerPhi();
	double outEta = track.track().outerEta();
	double outTheta = track.track().outerTheta();
	double outRadius = track.track().outerRadius();
	
	themyHistoClassDbe->fillTrackPlot(charge, chi2,normalizedChi2, 
					  numValidHits, numLostHits,
					  outP, outPt, outPhi,
					  outEta, outTheta);

	LogTrace("RPCEffTrackExtrapolation") << "Just before TSOSINNER!"<<std::endl;  
	
	if(chi2 > 180) continue;
	if(numValidHits < 26) continue;

	if(outPt < 10) continue;

	TrajectoryStateOnSurface tsosinner = track.innermostMeasurementState();
	if(!tsosinner.isValid()) continue;
	
	FreeTrajectoryState* fts = tsosinner.freeState(true);

	DetId idinner = track.recHit(track.recHitsSize()-1)->geographicalId();

        std::vector<const DetLayer*> vdetlayers;
        const DetLayer *initialLayer = theService->detLayerGeometry()->idToLayer( idinner );

        vdetlayers = this->compatibleLayers(initialLayer,*fts,alongMomentum);
        
	FreeTrajectoryState ftstemp = *fts;
		
	for(std::vector<const DetLayer*>::iterator itdet = vdetlayers.begin(); itdet != vdetlayers.end(); ++itdet){
	  
	  DetId idDetLay = (*itdet)->basicComponents().front()->geographicalId();
	  
	  if((idDetLay.det() == DetId::Muon && idDetLay.subdetId() == MuonSubdetId::DT) ||
	     (idDetLay.det() == DetId::Muon && idDetLay.subdetId() == MuonSubdetId::CSC)) continue;

	  const TrajectoryStateOnSurface& tsos = theService->propagator(thePropagatorName)->propagate(ftstemp,(*itdet)->surface());
	  if(!tsos.isValid()) continue;

	  const Propagator *p = &*(theService->propagator(thePropagatorName));
	  std::vector<DetWithState> resultDetWithState = (*itdet)->compatibleDets(tsos,*p,*theEstimator);
	  
	  if(resultDetWithState.size() > 0){

	    std::cout<<"Understanding --> Good events counter: "<<count_goodevent<<std::endl;
	    const GeomDet* geomd = resultDetWithState[0].first;
	    
	    RPCDetId rpcid(geomd->geographicalId().rawId());

	    const RPCRoll* rollasociated = rpcGeo->roll(rpcid);
	    const BoundPlane & RPCSurface = rollasociated->surface();
	    GlobalPoint rpRPC = RPCSurface.toGlobal(LocalPoint(0.,0.,0.));

	    FreeTrajectoryState* fts;
	    TrajectoryMeasurement  tMt;
	    float minDistance = 999.;
	    
	    for(Trajectories::const_iterator trajectory = trajectories.begin();  // added for trajectory from YoungKwon's code
		trajectory != trajectories.end(); ++trajectory){

	      TransientTrackingRecHit::ConstRecHitContainer rechits = trajectory->recHits();
	      for(TransientTrackingRecHit::ConstRecHitContainer::const_iterator recHit = rechits.begin(); 
		  recHit != rechits.end(); ++recHit){
		if((*recHit)->isValid()){
		  const GeomDet* geomDet = TrackingGeometry->idToDet((*recHit)->geographicalId());
		  if(DetId(geomDet->geographicalId().rawId()).det() == DetId::Muon && DetId(geomDet->geographicalId().rawId()).subdetId() == MuonSubdetId::DT){
		      DTChamberId chdtid(geomDet->geographicalId().rawId());
		      if (chdtid.wheel() == rpcid.ring() && 
			  chdtid.sector() == rpcid.sector() && 
			  chdtid.station() == rpcid.station()){
			
			GlobalPoint dtgp = geomDet->toGlobal((*recHit)->localPosition());
			if (sqrt(pow(dtgp.x()-rpRPC.x(),2)+pow(dtgp.y()-rpRPC.y(),2)+pow(dtgp.z()-rpRPC.z(),2)) < minDistance){
			  minDistance = sqrt(pow(dtgp.x()-rpRPC.x(),2)+pow(dtgp.y()-rpRPC.y(),2)+pow(dtgp.z()-rpRPC.z(),2));
			  tMt = trajectory->closestMeasurement(geomDet->toGlobal((*recHit)->localPosition()));
			}
		      }
		  }
		  else if(DetId(geomDet->geographicalId().rawId()).det() == DetId::Muon && DetId(geomDet->geographicalId().rawId()).subdetId() == MuonSubdetId::CSC){
		    CSCDetId chdtid(geomDet->geographicalId().rawId());
		    if (chdtid.ring() == rpcid.ring() && 
			chdtid.station() == rpcid.station()){
		      
		      GlobalPoint dtgp = geomDet->toGlobal((*recHit)->localPosition());
		      if (sqrt(pow(dtgp.x()-rpRPC.x(),2)+pow(dtgp.y()-rpRPC.y(),2)+pow(dtgp.z()-rpRPC.z(),2)) < minDistance){
			minDistance = sqrt(pow(dtgp.x()-rpRPC.x(),2)+pow(dtgp.y()-rpRPC.y(),2)+pow(dtgp.z()-rpRPC.z(),2));
			tMt = trajectory->closestMeasurement(geomDet->toGlobal((*recHit)->localPosition()));
		      }
		    }
		  }
		}
	      }
	    }
	  
	    TrajectoryStateOnSurface upd2 = (tMt).updatedState(); // added for trajectory from YoungKwon Jo's code 
	    if (upd2.isValid()) fts = upd2.freeState(true); 

	    TrajectoryStateOnSurface tsosAtRPC = theService->propagator(thePropagatorName)->propagate(*fts,RPCSurface);
	    
	    LocalPoint xmin; LocalPoint xmax; float rsize; float stripl; float stripw;
	    
	    if(rollasociated->isBarrel()){
	      const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
	      xmin = top_->localPosition(0.);
	      xmax = top_->localPosition((float)rollasociated->nstrips());
	      rsize = fabs( xmax.x()-xmin.x()); 
	      stripl = top_->stripLength();
	      stripw = top_->pitch();
	    }
	    else{
	      const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&(rollasociated->topology()));
	      //aggiungere endcap parameters
	    }
	    
	    const Bounds & rollbound = RPCSurface.bounds();
	    float boundlenght = rollbound.length();
	    float boundwidth = rollbound.width();
	    float boundthickness = rollbound.thickness();
	    
	    const float lxextrap = tsosAtRPC.localPosition().x();
	    const float lyextrap = tsosAtRPC.localPosition().y();
	    const float lzextrap = tsosAtRPC.localPosition().z();
	    
	    float stripPredicted;
	    
	    float gxextrap; float gyextrap; float gzextrap;
	    float gvX; float gvY; float gvZ; float zangle;
	    float gdX; float gdY; float gdZ;
	    double mfieldX; double mfieldY; double mfieldZ;
	    int rawrpcid; std::string nameRoll; RPCDetId rpcId;

	    rawrpcid = rpcid.rawId();
	    RPCGeomServ rpcsrv(rpcid);
	    nameRoll = rpcsrv.name();
	    
	    stripPredicted =rollasociated->strip(LocalPoint(tsosAtRPC.localPosition().x(),tsosAtRPC.localPosition().y(),0.));
	    
	    GlobalPoint PointRollGlobal = RPCSurface.toGlobal(tsosAtRPC.localPosition());
	    gxextrap = PointRollGlobal.x();
	    gyextrap = PointRollGlobal.y();
	    gzextrap = PointRollGlobal.z();
	    
	    std::cout<<"Understanding ---> Extrapolated hits coord : "<<PointRollGlobal<<std::endl;

	    //	    zangle = tsosAtRPC.localDirection().theta();
	    LocalVector lv = tsosAtRPC.localDirection();
	    zangle = atan(lv.z()/lv.x()); 

	    if(field->isDefined(PointRollGlobal)) mfield = field->inTesla(PointRollGlobal);
	    
	    mfieldX = mfield.x();
	    mfieldY = mfield.y();
	    mfieldZ = mfield.z();
	    
	    int efficiency1 = 0;
	    int efficiency2 = 0;
	    int efficiency3 = 0;
	    int efficiency4 = 0;
	    int efficiency5 = 0;
	    int efficiency6 = 0;
	    int efficiency7 = 0;
	    int efficiency8 = 0;
	    
	    MeasurementContainer result = theMeasurementExtractor->measurements((*itdet), tsos,*theService->propagator(thePropagatorName), *theEstimator, iEvent);
	    
	    std::cout<<"Understanding ---> Roll and Extrapolated point: "<<"Event num: "<<nevent<<"  "<<counttrack<<"  "<<nameRoll<<"  "<<PointRollGlobal<<std::endl;


	    if(result.size() != 0) {
	      TrajectoryMeasurement bestMeasurement = result[0];
	      
	      RPCDetId rpcid1 ;
	      if(idDetLay.det() == DetId::Muon && idDetLay.subdetId() == MuonSubdetId::RPC){
		rpcid1 = RPCDetId(bestMeasurement.recHit()->geographicalId());
	      }

	      if(rpcid1.rawId() == rpcid.rawId()){
		count_goodevent++;
		count_effevent++;
		std::cout<<"Understanding --> Eff events counter: "<<count_effevent<<std::endl;

		//--------------- RecHit and efficiency section ---------------------------------

		float res=0.,minres=999.;int count = 0;int minbx = -999;
		double rhxminres; 
		int clsminres;
		RPCRecHitCollection::const_iterator recIt;
		RPCRecHitCollection::range rpcRecHitRange = rpcHits->get(rpcid);
		for (recIt = rpcRecHitRange.first; recIt!=rpcRecHitRange.second; ++recIt){
		  
		  LocalPoint rhitlocal = (*recIt).localPosition();
		  
		  double rhitposX = rhitlocal.x();  
		  int cls = (*recIt).clusterSize();
		  LocalError RecError = (*recIt).localPositionError();
		  double sigmaRec = RecError.xx();
		  res = (double)(lxextrap - rhitposX);
		  int bx = (*recIt).BunchX();

		  // 		  std::cout<<"Understanging ---> RecHit on the RPC Roll ---> "<<"Event num: "<<nevent<<"  "<<nameRoll<<"  "<<lxextrap<<"  "<<rhitposX<<"  "<<res<<std::endl;
		  
		  themyHistoClassDbe->fillPlotRoll(rawrpcid, nameRoll, bx, cls, res);
		  themyHistoClassDbe->fillDetectorPlots2D(bx, cls, res, outPt);

		  if(count == 0) {
		    minres = res;
		    clsminres = cls;
		    minbx = bx;
		  }
		  else if(count >0 && fabs(res) < fabs(minres)){
		    minres = res;
		    clsminres = cls;
		    minbx = bx;
		  }
		  count++;
		}
		
		if(count > 0 && fabs(minres) <= 1.0) efficiency1 = 1;
		if(count > 0 && fabs(minres) <= 2.0) efficiency2 = 1;
		if(count > 0 && fabs(minres) <= 3.0) efficiency3 = 1;
		if(count > 0 && fabs(minres) <= 4.0) efficiency4 = 1;
		if(count > 0 && fabs(minres) <= 7.0) efficiency5 = 1;
		if(count > 0 && fabs(minres) <= 10.0) efficiency6 = 1;
		if(count > 0 && fabs(minres) <= 13.0) efficiency7 = 1;
		if(count > 0 && fabs(minres) <= 100) efficiency8 = 1;


		if(efficiency8 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff8"<<std::endl;
		}
		if(efficiency7 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff7"<<std::endl;
		}
		if(efficiency6 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff6"<<std::endl;
		}
		if(efficiency5 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff5"<<std::endl;
		}
		if(efficiency4 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff4"<<std::endl;
		}
		if(efficiency3 == 0) {
		  std::cout<<"Understanging --->  !!!!!!!!!!!!!!"<<  "Event num: "<<nevent<<"  "<<"Roll name: "<<nameRoll<<"  "<<"MinRes: "<<fabs(minres)<< "  "<<"Eff3"<<std::endl;
		}


		themyHistoClassDbe->fillPlotRollTrack(rawrpcid, nameRoll, stripPredicted,
						      clsminres, minres, minbx, 
						      efficiency1, efficiency2, efficiency3, efficiency4, 
						      efficiency5, efficiency6, efficiency7, efficiency8);
		
		themyHistoClassDbe->fillGeneralPlots2D(mfieldX, mfieldY, mfieldZ,outP, outPt, zangle,
						       outPhi, outEta, outTheta,count, clsminres, minres, minbx);
		
	      }
	      else{
 		std::cout<<"Understanding --> The roll expected is Inefficient because is different from the got one from the measurement!"<<"  "<<nameRoll<<std::endl;
	      }
	    }
	    else {
	      std::cout<<"Understanding --> The measurement container is empty, therefore the Roll is inefficient!"<<"  "<<"Event num: "<<nevent<<"  "<<nameRoll<<std::endl;
	      count_goodevent++;
	      efficiency1 = 0;
	      efficiency2 = 0;
	      efficiency3 = 0;
	      efficiency4 = 0;
	      efficiency5 = 0;
	      efficiency6 = 0;
	      efficiency7 = 0;
	      efficiency8 = 0;
	      
	      int clsminres = 0,minbx = -999;
	      float minres = -999; 
	      
	      themyHistoClassDbe->fillPlotRollTrack(rawrpcid, nameRoll, stripPredicted,
						    clsminres, minres, minbx, 
						    efficiency1, efficiency2, efficiency3, efficiency4, 
						    efficiency5, efficiency6, efficiency7, efficiency8);
	    }

	    ftstemp = *tsosAtRPC.freeState(true);
	  }
	}
      }
    }
  }catch ( cms::Exception& er ){
    LogTrace("RPCEff")<<"caught std::exception  "<< er.what()<< std::endl;
  }
}


void RPCEffTrackExtrapolationNew::endJob(){
  if(themyHistoClassDbe)
    delete themyHistoClassDbe;

}

RPCEffTrackExtrapolationNew::~RPCEffTrackExtrapolationNew(){

  float eff = ((float)count_effevent/count_goodevent);
  std::cout<<"Understanding --> Efficiency:  "<<count_goodevent<<"  "<<count_effevent<<"  "<<eff<<std::endl;
}

DEFINE_FWK_MODULE(RPCEffTrackExtrapolationNew);
