#ifndef RPCEffTrackExtraoplationNew_h
#define RPCEffTrackExtrapolationNew_h

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

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDAnalyzer.h>
#include <FWCore/Framework/interface/Event.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "RecoMuon/MeasurementDet/interface/MuonDetLayerMeasurements.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"

#include "RecoMuon/Navigation/interface/DirectMuonNavigation.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h" 
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h" 
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "RecoMuon/TransientTrackingRecHit/interface/MuonTransientTrackingRecHit.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerBase.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


#include "./MyHistoClassDbeNew.h"

#include<string>
#include<map>
#include<fstream>

#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"

class RPCDetId;
class Trajectory;
class Propagator;
class GeomDet;
class TrajectoryStateOnSurface;


typedef std::vector<TrajectoryMeasurement>          MeasurementContainer;
typedef std::pair<const GeomDet*,TrajectoryStateOnSurface> DetWithState;
typedef std::vector<Trajectory> Trajectories;
typedef MuonTransientTrackingRecHit::MuonRecHitPointer MuonRecHitPointer;
typedef MuonTransientTrackingRecHit::MuonRecHitContainer MuonRecHitContainer;

class RPCEffTrackExtrapolationNew : public edm::EDAnalyzer {
   public:
      explicit RPCEffTrackExtrapolationNew(const edm::ParameterSet&);
      ~RPCEffTrackExtrapolationNew();
      std::vector<const DetLayer*> compatibleLayers(const DetLayer *initialLayer,
					       FreeTrajectoryState& fts,
					       PropagationDirection propDir);

   private:

      int count_goodevent;
      int count_goodeventBarrel;
      int count_goodeventEndcap;
      int count_goodeventDiskm1;
      int count_goodeventDiskm2;
      int count_goodeventDiskm3;
      int count_goodeventDisk1;
      int count_goodeventDisk2;
      int count_goodeventDisk3;

      int count_effevent;
      int count_effeventBarrel;
      int count_effeventEndcap;
      int count_effeventDiskm1;
      int count_effeventDiskm2;
      int count_effeventDiskm3;
      int count_effeventDisk1;
      int count_effeventDisk2;
      int count_effeventDisk3;



      virtual void beginJob(const edm::EventSetup&) ;
      //virtual void beginRun( edm::Run&, const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      //virtual void endRun( edm::Run&, const edm::EventSetup& );
      virtual void endJob();

      TrackTransformerBase *theTrackTransformer;

      edm::InputTag dt4DSegments;
      edm::InputTag cscSegments;
      edm::InputTag RPCRecHits;
      std::string theNavigationType;

      std::string fileName;

      TFile* fOutputFile;
      double maxRes;
      bool EffSaveRootFile;

      int wh;
      bool cosmic;

      int Run;
      time_t aTime;

      std::string root_file_name;

      unsigned int bsec;
      unsigned int esec;

      ofstream* effres;
      std::string EffRootFileName;
      std::string TjInput;
      std::string thePropagatorName;
      mutable Propagator* thePropagator;

      bool saveRootFile;
      std::string RootFileName;

      MyHistoClassDbeNew* themyHistoClassDbe;
      DQMStore * dbe;
      MuonServiceProxy* theService;
      MuonDetLayerMeasurements* theMeasurementExtractor;
      Chi2MeasurementEstimator* theEstimator;

      typedef std::vector<Trajectory> Trajectories;


};
#endif
