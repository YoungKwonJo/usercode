#ifndef MYHISTOCLASSDBENEW_H
#define MYHISTOCLASSDBENEW_H

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Common/interface/Handle.h"
#include<string>
#include<map>
#include<cmath>

#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>

// Root stuff
#include "TROOT.h"
#include "TDirectory.h"
#include "TFolder.h"
#include "TObject.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TTree.h"
#include <cstdlib>
#include <stdio.h>
#include "stdlib.h"
#include <string>
#include <sstream>
#include <memory> // needed for the auto_ptr below

class MyHistoClassDbeNew{
 public: 
  MyHistoClassDbeNew(); 
  ~MyHistoClassDbeNew();

  void init(DQMStore * dbe,const RPCGeometry *pGeom,std::string RootFileName);
  void fillPlotRollTrack(int rawrpcid, std::string nameRoll, float stripPredicted,
			 int clsminres,float minres, int minbx, float x,float y,
			 int efficiency1, int efficiency2, int efficiency3, int efficiency4, 
			 int efficiency5, int efficiency6, int efficiency7, int efficiency8);

  void fillPlotRoll(int rawrpcid, std::string nameRoll, int bx, int cls, float res);

  void fillTrackPlot(int charge, float chi2,float normalizedChi2, 
		     int numberOfValidHits, int numberOfLostHits,
		     float outerP, float outerPt, float outerPhi,
		     float outerEta, float outerTheta);

  void fillGeneralPlots2D(uint32_t rawid, float mfieldX, float mfieldY, float mfieldZ,
			  float outerP, float outerPt, float zangle, float stripw,
			  float outerPhi,float outerEta, float outerTheta,
			  int count, int clsminres,float minres, int minbx);

  void fillDetectorPlots2D(int bx, int cls, float res,float outerPt);

 private:

  DQMStore * _dbe;
  std::map<uint32_t, std::map<std::string, MonitorElement*> >  meCollection;
  std::map<std::string, std::string> _mapBWPathAndNameRoll;

  const RPCGeometry *_pGeom;

  TFile* theFile;
  std::string _rootFileName;

  MonitorElement* _hglobal_chi2;
  MonitorElement* _hglobal_num_valid_hit;
  MonitorElement* _hglobal_momentum;
  MonitorElement* _hglobal_transverse_momentum;
  MonitorElement* _hglobal_eta;
  MonitorElement* _hglobal_phi;
  MonitorElement* _hglobal_theta;
  MonitorElement* _hglobal_zangle;

  MonitorElement* _hglobal_theta_vs_minres_2d;
  MonitorElement* _hglobal_phi_vs_minres_2d;
  MonitorElement* _hglobal_theta_vs_mincls_2d;
  MonitorElement* _hglobal_phi_vs_mincls_2d;
  MonitorElement* _hglobal_minres_vs_mincls_2d;
  MonitorElement* _hglobal_minres_vs_zangler_2d;
  MonitorElement* _hglobal_mincls_vs_zangler_2d;
  MonitorElement* _hglobal_minres_vs_p_2d;
  MonitorElement* _hglobal_mincls_vs_p_2d;
  MonitorElement* _hglobal_minres_vs_pt_2d;
  MonitorElement* _hglobal_mincls_vs_pt_2d;
  MonitorElement* _hglobal_mgf_vs_pt_2d;
  MonitorElement* _hglobal_res_vs_bx;
  MonitorElement* _hglobal_res_vs_cls;
  MonitorElement* _hglobal_cls_vs_bx;
  MonitorElement* _hglobal_bx_vs_pt;
  MonitorElement* _hglobal_cls_vs_pt;

  MonitorElement* _hglobal_pt_vs_zangle;
  MonitorElement* _hglobal_cls_verticaltracks;

  MonitorElement*  _hglobal_cls_layer1;
  MonitorElement*  _hglobal_cls_layer2;
  MonitorElement*  _hglobal_cls_layer3;
  MonitorElement*  _hglobal_cls_layer4;
  MonitorElement*  _hglobal_cls_layer5;
  MonitorElement*  _hglobal_cls_layer6;

  MonitorElement*  _hglobal_cls_layer1_verticaltracks;
  MonitorElement*  _hglobal_cls_layer2_verticaltracks;
  MonitorElement*  _hglobal_cls_layer3_verticaltracks;
  MonitorElement*  _hglobal_cls_layer4_verticaltracks;
  MonitorElement*  _hglobal_cls_layer5_verticaltracks;
  MonitorElement*  _hglobal_cls_layer6_verticaltracks;

  MonitorElement* _hglobal_residuals_barrel;
  MonitorElement* _hglobal_residuals_endcap;

  MonitorElement* _hglobal_residuals_diskm3;
  MonitorElement* _hglobal_residuals_diskm2;
  MonitorElement* _hglobal_residuals_diskm1;

  MonitorElement* _hglobal_residuals_disk3;
  MonitorElement* _hglobal_residuals_disk2;
  MonitorElement* _hglobal_residuals_disk1;

  MonitorElement* _hglobal_residuals_layer1;
  MonitorElement* _hglobal_residuals_layer2;
  MonitorElement* _hglobal_residuals_layer3;
  MonitorElement* _hglobal_residuals_layer4;
  MonitorElement* _hglobal_residuals_layer5;
  MonitorElement* _hglobal_residuals_layer6;

  MonitorElement* _hglobal_residuals_layer1_cls1;
  MonitorElement* _hglobal_residuals_layer2_cls1;
  MonitorElement* _hglobal_residuals_layer3_cls1;
  MonitorElement* _hglobal_residuals_layer4_cls1;
  MonitorElement* _hglobal_residuals_layer5_cls1;
  MonitorElement* _hglobal_residuals_layer6_cls1;

  MonitorElement* _hglobal_residuals_layer1_cls2;
  MonitorElement* _hglobal_residuals_layer2_cls2;
  MonitorElement* _hglobal_residuals_layer3_cls2;
  MonitorElement* _hglobal_residuals_layer4_cls2;
  MonitorElement* _hglobal_residuals_layer5_cls2;
  MonitorElement* _hglobal_residuals_layer6_cls2;

  MonitorElement* _hglobal_residuals_layer1_cls3;
  MonitorElement* _hglobal_residuals_layer2_cls3;
  MonitorElement* _hglobal_residuals_layer3_cls3;
  MonitorElement* _hglobal_residuals_layer4_cls3;
  MonitorElement* _hglobal_residuals_layer5_cls3;
  MonitorElement* _hglobal_residuals_layer6_cls3;

  MonitorElement* _DistBorderClu1La1;
  MonitorElement* _DistBorderClu1La2;
  MonitorElement* _DistBorderClu1La3;
  MonitorElement* _DistBorderClu1La4;
  MonitorElement* _DistBorderClu1La5;
  MonitorElement* _DistBorderClu1La6;

  MonitorElement*  _DistBorderClu2La1;
  MonitorElement*  _DistBorderClu2La2;
  MonitorElement*  _DistBorderClu2La3;
  MonitorElement*  _DistBorderClu2La4;
  MonitorElement*  _DistBorderClu2La5;
  MonitorElement*  _DistBorderClu2La6;

  MonitorElement*  _DistBorderClu3La1;
  MonitorElement*  _DistBorderClu3La2;
  MonitorElement*  _DistBorderClu3La3;
  MonitorElement*  _DistBorderClu3La4;
  MonitorElement*  _DistBorderClu3La5;
  MonitorElement*  _DistBorderClu3La6;

};

MyHistoClassDbeNew::MyHistoClassDbeNew() {}

MyHistoClassDbeNew::~MyHistoClassDbeNew() {

  _dbe->save(_rootFileName); 
  
}

void MyHistoClassDbeNew::init(DQMStore * dbe, const RPCGeometry *pGeom,std::string RootFileName){
  _pGeom = pGeom;
  _dbe = dbe;
  _rootFileName = RootFileName;


  _dbe->setCurrentFolder("STAMuonTrackEff/GeneralTrackPlots");  
  _hglobal_chi2 = _dbe->book1D("hglobal_chi2","hglobal_chi2",200,0,1000);
  _hglobal_num_valid_hit = _dbe->book1D("hglobal_num_valid_hit","hglobal_num_valid_hit",100,0,100);
  _hglobal_momentum = _dbe->book1D("hglobal_momentum","hglobal_momentum",100,0,100);
  _hglobal_transverse_momentum = _dbe->book1D("hglobal_transverse_momentum","hglobal_transverse_momentum",100,0,100);
  _hglobal_eta = _dbe->book1D("_hglobal_eta","_hglobal_eta",50,-2.4,2.4);
  _hglobal_phi = _dbe->book1D("_hglobal_phi","_hglobal_phi",50,-3.14,3.14);
  _hglobal_theta = _dbe->book1D("_hglobal_theta","_hglobal_theta",50,-3.14,3.14);
  _hglobal_zangle = _dbe->book1D("_hglobal_zangle","_hglobal_zangle",50,-3.14,3.14);

  _hglobal_cls_layer1 =  _dbe->book1D("_hglobal_cls_layer1","_hglobal_cls_layer1",30,0.5,30.5);
  _hglobal_cls_layer2 =  _dbe->book1D("_hglobal_cls_layer2","_hglobal_cls_layer2",30,0.5,30.5);
  _hglobal_cls_layer3 =  _dbe->book1D("_hglobal_cls_layer3","_hglobal_cls_layer3",30,0.5,30.5);
  _hglobal_cls_layer4 =  _dbe->book1D("_hglobal_cls_layer4","_hglobal_cls_layer4",30,0.5,30.5);
  _hglobal_cls_layer5 =  _dbe->book1D("_hglobal_cls_layer5","_hglobal_cls_layer5",30,0.5,30.5);
  _hglobal_cls_layer6 =  _dbe->book1D("_hglobal_cls_layer6","_hglobal_cls_layer6",30,0.5,30.5);

  _hglobal_cls_verticaltracks = _dbe->book1D("_hglobal_cls_verticaltracks","_hglobal_cls_verticaltracks",30,0.5,30.5);

  _hglobal_cls_layer1_verticaltracks =  _dbe->book1D("_hglobal_cls_layer1_verticaltracks","_hglobal_cls_layer1_verticaltracks",30,0.5,30.5);
  _hglobal_cls_layer2_verticaltracks =  _dbe->book1D("_hglobal_cls_layer2_verticaltracks","_hglobal_cls_layer2_verticaltracks",30,0.5,30.5);
  _hglobal_cls_layer3_verticaltracks =  _dbe->book1D("_hglobal_cls_layer3_verticaltracks","_hglobal_cls_layer3_verticaltracks",30,0.5,30.5);
  _hglobal_cls_layer4_verticaltracks =  _dbe->book1D("_hglobal_cls_layer4_verticaltracks","_hglobal_cls_layer4_verticaltracks",30,0.5,30.5);
  _hglobal_cls_layer5_verticaltracks =  _dbe->book1D("_hglobal_cls_layer5_verticaltracks","_hglobal_cls_layer5_verticaltracks",30,0.5,30.5);
  _hglobal_cls_layer6_verticaltracks =  _dbe->book1D("_hglobal_cls_layer6_verticaltracks","_hglobal_cls_layer6_verticaltracks",30,0.5,30.5);

  _hglobal_residuals_barrel =  _dbe->book1D("_hglobal_residuals_barrel","_hglobal_residuals_barrel",101,-50.,50.);
  _hglobal_residuals_endcap =  _dbe->book1D("_hglobal_residuals_endcap","_hglobal_residuals_endcap",101,-50.,50.);

  _hglobal_residuals_diskm3 =  _dbe->book1D("_hglobal_residuals_diskm3","_hglobal_residuals_diskm3",101,-50.,50.);
  _hglobal_residuals_diskm2 =  _dbe->book1D("_hglobal_residuals_diskm2","_hglobal_residuals_diskm2",101,-50.,50.);
  _hglobal_residuals_diskm1 =  _dbe->book1D("_hglobal_residuals_diskm1","_hglobal_residuals_diskm1",101,-50.,50.);

  _hglobal_residuals_disk3 =  _dbe->book1D("_hglobal_residuals_disk3","_hglobal_residuals_disk3",101,-50.,50.);
  _hglobal_residuals_disk2 =  _dbe->book1D("_hglobal_residuals_disk2","_hglobal_residuals_disk2",101,-50.,50.);
  _hglobal_residuals_disk1 =  _dbe->book1D("_hglobal_residuals_disk1","_hglobal_residuals_disk1",101,-50.,50.);


  _hglobal_residuals_layer1 =  _dbe->book1D("_hglobal_residuals_layer1","_hglobal_residuals_layer1",101,-50.,50.);
  _hglobal_residuals_layer2 =  _dbe->book1D("_hglobal_residuals_layer2","_hglobal_residuals_layer2",101,-50.,50.);
  _hglobal_residuals_layer3 =  _dbe->book1D("_hglobal_residuals_layer3","_hglobal_residuals_layer3",101,-50.,50.);
  _hglobal_residuals_layer4 =  _dbe->book1D("_hglobal_residuals_layer4","_hglobal_residuals_layer4",101,-50.,50.);
  _hglobal_residuals_layer5 =  _dbe->book1D("_hglobal_residuals_layer5","_hglobal_residuals_layer5",101,-50.,50.);
  _hglobal_residuals_layer6 =  _dbe->book1D("_hglobal_residuals_layer6","_hglobal_residuals_layer6",101,-50.,50.);

  _hglobal_residuals_layer1_cls1 =  _dbe->book1D("_hglobal_residuals_layer1_cls1","_hglobal_residuals_layer1_cls1",101,-50.,50.);
  _hglobal_residuals_layer2_cls1 =  _dbe->book1D("_hglobal_residuals_layer2_cls1","_hglobal_residuals_layer2_cls1",101,-50.,50.);
  _hglobal_residuals_layer3_cls1 =  _dbe->book1D("_hglobal_residuals_layer3_cls1","_hglobal_residuals_layer3_cls1",101,-50.,50.);
  _hglobal_residuals_layer4_cls1 =  _dbe->book1D("_hglobal_residuals_layer4_cls1","_hglobal_residuals_layer4_cls1",101,-50.,50.);
  _hglobal_residuals_layer5_cls1 =  _dbe->book1D("_hglobal_residuals_layer5_cls1","_hglobal_residuals_layer5_cls1",101,-50.,50.);
  _hglobal_residuals_layer6_cls1 =  _dbe->book1D("_hglobal_residuals_layer6_cls1","_hglobal_residuals_layer6_cls1",101,-50.,50.);

  _hglobal_residuals_layer1_cls2 =  _dbe->book1D("_hglobal_residuals_layer1_cls2","_hglobal_residuals_layer1_cls2",101,-50.,50.);
  _hglobal_residuals_layer2_cls2 =  _dbe->book1D("_hglobal_residuals_layer2_cls2","_hglobal_residuals_layer2_cls2",101,-50.,50.);
  _hglobal_residuals_layer3_cls2 =  _dbe->book1D("_hglobal_residuals_layer3_cls2","_hglobal_residuals_layer3_cls2",101,-50.,50.);
  _hglobal_residuals_layer4_cls2 =  _dbe->book1D("_hglobal_residuals_layer4_cls2","_hglobal_residuals_layer4_cls2",101,-50.,50.);
  _hglobal_residuals_layer5_cls2 =  _dbe->book1D("_hglobal_residuals_layer5_cls2","_hglobal_residuals_layer5_cls2",101,-50.,50.);
  _hglobal_residuals_layer6_cls2 =  _dbe->book1D("_hglobal_residuals_layer6_cls2","_hglobal_residuals_layer6_cls2",101,-50.,50.);

  _hglobal_residuals_layer1_cls3 =  _dbe->book1D("_hglobal_residuals_layer1_cls3","_hglobal_residuals_layer1_cls3",101,-50.,50.);
  _hglobal_residuals_layer2_cls3 =  _dbe->book1D("_hglobal_residuals_layer2_cls3","_hglobal_residuals_layer2_cls3",101,-50.,50.);
  _hglobal_residuals_layer3_cls3 =  _dbe->book1D("_hglobal_residuals_layer3_cls3","_hglobal_residuals_layer3_cls3",101,-50.,50.);
  _hglobal_residuals_layer4_cls3 =  _dbe->book1D("_hglobal_residuals_layer4_cls3","_hglobal_residuals_layer4_cls3",101,-50.,50.);
  _hglobal_residuals_layer5_cls3 =  _dbe->book1D("_hglobal_residuals_layer5_cls3","_hglobal_residuals_layer5_cls3",101,-50.,50.);
  _hglobal_residuals_layer6_cls3 =  _dbe->book1D("_hglobal_residuals_layer6_cls3","_hglobal_residuals_layer6_cls3",101,-50.,50.);

  _DistBorderClu1La1 = dbe->book1D("DistBorderClu1La1","Distance to the Border of the Strip Layer 1 Cluster Size 1",200,-2.,3.);
  _DistBorderClu1La2 = dbe->book1D("DistBorderClu1La2","Distance to the Border of the Strip Layer 2 Cluster Size 1",200,-2.,3.);
  _DistBorderClu1La3 = dbe->book1D("DistBorderClu1La3","Distance to the Border of the Strip Layer 3 Cluster Size 1",200,-2.,3.);
  _DistBorderClu1La4 = dbe->book1D("DistBorderClu1La4","Distance to the Border of the Strip Layer 4 Cluster Size 1",200,-2.,3.);
  _DistBorderClu1La5 = dbe->book1D("DistBorderClu1La5","Distance to the Border of the Strip Layer 5 Cluster Size 1",200,-2.,3.);
  _DistBorderClu1La6 = dbe->book1D("DistBorderClu1La6","Distance to the Border of the Strip Layer 6 Cluster Size 1",200,-2.,3.);

  _DistBorderClu2La1 = dbe->book1D("DistBorderClu2La1","Distance to the Border of the Strip Layer 1 Cluster Size 2",200,-2.,3.);
  _DistBorderClu2La2 = dbe->book1D("DistBorderClu2La2","Distance to the Border of the Strip Layer 2 Cluster Size 2",200,-2.,3.);
  _DistBorderClu2La3 = dbe->book1D("DistBorderClu2La3","Distance to the Border of the Strip Layer 3 Cluster Size 2",200,-2.,3.);
  _DistBorderClu2La4 = dbe->book1D("DistBorderClu2La4","Distance to the Border of the Strip Layer 4 Cluster Size 2",200,-2.,3.);
  _DistBorderClu2La5 = dbe->book1D("DistBorderClu2La5","Distance to the Border of the Strip Layer 5 Cluster Size 2",200,-2.,3.);
  _DistBorderClu2La6 = dbe->book1D("DistBorderClu2La6","Distance to the Border of the Strip Layer 6 Cluster Size 2",200,-2.,3.);

  _DistBorderClu3La1 = dbe->book1D("DistBorderClu3La1","Distance to the Border of the Strip Layer 1 Cluster Size 3",200,-2.,3.);
  _DistBorderClu3La2 = dbe->book1D("DistBorderClu3La2","Distance to the Border of the Strip Layer 2 Cluster Size 3",200,-2.,3.);
  _DistBorderClu3La3 = dbe->book1D("DistBorderClu3La3","Distance to the Border of the Strip Layer 3 Cluster Size 3",200,-2.,3.);
  _DistBorderClu3La4 = dbe->book1D("DistBorderClu3La4","Distance to the Border of the Strip Layer 4 Cluster Size 3",200,-2.,3.);
  _DistBorderClu3La5 = dbe->book1D("DistBorderClu3La5","Distance to the Border of the Strip Layer 5 Cluster Size 3",200,-2.,3.);
  _DistBorderClu3La6 = dbe->book1D("DistBorderClu3La6","Distance to the Border of the Strip Layer 6 Cluster Size 3",200,-2.,3.);

  _dbe->setCurrentFolder("STAMuonTrackEff/GeneralPlots2D");

  _hglobal_theta_vs_minres_2d = _dbe->book2D("hglobal_theta_vs_minres_2d","hglobal_theta_vs_minres_2d",200,-100,100,50,-3.14,3.14);
  _hglobal_phi_vs_minres_2d = _dbe->book2D("hglobal_phi_vs_minres_2d","hglobal_phi_vs_minres_2d",200,-100,100,50,-3.14,3.14);
  _hglobal_theta_vs_mincls_2d = _dbe->book2D("hglobal_theta_vs_mincls_2d","hglobal_theta_vs_mincls_2d",30,0.5,30.5,50,-3.14,3.14);
  _hglobal_phi_vs_mincls_2d = _dbe->book2D("hglobal_theta_vs_mincls_2d","hglobal_theta_vs_mincls_2d",30,0.5,30.5,50,-3.14,3.14);
  _hglobal_minres_vs_mincls_2d = _dbe->book2D("hglobal_minres_vs_mincls_2d","hglobal_minres_vs_mincls_2d",30,0.5,30.5,200,-100,100);
  _hglobal_minres_vs_zangler_2d = _dbe->book2D("hglobal_minres_vs_zangler_2d","hglobal_minres_vs_zangler_2d",50,-3.14,3.14,200,-100,100);
  _hglobal_mincls_vs_zangler_2d = _dbe->book2D("hglobal_mincls_vs_zangler_2d","hglobal_mincls_vs_zangler_2d",50,-3.14,3.14,30,0.5,30.5);
  _hglobal_minres_vs_p_2d = _dbe->book2D("hglobal_minres_vs_p_2d","hglobal_minres_vs_p_2d",100,0.,200.,200,-100,100);
  _hglobal_mincls_vs_p_2d = _dbe->book2D("hglobal_mincls_vs_p_2d","hglobal_mincls_vs_p_2d",100,0.,200.,30,0.5,30.5);
  _hglobal_minres_vs_pt_2d = _dbe->book2D("hglobal_minres_vs_pt_2d","hglobal_minres_vs_pt_2d",100,0.,200.,200,-100,100);
  _hglobal_mincls_vs_pt_2d = _dbe->book2D("hglobal_mincls_vs_pt_2d","hglobal_mincls_vs_pt_2d",100,0.,200.,30,0.5,30.5);
  _hglobal_mgf_vs_pt_2d = _dbe->book2D("hglobal_mgf_vs_pt_2d","hglobal_mgf_vs_pt_2d",100,0.,200.,100,0.,100.);

  _hglobal_res_vs_bx = _dbe->book2D("hglobal_res_vs_bx","hglobal_res_vs_bx",15,-7.5,7.5,200,-100,100);
  _hglobal_res_vs_cls = _dbe->book2D("hglobal_res_vs_cls","hglobal_res_vs_cls",30,0.5,30.5,200,-100,100);
  _hglobal_cls_vs_bx = _dbe->book2D("hglobal_cls_vs_bx","hglobal_cls_vs_bx",15,-7.5,7.5,30,0.5,30.5);
  _hglobal_bx_vs_pt = _dbe->book2D("hglobal_bx_vs_pt","hglobal_bx_vs_pt",100,0.,100.,15,-7.5,7.5);
  _hglobal_cls_vs_pt = _dbe->book2D("hglobal_cls_vs_pt","hglobal_cls_vs_pt",100,0.,100.,30,0.5,30.5);
  _hglobal_pt_vs_zangle = _dbe->book2D("hglobal_pt_vs_zangle","hglobal_pt_vs_zangle",50,-3.14,3.14,100,0.,100.);

  std::map<std::string, MonitorElement*> mapBWMeAndNameHisto;

  float scale2D = 0.6;
  float stripw = 0.;
  float stripl = 0.;
  float nstrips = 0.;

  for(TrackingGeometry::DetContainer::const_iterator it = pGeom->dets().begin(); it != pGeom->dets().end(); it++){
    
    if( dynamic_cast< RPCChamber* >( *it ) != 0 ){
      RPCChamber* ch = dynamic_cast< RPCChamber* >( *it ); 
      
      std::vector< const RPCRoll*> rollsRaf = (ch->rolls());
      for(std::vector<const RPCRoll*>::iterator r = rollsRaf.begin();
	  r != rollsRaf.end(); ++r){
	  
	RPCGeomServ rpcsrv((*r)->id());
	std::string nameRoll = rpcsrv.name();
	std::vector<std::string> vname;

	nstrips = (float)(*r)->nstrips();
	
	if((*r)->isBarrel()){
	  const RectangularStripTopology* top_= dynamic_cast<const RectangularStripTopology*> (&((*r)->topology()));
	  stripw = top_->pitch();
	  stripl = top_->stripLength();

	}
	else{
	  const TrapezoidalStripTopology* top_= dynamic_cast<const TrapezoidalStripTopology*> (&((*r)->topology()));
	  //stripw = top_->pitch();
	}

	std::string path = "STAMuonTrackEff/PlotsRolls/";

	if ( ((*r)->id()).region() == 0 ) {
	  path += Form( ( ((*r)->id()).ring() == 0 ? "W%d/" : "W%+d/" ), ((*r)->id()).ring() );
	  path += Form("Sec%d/Station%d",((*r)->id()).sector(),((*r)->id()).station());
	}
	else{
	  path += Form( "D%+d/" , ((*r)->id()).ring() );
	  path += Form("Sec%d/Station%d",((*r)->id()).sector(),((*r)->id()).station());
	}

	_dbe->setCurrentFolder(path);  

	std::string namehistoExpOcc = "ExpOcc_"+nameRoll;
	
	std::string namehistoRealOcc1 = "RealOcc1_"+nameRoll;
	std::string namehistoRealOcc2 = "RealOcc2_"+nameRoll;
	std::string namehistoRealOcc3 = "RealOcc3_"+nameRoll;
	std::string namehistoRealOcc4 = "RealOcc4_"+nameRoll;
	std::string namehistoRealOcc5 = "RealOcc5_"+nameRoll;
	std::string namehistoRealOcc6 = "RealOcc6_"+nameRoll;
	std::string namehistoRealOcc7 = "RealOcc7_"+nameRoll;
	std::string namehistoRealOcc8 = "RealOcc8_"+nameRoll;
	
	std::string namehistoExpOcc2D = "2DExpOcc_"+nameRoll;
	std::string namehistoRealOcc2D = "2DRealOcc_"+nameRoll;

	std::string namehistoCls = "Cls_"+nameRoll;
	std::string namehistoClsmin = "Clsmin_"+nameRoll;
	std::string namehistoBx = "Bx_"+nameRoll;
	std::string namehistoBxmin = "Bxmin_"+nameRoll;
	std::string namehistoResmin = "Residualsmin_"+nameRoll;
	std::string namehistoRes = "Residuals_"+nameRoll;
	std::string namehistoResCls1 = "Residuals_cls1_"+nameRoll;
	std::string namehistoResCls2 = "Residuals_cls2_"+nameRoll;
	std::string namehistoResCls3 = "Residuals_cls3_"+nameRoll;

	mapBWMeAndNameHisto.clear();
	
	mapBWMeAndNameHisto[namehistoExpOcc] = _dbe->book1D(namehistoExpOcc.c_str(),namehistoExpOcc.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc1] = _dbe->book1D(namehistoRealOcc1.c_str(),namehistoRealOcc1.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc2] = _dbe->book1D(namehistoRealOcc2.c_str(),namehistoRealOcc2.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc3] = _dbe->book1D(namehistoRealOcc3.c_str(),namehistoRealOcc3.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc4] = _dbe->book1D(namehistoRealOcc4.c_str(),namehistoRealOcc4.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc5] = _dbe->book1D(namehistoRealOcc5.c_str(),namehistoRealOcc5.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc6] = _dbe->book1D(namehistoRealOcc6.c_str(),namehistoRealOcc6.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc7] = _dbe->book1D(namehistoRealOcc7.c_str(),namehistoRealOcc7.c_str(),100,0.,100.);
	mapBWMeAndNameHisto[namehistoRealOcc8] = _dbe->book1D(namehistoRealOcc8.c_str(),namehistoRealOcc8.c_str(),100,0.,100.);

	mapBWMeAndNameHisto[namehistoExpOcc2D] = 
	  _dbe->book2D(namehistoExpOcc2D.c_str(),namehistoExpOcc2D.c_str(),2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);	
		
	mapBWMeAndNameHisto[namehistoRealOcc2D] = 
	  _dbe->book2D(namehistoRealOcc2D.c_str(),namehistoRealOcc2D.c_str(),2*nstrips,-scale2D*nstrips*stripw,scale2D*nstrips*stripw,2*nstrips,-scale2D*stripl,scale2D*stripl);

	mapBWMeAndNameHisto[namehistoCls] = _dbe->book1D(namehistoCls.c_str(),namehistoCls.c_str(),30,0.5,30.5);
	mapBWMeAndNameHisto[namehistoClsmin] = _dbe->book1D(namehistoClsmin.c_str(),namehistoClsmin.c_str(),30,0.5,30.5);
	mapBWMeAndNameHisto[namehistoBx] = _dbe->book1D(namehistoBx.c_str(),namehistoBx.c_str(),19,-9.5,9.5);
	mapBWMeAndNameHisto[namehistoBxmin] = _dbe->book1D(namehistoBxmin.c_str(),namehistoBxmin.c_str(),19,-9.5,9.5);
	mapBWMeAndNameHisto[namehistoRes] = _dbe->book1D(namehistoRes.c_str(),namehistoRes.c_str(),1000,-100,100);
	mapBWMeAndNameHisto[namehistoResmin] = _dbe->book1D(namehistoResmin.c_str(),namehistoResmin.c_str(),101,-7.*stripw,7*stripw);
	mapBWMeAndNameHisto[namehistoResCls1] = _dbe->book1D(namehistoResCls1.c_str(),namehistoResCls1.c_str(),101,-7.*stripw,7*stripw);
	mapBWMeAndNameHisto[namehistoResCls2] = _dbe->book1D(namehistoResCls2.c_str(),namehistoResCls2.c_str(),101,-7.*stripw,7*stripw);
	mapBWMeAndNameHisto[namehistoResCls3] = _dbe->book1D(namehistoResCls3.c_str(),namehistoResCls3.c_str(),101,-7.*stripw,7*stripw);

	meCollection[(*r)->id()] = mapBWMeAndNameHisto;
      }
    }
  }
}

void MyHistoClassDbeNew::fillGeneralPlots2D(uint32_t rawid,float mfieldX, float mfieldY, float mfieldZ,
					    float outerP, float outerPt, float zangle,float stripw,
					    float outerPhi,float outerEta, float outerTheta,
					    int count, int clsminres,float minres, int minbx){
    
  
  RPCDetId rpcid(rawid);

  if(rpcid.region() != 0) {
    _hglobal_residuals_endcap->Fill(minres);

    if(rpcid.region() == -1 && rpcid.station() == 1) _hglobal_residuals_diskm1->Fill(minres);
    if(rpcid.region() == -1 && rpcid.station() == 2) _hglobal_residuals_diskm2->Fill(minres);
    if(rpcid.region() == -1 && rpcid.station() == 3) _hglobal_residuals_diskm3->Fill(minres);

    if(rpcid.region() == 1 && rpcid.station() == 1) _hglobal_residuals_disk1->Fill(minres);
    if(rpcid.region() == 1 && rpcid.station() == 2) _hglobal_residuals_disk2->Fill(minres);
    if(rpcid.region() == 1 && rpcid.station() == 3) _hglobal_residuals_disk3->Fill(minres);

  }

  if(rpcid.region() == 0 && rpcid.sector() != 1 && rpcid.sector() != 7){

    _hglobal_theta_vs_minres_2d->Fill(minres,outerTheta);
    _hglobal_phi_vs_minres_2d->Fill(minres,outerPhi);
    _hglobal_theta_vs_mincls_2d->Fill(clsminres,outerTheta);
    _hglobal_phi_vs_mincls_2d->Fill(clsminres,outerPhi); 
    _hglobal_minres_vs_mincls_2d->Fill(clsminres,minres); 
    _hglobal_minres_vs_zangler_2d->Fill(zangle,minres);
    _hglobal_mincls_vs_zangler_2d->Fill(zangle,clsminres);
    _hglobal_minres_vs_p_2d->Fill(outerP,minres);
    _hglobal_mincls_vs_p_2d->Fill(outerP,clsminres);
    _hglobal_minres_vs_pt_2d->Fill(outerPt,minres);
    _hglobal_mincls_vs_pt_2d->Fill(outerPt,clsminres);
    _hglobal_mgf_vs_pt_2d->Fill(outerPt,mfieldZ);
    _hglobal_zangle->Fill(zangle);
    
    _hglobal_cls_verticaltracks->Fill(clsminres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1) _hglobal_cls_layer1->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2) _hglobal_cls_layer2->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1) _hglobal_cls_layer3->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2) _hglobal_cls_layer4->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1) _hglobal_cls_layer5->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1) _hglobal_cls_layer6->Fill(clsminres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer1_verticaltracks->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer2_verticaltracks->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer3_verticaltracks->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer4_verticaltracks->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer5_verticaltracks->Fill(clsminres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && fabs(cos(zangle)) < 0.34 ) _hglobal_cls_layer6_verticaltracks->Fill(clsminres);
    
    if(rpcid.region() == 0) _hglobal_residuals_barrel->Fill(minres);


    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1) _hglobal_residuals_layer1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2) _hglobal_residuals_layer2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1) _hglobal_residuals_layer3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2) _hglobal_residuals_layer4->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1) _hglobal_residuals_layer5->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1) _hglobal_residuals_layer6->Fill(minres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 1) _hglobal_residuals_layer1_cls1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 1) _hglobal_residuals_layer2_cls1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 1) _hglobal_residuals_layer3_cls1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 1) _hglobal_residuals_layer4_cls1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 1) _hglobal_residuals_layer5_cls1->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 1) _hglobal_residuals_layer6_cls1->Fill(minres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 2) _hglobal_residuals_layer1_cls2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 2) _hglobal_residuals_layer2_cls2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 2) _hglobal_residuals_layer3_cls2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 2) _hglobal_residuals_layer4_cls2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 2) _hglobal_residuals_layer5_cls2->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 2) _hglobal_residuals_layer6_cls2->Fill(minres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 3) _hglobal_residuals_layer1_cls3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 3) _hglobal_residuals_layer2_cls3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 3) _hglobal_residuals_layer3_cls3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 3) _hglobal_residuals_layer4_cls3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 3) _hglobal_residuals_layer5_cls3->Fill(minres);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 3) _hglobal_residuals_layer6_cls3->Fill(minres);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 1) _DistBorderClu1La1->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 1) _DistBorderClu1La2->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 1) _DistBorderClu1La3->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 1) _DistBorderClu1La4->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 1) _DistBorderClu1La5->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 1) _DistBorderClu1La6->Fill((minres/stripw)+0.5);
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 2) _DistBorderClu2La1->Fill((minres/stripw));
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 2) _DistBorderClu2La2->Fill((minres/stripw));
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 2) _DistBorderClu2La3->Fill((minres/stripw));
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 2) _DistBorderClu2La4->Fill((minres/stripw));
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 2) _DistBorderClu2La5->Fill((minres/stripw));
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 2) _DistBorderClu2La6->Fill((minres/stripw));
    
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 1 && clsminres == 3) _DistBorderClu3La1->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 1 && rpcid.layer() == 2 && clsminres == 3) _DistBorderClu3La2->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 1 && clsminres == 3) _DistBorderClu3La3->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 2 && rpcid.layer() == 2 && clsminres == 3) _DistBorderClu3La4->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 3 && rpcid.layer() == 1 && clsminres == 3) _DistBorderClu3La5->Fill((minres/stripw)+0.5);
    if(rpcid.region() == 0 && rpcid.station() == 4 && rpcid.layer() == 1 && clsminres == 3) _DistBorderClu3La6->Fill((minres/stripw)+0.5);
  }
}

void MyHistoClassDbeNew::fillDetectorPlots2D(int bx, int cls, float res,float outerPt){

  _hglobal_res_vs_bx->Fill(bx,res);
  _hglobal_res_vs_cls->Fill(cls,res);
  _hglobal_cls_vs_bx->Fill(bx,cls);
  _hglobal_bx_vs_pt->Fill(outerPt,bx);
  _hglobal_cls_vs_pt->Fill(outerPt,cls);

}


void MyHistoClassDbeNew::fillTrackPlot(int charge, float chi2,float normalizedChi2, 
				    int numberOfValidHits, int numberOfLostHits,
				    float outerP, float outerPt, float outerPhi,
				    float outerEta, float outerTheta){

  _hglobal_chi2->Fill(chi2);
  _hglobal_num_valid_hit->Fill(numberOfValidHits);
  _hglobal_momentum->Fill(outerP);
  _hglobal_transverse_momentum->Fill(outerPt);
  _hglobal_eta->Fill(outerEta);
  _hglobal_phi->Fill(outerPhi);
  _hglobal_theta->Fill(outerTheta);

}


void MyHistoClassDbeNew::fillPlotRollTrack(int rawrpcid, std::string nameRoll, float stripPredicted,
					   int clsminres,float minres, int minbx, float x, float y,
					   int efficiency1, int efficiency2, int efficiency3, int efficiency4, 
					   int efficiency5, int efficiency6, int efficiency7, int efficiency8){


  std::map<std::string, MonitorElement*> mapBWMeAndNameHisto = meCollection[rawrpcid];
  
  std::string namehistoExpOcc = "ExpOcc_"+nameRoll;
  
  std::string namehistoExpOcc2D = "2DExpOcc_"+nameRoll;
  std::string namehistoRealOcc2D = "2DRealOcc_"+nameRoll;

  std::string namehistoRealOcc1 = "RealOcc1_"+nameRoll;
  std::string namehistoRealOcc2 = "RealOcc2_"+nameRoll;
  std::string namehistoRealOcc3 = "RealOcc3_"+nameRoll;
  std::string namehistoRealOcc4 = "RealOcc4_"+nameRoll;
  std::string namehistoRealOcc5 = "RealOcc5_"+nameRoll;
  std::string namehistoRealOcc6 = "RealOcc6_"+nameRoll;
  std::string namehistoRealOcc7 = "RealOcc7_"+nameRoll;
  std::string namehistoRealOcc8 = "RealOcc8_"+nameRoll;
  
  std::string namehistoClsmin = "Clsmin_"+nameRoll;
  std::string namehistoBxmin = "Bxmin_"+nameRoll;
  std::string namehistoResmin = "Residualsmin_"+nameRoll;

  if(mapBWMeAndNameHisto.find(namehistoExpOcc) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoExpOcc]->Fill(stripPredicted);

  if(mapBWMeAndNameHisto.find(namehistoExpOcc2D) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoExpOcc2D]->Fill(x,y);

  if(mapBWMeAndNameHisto.find(namehistoRealOcc1) != mapBWMeAndNameHisto.end() && efficiency1 == 1) mapBWMeAndNameHisto[namehistoRealOcc1]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc2) != mapBWMeAndNameHisto.end() && efficiency2 == 1) mapBWMeAndNameHisto[namehistoRealOcc2]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc3) != mapBWMeAndNameHisto.end() && efficiency3 == 1) mapBWMeAndNameHisto[namehistoRealOcc3]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc4) != mapBWMeAndNameHisto.end() && efficiency4 == 1) mapBWMeAndNameHisto[namehistoRealOcc4]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc5) != mapBWMeAndNameHisto.end() && efficiency5 == 1) mapBWMeAndNameHisto[namehistoRealOcc5]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc6) != mapBWMeAndNameHisto.end() && efficiency6 == 1) mapBWMeAndNameHisto[namehistoRealOcc6]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc7) != mapBWMeAndNameHisto.end() && efficiency7 == 1) mapBWMeAndNameHisto[namehistoRealOcc7]->Fill(stripPredicted);
  if(mapBWMeAndNameHisto.find(namehistoRealOcc8) != mapBWMeAndNameHisto.end() && efficiency8 == 1) mapBWMeAndNameHisto[namehistoRealOcc8]->Fill(stripPredicted);

  if(mapBWMeAndNameHisto.find(namehistoRealOcc2D) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoRealOcc2D]->Fill(x,y);

  if(mapBWMeAndNameHisto.find(namehistoClsmin) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoClsmin]->Fill(clsminres);
  if(mapBWMeAndNameHisto.find(namehistoBxmin) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoBxmin]->Fill(minbx);
  if(mapBWMeAndNameHisto.find(namehistoResmin) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoResmin]->Fill(minres);




}

void MyHistoClassDbeNew::fillPlotRoll(int rawrpcid, std::string nameRoll, int bx, int cls, float res){
  std::string namehistoCls = "Cls_"+nameRoll;
  std::string namehistoBx = "Bx_"+nameRoll;
  std::string namehistoRes = "Residuals_"+nameRoll;
  std::string namehistoResCls1 = "Residuals_cls1_"+nameRoll;
  std::string namehistoResCls2 = "Residuals_cls2_"+nameRoll;
  std::string namehistoResCls3 = "Residuals_cls3_"+nameRoll;

  std::map<std::string, MonitorElement*> mapBWMeAndNameHisto = meCollection[rawrpcid];
  if(mapBWMeAndNameHisto.find(namehistoCls) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoCls]->Fill(cls);
  if(mapBWMeAndNameHisto.find(namehistoBx) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoBx]->Fill(bx);
  if(mapBWMeAndNameHisto.find(namehistoRes) != mapBWMeAndNameHisto.end()) mapBWMeAndNameHisto[namehistoRes]->Fill(res);
  if(mapBWMeAndNameHisto.find(namehistoResCls1) != mapBWMeAndNameHisto.end() && cls == 1) mapBWMeAndNameHisto[namehistoResCls1]->Fill(res);
  if(mapBWMeAndNameHisto.find(namehistoResCls2) != mapBWMeAndNameHisto.end() && cls == 2) mapBWMeAndNameHisto[namehistoResCls2]->Fill(res);
  if(mapBWMeAndNameHisto.find(namehistoResCls3) != mapBWMeAndNameHisto.end() && cls == 3) mapBWMeAndNameHisto[namehistoResCls3]->Fill(res);
}




#endif

