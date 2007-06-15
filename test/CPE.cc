// File: CPE.cc
// Description:  see CPE.h
// Author:  C.Genta
//
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

#include "RecoLocalTracker/SiStripRecHitConverter/test/CPE.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoLocalTracker/Records/interface/TkStripCPERecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"

CPE::CPE(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}


// Virtual destructor needed.
CPE::~CPE() { }  

void CPE::beginJob(const edm::EventSetup& c) {  
  myfile_= new TFile("resolution.root","RECREATE");
  tibreshisto=new TProfile("TIBresvsproj","TIB res vs trackproj",NBINS,0,8);
  tobreshisto=new TProfile("TOBresvsproj","TOB res vs trackproj",NBINS,0,8);

  mtibreshistox=new TProfile("mTIBresvsprojx","matched x TIB res vs trackproj",NBINS,0,8);
  mtobreshistox=new TProfile("mTOBresvsprojx","matched x TOB res vs trackproj",NBINS,0,8);
  mtibreshistoy=new TProfile("mTIBresvsprojy","matched y TIB res vs trackproj",NBINS,0,8);
  mtobreshistoy=new TProfile("mTOBresvsprojy","matched y TOB res vs trackproj",NBINS,0,8);
  tiberrhisto=new TProfile("TIBerrvsproj","TIB err vs trackproj",NBINS,0,8);
  toberrhisto=new TProfile("TOBerrvsproj","TOB err vs trackproj",NBINS,0,8);
  tibsqrterrhisto=new TProfile("TIBsqrterrvsproj","TIB sqrt err vs trackproj",NBINS,0,8);
  tobsqrterrhisto=new TProfile("TOBsqrterrvsproj","TOB sqrt err vs trackproj",NBINS,0,8);
  mtibsqrterrhistox=new TProfile("mTIBsqrterrvsprojx","matched x TIB sqrt err vs trackproj",NBINS,0,8);
  mtobsqrterrhistox=new TProfile("mTOBsqrterrvsprojx","matched x TOB sqrt err vs trackproj",NBINS,0,8);
  mtibsqrterrhistoy=new TProfile("mTIBsqrterrvsprojy","matched y TIB sqrt err vs trackproj",NBINS,0,8);
  mtobsqrterrhistoy=new TProfile("mTOBsqrterrvsprojy","matched y TOB sqrt err vs trackproj",NBINS,0,8);

  mtibtksqrterrhistox=new TProfile("mTIBtksqrterrvsprojx","matched x TIB tk sqrt err vs trackproj",NBINS,0,8);
  mtobtksqrterrhistox=new TProfile("mTOBtksqrterrvsprojx","matched x TOB tk sqrt err vs trackproj",NBINS,0,8);
  mtibtksqrterrhistoy=new TProfile("mTIBtksqrterrvsprojy","matched y TIB tk sqrt err vs trackproj",NBINS,0,8);
  mtobtksqrterrhistoy=new TProfile("mTOBtksqrterrvsprojy","matched y TOB tk sqrt err vs trackproj",NBINS,0,8);

  tidreshisto=new TProfile("TIDresvsproj","TID res vs trackproj",NBINS,0,8);
  tecreshisto=new TProfile("TECresvsproj","TEC res vs trackproj",NBINS,0,8);

  mtidreshistox=new TProfile("mTIDresvsprojx","matched x TID res vs trackproj",NBINS,0,8);
  mtecreshistox=new TProfile("mTECresvsprojx","matched x TEC res vs trackproj",NBINS,0,8);
  mtidreshistoy=new TProfile("mTIDresvsprojy","matched y TID res vs trackproj",NBINS,0,8);
  mtecreshistoy=new TProfile("mTECresvsprojy","matched y TEC res vs trackproj",NBINS,0,8);
  tiderrhisto=new TProfile("TIDerrvsproj","TID err vs trackproj",NBINS,0,8);
  tecerrhisto=new TProfile("TECerrvsproj","TEC err vs trackproj",NBINS,0,8);
  tidsqrterrhisto=new TProfile("TIDsqrterrvsproj","TID sqrt err vs trackproj",NBINS,0,8);
  tecsqrterrhisto=new TProfile("TECsqrterrvsproj","TEC sqrt err vs trackproj",NBINS,0,8);

  mtidsqrterrhistox=new TProfile("mTIDsqrterrvsprojx","matched x TID sqrt err vs trackproj",NBINS,0,8);
  mtecsqrterrhistox=new TProfile("mTECsqrterrvsprojx","matched x TEC sqrt err vs trackproj",NBINS,0,8);
  mtidsqrterrhistoy=new TProfile("mTIDsqrterrvsprojy","matched y TID sqrt err vs trackproj",NBINS,0,8);
  mtecsqrterrhistoy=new TProfile("mTECsqrterrvsprojy","matched y TEC sqrt err vs trackproj",NBINS,0,8);

  mtidtksqrterrhistox=new TProfile("mTIDtksqrterrvsprojx","matched x TID tk sqrt err vs trackproj",NBINS,0,8);
  mtectksqrterrhistox=new TProfile("mTECtksqrterrvsprojx","matched x TEC tk sqrt err vs trackproj",NBINS,0,8);
  mtidtksqrterrhistoy=new TProfile("mTIDtksqrterrvsprojy","matched y TID tk sqrt err vs trackproj",NBINS,0,8);
  mtectksqrterrhistoy=new TProfile("mTECtksqrterrvsprojy","matched y TEC tk sqrt err vs trackproj",NBINS,0,8);

  tobproj=new TH1F("TOBproj","TOB trackproj",12,0,12);
  tibproj=new TH1F("TIBproj","TIB trackproj",12,0,12);
  mtobproj=new TH1F("mTOBproj","matched TOB trackproj",12,0,12);
  mtibproj=new TH1F("mTIBproj","matched TIB trackproj",12,0,12);

  toberr=new TH1F("TOBerr","TOB error",20,0,0.5);
  tiberr=new TH1F("TIBerr","TIB error",20,0,0.5);

  for(int i=0;i<NBINS;i++){
    tobres[i]=new TH1F(Form("TOBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution",40,-1,1);
    tibres[i]=new TH1F(Form("TIBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution",40,-1,1);
    tidres[i]=new TH1F(Form("TIDres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution",40,-1,1);
    tecres[i]=new TH1F(Form("TECres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution",40,-1,1);
    mtobres[i]=new TH1F(Form("mTOBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched x TOB resolution",40,-1,1);
    mtibres[i]=new TH1F(Form("mTIBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched x TIB resolution",40,-1,1);
    mtidres[i]=new TH1F(Form("mTIDres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched x TID resolution",40,-1,1);
    mtecres[i]=new TH1F(Form("mTECres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched x TEC resolution",40,-1,1);
    mtobres_y[i]=new TH1F(Form("mTOBresy_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched y TOB resolution",40,-0.01,0.01);
    mtibres_y[i]=new TH1F(Form("mTIBresy_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched y TIB resolution",40,-0.01,0.01);
    mtidres_y[i]=new TH1F(Form("mTIDresy_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched y TID resolution",40,-0.01,0.01);
    mtecres_y[i]=new TH1F(Form("mTECresy_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"matched y TEC resolution",40,-0.01,0.01);
  }
  tibtkreshisto=new TProfile("TIBtkresvsproj","TIB tk res vs trackproj",12,0,8);
  tobtkreshisto=new TProfile("TOBtkresvsproj","TOB tk res vs trackproj",12,0,8);
  tibtkerrhisto=new TProfile("TIBtkerrvsproj","TIB tk err vs trackproj",12,0,8);
  tobtkerrhisto=new TProfile("TOBtkerrvsproj","TOB tk err vs trackproj",12,0,8);

  tobtkproj=new TH1F("TOBtkproj","TOB trackproj",12,0,12);
  tibtkproj=new TH1F("TIBtkproj","TIB trackproj",12,0,12);

  tobtkerr=new TH1F("TOBtkerr","TOB error",20,0,0.5);
  tibtkerr=new TH1F("TIBtkerr","TIB error",20,0,0.5);

  mtobresx=new TH1F("mTOBres_x","matched TOB resolution",40,-1,1);
  mtibresx=new TH1F("mTIBres_x","matched TIB resolution",40,-1,1);
  mtidresx=new TH1F("mTIDres_x","matched TID resolution",40,-1,1);
  mtecresx=new TH1F("mTECres_x","matched TEC resolution",40,-1,1);
  mtobresy=new TH1F("mTOBres_y","matched TOB resolution",40,-0.01,0.01);
  mtibresy=new TH1F("mTIBres_y","matched TIB resolution",40,-0.01,0.01);
  mtidresy=new TH1F("mTIDres_y","matched TID resolution",40,-0.01,0.01);
  mtecresy=new TH1F("mTECres_y","matched TEC resolution",40,-0.01,0.01);
  
  tobtkres=new TH1F("TOBtkres","TOB resolution",40,-1,1);
  tibtkres=new TH1F("TIBtkres","TIB tk resolution",40,-1,1);

  //    gFile=NULL;
}
void CPE::endJob() {  
  myfile_->Write();
  myfile_->Close();
  }
// Functions that gets called by framework every event
void CPE::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  using namespace edm;
  std::string rechitProducer = conf_.getParameter<std::string>("RecHitProducer");
  
  // Step A: Get Inputs 
  edm::Handle<SiStripRecHit2DCollection> rechitsrphi;
  edm::Handle<SiStripRecHit2DCollection> rechitsstereo;
  e.getByLabel(rechitProducer,"rphiRecHit", rechitsrphi);
  e.getByLabel(rechitProducer,"stereoRecHit", rechitsstereo);
  std::string cpe = conf_.getParameter<std::string>("StripCPE");
  edm::ESHandle<StripClusterParameterEstimator> parameterestimator;
  es.get<TkStripCPERecord>().get(cpe, parameterestimator); 

  edm::InputTag TkTag = conf_.getParameter<edm::InputTag>("TkRecHits");
  edm::Handle<TrackingRecHitCollection> trackrechitCollection;
  e.getByLabel(TkTag,trackrechitCollection);

  edm::Handle<std::vector<Trajectory> > TrajectoryCollection;
  e.getByLabel(TkTag,TrajectoryCollection);
  
  TrackerHitAssociator  associate(e,conf_);
  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord> ().get (pDD);
  const TrackerGeometry &tracker(*pDD);
  
  const std::vector<DetId> detIDs = rechitsrphi->ids();
  for ( std::vector<DetId>::const_iterator detunit_iterator = detIDs.begin(); detunit_iterator != detIDs.end(); detunit_iterator++ ) {//loop over detectors
    unsigned int id = (*detunit_iterator).rawId();
    edm::OwnVector<SiStripRecHit2D> collector; 
    if(id!=999999999){ //if is valid detector
      SiStripRecHit2DCollection::range rechitRange = rechitsrphi->get((*detunit_iterator));
      SiStripRecHit2DCollection::const_iterator rechitRangeIteratorBegin = rechitRange.first;
      SiStripRecHit2DCollection::const_iterator rechitRangeIteratorEnd   = rechitRange.second;
      SiStripRecHit2DCollection::const_iterator iter=rechitRangeIteratorBegin;
      for(iter=rechitRangeIteratorBegin;iter!=rechitRangeIteratorEnd;++iter){//loop on the rechit
	//	  SiStripRecHit2D rechit=*iter;
	const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > clust=iter->cluster();
	if (clust.isNonnull ()){
	  //edm::LogInfo("CPE")<<"The cluster is valid";
	  std::vector<PSimHit> matched=associate.associateHit(*iter);
	  //edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	  if(!matched.empty()){
	    float dist = 999999;
	    float mindist = 999999;
	    PSimHit closest;
	    for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	      dist = abs((iter)->localPosition().x() - (*m).localPosition().x());
	      if(dist<mindist){
		mindist = dist;
		closest = (*m);
	      }
	    }
	      
	      //	      edm::LogInfo("CPE")<<"Match performed";
	    LocalVector tkdir=closest.localDirection();
	    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(*detunit_iterator);
	    const StripTopology &topol=(StripTopology&)stripdet->topology();

	    float thickness=stripdet->specificSurface().bounds().thickness();
	    float pitch=topol.localPitch(iter->localPosition());
	    //  float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	    LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( closest.localPosition(),closest.localDirection(),0);
	    float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	    StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	    float resolution=topol.measurementPosition(parameters.first).x()- topol.measurementPosition(closest.localPosition()).x();
	    SiStripDetId detid=SiStripDetId(id);
	    //edm::LogInfo("CPE")<<detid;
	    if(fabs(resolution)<1){
	      if(detid.subdetId() == StripSubdetector::TIB ){
		tibreshisto->Fill(fabs(trackproj),resolution*resolution);
		tiberrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tibsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		tibproj->Fill(trackproj);
		tiberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tibres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TOB){
		tobreshisto->Fill(fabs(trackproj),resolution*resolution);
		toberrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tobsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		tobproj->Fill(trackproj);
		toberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tobres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TID){
		tidreshisto->Fill(fabs(trackproj),resolution*resolution);
		tiderrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tidsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		//tidproj->Fill(trackproj);
		//tiderr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tidres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TEC){
		tecreshisto->Fill(fabs(trackproj),resolution*resolution);
		tecerrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tecsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		//tecproj->Fill(trackproj);
		//tecerr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tecres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	  }
	}
	else{
	  edm::LogError("CPE")<<"The cluster is empty!";
	}
      }
    }
  }
  
  //   const std::vector<DetId> detIDs = rechitsstereo->ids();
  for ( std::vector<DetId>::const_iterator detunit_iterator = detIDs.begin(); detunit_iterator != detIDs.end(); detunit_iterator++ ) {//loop over detectors
    unsigned int id = (*detunit_iterator).rawId();
    edm::OwnVector<SiStripRecHit2D> collector; 
    if(id!=999999999){ //if is valid detector
      SiStripRecHit2DCollection::range rechitRange = rechitsstereo->get((*detunit_iterator));
      SiStripRecHit2DCollection::const_iterator rechitRangeIteratorBegin = rechitRange.first;
      SiStripRecHit2DCollection::const_iterator rechitRangeIteratorEnd   = rechitRange.second;
      SiStripRecHit2DCollection::const_iterator iter=rechitRangeIteratorBegin;
      for(iter=rechitRangeIteratorBegin;iter!=rechitRangeIteratorEnd;++iter){//loop on the rechit
	//	  SiStripRecHit2D rechit=*iter;
	const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > clust=iter->cluster();
	if (clust.isNonnull ()){
	  std::vector<PSimHit> matched=associate.associateHit(*iter);
	  edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	  if(!matched.empty()){
	    float dist = 999999;
	    float mindist = 999999;
	    PSimHit closest;
	    for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	      dist = abs((iter)->localPosition().x() - (*m).localPosition().x());
	      if(dist<mindist){
		mindist = dist;
		closest = (*m);
	      }
	    }
	    //	    float resolution=iter->localPosition().x()- matched[0].localPosition().x();
	    LocalVector tkdir=closest.localDirection();
	    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(*detunit_iterator);
	    const StripTopology &topol=(StripTopology&)stripdet->topology();
	    //	    float resolution=topol.measurementPosition(iter->localPosition()).x()- topol.measurementPosition(closest.localPosition()).x();
	    float thickness=stripdet->specificSurface().bounds().thickness();
	    float pitch=topol.localPitch(iter->localPosition());
	    //	    float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	    LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( closest.localPosition(),closest.localDirection(),0);
	    float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	    StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	    float resolution=topol.measurementPosition(parameters.first).x()- topol.measurementPosition(closest.localPosition()).x();
	    SiStripDetId detid=SiStripDetId(id);
	    if(fabs(resolution)<1){
	      if(detid.subdetId() == StripSubdetector::TIB ){
		tibreshisto->Fill(fabs(trackproj),resolution*resolution);
		tiberrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tibsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		tibproj->Fill(trackproj);
		tiberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tibres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TOB){
		tobreshisto->Fill(fabs(trackproj),resolution*resolution);
		toberrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tobsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		tobproj->Fill(trackproj);
		toberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tobres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TID){
		tidreshisto->Fill(fabs(trackproj),resolution*resolution);
		tiderrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tidsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		//tidproj->Fill(trackproj);
		//tiderr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tidres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	      else if(detid.subdetId() == StripSubdetector::TEC){
		tecreshisto->Fill(fabs(trackproj),resolution*resolution);
		tecerrhisto->Fill(fabs(trackproj),topol.measurementError(parameters.first,parameters.second).uu());
		tecsqrterrhisto->Fill(fabs(trackproj),sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		//tecproj->Fill(trackproj);
		//tecerr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
		if(trackproj<8)tecres[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	  }
	}
	else{
	  edm::LogError("ReadRecHit")<<"The cluster is empty!";
	}
      }
    }
  }
  TrackingRecHitCollection::const_iterator tkrechit;
  for(tkrechit=trackrechitCollection->begin();tkrechit!=trackrechitCollection->end();++tkrechit){
    //    for ( std::vector<DetId>::const_iterator detunit_iterator = detIDs.begin(); detunit_iterator != detIDs.end(); detunit_iterator++ ) {//loop over detectors
    unsigned int id = tkrechit->geographicalId().rawId();
    edm::OwnVector<SiStripRecHit2D> collector; 
    if(id!=999999999){ //if is valid detector
      //	SiStripRecHit2DCollection::range rechitRange = rechitsstereo->get((*detunit_iterator));
      //	SiStripRecHit2DCollection::const_iterator rechitRangeIteratorBegin = rechitRange.first;
      //	SiStripRecHit2DCollection::const_iterator rechitRangeIteratorEnd   = rechitRange.second;
      //	SiStripRecHit2DCollection::const_iterator iter=rechitRangeIteratorBegin;
      //	for(iter=rechitRangeIteratorBegin;iter!=rechitRangeIteratorEnd;++iter){//loop on the rechit
      //	  SiStripRecHit2D rechit=*iter;
      //      const SiStripRecHit2D *hit=dynamic_cast<const SiStripRecHit2D*>( &(*tkrechit));
      std::vector<SiStripRecHit2D*> rechits=getRecHitComponents(&(*tkrechit));
	for(std::vector<SiStripRecHit2D*>::iterator ihit=rechits.begin();ihit!=rechits.end();++ihit){
	  const SiStripRecHit2D* hit=*ihit;
	  const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > clust=hit->cluster();
	  if (clust.isNonnull ()){
	    std::vector<PSimHit> matched=associate.associateHit((*hit));
	    edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	    if(!matched.empty()){
	    float dist = 999999;
	    float mindist = 999999;
	    PSimHit closest;
	    for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	      dist = fabs((hit)->localPosition().x() - (*m).localPosition().x());
	      if(dist<mindist){
		mindist = dist;
		closest = (*m);
	      }
	    }
	    float resolution=hit->localPosition().x()- closest.localPosition().x();
	    LocalVector tkdir=closest.localDirection();
	    SiStripDetId detid=(SiStripDetId) hit->geographicalId();
	    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(detid);
	    if(stripdet==0){
	      const GeomDet* geomdet=(const GeomDet*)tracker.idToDet(detid);
	      stripdet=(StripGeomDetUnit*)geomdet->components()[0];
	    }
	    const StripTopology &topol=(StripTopology&)stripdet->topology();
	    float thickness=stripdet->specificSurface().bounds().thickness();
	    float pitch=topol.localPitch(hit->localPosition());
	    //float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	    LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( closest.localPosition(),closest.localDirection(),0);
	    float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	    StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	    if(fabs(resolution)/pitch<0.8){
	      if(detid.subdetId() == StripSubdetector::TIB ){
		tibtkreshisto->Fill(fabs(trackproj),resolution*resolution/(pitch*pitch));
		tibtkerrhisto->Fill(fabs(trackproj),parameters.second.xx()/(pitch*pitch));
		tibtkproj->Fill(trackproj);
		  tibtkerr->Fill(sqrt(parameters.second.xx())/pitch);
		  tibtkres->Fill(resolution/pitch);
	      }
	      else if(detid.subdetId() == StripSubdetector::TOB){
		tobtkreshisto->Fill(fabs(trackproj),resolution*resolution/(pitch*pitch));
		tobtkerrhisto->Fill(fabs(trackproj),parameters.second.xx()/(pitch*pitch));
		tobtkproj->Fill(trackproj);
		tobtkerr->Fill(sqrt(parameters.second.xx())/pitch);
		tobtkres->Fill(resolution/pitch);
	      }
	    }
	    //std::cout<<sqrt(parameters.second.xx())<<std::endl;	
	  }
	}
	else{
	  edm::LogError("ReadRecHit")<<"The cluster is empty!";
	}
      }
    }
  }
  
  std::vector<Trajectory>::const_iterator theTraj;
  for(theTraj = TrajectoryCollection->begin(); theTraj!= TrajectoryCollection->end();theTraj++){
    std::cout<<"Loop on traj"<<std::endl;
    std::vector<TrajectoryMeasurement> TMeas=theTraj->measurements();
    std::vector<TrajectoryMeasurement>::iterator itm;
    for (itm=TMeas.begin();itm!=TMeas.end();itm++){
      TrajectoryStateOnSurface tsos=itm->updatedState();
      const TransientTrackingRecHit::ConstRecHitPointer thit=itm->recHit();
      const SiStripMatchedRecHit2D* hit=dynamic_cast<const SiStripMatchedRecHit2D*>((*thit).hit());
      if(hit){//if matched hit...
	SiStripDetId detid=(SiStripDetId)hit->geographicalId();
	LocalVector trackdirection=tsos.localDirection();
	std::vector<PSimHit> matched=associate.associateHit((*hit));
	edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	if(!matched.empty()){
	  //project simhit;
	  const GluedGeomDet* gluedDet = (const GluedGeomDet*)tracker.idToDet(hit->geographicalId());
	  const StripGeomDetUnit* partnerstripdet =(StripGeomDetUnit*) gluedDet->stereoDet();
	  std::pair<LocalPoint,LocalVector> hitPair;
	  float dist = 999999;
	  float mindist = 999999;
	  float distx = 999999;
	  float disty = 999999;
	  std::pair<LocalPoint,LocalVector> closestPair;
	  //std::cout << " RECHIT position = " << position << std::endl;          
	  std::vector<PSimHit>::const_iterator mfound;
	  for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	    //project simhit;
	    hitPair= projectHit((*m),partnerstripdet,gluedDet->surface());
	    distx = fabs(hit->localPosition().x() - hitPair.first.x());
	    disty = fabs(hit->localPosition().y() - hitPair.first.y());
	    dist = sqrt(distx*distx+disty*disty);
	    // std::cout << " Simhit position x = " << hitPair.first.x() 
	     //      << " y = " << hitPair.first.y() << " dist = " << dist << std::endl;   
	    if(dist<mindist){
	      mindist = dist;
	      closestPair = hitPair;
	      mfound=m;
	    }
	  }
	  //std::cout << " Closest position x = " << closestPair.first.x() 
	  //      << " y = " << closestPair.first.y() << " dist = " << dist << std::endl;         
	  const StripGeomDetUnit * stripdet=(StripGeomDetUnit*) gluedDet->monoDet();
	  const StripTopology &topol=(StripTopology&)stripdet->topology();
	    //	    float resolution=topol.measurementPosition(iter->localPosition()).x()- topol.measurementPosition(closest.localPosition()).x();
	  //float thickness=stripdet->specificSurface().bounds().thickness();
	  //float pitch=topol.localPitch(iter->localPosition());
	    //	    float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	  LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( mfound->localPosition(),mfound->localDirection(),0);
	  float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	    //	    StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	    float resolution=topol.measurementPosition(hit->localPosition()).x()- topol.measurementPosition(closestPair.first).x();
	    float resolutiony=topol.measurementPosition(hit->localPosition()).y()- topol.measurementPosition(closestPair.first).y();

	    //	  float       resolution = hit->localPosition().x() - closestPair.first.x();
	    // float      resolutiony = hit->localPosition().y() - closestPair.first.y();
	  

	  //	      float resolution=hit->localPosition().x()- closest.localPosition().x();
	  //     float resolutiony=hit->localPosition().y()- closest.localPosition().y();
	    // LocalVector tkdir=closestPair.second;
	    // float trackproj=tkdir.x();
	      std::cout<<"Filling histograms:"<<std::endl;
	      if(detid.subdetId() == StripSubdetector::TIB ){
		mtibreshistox->Fill(fabs(trackproj),resolution);
		mtibsqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).uu()));
		mtibtksqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).uu()));
		mtibreshistoy->Fill(fabs(trackproj),resolutiony);
		mtibsqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).vv()));
		mtibtksqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).vv()));
		mtibproj->Fill(trackproj);
		if(trackproj<8){
		  mtibres[int(trackproj/(8./NBINS))]->Fill(resolution);
		  mtibres_y[int(trackproj/(8./NBINS))]->Fill(resolutiony);
		}
		mtibresx->Fill(resolution);
		mtibresy->Fill(resolutiony);
	      std::cout<<"Filled TIB"<<std::endl;
	      }
	      else if(detid.subdetId() == StripSubdetector::TOB){
		mtobreshistox->Fill(fabs(trackproj),resolution);
		mtobsqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).uu()));
		mtobtksqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).uu()));
		mtobreshistoy->Fill(fabs(trackproj),resolutiony);
		mtobsqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).vv()));
		mtobtksqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).vv()));
		mtobproj->Fill(trackproj);
		if(trackproj<8){
		  mtobres[int(trackproj/(8./NBINS))]->Fill(resolution);
		  mtobres_y[int(trackproj/(8./NBINS))]->Fill(resolutiony);
		}
		mtobresx->Fill(resolution);
		mtobresy->Fill(resolutiony);
	      std::cout<<"Filled TOB"<<std::endl;
	      }
	      else if(detid.subdetId() == StripSubdetector::TID){
		mtidreshistox->Fill(fabs(trackproj),resolution);
		mtidsqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).uu()));
		mtidtksqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).uu()));
		mtidreshistoy->Fill(fabs(trackproj),resolutiony);
		mtidsqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).vv()));
		mtidtksqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).vv()));
		//mtidproj->Fill(trackproj);
		if(trackproj<8){
		  mtidres[int(trackproj/(8./NBINS))]->Fill(resolution);
		  mtidres_y[int(trackproj/(8./NBINS))]->Fill(resolutiony);
		}
		mtidresx->Fill(resolution);
		mtidresy->Fill(resolutiony);
	      std::cout<<"Filled TID"<<std::endl;
	      }
	      else if(detid.subdetId() == StripSubdetector::TEC){
		mtecreshistox->Fill(fabs(trackproj),resolution);
		mtecsqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).uu()));
		mtectksqrterrhistox->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).uu()));
		mtecreshistoy->Fill(fabs(trackproj),resolutiony);
		mtecsqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),hit->localPositionError()).vv()));
		mtectksqrterrhistoy->Fill(fabs(trackproj),sqrt(topol.measurementError(hit->localPosition(),tsos.localError().positionError()).vv()));
		//mtecproj->Fill(trackproj);
		if(trackproj<8){
		  mtecres[int(trackproj/(8./NBINS))]->Fill(resolution);
		  mtecres_y[int(trackproj/(8./NBINS))]->Fill(resolutiony);
		}
		mtecresx->Fill(resolution);
		mtecresy->Fill(resolutiony);
	      std::cout<<"Filled TEC"<<std::endl;
	      }
	    }
	 }
     }
 }

}
  
std::vector<SiStripRecHit2D*> CPE::getRecHitComponents(const TrackingRecHit* rechit){
	std::vector<SiStripRecHit2D*> output;
	const ProjectedSiStripRecHit2D* phit=dynamic_cast<const ProjectedSiStripRecHit2D*>(rechit);
	const SiStripMatchedRecHit2D* matchedhit=dynamic_cast<const SiStripMatchedRecHit2D*>(rechit);
        const SiStripRecHit2D* hit=dynamic_cast<const SiStripRecHit2D*>(rechit);
	if(phit) hit=&(phit->originalHit());
        if(matchedhit){
                const SiStripRecHit2D* monohit   =  matchedhit->monoHit();
                const SiStripRecHit2D* stereohit =  matchedhit->stereoHit();
		output.push_back(monohit->clone());
		output.push_back(stereohit->clone());
        }
        else if (hit){
		output.push_back(hit->clone());
        }
	return output;
}
  
float CPE::uProj(const StripGeomDetUnit * stripdet,LocalTrajectoryParameters ltp,StripCPE *stripcpe){

  LocalPoint middlepoint = ltp.position();
  LocalVector atrackUnit = ltp.momentum()/ltp.momentum().mag();

  LocalError eresult;
  LocalVector drift=LocalVector(0,0,1);
  //  DetId detId(det.geographicalId());
  const StripTopology &topol=(StripTopology&)stripdet->topology();
  
  drift= stripcpe->driftDirection(stripdet);
  float thickness=stripdet->surface().bounds().thickness();

  LocalVector trackDir = atrackUnit;
      
  if(trackDir.z()*drift.z() > 0.) trackDir *= -1.;

  const Bounds& bounds = stripdet->surface().bounds();

  float maxLength = sqrt( bounds.length()*bounds.length()+bounds.width()*bounds.width());
  drift *= fabs(thickness/drift.z());       
  if(trackDir.z() !=0.) {
    trackDir *= fabs(thickness/trackDir.z());
  } else {
    trackDir *= maxLength/trackDir.mag();
  }

  LocalVector middleOfProjection = 0.5*(trackDir + drift);
  
  LocalPoint middlePointOnStrips = middlepoint + 0.5*drift;
  
  LocalPoint p1 = LocalPoint(middlePointOnStrips.x() + middleOfProjection.x()
			     ,middlePointOnStrips.y() + middleOfProjection.y());
  LocalPoint p2 = LocalPoint(middlePointOnStrips.x() - middleOfProjection.x()
			     ,middlePointOnStrips.y() - middleOfProjection.y());
  MeasurementPoint m1 = topol.measurementPosition(p1);
  MeasurementPoint m2 = topol.measurementPosition(p2);
  float u1 = m1.x();
  float u2 = m2.x();
  int nstrips = topol.nstrips(); 
  float UProj = std::min( float(fabs( u1 - u2)), float(nstrips));
  return UProj;

}

std::pair<LocalPoint,LocalVector> CPE::projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
                                                                    const BoundPlane& plane) 
 {
   //  const StripGeomDetUnit* stripDet = dynamic_cast<const StripGeomDetUnit*>(hit.det());
   //if (stripDet == 0) throw MeasurementDetException("HitMatcher hit is not on StripGeomDetUnit");
   
   const StripTopology& topol = stripDet->specificTopology();
   GlobalPoint globalpos= stripDet->surface().toGlobal(hit.localPosition());
   LocalPoint localHit = plane.toLocal(globalpos);
   //track direction
   LocalVector locdir=hit.localDirection();
   //rotate track in new frame
   
   GlobalVector globaldir= stripDet->surface().toGlobal(locdir);
   LocalVector dir=plane.toLocal(globaldir);
   float scale = -localHit.z() / dir.z();
   
   LocalPoint projectedPos = localHit + scale*dir;
   
   //  std::cout << "projectedPos " << projectedPos << std::endl;
   
   float selfAngle = topol.stripAngle( topol.strip( hit.localPosition()));
   
   LocalVector stripDir( sin(selfAngle), cos(selfAngle), 0); // vector along strip in hit frame
   
   LocalVector localStripDir( plane.toLocal(stripDet->surface().toGlobal( stripDir)));
   
   return std::pair<LocalPoint,LocalVector>( projectedPos, localStripDir);
 }
