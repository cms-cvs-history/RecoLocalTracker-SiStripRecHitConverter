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
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/SiStripRecHitMatcher.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

CPE::CPE(edm::ParameterSet const& conf) : 
  conf_(conf)
{
}


// Virtual destructor needed.
CPE::~CPE() { }  

void CPE::beginJob(const edm::EventSetup& c) {  
  //  myfile_= new TFile("resolution_QCD.root","RECREATE");
  edm::Service<TFileService> fs;
  tobproj=fs->make<TH1F>("TOBproj","TOB trackproj",12,0,12);
  tibproj=fs->make<TH1F>("TIBproj","TIB trackproj",12,0,12);


  toberr=fs->make<TH1F>("TOBerr","TOB error",20,0,0.5);
  tiberr=fs->make<TH1F>("TIBerr","TIB error",20,0,0.5);


  for(int i=0;i<NBINS;i++){

    tobres[i]=fs->make<TH1F>(Form("TOBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution",1000,-10,10);
    tibres[i]=fs->make<TH1F>(Form("TIBres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution",1000,-10,10);
    tidres[i]=fs->make<TH1F>(Form("TIDres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution",1000,-10,10);
    tecres[i]=fs->make<TH1F>(Form("TECres_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution",1000,-10,10);

    tobres1[i]=fs->make<TH1F>(Form("TOBres1_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 1 strip",1000,-10,10);
    tibres1[i]=fs->make<TH1F>(Form("TIBres1_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 1 strip",1000,-10,10);
    tidres1[i]=fs->make<TH1F>(Form("TIDres1_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 1 strip",1000,-10,10);
    tecres1[i]=fs->make<TH1F>(Form("TECres1_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 1 strip",1000,-10,10);

    tobres2[i]=fs->make<TH1F>(Form("TOBres2_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 2 strip",1000,-10,10);
    tibres2[i]=fs->make<TH1F>(Form("TIBres2_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 2 strip",1000,-10,10);
    tidres2[i]=fs->make<TH1F>(Form("TIDres2_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 2 strip",1000,-10,10);
    tecres2[i]=fs->make<TH1F>(Form("TECres2_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 2 strip",1000,-10,10);

    tobres3[i]=fs->make<TH1F>(Form("TOBres3_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 3 strip",1000,-10,10);
    tibres3[i]=fs->make<TH1F>(Form("TIBres3_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 3 strip",1000,-10,10);
    tidres3[i]=fs->make<TH1F>(Form("TIDres3_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 3 strip",1000,-10,10);
    tecres3[i]=fs->make<TH1F>(Form("TECres3_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 3 strip",1000,-10,10);

    tobres4[i]=fs->make<TH1F>(Form("TOBres4_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 4 strip",1000,-10,10);
    tibres4[i]=fs->make<TH1F>(Form("TIBres4_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 4 strip",1000,-10,10);
    tidres4[i]=fs->make<TH1F>(Form("TIDres4_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 4 strip",1000,-10,10);
    tecres4[i]=fs->make<TH1F>(Form("TECres4_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 4 strip",1000,-10,10);

    tobres5[i]=fs->make<TH1F>(Form("TOBres5_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 5 strip",1000,-10,10);
    tibres5[i]=fs->make<TH1F>(Form("TIBres5_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 5 strip",1000,-10,10);
    tidres5[i]=fs->make<TH1F>(Form("TIDres5_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 5 strip",1000,-10,10);
    tecres5[i]=fs->make<TH1F>(Form("TECres5_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 5 strip",1000,-10,10);

    tobres6[i]=fs->make<TH1F>(Form("TOBres6_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TOB resolution 6 strip",1000,-10,10);
    tibres6[i]=fs->make<TH1F>(Form("TIBres6_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TIB resolution 6 strip",1000,-10,10);
    tidres6[i]=fs->make<TH1F>(Form("TIDres6_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TID resolution 6 strip",1000,-10,10);
    tecres6[i]=fs->make<TH1F>(Form("TECres6_%f-%f",i*8./NBINS,i*8./NBINS+8./NBINS),"TEC resolution 6 strip",1000,-10,10);
 
    tobreswc[i]=fs->make<TH1F>(Form("TOBreswc_%d_strip",i),"TOB resolution wide cluster",2000,-20,20);
    tibreswc[i]=fs->make<TH1F>(Form("TIBreswc_%d_strip",i),"TIB resolution wide cluster",2000,-20,20);
    tidreswc[i]=fs->make<TH1F>(Form("TIDreswc_%d_strip",i),"TID resolution wide cluster",2000,-20,20);
    tecreswc[i]=fs->make<TH1F>(Form("TECreswc_%d_strip",i),"TEC resolution wide cluster",2000,-20,20);

  }

  tobtkproj=fs->make<TH1F>("TOBtkproj","TOB trackproj",12,0,12);
  tibtkproj=fs->make<TH1F>("TIBtkproj","TIB trackproj",12,0,12);

  tobtkerr=fs->make<TH1F>("TOBtkerr","TOB error",20,0,0.5);
  tibtkerr=fs->make<TH1F>("TIBtkerr","TIB error",20,0,0.5);



}
void CPE::endJob() {  
  //  myfile_->Write();
  //myfile_->Close();
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

  
  TrackerHitAssociator  associate(e,conf_);
  edm::ESHandle<TrackerGeometry> pDD;
  es.get<TrackerDigiGeometryRecord> ().get (pDD);
  const TrackerGeometry &tracker(*pDD);

  edm::Handle< edm::DetSetVector<StripDigiSimLink> >  stripdigisimlink;
  e.getByLabel("siStripDigis", stripdigisimlink);

  

  for (SiStripRecHit2DCollection::const_iterator detunit_iterator = rechitsrphi->begin(); detunit_iterator != rechitsrphi->end(); ++detunit_iterator) {
    SiStripRecHit2DCollection::DetSet rechits = *detunit_iterator; 
    unsigned int id = rechits.detId();
    edm::OwnVector<SiStripRecHit2D> collector; 
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rechits.begin();
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rechits.end();
    SiStripRecHit2DCollection::DetSet::const_iterator iter=rechitRangeIteratorBegin;
    //SiStripRecHit2DCollection::DetSet::const_iterator iter2=rechitRangeIteratorBegin;
    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(rechits.detId());
    const StripTopology &topol=(StripTopology&)stripdet->topology();
    //    std::cout<<"DetId= "<<id<<std::endl;
    for(iter=rechitRangeIteratorBegin;iter!=rechitRangeIteratorEnd;++iter){//loop on the rechit
      //	  SiStripRecHit2D rechit=*iter;
      SiStripRecHit2D::ClusterRef clust=iter->cluster();
      //    const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > clust=iter->cluster();
      if (clust.isNonnull ()){ // cluster
 	//edm::LogInfo("CPE")<<"The cluster is valid";
	std::vector<PSimHit> matched=associate.associateHit(*iter);
	//edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	  float dist = 999999;
	  float mindist = 999999;
	  PSimHit closest;
	  bool simfound=false;
	  //	  std::cout<<std::endl;
	  for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	    //	    if(m->pabs()>0.9){
	    dist = fabs(topol.measurementPosition((iter)->localPosition()).x() - topol.measurementPosition((*m).localPosition()).x());
	      if(dist<mindist){
		mindist = dist;
		closest = (*m);
		simfound=true;
	      }
	    }
	  //	  }
	  if(simfound){
	  
	  //	      edm::LogInfo("CPE")<<"Match performed";
	  LocalVector tkdir=closest.localDirection();
	  
	  //float thickness=stripdet->specificSurface().bounds().thickness();
	  //float pitch=topol.localPitch(iter->localPosition());
	  //  float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	  LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( closest.localPosition(),closest.localDirection(),0);
	  float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	  StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	  float resolution=topol.measurementPosition(parameters.first).x()- topol.measurementPosition(closest.localPosition()).x();
	  SiStripDetId detid=SiStripDetId(id);
	  //edm::LogInfo("CPE")<<detid;
	  int size=clust->amplitudes().size();
	  //	  if(closest.pabs()>0.9){
	    if(detid.subdetId() == StripSubdetector::TIB ){
	      tibproj->Fill(trackproj);
	      tiberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tibreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tibres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tibres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tibres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tibres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tibres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tibres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tibres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	  }
	    else if(detid.subdetId() == StripSubdetector::TOB){
	      tobproj->Fill(trackproj);
	      toberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tobreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tobres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tobres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tobres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tobres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tobres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tobres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tobres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	    else if(detid.subdetId() == StripSubdetector::TID){
	      //tidproj->Fill(trackproj);
	      //tiderr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tidreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tidres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tidres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tidres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tidres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tidres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tidres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tidres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	    else if(detid.subdetId() == StripSubdetector::TEC){
	      //tecproj->Fill(trackproj);
	      //tecerr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tecreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tecres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tecres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tecres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tecres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tecres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tecres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tecres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	  }
	  //	}
      }
      else{
	edm::LogError("CPE")<<"The cluster is empty!";
      }
    }
  }
  
  for (SiStripRecHit2DCollection::const_iterator detunit_iterator = rechitsstereo->begin(); detunit_iterator != rechitsstereo->end(); ++detunit_iterator) {
    SiStripRecHit2DCollection::DetSet rechits = *detunit_iterator; 
    unsigned int id = rechits.detId();
    edm::OwnVector<SiStripRecHit2D> collector; 
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rechits.begin();
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rechits.end();
    SiStripRecHit2DCollection::DetSet::const_iterator iter=rechitRangeIteratorBegin;
    const StripGeomDetUnit * stripdet=(const StripGeomDetUnit*)tracker.idToDetUnit(rechits.detId());
    const StripTopology &topol=(StripTopology&)stripdet->topology();
    for(iter=rechitRangeIteratorBegin;iter!=rechitRangeIteratorEnd;++iter){//loop on the rechit
      //	  SiStripRecHit2D rechit=*iter;
      SiStripRecHit2D::ClusterRef clust=iter->cluster();
      //    const edm::Ref<edm::DetSetVector<SiStripCluster>, SiStripCluster, edm::refhelper::FindForDetSetVector<SiStripCluster> > clust=iter->cluster();
      if (clust.isNonnull ()){
	std::vector<PSimHit> matched=associate.associateHit(*iter);
	//	edm::LogInfo("CPE")<<"matched size= "<<matched.size();
	  float dist = 999999;
	  float mindist = 999999;
	  PSimHit closest;
	  bool simfound=false;
	  for(std::vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	    //	    if(m->pabs()>0.9){
	    dist = fabs(topol.measurementPosition((iter)->localPosition()).x() - topol.measurementPosition((*m).localPosition()).x());
	      if(dist<mindist){
		mindist = dist;
		closest = (*m);
		simfound=true;
	      }
	      //	    }
	  }
	  if(simfound){
	  //	    float resolution=iter->localPosition().x()- matched[0].localPosition().x();
	  LocalVector tkdir=closest.localDirection();
	  //	    float resolution=topol.measurementPosition(iter->localPosition()).x()- topol.measurementPosition(closest.localPosition()).x();
	  //float thickness=stripdet->specificSurface().bounds().thickness();
	  //float pitch=topol.localPitch(iter->localPosition());
	  //	    float trackproj=tkdir.x()/tkdir.z()*thickness/pitch;
	  LocalTrajectoryParameters tkparam= LocalTrajectoryParameters( closest.localPosition(),closest.localDirection(),0);
	  float trackproj=uProj(stripdet,tkparam,(StripCPE *)&(*parameterestimator));
	  StripClusterParameterEstimator::LocalValues parameters=parameterestimator->localParameters(*clust,*stripdet,tkparam);
	  float resolution=topol.measurementPosition(parameters.first).x()- topol.measurementPosition(closest.localPosition()).x();
	  SiStripDetId detid=SiStripDetId(id);
	  int size=clust->amplitudes().size();
	  //	  if(closest.pabs()>0.9){
	    if(detid.subdetId() == StripSubdetector::TIB ){
	      tibproj->Fill(trackproj);
	      tiberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tibreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tibres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tibres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tibres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tibres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tibres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tibres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tibres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	    else if(detid.subdetId() == StripSubdetector::TOB){
	      tobproj->Fill(trackproj);
	      toberr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tobreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tobres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tobres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tobres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tobres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tobres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tobres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tobres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	    else if(detid.subdetId() == StripSubdetector::TID){
	      //tidproj->Fill(trackproj);
	      //tiderr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tidreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tidres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tidres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tidres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tidres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tidres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tidres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tidres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	    else if(detid.subdetId() == StripSubdetector::TEC){
	      //tecproj->Fill(trackproj);
		//tecerr->Fill(sqrt(topol.measurementError(parameters.first,parameters.second).uu()));
	      if(((size - (trackproj+2.5)) > 1)&&size<NBINS)tecreswc[size]->Fill(resolution);
	      else if(!((size - (trackproj+2.5))>1)&&trackproj<8){
		tecres[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==1)tecres1[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==2)tecres2[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==3)tecres3[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==4)tecres4[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==5)tecres5[int(trackproj/(8./NBINS))]->Fill(resolution);
		if(size==6)tecres6[int(trackproj/(8./NBINS))]->Fill(resolution);
	      }
	    }
	  }
	  //	}
      }
      else{
	edm::LogError("ReadRecHit")<<"The cluster is empty!";
      }
    }
  }

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
   //   LocalVector locdir=hit.localDirection();
   LocalVector locdir=hit.entryPoint()-hit.exitPoint();
   //rotate track in new frame
   
   GlobalVector globaldir= stripDet->surface().toGlobal(locdir);
   LocalVector dir=plane.toLocal(globaldir);
   float scale = -localHit.z() / dir.z();
   
   LocalPoint projectedPos = localHit + scale*dir;
   
   //  std::cout << "projectedPos " << projectedPos << std::endl;
   
   float selfAngle = topol.stripAngle( topol.strip( hit.localPosition()));
   
   LocalVector stripDir( sin(selfAngle), cos(selfAngle), 0); // vector along strip in hit frame
   
   LocalVector localStripDir( plane.toLocal(stripDet->surface().toGlobal( stripDir)));
   
   return std::pair<LocalPoint,LocalVector>( projectedPos, dir);
 }
