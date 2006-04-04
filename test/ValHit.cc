// File: ValHit.cc
// Description:  see ValHit.h
// Author:  P. Azzi
// Creation Date:  PA Feb 2006 Initial version.
//
//--------------------------------------------
#include <memory>
#include <string>
#include <iostream>

#include "RecoLocalTracker/SiStripRecHitConverter/test/ValHit.h"

//--- for SimHit
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

//--- for RecHit
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DLocalPosCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DMatchedLocalPosCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"

//--- for StripDigiSimLink
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLinkCollection.h"

//--- framework stuff
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//--- for Geometry:
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "Geometry/Vector/interface/LocalPoint.h"
#include "Geometry/Vector/interface/GlobalPoint.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
                                                                                                                          


//for ntuple tracking
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

class TFile;

//---------------
// Constructor --
//---------------

using namespace std;
using namespace edm;

namespace cms
{
  
  // Functions that gets called by framework every event
  void ValHit::analyze(const edm::Event& e, const edm::EventSetup& es) {
    
    using namespace edm;
    
    //get event and run number
    
    myRun       = e.id().run();
    myEvent     = e.id().event();

    //--- get RecHits
    
    std::string rechitProducer = conf_.getParameter<std::string>("RecHitProducer");
    
    // Step A: Get Inputs 
    edm::Handle<SiStripRecHit2DMatchedLocalPosCollection> rechitsmatched;
    edm::Handle<SiStripRecHit2DLocalPosCollection> rechitsrphi;
    edm::Handle<SiStripRecHit2DLocalPosCollection> rechitsstereo;
    e.getByLabel(rechitProducer,"matchedRecHit", rechitsmatched);
    e.getByLabel(rechitProducer,"rphiRecHit", rechitsrphi);
    e.getByLabel(rechitProducer,"stereoRecHit", rechitsstereo);

    //--- get SimHits

   // Step A: Get Inputs
    theStripHits.clear();
    edm::Handle<edm::PSimHitContainer> TIBHitsLowTof;
    edm::Handle<edm::PSimHitContainer> TIBHitsHighTof;
    edm::Handle<edm::PSimHitContainer> TIDHitsLowTof;
    edm::Handle<edm::PSimHitContainer> TIDHitsHighTof;
    edm::Handle<edm::PSimHitContainer> TOBHitsLowTof;
    edm::Handle<edm::PSimHitContainer> TOBHitsHighTof;
    edm::Handle<edm::PSimHitContainer> TECHitsLowTof;
    edm::Handle<edm::PSimHitContainer> TECHitsHighTof;

    e.getByLabel("SimG4Object","TrackerHitsTIBLowTof", TIBHitsLowTof);
    e.getByLabel("SimG4Object","TrackerHitsTIBHighTof", TIBHitsHighTof);
    e.getByLabel("SimG4Object","TrackerHitsTIDLowTof", TIDHitsLowTof);
    e.getByLabel("SimG4Object","TrackerHitsTIDHighTof", TIDHitsHighTof);
    e.getByLabel("SimG4Object","TrackerHitsTOBLowTof", TOBHitsLowTof);
    e.getByLabel("SimG4Object","TrackerHitsTOBHighTof", TOBHitsHighTof);
    e.getByLabel("SimG4Object","TrackerHitsTECLowTof", TECHitsLowTof);
    e.getByLabel("SimG4Object","TrackerHitsTECHighTof", TECHitsHighTof);
    
    theStripHits.insert(theStripHits.end(), TIBHitsLowTof->begin(), TIBHitsLowTof->end()); 
    theStripHits.insert(theStripHits.end(), TIBHitsHighTof->begin(), TIBHitsHighTof->end());
    theStripHits.insert(theStripHits.end(), TIDHitsLowTof->begin(), TIDHitsLowTof->end()); 
    theStripHits.insert(theStripHits.end(), TIDHitsHighTof->begin(), TIDHitsHighTof->end());
    theStripHits.insert(theStripHits.end(), TOBHitsLowTof->begin(), TOBHitsLowTof->end()); 
    theStripHits.insert(theStripHits.end(), TOBHitsHighTof->begin(), TOBHitsHighTof->end());
    theStripHits.insert(theStripHits.end(), TECHitsLowTof->begin(), TECHitsLowTof->end()); 
    theStripHits.insert(theStripHits.end(), TECHitsHighTof->begin(), TECHitsHighTof->end());


    int numrechitrphi =0;
    int numrechitsas =0;
    int numsimhit =0;
      
    SimHitMap.clear();
    for (std::vector<PSimHit>::iterator isim = theStripHits.begin();
	 isim != theStripHits.end(); ++isim){
      SimHitMap[(*isim).detUnitId()].push_back((*isim));
    }
    
    //first instance tracking geometry
    edm::ESHandle<TrackerGeometry> pDD;
    es.get<TrackerDigiGeometryRecord> ().get (pDD);
    
    // loop over detunits
    for(TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++){
      uint32_t myid=((*it)->geographicalId()).rawId();       
      DetId detid = ((*it)->geographicalId());

      numrechitrphi =0;
      numrechitsas =0;
      numsimhit =0;
      Tiblayer =-99;
      Tibstereo=-99; 
      Tibnumsimhit=0;
      Tibnumrechitrphi=0;
      Tibnumrechitsas=0;
      Toblayer=-99;
      Tobstereo=-99; 
      Tobnumsimhit=0;
      Tobnumrechitrphi=0;
      Tobnumrechitsas=0;
      Tidlayer =-99;
      Tidstereo=-99; 
      Tidnumsimhit=0;
      Tidnumrechitrphi=0;
      Tidnumrechitsas=0;
      Teclayer=-99;
      Tecstereo=-99; 
      Tecnumsimhit=0;
      Tecnumrechitrphi=0;
      Tecnumrechitsas=0;
      
      // initialize here
      for(int i=0; i<MAXHIT; i++){
	simhitx[i] =0.;
	simhity[i] =0.;
	simhitz[i] =0.;
	simhitphi[i]=0;
	simhiteta[i]=0;
	rechitrphix[i] =0;
	rechitrphiy[i] =0;
	rechitrphiz[i] =0;
	rechitrphiphi[i] =0;
	rechitsasx[i] =0;
	rechitsasy[i] =0;
	rechitsasz[i] =0;
	rechitsasphi[i] =0;
	clusizrphi[i] =0;
	clusizsas[i] =0;
	cluchgrphi[i] =0;
	cluchgsas[i] =0;
	//tree vars	
	Tibsimhitx[i]=0;
	Tibsimhity[i]=0;
	Tibsimhitz[i]=0;
	Tibsimhitphi[i]=0;
	Tibsimhiteta[i]=0;
	Tibrphix[i]=0;
	Tibrphiy[i]=0;
	Tibrphiz[i]=0;
	Tibrphiphi[i]=0;
	Tibrphires[i]=0;
	Tibrphisiz[i]=0;
	Tibrphichg[i]=0;
	Tibsasx[i]=0;
	Tibsasy[i]=0;
	Tibsasz[i]=0;
	Tibsasphi[i]=0;
	Tibsasres[i]=0;
	Tibsassiz[i]=0;
	Tibsaschg[i]=0;
	
	Tobsimhitx[i]=0;
	Tobsimhity[i]=0;
	Tobsimhitz[i]=0;
	Tobsimhitphi[i]=0;
	Tobsimhiteta[i]=0;
	Tobrphix[i]=0;
	Tobrphiy[i]=0;
	Tobrphiz[i]=0;
	Tobrphiphi[i]=0;
	Tobrphires[i]=0;
	Tobrphisiz[i]=0;
	Tobrphichg[i]=0;
	Tobsasx[i]=0;
	Tobsasy[i]=0;
	Tobsasz[i]=0;
	Tobsasphi[i]=0;
	Tobsasres[i]=0;
	Tobsassiz[i]=0;
	Tobsaschg[i]=0;	

	Tidsimhitx[i]=0;
	Tidsimhity[i]=0;
	Tidsimhitz[i]=0;
	Tidsimhitphi[i]=0;
	Tidsimhiteta[i]=0;
	Tidrphix[i]=0;
	Tidrphiy[i]=0;
	Tidrphiz[i]=0;
	Tidrphiphi[i]=0;
	Tidrphires[i]=0;
	Tidrphisiz[i]=0;
	Tidrphichg[i]=0;
	Tidsasx[i]=0;
	Tidsasy[i]=0;
	Tidsasz[i]=0;
	Tidsasphi[i]=0;
	Tidsasres[i]=0;
	Tidsassiz[i]=0;
	Tidsaschg[i]=0;	

	Tecsimhitx[i]=0;
	Tecsimhity[i]=0;
	Tecsimhitz[i]=0;
	Tecsimhitphi[i]=0;
	Tecsimhiteta[i]=0;
	Tecrphix[i]=0;
	Tecrphiy[i]=0;
	Tecrphiz[i]=0;
	Tecrphiphi[i]=0;
	Tecrphires[i]=0;
	Tecrphisiz[i]=0;
	Tecrphichg[i]=0;
	Tecsasx[i]=0;
	Tecsasy[i]=0;
	Tecsasz[i]=0;
	Tecsasphi[i]=0;
	Tecsasres[i]=0;
	Tecsassiz[i]=0;
	Tecsaschg[i]=0;	
      }


      //---get simhit
      std::map<unsigned int, std::vector<PSimHit> >::const_iterator it = SimHitMap.find(myid);
      vector<PSimHit> simHit; 
      simHit.clear();
      if (it!= SimHitMap.end()){
	simHit = it->second;
	vector<PSimHit>::const_iterator simHitIter = simHit.begin();
	vector<PSimHit>::const_iterator simHitIterEnd = simHit.end();
        numsimhit = simHit.size();
	int i=0;
	for (;simHitIter != simHitIterEnd; ++simHitIter) {
	  const PSimHit ihit = *simHitIter;
	  LocalPoint simhitpos = ihit.localPosition();
	  simhitx[i] = simhitpos.x();
	  simhity[i] = simhitpos.y();
	  simhitz[i] = simhitpos.z();
	  simhitphi[i] = simhitpos.phi();
	  simhiteta[i] = simhitpos.eta();
	  i++;
	}
      }
      numrechitrphi =0;
      //loop over rechits-rphi in the same subdetector
      SiStripRecHit2DLocalPosCollection::range          rechitrphiRange = rechitsrphi->get(detid);
      SiStripRecHit2DLocalPosCollection::const_iterator rechitrphiRangeIteratorBegin = rechitrphiRange.first;
      SiStripRecHit2DLocalPosCollection::const_iterator rechitrphiRangeIteratorEnd   = rechitrphiRange.second;
      SiStripRecHit2DLocalPosCollection::const_iterator iterrphi=rechitrphiRangeIteratorBegin;
      
      numrechitrphi = rechitrphiRangeIteratorEnd - rechitrphiRangeIteratorBegin;   
      //cout << "TEST numrechitrphi = " << numrechitrphi << endl;
      if(numrechitrphi > 0 ){
	int i=0;
	for(iterrphi=rechitrphiRangeIteratorBegin; iterrphi!=rechitrphiRangeIteratorEnd;++iterrphi){
	  SiStripRecHit2DLocalPos const rechit=*iterrphi;
	  LocalPoint position=rechit.localPosition();
	  //LocalError error=rechit.localPositionError();
	  const std::vector<const SiStripCluster*> clust=rechit.cluster();
	  int clusiz=0;
	  int totcharge=0;
	  for(vector<const SiStripCluster*>::const_iterator ic = clust.begin(); ic!=clust.end(); ic++) {
	    clusiz = (*ic)->amplitudes().size();
	    const std::vector<short> amplitudes=(*ic)->amplitudes();
	    for(size_t i=0; i<amplitudes.size();i++){
	      totcharge+=amplitudes[i];
	    }
	  }
	  rechitrphix[i] = position.x();
	  rechitrphiy[i] = position.y();
	  rechitrphiz[i] = position.z();
	  rechitrphiphi[i] = position.phi();
	  clusizrphi[i] = clusiz;
	  cluchgrphi[i] = totcharge;
	  i++;
	}
      }
      
      //loop over rechits-sas in the same subdetector
      numrechitsas=0;
      SiStripRecHit2DLocalPosCollection::range rechitsasRange = rechitsstereo->get(detid);
      SiStripRecHit2DLocalPosCollection::const_iterator rechitsasRangeIteratorBegin = rechitsasRange.first;
      SiStripRecHit2DLocalPosCollection::const_iterator rechitsasRangeIteratorEnd   = rechitsasRange.second;
      SiStripRecHit2DLocalPosCollection::const_iterator itersas=rechitsasRangeIteratorBegin;
      numrechitsas = rechitsasRangeIteratorEnd - rechitsasRangeIteratorBegin;   
      if(numrechitsas > 0){
	int j=0;
	for(itersas=rechitsasRangeIteratorBegin; itersas!=rechitsasRangeIteratorEnd;++itersas){
	  SiStripRecHit2DLocalPos const rechit=*itersas;
	  LocalPoint position=rechit.localPosition();
	  //	LocalError error=rechit.localPositionError();
	  std::vector<const SiStripCluster*> clust=rechit.cluster();
	  int clusiz=0;
	  int totcharge=0;
	  for(vector<const SiStripCluster*>::const_iterator ic = clust.begin(); ic!=clust.end(); ic++) {
	    clusiz = (*ic)->amplitudes().size();
	    const std::vector<short> amplitudes=(*ic)->amplitudes();
	    for(size_t i=0; i<amplitudes.size();i++){
	      totcharge+=amplitudes[i];
	    }
	  }
	  rechitsasx[j] = position.x();
	  rechitsasy[j] = position.y();
	  rechitsasz[j] = position.z();
	  rechitsasphi[j] = position.phi();
	  clusizsas[j] = clusiz;
	  cluchgsas[j] = totcharge;
	  j++;
	}
      }
      
      //--- fill histograms
      
      // 	for(int k=0; k<numsimhit; k++){
      // 	  std::cout << " simhit (x,y,z) = "  <<  simhitx[k] << ", "  <<  simhity[k] << ", " <<  simhitz[k] << std::endl;
      // 	}
      // 	for(int k=0; k<numrechitrphi; k++){
      // 	  std::cout << " rechitphi (x,y,z) = "  <<  rechitrphix[k] << ", " <<  rechitrphiy[k] << ", " << rechitrphiz[k] << std::endl;
      // 	}
      // 	for(int k=0; k<numrechitsas; k++){
      // 	  std::cout << " rechitsas (x,y,z) = "  <<  rechitsasx[k] << ", " <<  rechitsasy[k] << ", " << rechitsasz[k] << std::endl;
      // 	}
      
      
      //std::cout << " det id= " << myid << " N(simhit) = " << numsimhit 
      //<< " N(rechitrphi) = " << numrechitrphi << " N(rechitsas)= " << numrechitsas << std::endl;      
      
      if(numsimhit>0 || numrechitrphi>0 || numrechitsas>0 ){
	if( detid.subdetId() == int(StripSubdetector::TIB)){
	  TIBDetId tibid(myid); 
	  tiblayer   = tibid.layer();
	  // 	  tibfw_bw   = tibid.string()[0];
	  // 	  tibext_int = tibid.string()[1];
	  // 	  tibstring  = tibid.string()[2];
	  // 	  tibmodule  = tibid.module();
	  tibstereo  = tibid.stereo();
	  
	  Tiblayer = tiblayer;
	  Tibstereo = tibstereo;
	  Tibnumsimhit = numsimhit;
	  Tibnumrechitrphi = numrechitrphi;
	  Tibnumrechitsas = numrechitsas;
	  
	  if(tibstereo == 0) {
	    for(int k=0; k<numsimhit; k++){
	      Tibsimhitx[k] = simhitx[k];
	      Tibsimhity[k] = simhity[k];
	      Tibsimhitz[k] = simhitz[k];
	      Tibsimhitphi[k] = simhitphi[k];
	      Tibsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitrphi; kk++){
	      Tibrphix[kk] = rechitrphix[kk];
	      Tibrphiy[kk] = rechitrphiy[kk];
	      Tibrphiz[kk] = rechitrphiz[kk];
	      Tibrphisiz[kk] = clusizrphi[kk];
	      Tibrphichg[kk] = cluchgrphi[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitrphi>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitrphi; kk++){
	    // 		Tibrphires[kk] = rechitrphix[kk]-simhitx[k];
	    // 		if(tiblayer == 1){
	    // 		  tibclu1rphi->Fill(clusizrphi[kk]);
	    // 		  tibres1rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(tiblayer == 2){
	    // 		  tibclu2rphi->Fill(clusizrphi[kk]);
	    // 		  tibres2rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(tiblayer == 3){
	    // 		  tibclu3rphi->Fill(clusizrphi[kk]);
	    // 		  tibres3rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		  }
	    // 		if(tiblayer == 4){
	    // 		  tibclu4rphi->Fill(clusizrphi[kk]);
	    // 		  tibres4rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }
	    // 	  }
	  } else if (tibstereo == 1) {
	    for(int k=0; k<numsimhit; k++){
	      Tibsimhitx[k] = simhitx[k];
	      Tibsimhity[k] = simhity[k];
	      Tibsimhitz[k] = simhitz[k];
	      Tibsimhitphi[k] = simhitphi[k];
	      Tibsimhiteta[k] = simhiteta[k];
	    }	    
	    for(int kk=0;kk<numrechitsas; kk++){
	      Tibsasx[kk] = rechitsasx[kk];
	      Tibsasy[kk] = rechitsasy[kk];
	      Tibsasz[kk] = rechitsasz[kk];
	      Tibsassiz[kk] = clusizsas[kk];
	      Tibsaschg[kk] = cluchgsas[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitsas>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitsas; kk++){
	    // 		Tibsasres[kk] = rechitsasx[kk]-simhitx[k];
	    // 		if(tiblayer == 1){
	    // 		    tibclu1sas->Fill(clusizsas[kk]);
	    // 		    tibres1sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		  }
	    // 		if(tiblayer == 2){
	    // 		  tibclu2sas->Fill(clusizsas[kk]);
	    // 		  tibres2sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }   
	    // 	  }
	  }
	}
	if( detid.subdetId() == int( StripSubdetector::TOB)){
	  TOBDetId tobid(myid); 
	  toblayer   = tobid.layer();
	  tobstereo  = tobid.stereo();
	  Toblayer = toblayer;
	  Tobstereo = tobstereo;
	  Tobnumsimhit = numsimhit;
	  Tobnumrechitrphi = numrechitrphi;
	  Tobnumrechitsas = numrechitsas;
	  if(tobstereo == 0) {
	    for(int k=0; k<numsimhit; k++){
	      Tobsimhitx[k] = simhitx[k];
	      Tobsimhity[k] = simhity[k];
	      Tobsimhitz[k] = simhitz[k];
	      Tobsimhitphi[k] = simhitphi[k];
	      Tobsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitrphi; kk++){
	      Tobrphix[kk] = rechitrphix[kk];
	      Tobrphiy[kk] = rechitrphiy[kk];
	      Tobrphiz[kk] = rechitrphiz[kk];
	      Tobrphisiz[kk] = clusizrphi[kk];
	      Tobrphichg[kk] = cluchgrphi[kk];
	    }	
	    // 	  if(numsimhit>0 && numrechitrphi>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitrphi; kk++){
	    // 		Tobrphires[kk] = rechitrphix[kk]-simhitx[k];
	    // 		if(toblayer == 1){
	    // 		  tobclu1rphi->Fill(clusizrphi[kk]);
	    // 		  tobres1rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 2){
	    // 		  tobclu2rphi->Fill(clusizrphi[kk]);
	    // 		  tobres2rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 3){
	    // 		  tobclu3rphi->Fill(clusizrphi[kk]);
	    // 		  tobres3rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 4){
	    // 		  tobclu4rphi->Fill(clusizrphi[kk]);
	    // 		  tobres4rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 5){
	    // 		  tobclu5rphi->Fill(clusizrphi[kk]);
	    // 		  tobres5rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 6){
	    // 		  tobclu6rphi->Fill(clusizrphi[kk]);
	    // 		  tobres6rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }
	    // 	  }
	  } else if (tobstereo == 1) {
	    for(int k=0; k<numsimhit; k++){
	      Tobsimhitx[k] = simhitx[k];
	      Tobsimhity[k] = simhity[k];
	      Tobsimhitz[k] = simhitz[k];
	      Tobsimhitphi[k] = simhitphi[k];
	      Tobsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitsas; kk++){
	      Tobsasx[kk] = rechitsasx[kk];
	      Tobsasy[kk] = rechitsasy[kk];
	      Tobsasz[kk] = rechitsasz[kk];
	      Tobsassiz[kk] = clusizsas[kk];
	      Tobsaschg[kk] = cluchgsas[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitsas>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitsas; kk++){
	    // 		Tobsasres[kk] = rechitsasx[kk]-simhitx[k];
	    // 		if(toblayer == 1){
	    // 		  tobclu1sas->Fill(clusizsas[kk]);
	    // 		  tobres1sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 		if(toblayer == 2){
	    // 		  tobclu2sas->Fill(clusizsas[kk]);
	    // 		  tobres2sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }   
	  }
	}
	
	if( detid.subdetId() == int (StripSubdetector::TID)){
	  
	  TIDDetId tidid(myid); 
	  tidlayer   = tidid.wheel();
	  tidstereo  = tidid.stereo();
	  Tidlayer   = tidlayer;
	  Tidring    = tidid.ring();
	  Tidstereo = tidstereo;
	  Tidnumsimhit = numsimhit;
	  Tidnumrechitrphi = numrechitrphi;
	  Tidnumrechitsas = numrechitsas;
	  if(tidstereo == 0) {
	    for(int k=0; k<numsimhit; k++){
	      Tidsimhitx[k] = simhitx[k];
	      Tidsimhity[k] = simhity[k];
	      Tidsimhitz[k] = simhitz[k];
	      Tidsimhitphi[k] = simhitphi[k];
	      Tidsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitrphi; kk++){
	      Tidrphix[kk] = rechitrphix[kk];
	      Tidrphiy[kk] = rechitrphiy[kk];
	      Tidrphiz[kk] = rechitrphiz[kk];
	      Tidrphisiz[kk] = clusizrphi[kk];
	      Tidrphichg[kk] = cluchgrphi[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitrphi>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitrphi; kk++){
	    // 		Tidrphires[kk] = rechitrphix[kk]-simhitx[k];
	    // 		if(tidlayer == 1){
	    // 		  Tidrphires[kk] = rechitrphix[kk]-simhitx[k];
	    // 		  tidclu1rphi->Fill(clusizrphi[kk]);
	    // 		  tidres1rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(tidlayer == 2){
	    // 		  tidclu2rphi->Fill(clusizrphi[kk]);
	    // 		  tidres2rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(tidlayer == 3){
	    // 		  tidclu3rphi->Fill(clusizrphi[kk]);
	    // 		  tidres3rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }
	    // 	  }
	  } else if (tidstereo == 1) {
	    for(int k=0; k<numsimhit; k++){
	      Tidsimhitx[k] = simhitx[k];
	      Tidsimhity[k] = simhity[k];
	      Tidsimhitz[k] = simhitz[k];
	      Tidsimhitphi[k] = simhitphi[k];
	      Tidsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitsas; kk++){
	      Tidsasx[kk] = rechitsasx[kk];
	      Tidsasy[kk] = rechitsasy[kk];
	      Tidsasz[kk] = rechitsasz[kk];
	      Tidsassiz[kk] = clusizsas[kk];
	      Tidsaschg[kk] = cluchgsas[kk];
	    }
	    // if(numsimhit>0 && numrechitsas>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitsas; kk++){
	    // 		Tidsasres[kk] = rechitsasx[kk]-simhitx[k];
	    // 		if(tidlayer == 1){
	    // 		  tidclu1sas->Fill(clusizsas[kk]);
	    // 		  tidres1sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 		if(tidlayer == 2){
	    // 		  tidclu2sas->Fill(clusizsas[kk]);
	    // 		  tidres2sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }   
	    // 	  }
	  }
	}
	if( detid.subdetId() == int(StripSubdetector::TEC)){
	  
	  TECDetId tecid(myid); 
	  teclayer   = tecid.wheel();
	  tecstereo  = tecid.stereo();
	  Teclayer = teclayer;
	  Tecstereo = tecstereo;
	  Tecring    = tecid.ring();
	  Tecnumsimhit = numsimhit;
	  Tecnumrechitrphi = numrechitrphi;
	  Tecnumrechitsas = numrechitsas;
	  if(tecstereo == 0) {
	    for(int k=0; k<numsimhit; k++){
	      Tecsimhitx[k] = simhitx[k];
	      Tecsimhity[k] = simhity[k];
	      Tecsimhitz[k] = simhitz[k];
	      Tecsimhitphi[k] = simhitphi[k];
	      Tecsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitrphi; kk++){
	      Tecrphix[kk] = rechitrphix[kk];
	      Tecrphiy[kk] = rechitrphiy[kk];
	      Tecrphiz[kk] = rechitrphiz[kk];
	      Tecrphisiz[kk] = clusizrphi[kk];
	      Tecrphichg[kk] = cluchgrphi[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitrphi>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitrphi; kk++){
	    // 		Tecrphires[kk] = rechitrphix[kk]-simhitx[k];
	    // 		if(teclayer == 1){
	    // 		  tecclu1rphi->Fill(clusizrphi[kk]);
	    // 		  tecres1rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 2){
	    // 		  tecclu2rphi->Fill(clusizrphi[kk]);
	    // 		  tecres2rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 3){
	    // 		  tecclu3rphi->Fill(clusizrphi[kk]);
	    // 		  tecres3rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 4){
	    // 		  tecclu4rphi->Fill(clusizrphi[kk]);
	    // 		  tecres4rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 5){
	    // 		  tecclu5rphi->Fill(clusizrphi[kk]);
	    // 		  tecres5rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 6){
	    // 		  tecclu6rphi->Fill(clusizrphi[kk]);
	    // 		  tecres6rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 7){
	    // 		  tecclu7rphi->Fill(clusizrphi[kk]);
	    // 		  tecres7rphi->Fill(rechitrphix[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }
	    // 	  }
	  } else if (tecstereo == 1) {
	    for(int k=0; k<numsimhit; k++){
	      Tecsimhitx[k] = simhitx[k];
	      Tecsimhity[k] = simhity[k];
	      Tecsimhitz[k] = simhitz[k];
	      Tecsimhitphi[k] = simhitphi[k];
	      Tecsimhiteta[k] = simhiteta[k];
	    }
	    for(int kk=0;kk<numrechitsas; kk++){
	      Tecsasx[kk] = rechitsasx[kk];
	      Tecsasy[kk] = rechitsasy[kk];
	      Tecsasz[kk] = rechitsasz[kk];
	      Tecsassiz[kk] = clusizsas[kk];
	      Tecsaschg[kk] = cluchgsas[kk];
	    }
	    // 	  if(numsimhit>0 && numrechitsas>0){
	    // 	    for(int k=0; k<numsimhit; k++){
	    // 	      for(int kk=0;kk<numrechitsas; kk++){
	    // 		Tecsasres[kk] = rechitsasx[kk]-simhitx[k];
	    // 		if(teclayer == 1){
	    // 		  tecclu1sas->Fill(clusizsas[kk]);
	    // 		  tecres1sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 2){
	    // 		  tecclu2sas->Fill(clusizsas[kk]);
	    // 		  tecres2sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 		if(teclayer == 5){
	    // 		  tecclu5sas->Fill(clusizsas[kk]);
	    // 		  tecres5sas->Fill(rechitsasx[kk]-simhitx[k]);
	    // 		}
	    // 	      }
	    // 	    }   
	    // 	  }
	  }
	}
	myTree->Fill();     
      }
    }
    cout << " === calling end job " << endl;  
  }
  
  
  ValHit::ValHit(edm::ParameterSet const& conf) : 
    conf_(conf),filename_(conf.getParameter<std::string>("fileName")) 
    //histograms here
  {
    //TIB
    tibres1rphi = new TH1F("tibres1rphi", "res TIB lay 1 rphi", 100, -0.03, +0.03);   
    tibres2rphi = new TH1F("tibres2rphi", "res TIB lay 2 rphi", 100, -0.03, +0.03); 
    tibres3rphi = new TH1F("tibres3rphi", "res TIB lay 3 rphi", 100, -0.03, +0.03); 
    tibres4rphi = new TH1F("tibres4rphi", "res TIB lay 4 rphi", 100, -0.03, +0.03); 
    tibres1sas  = new TH1F("tibres1sas", "res TIB lay 1 sas", 100, -0.03, +0.03);
    tibres2sas  = new TH1F("tibres2sas", "res TIB lay 2 sas", 100, -0.03, +0.03); 
    tibclu1rphi = new TH1F("tibclu1rphi", "cluster size TIB lay 1 rphi", 5, -0.5, 4.5);   
    tibclu2rphi = new TH1F("tibclu2rphi", "cluster size TIB lay 2 rphi", 5, -0.5, 4.5);   
    tibclu3rphi = new TH1F("tibclu3rphi", "cluster size TIB lay 3 rphi", 5, -0.5, 4.5);   
    tibclu4rphi = new TH1F("tibclu4rphi", "cluster size TIB lay 4 rphi", 5, -0.5, 4.5);   
    tibclu1sas  = new TH1F("tibclu1sas", "cluster size TIB lay 1 sas", 5, -0.5, 4.5);   
    tibclu2sas  = new TH1F("tibclu2sas", "cluster size TIB lay 2 sas", 5, -0.5, 4.5);   
    tibhits     = new TH1F("tibhits", "# rechits TIB",6, -0.5, 6.5);

    tib_rphi_digi_chg  = new TH1F("tib_rphi_digi_chg", "charge digi",300,0, 300);
    tib_rphi_tot_nst  = new TH2F("tib_rphi_tot_nst", "charge vs nstp", 5, 0, 5, 300, 0, 300);
    tib_rphi_pos_nst  = new TH2F("tib_rphi_pos_nst", "bary-1st vs nstp",5, 0, 5, 100,0, 4);
 
    //TOB
    tobres1rphi = new TH1F("tobres1rphi", "res TOB lay 1 rphi", 100, -0.03, +0.03);   
    tobres2rphi = new TH1F("tobres2rphi", "res TOB lay 2 rphi", 100, -0.03, +0.03); 
    tobres3rphi = new TH1F("tobres3rphi", "res TOB lay 3 rphi", 100, -0.03, +0.03); 
    tobres4rphi = new TH1F("tobres4rphi", "res TOB lay 4 rphi", 100, -0.03, +0.03); 
    tobres5rphi = new TH1F("tobres5rphi", "res TOB lay 5 rphi", 100, -0.03, +0.03); 
    tobres6rphi = new TH1F("tobres6rphi", "res TOB lay 6 rphi", 100, -0.03, +0.03); 
    tobres1sas  = new TH1F("tobres1sas", "res TOB lay 1 sas", 100, -0.03, +0.03);
    tobres2sas  = new TH1F("tobres2sas", "res TOB lay 2 sas", 100, -0.03, +0.03); 
    tobclu1rphi = new TH1F("tobclu1rphi", "cluster size TOB lay 1 rphi", 5, -0.5, 4.5);   
    tobclu2rphi = new TH1F("tobclu2rphi", "cluster size TOB lay 2 rphi", 5, -0.5, 4.5);   
    tobclu3rphi = new TH1F("tobclu3rphi", "cluster size TOB lay 3 rphi", 5, -0.5, 4.5);   
    tobclu4rphi = new TH1F("tobclu4rphi", "cluster size TOB lay 4 rphi", 5, -0.5, 4.5);   
    tobclu5rphi = new TH1F("tobclu5rphi", "cluster size TOB lay 5 rphi", 5, -0.5, 4.5);   
    tobclu6rphi = new TH1F("tobclu6rphi", "cluster size TOB lay 6 rphi", 5, -0.5, 4.5);   
    tobclu1sas  = new TH1F("tobclu1sas", "cluster size TOB lay 1 sas", 5, -0.5, 4.5);   
    tobclu2sas  = new TH1F("tobclu2sas", "cluster size TOB lay 2 sas", 5, -0.5, 4.5);   

    //TID
    tidres1rphi = new TH1F("tidres1rphi", "res TID lay 1 rphi", 100, -1.0, +1.0);   
    tidres2rphi = new TH1F("tidres2rphi", "res TID lay 2 rphi", 100, -1.0, +1.0);   
    tidres3rphi = new TH1F("tidres3rphi", "res TID lay 3 rphi", 100, -1.0, +1.0);   
    tidres1sas = new TH1F("tidres1sas", "res TID lay 1 sas", 100, -1.0, +1.0);   
    tidres2sas = new TH1F("tidres2sas", "res TID lay 2 sas", 100, -1.0, +1.0);   
    tidclu1rphi = new TH1F("tidclu1rphi", "cluster size TID lay 1 rphi", 5, -0.5, 4.5);   
    tidclu2rphi = new TH1F("tidclu2rphi", "cluster size TID lay 1 rphi", 5, -0.5, 4.5);   
    tidclu3rphi = new TH1F("tidclu3rphi", "cluster size TID lay 1 rphi", 5, -0.5, 4.5);   
    tidclu1sas = new TH1F("tidclu1sas", "cluster size TID lay 1 sas", 5, -0.5, 4.5);   
    tidclu2sas = new TH1F("tidclu2sas", "cluster size TID lay 1 sas", 5, -0.5, 4.5);   

    //TEC
    tecres1rphi = new TH1F("tecres1rphi", "res TEC lay 1 rphi", 100, -1.0, +1.0);   
    tecres2rphi = new TH1F("tecres2rphi", "res TEC lay 2 rphi", 100, -1.0, +1.0);   
    tecres3rphi = new TH1F("tecres3rphi", "res TEC lay 3 rphi", 100, -1.0, +1.0);   
    tecres4rphi = new TH1F("tecres4rphi", "res TEC lay 4 rphi", 100, -1.0, +1.0);   
    tecres5rphi = new TH1F("tecres5rphi", "res TEC lay 5 rphi", 100, -1.0, +1.0);   
    tecres6rphi = new TH1F("tecres6rphi", "res TEC lay 6 rphi", 100, -1.0, +1.0);   
    tecres7rphi = new TH1F("tecres7rphi", "res TEC lay 7 rphi", 100, -1.0, +1.0);   
    tecres1sas = new TH1F("tecres1sas", "res TEC lay 1 sas", 100, -1.0, +1.0);   
    tecres2sas = new TH1F("tecres2sas", "res TEC lay 2 sas", 100, -1.0, +1.0);   
    tecres5sas = new TH1F("tecres5sas", "res TEC lay 5 sas", 100, -1.0, +1.0);   
    tecclu1rphi = new TH1F("tecclu1rphi", "cluster size TEC lay 1 rphi", 5, -0.5, 4.5);   
    tecclu2rphi = new TH1F("tecclu2rphi", "cluster size TEC lay 2 rphi", 5, -0.5, 4.5);   
    tecclu3rphi = new TH1F("tecclu3rphi", "cluster size TEC lay 3 rphi", 5, -0.5, 4.5);   
    tecclu4rphi = new TH1F("tecclu4rphi", "cluster size TEC lay 4 rphi", 5, -0.5, 4.5);   
    tecclu5rphi = new TH1F("tecclu5rphi", "cluster size TEC lay 5 rphi", 5, -0.5, 4.5);   
    tecclu6rphi = new TH1F("tecclu6rphi", "cluster size TEC lay 6 rphi", 5, -0.5, 4.5);   
    tecclu7rphi = new TH1F("tecclu7rphi", "cluster size TEC lay 7 rphi", 5, -0.5, 4.5);   
    tecclu1sas = new TH1F("tecclu1sas", "cluster size TEC lay 1 sas", 5, -0.5, 4.5);   
    tecclu2sas = new TH1F("tecclu2sas", "cluster size TEC lay 2 sas", 5, -0.5, 4.5);   
    tecclu5sas = new TH1F("tecclu5sas", "cluster size TEC lay 5 sas", 5, -0.5, 4.5);   


    //======================================================
    myFile = new TFile(filename_.c_str(),"RECREATE");
    myTree = new TTree("HitTree","Tracker Validation tree");
    // GENERAL block
    myTree->Branch("Run", &myRun, "Run/I");
    myTree->Branch("Event", &myEvent, "Event/I");

    //TIB info
    myTree->Branch("Tiblayer", &Tiblayer, "Tiblayer/I");
    myTree->Branch("Tibstereo", &Tibstereo, "Tibstereo/I");
 
    Tibnumsimhit=0;
    myTree->Branch("Tibnumsimhit", &Tibnumsimhit, "Tibnumsimhit/I");
    myTree->Branch("Tibsimx", &Tibsimhitx, "Tibsimx[Tibnumsimhit]/F");
    myTree->Branch("Tibsimy", &Tibsimhity, "Tibsimy[Tibnumsimhit]/F");
    myTree->Branch("Tibsimz", &Tibsimhitz, "Tibsimz[Tibnumsimhit]/F");
    myTree->Branch("Tibsimphi", &Tibsimhitphi, "Tibsimphi[Tibnumsimhit]/F");
    myTree->Branch("Tibsimeta", &Tibsimhiteta, "Tibsimeta[Tibnumsimhit]/F");
    Tibnumrechitrphi=0;
    myTree->Branch("Tibnumrechitrphi", &Tibnumrechitrphi, "Tibnumrechitrphi/I");
    myTree->Branch("Tibrphix", &Tibrphix, "Tibrphix[Tibnumrechitrphi]/F"); 
    myTree->Branch("Tibrphiy", &Tibrphiy, "Tibrphiy[Tibnumrechitrphi]/F");
    myTree->Branch("Tibrphiz", &Tibrphiz, "Tibrphiz[Tibnumrechitrphi]/F");
    myTree->Branch("Tibrphires", &Tibrphires, "Tibrphires[Tibnumrechitrphi]/I");
    myTree->Branch("Tibrphisiz", &Tibrphisiz, "Tibrphisiz[Tibnumrechitrphi]/I");
    myTree->Branch("Tibrphichg", &Tibrphichg, "Tibrphichg[Tibnumrechitrphi]/F");
    Tibnumrechitsas=0;
    myTree->Branch("Tibnumrechitsas", &Tibnumrechitsas, "Tibnumrechitsas/I");
    myTree->Branch("Tibsasx", &Tibsasx, "Tibsasx[Tibnumrechitsas]/F"); 
    myTree->Branch("Tibsasy", &Tibsasy, "Tibsasy[Tibnumrechitsas]/F");
    myTree->Branch("Tibsasz", &Tibsasz, "Tibsasz[Tibnumrechitsas]/F");        
    myTree->Branch("Tibsasres", &Tibsasres, "Tibsasres[Tibnumrechitsas]/I");
    myTree->Branch("Tibsassiz", &Tibsassiz, "Tibsassiz[Tibnumrechitsas]/I");
    myTree->Branch("Tibsaschg", &Tibsaschg, "Tibsaschg[Tibnumrechitsas]/F");

    //TOB info
    myTree->Branch("Toblayer", &Toblayer, "Toblayer/I");
    myTree->Branch("Tobstereo", &Tobstereo, "Tobstereo/I");
    Tobnumsimhit=0;
    myTree->Branch("Tobnumsimhit", &Tobnumsimhit, "Tobnumsimhit/I");
    myTree->Branch("Tobsimx", &Tobsimhitx, "Tobsimx[Tobnumsimhit]/F");
    myTree->Branch("Tobsimy", &Tobsimhity, "Tobsimy[Tobnumsimhit]/F");
    myTree->Branch("Tobsimz", &Tobsimhitz, "Tobsimz[Tobnumsimhit]/F");
    myTree->Branch("Tobsimphi", &Tobsimhitphi, "Tobsimphi[Tobnumsimhit]/F");
    myTree->Branch("Tobsimeta", &Tobsimhiteta, "Tobsimeta[Tobnumsimhit]/F");
    Tobnumrechitrphi=0;
    myTree->Branch("Tobnumrechitrphi", &Tobnumrechitrphi, "Tobnumrechitrphi/I");
    myTree->Branch("Tobrphix", &Tobrphix, "Tobrphix[Tobnumrechitrphi]/F"); 
    myTree->Branch("Tobrphiy", &Tobrphiy, "Tobrphiy[Tobnumrechitrphi]/F");
    myTree->Branch("Tobrphiz", &Tobrphiz, "Tobrphiz[Tobnumrechitrphi]/F");
    myTree->Branch("Tobrphires", &Tobrphires, "Tobrphires[Tobnumrechitrphi]/I");
    myTree->Branch("Tobrphisiz", &Tobrphisiz, "Tobrphisiz[Tobnumrechitrphi]/I");
    myTree->Branch("Tobrphichg", &Tobrphichg, "Tobrphichg[Tobnumrechitrphi]/F");
    Tobnumrechitsas=0;
    myTree->Branch("Tobnumrechitsas", &Tobnumrechitsas, "Tobnumrechitsas/I");
    myTree->Branch("Tobsasx", &Tobsasx, "Tobsasx[Tobnumrechitsas]/F"); 
    myTree->Branch("Tobsasy", &Tobsasy, "Tobsasy[Tobnumrechitsas]/F");
    myTree->Branch("Tobsasz", &Tobsasz, "Tobsasz[Tobnumrechitsas]/F");        
    myTree->Branch("Tobsasres", &Tobsasres, "Tobsasres[Tobnumrechitsas]/I");
    myTree->Branch("Tobsassiz", &Tobsassiz, "Tobsassiz[Tobnumrechitsas]/I");
    myTree->Branch("Tobsaschg", &Tobsaschg, "Tobsaschg[Tobnumrechitsas]/F");

    //TID info
    myTree->Branch("Tidlayer", &Tidlayer, "Tidlayer/I");
    myTree->Branch("Tidstereo", &Tidstereo, "Tidstereo/I");
    myTree->Branch("Tidring", &Tidring, "Tidring/I");
    Tidnumsimhit=0;
    myTree->Branch("Tidnumsimhit", &Tidnumsimhit, "Tidnumsimhit/I");
    myTree->Branch("Tidsimx", &Tidsimhitx, "Tidsimx[Tidnumsimhit]/F");
    myTree->Branch("Tidsimy", &Tidsimhity, "Tidsimy[Tidnumsimhit]/F");
    myTree->Branch("Tidsimz", &Tidsimhitz, "Tidsimz[Tidnumsimhit]/F");
    myTree->Branch("Tidsimphi", &Tidsimhitphi, "Tidsimphi[Tidnumsimhit]/F");
    myTree->Branch("Tidsimeta", &Tidsimhiteta, "Tidsimeta[Tidnumsimhit]/F");
    Tidnumrechitrphi=0;
    myTree->Branch("Tidnumrechitrphi", &Tidnumrechitrphi, "Tidnumrechitrphi/I");
    myTree->Branch("Tidrphix", &Tidrphix, "Tidrphix[Tidnumrechitrphi]/F"); 
    myTree->Branch("Tidrphiy", &Tidrphiy, "Tidrphiy[Tidnumrechitrphi]/F");
    myTree->Branch("Tidrphiz", &Tidrphiz, "Tidrphiz[Tidnumrechitrphi]/F");
    myTree->Branch("Tidrphires", &Tidrphires, "Tidrphires[Tidnumrechitrphi]/I");
    myTree->Branch("Tidrphisiz", &Tidrphisiz, "Tidrphisiz[Tidnumrechitrphi]/I");
    myTree->Branch("Tidrphichg", &Tidrphichg, "Tidrphichg[Tidnumrechitrphi]/F");
    Tidnumrechitsas=0;
    myTree->Branch("Tidnumrechitsas", &Tidnumrechitsas, "Tidnumrechitsas/I");
    myTree->Branch("Tidsasx", &Tidsasx, "Tidsasx[Tidnumrechitsas]/F"); 
    myTree->Branch("Tidsasy", &Tidsasy, "Tidsasy[Tidnumrechitsas]/F");
    myTree->Branch("Tidsasz", &Tidsasz, "Tidsasz[Tidnumrechitsas]/F");        
    myTree->Branch("Tidsasres", &Tidsasres, "Tidsasres[Tidnumrechitsas]/I");
    myTree->Branch("Tidsassiz", &Tidsassiz, "Tidsassiz[Tidnumrechitsas]/I");
    myTree->Branch("Tidsaschg", &Tidsaschg, "Tidsaschg[Tidnumrechitsas]/F");

    //TEC info
    myTree->Branch("Teclayer", &Teclayer, "Teclayer/I");
    myTree->Branch("Tecstereo", &Tecstereo, "Tecstereo/I");
    myTree->Branch("Tecring", &Tecring, "Tecring/I");
    Tecnumsimhit=0;
    myTree->Branch("Tecnumsimhit", &Tecnumsimhit, "Tecnumsimhit/I");
    myTree->Branch("Tecsimx", &Tecsimhitx, "Tecsimx[Tecnumsimhit]/F");
    myTree->Branch("Tecsimy", &Tecsimhity, "Tecsimy[Tecnumsimhit]/F");
    myTree->Branch("Tecsimz", &Tecsimhitz, "Tecsimz[Tecnumsimhit]/F");
    myTree->Branch("Tecsimphi", &Tecsimhitphi, "Tecsimphi[Tecnumsimhit]/F");
    myTree->Branch("Tecsimeta", &Tecsimhiteta, "Tecsimeta[Tecnumsimhit]/F");
    Tecnumrechitrphi=0;
    myTree->Branch("Tecnumrechitrphi", &Tecnumrechitrphi, "Tecnumrechitrphi/I");
    myTree->Branch("Tecrphix", &Tecrphix, "Tecrphix[Tecnumrechitrphi]/F"); 
    myTree->Branch("Tecrphiy", &Tecrphiy, "Tecrphiy[Tecnumrechitrphi]/F");
    myTree->Branch("Tecrphiz", &Tecrphiz, "Tecrphiz[Tecnumrechitrphi]/F");
    myTree->Branch("Tecrphires", &Tecrphires, "Tecrphires[Tecnumrechitrphi]/I");
    myTree->Branch("Tecrphisiz", &Tecrphisiz, "Tecrphisiz[Tecnumrechitrphi]/I");
    myTree->Branch("Tecrphichg", &Tecrphichg, "Tecrphichg[Tecnumrechitrphi]/F");
    Tecnumrechitsas=0;
    myTree->Branch("Tecnumrechitsas", &Tecnumrechitsas, "Tecnumrechitsas/I");
    myTree->Branch("Tecsasx", &Tecsasx, "Tecsasx[Tecnumrechitsas]/F"); 
    myTree->Branch("Tecsasy", &Tecsasy, "Tecsasy[Tecnumrechitsas]/F");
    myTree->Branch("Tecsasz", &Tecsasz, "Tecsasz[Tecnumrechitsas]/F");        
    myTree->Branch("Tecsasres", &Tecsasres, "Tecsasres[Tecnumrechitsas]/I");
    myTree->Branch("Tecsassiz", &Tecsassiz, "Tecsassiz[Tecnumrechitsas]/I");
    myTree->Branch("Tecsaschg", &Tecsaschg, "Tecsaschg[Tecnumrechitsas]/F");

  }
  
  
  // Virtual destructor needed.
  ValHit::~ValHit(){ 
    delete myFile;
  }
  


  void ValHit::endJob() {  
    cout << ">>> ending histograms" << endl;
    myFile->cd();
    myTree->Write();
    tibres1rphi->Write(); 
    tibres2rphi->Write(); 
    tibres3rphi->Write(); 
    tibres4rphi->Write();
    tibres1sas->Write();  
    tibres2sas->Write();
    tibclu1rphi->Write(); 
    tibclu2rphi->Write(); 
    tibclu3rphi->Write(); 
    tibclu4rphi->Write();
    tibclu1sas->Write();
    tibclu2sas->Write();
    tib_rphi_digi_chg->Write();
    tib_rphi_tot_nst->Write();
    tib_rphi_pos_nst->Write();
    //TOB
    tobres1rphi->Write(); 
    tobres2rphi->Write(); 
    tobres3rphi->Write(); 
    tobres4rphi->Write();
    tobres5rphi->Write();
    tobres6rphi->Write();
    tobres1sas->Write();  
    tobres2sas->Write();
    tobclu1rphi->Write(); 
    tobclu2rphi->Write(); 
    tobclu3rphi->Write(); 
    tobclu4rphi->Write();
    tobclu5rphi->Write();
    tobclu6rphi->Write();
    tobclu1sas->Write();
    tobclu2sas->Write();
    //TID
    tidres1rphi->Write(); 
    tidres2rphi->Write(); 
    tidres3rphi->Write(); 
    tidres1sas->Write();  
    tidres2sas->Write();
    tidclu1rphi->Write(); 
    tidclu2rphi->Write(); 
    tidclu3rphi->Write(); 
    tidclu1sas->Write();
    tidclu2sas->Write();
    //TEC
    tecres1rphi->Write(); 
    tecres2rphi->Write(); 
    tecres3rphi->Write(); 
    tecres4rphi->Write();
    tecres5rphi->Write();
    tecres6rphi->Write();
    tecres7rphi->Write();
    tecres1sas->Write();  
    tecres2sas->Write();
    tecres5sas->Write();
    tecclu1rphi->Write(); 
    tecclu2rphi->Write(); 
    tecclu3rphi->Write(); 
    tecclu4rphi->Write();
    tecclu5rphi->Write();
    tecclu6rphi->Write();
    tecclu7rphi->Write();
    tecclu1sas->Write();
    tecclu2sas->Write();
    tecclu5sas->Write();
    //------------------------
    tibhits->Write();
    myFile->Close();
    cout << ">>> File closed " << endl;
  }


}




