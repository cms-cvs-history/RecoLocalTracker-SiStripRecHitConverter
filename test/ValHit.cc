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

#include "Geometry/Vector/interface/LocalPoint.h"
#include "Geometry/Vector/interface/GlobalPoint.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/GluedGeomDet.h"
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

    std::cout << " =========== Event =  " << myEvent << " ================== " << std::endl;

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
    int numrechitmatched=0;
      
    SimHitMap.clear();
    for (std::vector<PSimHit>::iterator isim = theStripHits.begin();
	 isim != theStripHits.end(); ++isim){
      SimHitMap[(*isim).detUnitId()].push_back((*isim));
    }
    
    //first instance tracking geometry
    edm::ESHandle<TrackerGeometry> pDD;
    es.get<TrackerDigiGeometryRecord> ().get (pDD);
    const TrackerGeometry &tracker(*pDD);

    // loop over detunits
    for(TrackerGeometry::DetContainer::const_iterator it = pDD->dets().begin(); it != pDD->dets().end(); it++){
      uint32_t myid=((*it)->geographicalId()).rawId();       
      DetId detid = ((*it)->geographicalId());

      numrechitrphi =0;
      numrechitsas =0;
      numrechitmatched=0;
      numsimhit =0;
      Tiblayer =-99;
      Tibstereo=-99; 
      Tibnumsimhit=0;
      Tibnumrechitrphi=0;
      Tibnumrechitsas=0;
      Tibnumrechitmatched=0;
      Toblayer=-99;
      Tobstereo=-99; 
      Tobnumsimhit=0;
      Tobnumrechitrphi=0;
      Tobnumrechitsas=0;
      Tobnumrechitmatched=0;
      Tidlayer =-99;
      Tidstereo=-99; 
      Tidnumsimhit=0;
      Tidnumrechitrphi=0;
      Tidnumrechitsas=0;
      Tidnumrechitmatched=0;
      Teclayer=-99;
      Tecstereo=-99; 
      Tecnumsimhit=0;
      Tecnumrechitrphi=0;
      Tecnumrechitsas=0;
      Tecnumrechitmatched=0;
      
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
	rechitrphires[i]=-999.;
	rechitsasres[i]=-999.;
	rechitrphimatch[i]=0;
	rechitrphimatch[i]=0;

	rechitmatchedx[i] =0;
	rechitmatchedy[i] =0;
	rechitmatchedz[i] =0;
	rechitmatchederrxx[i] =0;
	rechitmatchederrxy[i] =0;
	rechitmatchederryy[i] =0;
	rechitmatchedmatch[i]=0;
	rechitmatchedresx[i]=0;
	rechitmatchedresy[i]=0;

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
	Tibrphires[i]=-999.;
	Tibrphimatch[i]=0;
	Tibrphisiz[i]=0;
	Tibrphichg[i]=0;
	Tibsasx[i]=0;
	Tibsasy[i]=0;
	Tibsasz[i]=0;
	Tibsasphi[i]=0;
	Tibsasres[i]=-999;
	Tibsasmatch[i]=0;
	Tibsassiz[i]=0;
	Tibsaschg[i]=0;
	Tibmatchedx[i]=0;
	Tibmatchedy[i]=0;
	Tibmatchedz[i]=0;
	Tibmatchederrxx[i]=0;
	Tibmatchederrxy[i]=0;
	Tibmatchederryy[i]=0;
	Tibmatchedresx[i]=-999;	
	Tibmatchedresy[i]=-999;
	Tibmatchedmatch[i]=0;
	
	Tobsimhitx[i]=0;
	Tobsimhity[i]=0;
	Tobsimhitz[i]=0;
	Tobsimhitphi[i]=0;
	Tobsimhiteta[i]=0;
	Tobrphix[i]=0;
	Tobrphiy[i]=0;
	Tobrphiz[i]=0;
	Tobrphiphi[i]=0;
	Tobrphires[i]=-999;
	Tobrphimatch[i]=0;
	Tobrphisiz[i]=0;
	Tobrphichg[i]=0;
	Tobsasx[i]=0;
	Tobsasy[i]=0;
	Tobsasz[i]=0;
	Tobsasphi[i]=0;
	Tobsasres[i]=-999;
	Tobsasmatch[i]=0;
	Tobsassiz[i]=0;
	Tobsaschg[i]=0;	
	Tobmatchedx[i]=0;
	Tobmatchedy[i]=0;
	Tobmatchedz[i]=0;
	Tobmatchederrxx[i]=0;
	Tobmatchederrxy[i]=0;
	Tobmatchederryy[i]=0;
	Tobmatchedresx[i]=-999;	
	Tobmatchedresy[i]=-999;
	Tobmatchedmatch[i]=0;

	Tidsimhitx[i]=0;
	Tidsimhity[i]=0;
	Tidsimhitz[i]=0;
	Tidsimhitphi[i]=0;
	Tidsimhiteta[i]=0;
	Tidrphix[i]=0;
	Tidrphiy[i]=0;
	Tidrphiz[i]=0;
	Tidrphiphi[i]=0;
	Tidrphires[i]=-999;
	Tidrphimatch[i]=0;
	Tidrphisiz[i]=0;
	Tidrphichg[i]=0;
	Tidsasx[i]=0;
	Tidsasy[i]=0;
	Tidsasz[i]=0;
	Tidsasphi[i]=0;
	Tidsasres[i]=-999;
	Tidsasmatch[i]=0;
	Tidsassiz[i]=0;
	Tidsaschg[i]=0;	
	Tidmatchedx[i]=0;
	Tidmatchedy[i]=0;
	Tidmatchedz[i]=0;
	Tidmatchederrxx[i]=0;
	Tidmatchederrxy[i]=0;
	Tidmatchederryy[i]=0;
	Tidmatchedresx[i]=-999;	
	Tidmatchedresy[i]=-999;
	Tidmatchedmatch[i]=0;

	Tecsimhitx[i]=0;
	Tecsimhity[i]=0;
	Tecsimhitz[i]=0;
	Tecsimhitphi[i]=0;
	Tecsimhiteta[i]=0;
	Tecrphix[i]=0;
	Tecrphiy[i]=0;
	Tecrphiz[i]=0;
	Tecrphiphi[i]=0;
	Tecrphimatch[i]=0;
	Tecrphires[i]=-999;
	Tecrphisiz[i]=0;
	Tecrphichg[i]=0;
	Tecsasx[i]=0;
	Tecsasy[i]=0;
	Tecsasz[i]=0;
	Tecsasphi[i]=0;
	Tecsasres[i]=-999;
	Tecsasmatch[i]=0;
	Tecsassiz[i]=0;
	Tecsaschg[i]=0;	
	Tecmatchedx[i]=0;
	Tecmatchedy[i]=0;
	Tecmatchedz[i]=0;
	Tecmatchederrxx[i]=0;
	Tecmatchederrxy[i]=0;
	Tecmatchederryy[i]=0;
	Tecmatchedresx[i]=-999;	
	Tecmatchedresy[i]=-999;
	Tecmatchedmatch[i]=0;
      }



      TrackerHitAssociator  associate(e);

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
	  //try association here
	  //--- add matching code here for testing

	  //	  cout << "ValHit ---> try association RPHI! " << endl;
	  matched.clear();
	  matched = associate.associateHit(rechit);
	  if(!matched.empty()){
// 	    cout << " detector =  " << myid << " clusize = " << clusiz << " Rechit x= " << position.x() 
// 		 << " y = "<< position.y() << " z = " << position.z() << endl;
// 	    cout << " matched = " << matched.size() << endl;
// 	    for(vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
// 	      cout << " hit  ID = " << (*m).trackId() << " Simhit x = " << (*m).localPosition().x() 
// 		   << " y = " <<  (*m).localPosition().y() << " z = " <<  (*m).localPosition().x() << endl;
// 	    }
	    rechitrphimatch[i] = 1;
	    rechitrphires[i] = rechitrphix[i] - matched[0].localPosition().x();
	  }
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

	  //try association here
	  //--- add matching code here for testing

	  //	  cout << "ValHit ---> try association SAS! " << endl;
	  matched.clear();
	  matched = associate.associateHit(rechit);
	  if(!matched.empty()){
// 	    cout << " detector = " << myid << " clusize = " << clusiz << " Rechit x = " << position.x() 
// 		 << " y = "<< position.y() << " z = " << position.z() << endl;
// 	    cout << " matched = " << matched.size() << endl;
// 	    for(vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
// 	      cout << " hit  ID = " << (*m).trackId() << " Simhit x = " << (*m).localPosition().x() 
// 		   << " y = " <<  (*m).localPosition().y() << " z = " <<  (*m).localPosition().x() << endl;
// 	    }
	    rechitsasmatch[j] = 1;
	    rechitsasres[j] = rechitsasx[j] - matched[0].localPosition().x();
	  }
	  j++;
	}
      }

      //now matched hits

      //loop over rechits-matched in the same subdetector
      numrechitmatched=0;
      SiStripRecHit2DMatchedLocalPosCollection::range rechitmatchedRange = rechitsmatched->get(detid);
      SiStripRecHit2DMatchedLocalPosCollection::const_iterator rechitmatchedRangeIteratorBegin = rechitmatchedRange.first;
      SiStripRecHit2DMatchedLocalPosCollection::const_iterator rechitmatchedRangeIteratorEnd   = rechitmatchedRange.second;
      SiStripRecHit2DMatchedLocalPosCollection::const_iterator itermatched=rechitmatchedRangeIteratorBegin;
      numrechitmatched = rechitmatchedRangeIteratorEnd - rechitmatchedRangeIteratorBegin;   
      if(numrechitmatched > 0){
	int j=0;
	for(itermatched=rechitmatchedRangeIteratorBegin; itermatched!=rechitmatchedRangeIteratorEnd;++itermatched){
	  SiStripRecHit2DMatchedLocalPos const rechit=*itermatched;
	  LocalPoint position=rechit.localPosition();
	  LocalError error=rechit.localPositionError();

	  //try association here
	  //--- add matching code here for testing 

	  cout << "ValHit ---> try association matched! " << endl;
	  matched.clear();
	  const SiStripRecHit2DLocalPos *mono = rechit.monoHit();
	  const SiStripRecHit2DLocalPos *st = rechit.stereoHit();
	  LocalPoint monopos = mono->localPosition();
	  LocalPoint stpos   = st->localPosition();

	  rechitmatchedx[j] = position.x();
	  rechitmatchedy[j] = position.y();
	  rechitmatchedz[j] = position.z();
	  //rechitmatchedphi[j] = position.phi();
	  rechitmatchederrxx[j] = error.xx();
	  rechitmatchederrxy[j] = error.xy();
	  rechitmatchederryy[j] = error.yy();
	  matched = associate.associateHit(*st);
	  if(!matched.empty()){
	    cout << " detector = " << myid << " #match = " << matched.size() << endl;
	    cout << " Matched x = " << position.x() << " y = "<< position.y() << " z = " << position.z() << endl;
	    cout << " Mono    x = " << monopos.x() << " y = "<< monopos.y() << " z = " << monopos.z() << endl;
	    cout << " Stereo  x = " << stpos.x() << " y = "<< stpos.y() << " z = " << stpos.z() << endl;

	    for(vector<PSimHit>::const_iterator m=matched.begin(); m<matched.end(); m++){
	      cout << " hit  ID = " << (*m).trackId() << " Simhit x = " << (*m).localPosition().x() 
		   << " y = " <<  (*m).localPosition().y() << " z = " <<  (*m).localPosition().x() << endl;
	    }
	    rechitmatchedmatch[j] = 1;
	    //project simhit;
	    const GluedGeomDet* gluedDet = (const GluedGeomDet*)tracker.idToDet(rechit.geographicalId());
	    const StripGeomDetUnit* partnerstripdet =(StripGeomDetUnit*) gluedDet->stereoDet();
	    std::pair<LocalPoint,LocalVector> hitPair= projectHit(matched[0],partnerstripdet,gluedDet->surface());
	    //	    rechitmatchedresx[j] = rechitmatchedx[j] - matched[0].localPosition().x();
	    // rechitmatchedresy[j] = rechitmatchedy[j] - matched[0].localPosition().y();
	    rechitmatchedresx[j] = rechitmatchedx[j] - hitPair.first.x();
	    rechitmatchedresy[j] = rechitmatchedy[j] - hitPair.first.y();
	    
	    cout << " res x = " << rechitmatchedresx[j] << " rec(x) = " <<  rechitmatchedx[j] 
	      //		 << " sim(x) = " << matched[0].localPosition().x() << endl;
	      		 << " sim(x) = " << hitPair.first.x() << endl;
	    cout << " res y = " << rechitmatchedresy[j] << " rec(y) = " <<  rechitmatchedy[j] 
	      //		 << " sim(x) = " << matched[0].localPosition().y() << endl;
	      		 << " sim(y) = " <<  hitPair.first.y()<< endl;
	  }

	  j++;
	}
      }
      
      //--- fill histograms
      
      
      if(numsimhit>0 || numrechitrphi>0 || numrechitsas>0 ){
// 	if(numrechitmatched>0){ 
// 	  std::cout << " det id= " << myid << " N(simhit) = " << numsimhit 
// 		    << " N(rechitrphi) = " << numrechitrphi << " N(rechitsas)= " << numrechitsas 
// 		    << " N(matched) = " << numrechitmatched << std::endl;      
// 	}
	if ( detid.subdetId() == int(StripSubdetector::TIB)){
	  TIBDetId tibid(myid); 
	  Tiblayer   = tibid.layer();
	  Tibfw_bw   = tibid.string()[0];
	  Tibext_int = tibid.string()[1];
	  Tibstring  = tibid.string()[2];
	  Tibmodule  = tibid.module();
	  Tibstereo  = tibid.stereo();
	  
	  //	  Tiblayer = tiblayer;
	  //Tibstereo = tibstereo;
	  Tibnumsimhit = numsimhit;
	  Tibnumrechitrphi = numrechitrphi;
	  Tibnumrechitsas = numrechitsas;
	  Tibnumrechitmatched = numrechitmatched;
	  
	  for(int kk=0;kk<numrechitmatched; kk++){
	    Tibmatchedx[kk] = rechitmatchedx[kk];
	    Tibmatchedy[kk] = rechitmatchedy[kk];
	    Tibmatchedz[kk] = rechitmatchedz[kk];
	    Tibmatchederrxx[kk] = rechitmatchederrxx[kk];
	    Tibmatchederrxy[kk] = rechitmatchederrxy[kk];
	    Tibmatchederryy[kk] = rechitmatchederryy[kk];
	    Tibmatchedresx[kk]= rechitmatchedresx[kk];
	    Tibmatchedresy[kk]= rechitmatchedresy[kk];
	    Tibmatchedmatch[kk]= rechitmatchedmatch[kk];
	  }
	  
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
	      Tibrphires[kk]= rechitrphires[kk];
	      Tibrphimatch[kk]= rechitrphimatch[kk];
	      Tibrphisiz[kk] = clusizrphi[kk];
	      Tibrphichg[kk] = cluchgrphi[kk];
	    }
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
	      Tibsasres[kk]= rechitsasres[kk];
	      Tibsasmatch[kk]= rechitsasmatch[kk];
	      Tibsassiz[kk] = clusizsas[kk];
	      Tibsaschg[kk] = cluchgsas[kk];
	    }
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
	  Tobnumrechitmatched = numrechitmatched;

	  for(int kk=0;kk<numrechitmatched; kk++){
	    Tobmatchedx[kk] = rechitmatchedx[kk];
	    Tobmatchedy[kk] = rechitmatchedy[kk];
	    Tobmatchedz[kk] = rechitmatchedz[kk];
	    Tobmatchederrxx[kk] = rechitmatchederrxx[kk];
	    Tobmatchederrxy[kk] = rechitmatchederrxy[kk];
	    Tobmatchederryy[kk] = rechitmatchederryy[kk];
	    Tobmatchedresx[kk]= rechitmatchedresx[kk];
	    Tobmatchedresy[kk]= rechitmatchedresy[kk];
	    Tobmatchedmatch[kk]= rechitmatchedmatch[kk];
	  }
	  
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
	      Tobrphires[kk]= rechitrphires[kk];
	      Tobrphimatch[kk]= rechitrphimatch[kk];
	      Tobrphisiz[kk] = clusizrphi[kk];
	      Tobrphichg[kk] = cluchgrphi[kk];
	    }	
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
	      Tobsasres[kk]= rechitsasres[kk];
	      Tobsasmatch[kk]= rechitsasmatch[kk];
	      Tobsassiz[kk] = clusizsas[kk];
	      Tobsaschg[kk] = cluchgsas[kk];
	    }
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
	  Tidnumrechitmatched = numrechitmatched;
	  for(int kk=0;kk<numrechitmatched; kk++){
	    Tidmatchedx[kk] = rechitmatchedx[kk];
	    Tidmatchedy[kk] = rechitmatchedy[kk];
	    Tidmatchedz[kk] = rechitmatchedz[kk];
	    Tidmatchederrxx[kk] = rechitmatchederrxx[kk];
	    Tidmatchederrxy[kk] = rechitmatchederrxy[kk];
	    Tidmatchederryy[kk] = rechitmatchederryy[kk];
	    Tidmatchedresx[kk]= rechitmatchedresx[kk];
	    Tidmatchedresy[kk]= rechitmatchedresy[kk];
	    Tidmatchedmatch[kk]= rechitmatchedmatch[kk];
	  }
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
	      Tidrphires[kk]= rechitrphires[kk];
	      Tidrphimatch[kk]= rechitrphimatch[kk];
	      Tidrphisiz[kk] = clusizrphi[kk];
	      Tidrphichg[kk] = cluchgrphi[kk];
	    }

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
	      Tidsasres[kk]= rechitsasres[kk];
	      Tidsasmatch[kk]= rechitsasmatch[kk];
	      Tidsassiz[kk] = clusizsas[kk];
	      Tidsaschg[kk] = cluchgsas[kk];
	    }
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
	  Tecnumrechitmatched = numrechitmatched;
	  for(int kk=0;kk<numrechitmatched; kk++){
	    Tecmatchedx[kk] = rechitmatchedx[kk];
	    Tecmatchedy[kk] = rechitmatchedy[kk];
	    Tecmatchedz[kk] = rechitmatchedz[kk];
	    Tecmatchederrxx[kk] = rechitmatchederrxx[kk];
	    Tecmatchederrxy[kk] = rechitmatchederrxy[kk];
	    Tecmatchederryy[kk] = rechitmatchederryy[kk];
	    Tecmatchedresx[kk]= rechitmatchedresx[kk];
	    Tecmatchedresy[kk]= rechitmatchedresy[kk];
	    Tecmatchedmatch[kk]= rechitmatchedmatch[kk];
	  }
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
	      Tecrphires[kk]= rechitrphires[kk];
	      Tecrphimatch[kk]= rechitrphimatch[kk];
	      Tecrphisiz[kk] = clusizrphi[kk];
	      Tecrphichg[kk] = cluchgrphi[kk];
	    }
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
	      Tecsasres[kk]= rechitsasres[kk];
	      Tecsasmatch[kk]= rechitsasmatch[kk];
	      Tecsassiz[kk] = clusizsas[kk];
	      Tecsaschg[kk] = cluchgsas[kk];
	    }
	  }
	}
	myTree->Fill();     
      }
    }
    cout << " === calling end job " << endl;  
  }
  
  
  ValHit::ValHit(edm::ParameterSet const& conf) : 
    conf_(conf),filename_(conf.getParameter<std::string>("fileName")) 

  {
    myFile = new TFile(filename_.c_str(),"RECREATE");
    myTree = new TTree("HitTree","Tracker Validation tree");
    // GENERAL block
    myTree->Branch("Run", &myRun, "Run/I");
    myTree->Branch("Event", &myEvent, "Event/I");

    //TIB info
    myTree->Branch("Tiblayer", &Tiblayer, "Tiblayer/I");
    myTree->Branch("Tibstereo", &Tibstereo, "Tibstereo/I");
    myTree->Branch("Tibfw_bw", &Tibfw_bw, "Tibfw_bw/I");
    myTree->Branch("Tibext_int", &Tibext_int, "Tibext_int/I");
    myTree->Branch("Tibstring", &Tibstring, "Tibstring/I");
    myTree->Branch("Tibmodule", &Tibmodule, "Tibmodule/I");
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
    myTree->Branch("Tibrphires", &Tibrphires, "Tibrphires[Tibnumrechitrphi]/F");
    myTree->Branch("Tibrphimatch", &Tibrphimatch, "Tibrphimatch[Tibnumrechitrphi]/I");
    myTree->Branch("Tibrphisiz", &Tibrphisiz, "Tibrphisiz[Tibnumrechitrphi]/I");
    myTree->Branch("Tibrphichg", &Tibrphichg, "Tibrphichg[Tibnumrechitrphi]/F");
    Tibnumrechitsas=0;
    myTree->Branch("Tibnumrechitsas", &Tibnumrechitsas, "Tibnumrechitsas/I");
    myTree->Branch("Tibsasx", &Tibsasx, "Tibsasx[Tibnumrechitsas]/F"); 
    myTree->Branch("Tibsasy", &Tibsasy, "Tibsasy[Tibnumrechitsas]/F");
    myTree->Branch("Tibsasz", &Tibsasz, "Tibsasz[Tibnumrechitsas]/F");        
    myTree->Branch("Tibsasres", &Tibsasres, "Tibsasres[Tibnumrechitsas]/F");
    myTree->Branch("Tibsasmatch", &Tibsasmatch, "Tibsasmatch[Tibnumrechitsas]/I");
    myTree->Branch("Tibsassiz", &Tibsassiz, "Tibsassiz[Tibnumrechitsas]/I");
    myTree->Branch("Tibsaschg", &Tibsaschg, "Tibsaschg[Tibnumrechitsas]/F");
    Tibnumrechitmatched=0;
    myTree->Branch("Tibnumrechitmatched", &Tibnumrechitmatched, "Tibnumrechitmatched/I");
    myTree->Branch("Tibmatchedx", &Tibmatchedx, "Tibmatchedx[Tibnumrechitmatched]/F"); 
    myTree->Branch("Tibmatchedy", &Tibmatchedy, "Tibmatchedy[Tibnumrechitmatched]/F");
    myTree->Branch("Tibmatchedz", &Tibmatchedz, "Tibmatchedz[Tibnumrechitmatched]/F");        
    myTree->Branch("Tibmatchederrxx", &Tibmatchederrxx, "Tibmatchederrxx[Tibnumrechitmatched]/F"); 
    myTree->Branch("Tibmatchederrxy", &Tibmatchederrxy, "Tibmatchederrxy[Tibnumrechitmatched]/F");
    myTree->Branch("Tibmatchederryy", &Tibmatchederryy, "Tibmatchederryy[Tibnumrechitmatched]/F");        
    myTree->Branch("Tibmatchedresx", &Tibmatchedresx, "Tibmatchedresx[Tibnumrechitmatched]/F");
    myTree->Branch("Tibmatchedresy", &Tibmatchedresy, "Tibmatchedresy[Tibnumrechitmatched]/F");
    myTree->Branch("Tibmatchedmatch", &Tibmatchedmatch, "Tibmatchedmatch[Tibnumrechitmatched]/I");
 
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
    myTree->Branch("Tobrphires", &Tobrphires, "Tobrphires[Tobnumrechitrphi]/F");
    myTree->Branch("Tobrphimatch", &Tobrphimatch, "Tobrphimatch[Tobnumrechitrphi]/I");
    myTree->Branch("Tobrphisiz", &Tobrphisiz, "Tobrphisiz[Tobnumrechitrphi]/I");
    myTree->Branch("Tobrphichg", &Tobrphichg, "Tobrphichg[Tobnumrechitrphi]/F");
     Tobnumrechitsas=0;
    myTree->Branch("Tobnumrechitsas", &Tobnumrechitsas, "Tobnumrechitsas/I");
    myTree->Branch("Tobsasx", &Tobsasx, "Tobsasx[Tobnumrechitsas]/F"); 
    myTree->Branch("Tobsasy", &Tobsasy, "Tobsasy[Tobnumrechitsas]/F");
    myTree->Branch("Tobsasz", &Tobsasz, "Tobsasz[Tobnumrechitsas]/F");        
    myTree->Branch("Tobsasres", &Tobsasres, "Tobsasres[Tobnumrechitsas]/F");
    myTree->Branch("Tobsasmatch", &Tobsasmatch, "Tobsasmatch[Tobnumrechitsas]/I");
    myTree->Branch("Tobsassiz", &Tobsassiz, "Tobsassiz[Tobnumrechitsas]/I");
    myTree->Branch("Tobsaschg", &Tobsaschg, "Tobsaschg[Tobnumrechitsas]/F");
    Tobnumrechitmatched=0;
    myTree->Branch("Tobnumrechitmatched", &Tobnumrechitmatched, "Tobnumrechitmatched/I");
    myTree->Branch("Tobmatchedx", &Tobmatchedx, "Tobmatchedx[Tobnumrechitmatched]/F"); 
    myTree->Branch("Tobmatchedy", &Tobmatchedy, "Tobmatchedy[Tobnumrechitmatched]/F");
    myTree->Branch("Tobmatchedz", &Tobmatchedz, "Tobmatchedz[Tobnumrechitmatched]/F");        
    myTree->Branch("Tobmatchederrxx", &Tobmatchederrxx, "Tobmatchederrxx[Tobnumrechitmatched]/F"); 
    myTree->Branch("Tobmatchederrxy", &Tobmatchederrxy, "Tobmatchederrxy[Tobnumrechitmatched]/F");
    myTree->Branch("Tobmatchederryy", &Tobmatchederryy, "Tobmatchederryy[Tobnumrechitmatched]/F");        
    myTree->Branch("Tobmatchedresx", &Tobmatchedresx, "Tobmatchedresx[Tobnumrechitmatched]/F");
    myTree->Branch("Tobmatchedresy", &Tobmatchedresy, "Tobmatchedresy[Tobnumrechitmatched]/F");
    myTree->Branch("Tobmatchedmatch", &Tobmatchedmatch, "Tobmatchedmatch[Tobnumrechitmatched]/I");
 
 
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
    myTree->Branch("Tidrphires", &Tidrphires, "Tidrphires[Tidnumrechitrphi]/F");
    myTree->Branch("Tidrphimatch", &Tidrphimatch, "Tidrphimatch[Tidnumrechitrphi]/I");
    myTree->Branch("Tidrphisiz", &Tidrphisiz, "Tidrphisiz[Tidnumrechitrphi]/I");
    myTree->Branch("Tidrphichg", &Tidrphichg, "Tidrphichg[Tidnumrechitrphi]/F");
     Tidnumrechitsas=0;
    myTree->Branch("Tidnumrechitsas", &Tidnumrechitsas, "Tidnumrechitsas/I");
    myTree->Branch("Tidsasx", &Tidsasx, "Tidsasx[Tidnumrechitsas]/F"); 
    myTree->Branch("Tidsasy", &Tidsasy, "Tidsasy[Tidnumrechitsas]/F");
    myTree->Branch("Tidsasz", &Tidsasz, "Tidsasz[Tidnumrechitsas]/F");        
    myTree->Branch("Tidsasres", &Tidsasres, "Tidsasres[Tidnumrechitsas]/F");
    myTree->Branch("Tidsasmatch", &Tidsasmatch, "Tidsasmatch[Tidnumrechitsas]/I");
    myTree->Branch("Tidsassiz", &Tidsassiz, "Tidsassiz[Tidnumrechitsas]/I");
    myTree->Branch("Tidsaschg", &Tidsaschg, "Tidsaschg[Tidnumrechitsas]/F");
     Tidnumrechitmatched=0;
    myTree->Branch("Tidnumrechitmatched", &Tidnumrechitmatched, "Tidnumrechitmatched/I");
    myTree->Branch("Tidmatchedx", &Tidmatchedx, "Tidmatchedx[Tidnumrechitmatched]/F"); 
    myTree->Branch("Tidmatchedy", &Tidmatchedy, "Tidmatchedy[Tidnumrechitmatched]/F");
    myTree->Branch("Tidmatchedz", &Tidmatchedz, "Tidmatchedz[Tidnumrechitmatched]/F");        
    myTree->Branch("Tidmatchederrxx", &Tidmatchederrxx, "Tidmatchederrxx[Tidnumrechitmatched]/F"); 
    myTree->Branch("Tidmatchederrxy", &Tidmatchederrxy, "Tidmatchederrxy[Tidnumrechitmatched]/F");
    myTree->Branch("Tidmatchederryy", &Tidmatchederryy, "Tidmatchederryy[Tidnumrechitmatched]/F");        
    myTree->Branch("Tidmatchedresx", &Tidmatchedresx, "Tidmatchedresx[Tidnumrechitmatched]/F");
    myTree->Branch("Tidmatchedresy", &Tidmatchedresy, "Tidmatchedresy[Tidnumrechitmatched]/F");
    myTree->Branch("Tidmatchedmatch", &Tidmatchedmatch, "Tidmatchedmatch[Tidnumrechitmatched]/I");
 

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
    myTree->Branch("Tecrphires", &Tecrphires, "Tecrphires[Tecnumrechitrphi]/F");
    myTree->Branch("Tecrphimatch", &Tecrphimatch, "Tecrphimatch[Tecnumrechitrphi]/I");
    myTree->Branch("Tecrphisiz", &Tecrphisiz, "Tecrphisiz[Tecnumrechitrphi]/I");
    myTree->Branch("Tecrphichg", &Tecrphichg, "Tecrphichg[Tecnumrechitrphi]/F");
     Tecnumrechitsas=0;
    myTree->Branch("Tecnumrechitsas", &Tecnumrechitsas, "Tecnumrechitsas/I");
    myTree->Branch("Tecsasx", &Tecsasx, "Tecsasx[Tecnumrechitsas]/F"); 
    myTree->Branch("Tecsasy", &Tecsasy, "Tecsasy[Tecnumrechitsas]/F");
    myTree->Branch("Tecsasz", &Tecsasz, "Tecsasz[Tecnumrechitsas]/F");        
    myTree->Branch("Tecsasres", &Tecsasres, "Tecsasres[Tecnumrechitsas]/F");
    myTree->Branch("Tecsasmatch", &Tecsasmatch, "Tecsasmatch[Tecnumrechitsas]/I");
    myTree->Branch("Tecsassiz", &Tecsassiz, "Tecsassiz[Tecnumrechitsas]/I");
    myTree->Branch("Tecsaschg", &Tecsaschg, "Tecsaschg[Tecnumrechitsas]/F");
    Tecnumrechitmatched=0;
    myTree->Branch("Tecnumrechitmatched", &Tecnumrechitmatched, "Tecnumrechitmatched/I");
    myTree->Branch("Tecmatchedx", &Tecmatchedx, "Tecmatchedx[Tecnumrechitmatched]/F"); 
    myTree->Branch("Tecmatchedy", &Tecmatchedy, "Tecmatchedy[Tecnumrechitmatched]/F");
    myTree->Branch("Tecmatchedz", &Tecmatchedz, "Tecmatchedz[Tecnumrechitmatched]/F");        
    myTree->Branch("Tecmatchederrxx", &Tecmatchederrxx, "Tecmatchederrxx[Tecnumrechitmatched]/F"); 
    myTree->Branch("Tecmatchederrxy", &Tecmatchederrxy, "Tecmatchederrxy[Tecnumrechitmatched]/F");
    myTree->Branch("Tecmatchederryy", &Tecmatchederryy, "Tecmatchederryy[Tecnumrechitmatched]/F");        
    myTree->Branch("Tecmatchedresx", &Tecmatchedresx, "Tecmatchedresx[Tecnumrechitmatched]/F");
    myTree->Branch("Tecmatchedresy", &Tecmatchedresy, "Tecmatchedresy[Tecnumrechitmatched]/F");
    myTree->Branch("Tecmatchedmatch", &Tecmatchedmatch, "Tecmatchedmatch[Tecnumrechitmatched]/I");
 
 
  }
  
  
  // Virtual destructor needed.
  ValHit::~ValHit(){ 
    delete myFile;
  }
  
  std::pair<LocalPoint,LocalVector> ValHit::projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
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


  void ValHit::endJob() {  
    cout << ">>> ending histograms" << endl;
    myFile->cd();
    myTree->Write();
    myFile->Close();
    cout << ">>> File closed " << endl;
  }


}




