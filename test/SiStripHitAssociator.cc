// File: SiStripHitAssociator.cc

#include <memory>
#include <string>
#include <vector>

#include "RecoLocalTracker/SiStripRecHitConverter/test/SiStripHitAssociator.h"

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
//#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//--- for Geometry:
#include "DataFormats/DetId/interface/DetId.h"


using namespace std;
using namespace edm;

namespace cms{
  
  std::vector<PSimHit> 
  SiStripHitAssociator::associateSimpleRecHit(const SiStripRecHit2DLocalPos & rechit) {
    
    
    //vector with the matched SimHit
    std::vector<PSimHit> result; 
    
    //vector with the trackIds
    std::vector<unsigned int> simtrackid; 
    
    //get the Detector of the rechit
    DetId detid=  rechit.geographicalId();
    unsigned int detID = detid.rawId();
    //    cout << "Associator ---> get Detid " << detID << endl;
    
    //get the range of digis in detid
    linkrange =  stripdigisimlink->get(detID);
    int numlink=0;
    StripDigiSimLinkCollection::ContainerIterator linkbegin=linkrange.first;
    StripDigiSimLinkCollection::ContainerIterator linkend=linkrange.second;
    numlink = linkend-linkbegin;
    if(numlink>0){

      //      cout << "Associator ---> get digilink! in Detid n = " << numlink << endl;
 
      //get the cluster info
      const std::vector<const SiStripCluster*> clust=rechit.cluster();
      //   cout << "Associator ---> get cluster info " << endl;
      
      for(vector<const SiStripCluster*>::const_iterator ic = clust.begin(); ic!=clust.end(); ic++) {
	unsigned int clusiz = (*ic)->amplitudes().size();
	unsigned int first  = (*ic)->firstStrip();     
	unsigned int last   = first + clusiz;
	//	cout << "Associator ---> clus size = " << clusiz << " first = " << first << " last = " << last << endl;
	for(StripDigiSimLinkCollection::ContainerIterator linkiter=linkrange.first; linkiter!=linkrange.second;++linkiter){
	  StripDigiSimLink link = *linkiter;
	  if( link.channel() >= first  && link.channel() < last ){
	    simtrackid.push_back(link.SimTrackId());
	    //  cout << "Associator --> digi list first= " << first << " last = " << last << endl;
	    //cout << "Associator link--> channel= " << link.channel() << "  trackid = " << link.SimTrackId() << endl;
	  }
	}
      }
      
      //now get the SimHit from the trackid
      vector<PSimHit> simHit; 
      std::map<unsigned int, std::vector<PSimHit> >::const_iterator it = SimHitMap.find(detID);
      simHit.clear();
      if (it!= SimHitMap.end()){
	simHit = it->second;
	vector<PSimHit>::const_iterator simHitIter = simHit.begin();
	vector<PSimHit>::const_iterator simHitIterEnd = simHit.end();
	for (;simHitIter != simHitIterEnd; ++simHitIter) {
	  const PSimHit ihit = *simHitIter;
	  unsigned int simHitid = ihit.trackId();
	  for(size_t i=0; i<simtrackid.size();i++){
	    //cout << " Associator -->  check sihit id's = " << simHitid << endl;
	    if(simHitid == simtrackid[i] && simtrackid[i]!= 65535){ //exclude the geant particles. they all have the same id
	      // cout << "Associator ---> ID" << ihit.trackId() << " Simhit x= " << ihit.localPosition().x() 
	      //	   << " y= " <<  ihit.localPosition().y() << " z= " <<  ihit.localPosition().x() << endl;	    
	      result.push_back(ihit);
	    }
	  }
	}
      } 
    }
    return result; 
  }
  SiStripHitAssociator::SiStripHitAssociator(const edm::Event& e)  : myEvent_(e)  {
    //  using namespace edm;
    
    //get stuff from the event
    e.getByLabel("stripdigi", stripdigisimlink);
        
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
    
    SimHitMap.clear();
    for (std::vector<PSimHit>::iterator isim = theStripHits.begin();
	 isim != theStripHits.end(); ++isim){
      SimHitMap[(*isim).detUnitId()].push_back((*isim));
    }
  }
}
