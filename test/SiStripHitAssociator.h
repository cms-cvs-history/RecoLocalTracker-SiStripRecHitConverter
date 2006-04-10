#ifndef SiStripHitAssociator_h
#define SiStripHitAssociator_h

/* \class SiStripAssociator
 *
 ** Associates SimHits and RecHits based on information produced during
 *  digitisation (StripDigiSimLinks).
 *  The association works in both ways: from a SimHit to RecHits and
 *  from a RecHit to SimHits.
 *
 * \author Patrizia Azzi, INFN PD
 *
 * \version   1st version April 2006  
 *
 *
 ************************************************************/

//#include <vector>
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/EDProduct.h"

//--- for SimHit
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLinkCollection.h"

//--- for RecHit
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DLocalPos.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DMatchedLocalPos.h"


namespace cms{

   class SiStripHitAssociator {

  public:
     
     SiStripHitAssociator(const edm::Event& e);
     virtual ~SiStripHitAssociator(){}
     
     std::vector<PSimHit> associateSimpleRecHit(const SiStripRecHit2DLocalPos & rechit);
     
     //will do next
     //  vector<const SimHit*> associateMatchedRecHit( const RecHit&) const;
     
     std::vector<PSimHit> theStripHits;
     typedef std::map<unsigned int, std::vector<PSimHit> > simhit_map;
     typedef simhit_map::iterator simhit_map_iterator;
     simhit_map SimHitMap;
     StripDigiSimLinkCollection::Range linkrange;

  private:
     const edm::Event& myEvent_;
     edm::Handle<StripDigiSimLinkCollection> stripdigisimlink;
   };  
   
}

#endif

