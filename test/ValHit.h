#ifndef ValHit_h
#define ValHit_h

/* \class ValHit
 *
 * ValHit is the analyzer which validates hits in the tracker: 
 * SimHit, Digi, Clusters, RecHit
 *
 * \author Patrizia Azzi, INFN PD
 *
 * \version   1st version Feb 2006  

 *
 ************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
 
//--- for SimHit
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

//--- for RecHit
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/Common/interface/OwnVector.h"

//--- for StripDigiSimLink
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

//needed for the geometry:
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
 
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH2.h"

namespace cms {

class ValHit : public edm::EDAnalyzer
{
 public:
  
  explicit ValHit(const edm::ParameterSet& conf);
  
  virtual ~ValHit();

  void endJob();
  std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
							const BoundPlane& plane);
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
  std::vector<PSimHit> matched;

 private:
  
  edm::ParameterSet conf_;
  const StripTopology* topol;
  int numStrips;    // number of strips in the module
  std::vector<PSimHit> theStripHits;
  typedef std::map<unsigned int, std::vector<PSimHit> > simhit_map;
  typedef simhit_map::iterator simhit_map_iterator;
  simhit_map SimHitMap;
  
  //for the output NTUPLE
  std::string filename_;
  static const int MAXHIT = 100;

  int tiblayer, tibstereo, tibfw_bw, tibext_int, tibstring, tibmodule;
  float simhitx[MAXHIT];
  float simhity[MAXHIT];
  float simhitz[MAXHIT];
  float simhitphi[MAXHIT];
  float simhiteta[MAXHIT];
  float rechitrphix[MAXHIT];
  float rechitrphiy[MAXHIT];
  float rechitrphiz[MAXHIT];
  float rechitrphiphi[MAXHIT];
  float rechitrphires[MAXHIT];
  int rechitrphimatch[MAXHIT];
  int clusizrphi[MAXHIT];
  float cluchgrphi[MAXHIT];
  float rechitsasx[MAXHIT];
  float rechitsasy[MAXHIT];
  float rechitsasz[MAXHIT];
  float rechitsasphi[MAXHIT];
  float rechitsasres[MAXHIT];
  int rechitsasmatch[MAXHIT];
  int clusizsas[MAXHIT];
  float cluchgsas[MAXHIT];
  float rechitmatchedx[MAXHIT];
  float rechitmatchedy[MAXHIT];
  float rechitmatchedz[MAXHIT];
  float rechitmatchederrxx[MAXHIT];
  float rechitmatchederrxy[MAXHIT];
  float rechitmatchederryy[MAXHIT];
  float rechitmatchedphi[MAXHIT];
  float rechitmatchedresx[MAXHIT];
  float rechitmatchedresy[MAXHIT];
  int rechitmatchedmatch[MAXHIT];

  //tree variable declaration
  TFile* myFile;
  TTree* myTree;
  unsigned int myRun;
  unsigned int myEvent;
  int Tiblayer;
  int Tibstereo; 
  int Tibfw_bw; 
  int Tibext_int; 
  int Tibstring; 
  int Tibmodule; 
  int Tibnumsimhit;
  int Tibnumrechitrphi;
  int Tibnumrechitsas;
  int Tibnumrechitmatched;
  float Tibsimhitx[MAXHIT];
  float Tibsimhity[MAXHIT];
  float Tibsimhitz[MAXHIT];
  float Tibsimhitphi[MAXHIT];
  float Tibsimhiteta[MAXHIT];
  float Tibrphix[MAXHIT];
  float Tibrphiy[MAXHIT];
  float Tibrphiz[MAXHIT];
  float Tibrphiphi[MAXHIT];
  float Tibrphires[MAXHIT]; 
  int Tibrphimatch[MAXHIT];
  int Tibrphisiz[MAXHIT];
  float Tibrphichg[MAXHIT];
  float Tibsasx[MAXHIT];
  float Tibsasy[MAXHIT];
  float Tibsasz[MAXHIT];
  float Tibsasphi[MAXHIT];
  float Tibsasres[MAXHIT];
  int Tibsasmatch[MAXHIT];
  int Tibsassiz[MAXHIT];
  float Tibsaschg[MAXHIT];
  float Tibmatchedx[MAXHIT];
  float Tibmatchedy[MAXHIT];
  float Tibmatchedz[MAXHIT];
  float Tibmatchederrxx[MAXHIT];
  float Tibmatchederrxy[MAXHIT];
  float Tibmatchederryy[MAXHIT];
  float Tibmatchedresx[MAXHIT];
  float Tibmatchedresy[MAXHIT];
  int Tibmatchedmatch[MAXHIT];

  unsigned int toblayer;
  unsigned int tobstereo;
  int Toblayer;
  int Tobstereo; 
  int Tobnumsimhit;
  int Tobnumrechitrphi;
  int Tobnumrechitsas;
  int Tobnumrechitmatched;
  float Tobsimhitx[MAXHIT];
  float Tobsimhity[MAXHIT];
  float Tobsimhitz[MAXHIT];
  float Tobsimhiteta[MAXHIT];
  float Tobsimhitphi[MAXHIT];
  float Tobrphix[MAXHIT];
  float Tobrphiy[MAXHIT];
  float Tobrphiz[MAXHIT];
  float Tobrphiphi[MAXHIT];
  float Tobrphires[MAXHIT];
  int Tobrphimatch[MAXHIT];
  int Tobrphisiz[MAXHIT];
  float Tobrphichg[MAXHIT];
  float Tobsasx[MAXHIT];
  float Tobsasy[MAXHIT];
  float Tobsasz[MAXHIT];
  float Tobsasphi[MAXHIT];
  float Tobsasres[MAXHIT];
  int Tobsasmatch[MAXHIT];
  int Tobsassiz[MAXHIT];
  float Tobsaschg[MAXHIT];
  float Tobmatchedx[MAXHIT];
  float Tobmatchedy[MAXHIT];
  float Tobmatchedz[MAXHIT];
  float Tobmatchederrxx[MAXHIT];
  float Tobmatchederrxy[MAXHIT];
  float Tobmatchederryy[MAXHIT];
  float Tobmatchedresx[MAXHIT];
  float Tobmatchedresy[MAXHIT];
  int Tobmatchedmatch[MAXHIT];


  unsigned int tidlayer;
  unsigned int tidstereo;
  int Tidlayer;
  int Tidstereo; 
  int Tidring; 
  int Tidnumsimhit;
  int Tidnumrechitrphi;
  int Tidnumrechitsas;
  int Tidnumrechitmatched;
  float Tidsimhitx[MAXHIT];
  float Tidsimhity[MAXHIT];
  float Tidsimhitz[MAXHIT];
  float Tidsimhiteta[MAXHIT];
  float Tidsimhitphi[MAXHIT];
  float Tidrphix[MAXHIT];
  float Tidrphiy[MAXHIT];
  float Tidrphiz[MAXHIT];
  float Tidrphiphi[MAXHIT];
  float Tidrphires[MAXHIT];
  int Tidrphimatch[MAXHIT];
  int Tidrphisiz[MAXHIT];
  float Tidrphichg[MAXHIT];
  float Tidsasx[MAXHIT];
  float Tidsasy[MAXHIT];
  float Tidsasz[MAXHIT];
  float Tidsasphi[MAXHIT];
  float Tidsasres[MAXHIT];
  int Tidsasmatch[MAXHIT];
  int Tidsassiz[MAXHIT];
  float Tidsaschg[MAXHIT];
  float Tidmatchedx[MAXHIT];
  float Tidmatchedy[MAXHIT];
  float Tidmatchedz[MAXHIT];
  float Tidmatchederrxx[MAXHIT];
  float Tidmatchederrxy[MAXHIT];
  float Tidmatchederryy[MAXHIT];
  float Tidmatchedresx[MAXHIT];
  float Tidmatchedresy[MAXHIT];
  int Tidmatchedmatch[MAXHIT];


  unsigned int teclayer;
  unsigned int tecstereo;
  int Teclayer;
  int Tecstereo; 
  int Tecring;
  int Tecnumsimhit;
  int Tecnumrechitrphi;
  int Tecnumrechitsas;
  int Tecnumrechitmatched;
  float Tecsimhitx[MAXHIT];
  float Tecsimhity[MAXHIT];
  float Tecsimhitz[MAXHIT];
  float Tecsimhitphi[MAXHIT];
  float Tecsimhiteta[MAXHIT];
  float Tecrphix[MAXHIT];
  float Tecrphiy[MAXHIT];
  float Tecrphiz[MAXHIT];
  float Tecrphiphi[MAXHIT];
  float Tecrphires[MAXHIT];
  int Tecrphimatch[MAXHIT];
  int Tecrphisiz[MAXHIT];
  float Tecrphichg[MAXHIT];
  float Tecsasx[MAXHIT];
  float Tecsasy[MAXHIT];
  float Tecsasz[MAXHIT];
  float Tecsasphi[MAXHIT];
  float Tecsasres[MAXHIT];
  int Tecsassiz[MAXHIT];
  int Tecsasmatch[MAXHIT];
  float Tecsaschg[MAXHIT];
  float Tecmatchedx[MAXHIT];
  float Tecmatchedy[MAXHIT];
  float Tecmatchedz[MAXHIT];
  float Tecmatchederrxx[MAXHIT];
  float Tecmatchederrxy[MAXHIT];
  float Tecmatchederryy[MAXHIT];
  float Tecmatchedresx[MAXHIT];
  float Tecmatchedresy[MAXHIT];
  int Tecmatchedmatch[MAXHIT];
    
};

}

#endif
