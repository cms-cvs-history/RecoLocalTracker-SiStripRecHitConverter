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
 
//#include "RecoLocalTracker/SiStripRecHitConverter/test/ReadRecHitAlgorithm.h"
//--- for SimHit
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
                                                                                                                          
//needed for the geometry:
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

 
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
  
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);
  
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

  int tiblayer, tibstereo;
  float simhitx[MAXHIT];
  float simhity[MAXHIT];
  float simhitz[MAXHIT];
  float simhitphi[MAXHIT];
  float simhiteta[MAXHIT];
  float rechitrphix[MAXHIT];
  float rechitrphiy[MAXHIT];
  float rechitrphiz[MAXHIT];
  float rechitrphiphi[MAXHIT];
  float rechitsasx[MAXHIT];
  float rechitsasy[MAXHIT];
  float rechitsasz[MAXHIT];
  float rechitsasphi[MAXHIT];
  int clusizrphi[MAXHIT];
  int clusizsas[MAXHIT];
  float cluchgrphi[MAXHIT];
  float cluchgsas[MAXHIT];


  //tree variable declaration
  TFile* myFile;
  TTree* myTree;
  unsigned int myRun;
  unsigned int myEvent;
  int Tiblayer;
  int Tibstereo; 
  int Tibnumsimhit;
  int Tibnumrechitrphi;
  int Tibnumrechitsas;

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
  int Tibrphisiz[MAXHIT];
  float Tibrphichg[MAXHIT];
  float Tibsasx[MAXHIT];
  float Tibsasy[MAXHIT];
  float Tibsasz[MAXHIT];
  float Tibsasphi[MAXHIT];
  float Tibsasres[MAXHIT];
  int Tibsassiz[MAXHIT];
  float Tibsaschg[MAXHIT];

  int Toblayer;
  int Tobstereo; 
  int Tobnumsimhit;
  int Tobnumrechitrphi;
  int Tobnumrechitsas;
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
  int Tobrphisiz[MAXHIT];
  float Tobrphichg[MAXHIT];
  float Tobsasx[MAXHIT];
  float Tobsasy[MAXHIT];
  float Tobsasz[MAXHIT];
  float Tobsasphi[MAXHIT];
  float Tobsasres[MAXHIT];
  int Tobsassiz[MAXHIT];
  float Tobsaschg[MAXHIT];


  int Tidlayer;
  int Tidstereo; 
  int Tidring; 
  int Tidnumsimhit;
  int Tidnumrechitrphi;
  int Tidnumrechitsas;
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
  int Tidrphisiz[MAXHIT];
  float Tidrphichg[MAXHIT];
  float Tidsasx[MAXHIT];
  float Tidsasy[MAXHIT];
  float Tidsasz[MAXHIT];
  float Tidsasphi[MAXHIT];
  float Tidsasres[MAXHIT];
  int Tidsassiz[MAXHIT];
  float Tidsaschg[MAXHIT];


  int Teclayer;
  int Tecstereo; 
  int Tecring;
  int Tecnumsimhit;
  int Tecnumrechitrphi;
  int Tecnumrechitsas;
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
  int Tecrphisiz[MAXHIT];
  float Tecrphichg[MAXHIT];
  float Tecsasx[MAXHIT];
  float Tecsasy[MAXHIT];
  float Tecsasz[MAXHIT];
  float Tecsasphi[MAXHIT];
  float Tecsasres[MAXHIT];
  int Tecsassiz[MAXHIT];
  float Tecsaschg[MAXHIT];
    
  TH1F * tibres1rphi, * tibres2rphi, * tibres3rphi, * tibres4rphi, * tibres1sas, * tibres2sas;
  TH1F * tibclu1rphi, * tibclu2rphi, * tibclu3rphi, * tibclu4rphi, * tibclu1sas, * tibclu2sas;
  TH1F * tibhits;    
  TH1F * tib_rphi_digi_chg;
  TH2F * tib_rphi_tot_nst;
  TH2F * tib_rphi_pos_nst;


  unsigned int toblayer;
  unsigned int tobstereo;
  TH1F * tobres1rphi, * tobres2rphi, * tobres3rphi, * tobres4rphi, * tobres5rphi, * tobres6rphi, * tobres1sas, * tobres2sas;
  TH1F * tobclu1rphi, * tobclu2rphi, * tobclu3rphi, * tobclu4rphi, * tobclu5rphi, * tobclu6rphi, * tobclu1sas, * tobclu2sas;

  unsigned int tidlayer;
  unsigned int tidstereo;
  TH1F * tidres1rphi, * tidres2rphi, * tidres3rphi, * tidres1sas, * tidres2sas;
  TH1F * tidclu1rphi, * tidclu2rphi, * tidclu3rphi, * tidclu1sas, * tidclu2sas;

  unsigned int teclayer;
  unsigned int tecstereo;
  TH1F * tecres1rphi, * tecres2rphi, * tecres3rphi, * tecres4rphi,* tecres5rphi,* tecres6rphi,* tecres7rphi, * tecres1sas, * tecres2sas, * tecres5sas;
  TH1F * tecclu1rphi, * tecclu2rphi, * tecclu3rphi, * tecclu4rphi,* tecclu5rphi,* tecclu6rphi,* tecclu7rphi, * tecclu1sas, * tecclu2sas, * tecclu5sas;

};

}

#endif
