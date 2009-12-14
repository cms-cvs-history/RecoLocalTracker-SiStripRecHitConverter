#ifndef CPE_h
#define CPE_h

/** \class CPE
 *
 * CPE is a analyzer which reads rechits
 *
 * \author C. Genta
 *
 *
 ************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoLocalTracker/SiStripRecHitConverter/interface/StripCPE.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "TFile.h"
#include "TProfile.h"
#include "TH2F.h"
#define NBINS 32

class CPE : public edm::EDAnalyzer
{
 public:
  
  CPE(const edm::ParameterSet& conf);
  
  virtual ~CPE();
  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

  float uProj(const StripGeomDetUnit * stripdet,LocalTrajectoryParameters ltp,StripCPE *stripcpe);
  std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
						const BoundPlane& plane);

  std::vector<SiStripRecHit2D*> getRecHitComponents(const TrackingRecHit* rechit);

  //  TFile *myfile_;
  TH1F *tibproj, *tobproj, *tibres[NBINS], *tobres[NBINS], *tiberr, *toberr;
  TH1F *tidres[NBINS], *tecres[NBINS];
  TH1F *tibres1[NBINS], *tobres1[NBINS], *tidres1[NBINS], *tecres1[NBINS];
  TH1F *tibres2[NBINS], *tobres2[NBINS], *tidres2[NBINS], *tecres2[NBINS];
  TH1F *tibres3[NBINS], *tobres3[NBINS], *tidres3[NBINS], *tecres3[NBINS];
  TH1F *tibres4[NBINS], *tobres4[NBINS], *tidres4[NBINS], *tecres4[NBINS];
  TH1F *tibres5[NBINS], *tobres5[NBINS], *tidres5[NBINS], *tecres5[NBINS];
  TH1F *tibres6[NBINS], *tobres6[NBINS], *tidres6[NBINS], *tecres6[NBINS];
  TH1F *tibreswc[NBINS], *tobreswc[NBINS],*tidreswc[NBINS], *tecreswc[NBINS];
  TH1F *tibtkproj, *tobtkproj, *tibtkres, *tobtkres, *tibtkerr, *tobtkerr;


 private:
  edm::ParameterSet conf_;
};


#endif
