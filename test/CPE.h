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
#define NBINS 32

class CPE : public edm::EDAnalyzer
{
 public:
  
  CPE(const edm::ParameterSet& conf);
  
  virtual ~CPE();
  virtual void beginJob(const edm::EventSetup& c);
  virtual void endJob();
  virtual void analyze(const edm::Event& e, const edm::EventSetup& c);

  float uProj(const StripGeomDetUnit * stripdet,LocalTrajectoryParameters ltp,StripCPE *stripcpe);
  std::pair<LocalPoint,LocalVector> projectHit( const PSimHit& hit, const StripGeomDetUnit* stripDet,
						const BoundPlane& plane);

  std::vector<SiStripRecHit2D*> getRecHitComponents(const TrackingRecHit* rechit);
  TProfile *tibreshisto, *tobreshisto, *tiberrhisto, *toberrhisto;
  TProfile *tidreshisto, *tecreshisto, *tiderrhisto, *tecerrhisto;
  TProfile *mtibreshistox, *mtobreshistox, *mtidreshistox, *mtecreshistox;
  TProfile *mtibreshistoy, *mtobreshistoy, *mtidreshistoy, *mtecreshistoy;
  TProfile *tibsqrterrhisto, *tobsqrterrhisto, *tidsqrterrhisto, *tecsqrterrhisto;
  TProfile *mtibsqrterrhistox, *mtobsqrterrhistox, *mtidsqrterrhistox, *mtecsqrterrhistox;
  TProfile *mtibsqrterrhistoy, *mtobsqrterrhistoy, *mtidsqrterrhistoy, *mtecsqrterrhistoy;
  TProfile *mtibtksqrterrhistox, *mtobtksqrterrhistox, *mtidtksqrterrhistox, *mtectksqrterrhistox;
  TProfile *mtibtksqrterrhistoy, *mtobtksqrterrhistoy, *mtidtksqrterrhistoy, *mtectksqrterrhistoy;
  TProfile *tibtkreshisto, *tobtkreshisto, *tibtkerrhisto, *tobtkerrhisto;
  TFile *myfile_;
  TH1F *tibproj, *tobproj, *tibres[NBINS], *tobres[NBINS], *tiberr, *toberr;
  TH1F *tidres[NBINS], *tecres[NBINS];
  TH1F *mtibproj, *mtobproj, *mtibres[NBINS], *mtobres[NBINS];
  TH1F *mtibresx, *mtobresx, *mtidresx, *mtecresx;
  TH1F *mtibresy, *mtobresy, *mtidresy, *mtecresy;
  TH1F *mtidres[NBINS], *mtecres[NBINS];
  TH1F *mtidres_y[NBINS], *mtecres_y[NBINS], *mtibres_y[NBINS], *mtobres_y[NBINS];
  TH1F *tibtkproj, *tobtkproj, *tibtkres, *tobtkres, *tibtkerr, *tobtkerr;
 private:
  edm::ParameterSet conf_;
};


#endif
