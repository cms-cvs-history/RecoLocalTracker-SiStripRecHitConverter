
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoLocalTracker/SiStripRecHitConverter/test/ReadRecHit.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/ClusterFilter.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/CPE.h"

using cms::ReadRecHit;


DEFINE_FWK_MODULE(ReadRecHit);
DEFINE_FWK_MODULE(ClusterFilter);
DEFINE_FWK_MODULE(CPE);


