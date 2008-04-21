
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoLocalTracker/SiStripRecHitConverter/test/ReadRecHit.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/ClusterFilter.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/CPE.h"

using cms::ReadRecHit;

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(ReadRecHit);
DEFINE_ANOTHER_FWK_MODULE(ClusterFilter);
DEFINE_ANOTHER_FWK_MODULE(CPE);


