
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoLocalTracker/SiStripRecHitConverter/test/ReadRecHit.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/ClusterFilter.h"
#include "RecoLocalTracker/SiStripRecHitConverter/test/ValHit.h"

using cms::ReadRecHit;
using cms::ValHit;

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(ValHit);
DEFINE_ANOTHER_FWK_MODULE(ReadRecHit);
DEFINE_ANOTHER_FWK_MODULE(ClusterFilter);


