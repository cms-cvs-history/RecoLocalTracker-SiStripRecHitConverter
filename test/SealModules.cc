
#include "PluginManager/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoLocalTracker/ReadRecHit/interface/ReadRecHit.h"

using cms::ReadRecHit;
DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(ReadRecHit)

