import FWCore.ParameterSet.Config as cms

siStripMatchedRecHits = cms.EDFilter("SiStripRecHitConverter",
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    Regional = cms.bool(False),
    stereoRecHits = cms.string('stereoRecHit'),
    Matcher = cms.string('StandardMatcher'),
    matchedRecHits = cms.string('matchedRecHit'),
    maximumHits2BeforeMatching = cms.uint32(1000), # skip modules where #hits(mono)*#hits(stereo) > 1000 (sanity check to avoid crashing tier0)
    # next label (LazyGetterProducer) is only used if Regional is true
    LazyGetterProducer = cms.string('SiStripRawToClustersFacility'),
    ClusterProducer = cms.string('siStripClusters'),
    VerbosityLevel = cms.untracked.int32(1),
    rphiRecHits = cms.string('rphiRecHit'),
    useSiStripQuality = cms.bool(False),
    siStripQualityLabel = cms.string(''),
    MaskBadAPVFibers = cms.bool(False),
)


