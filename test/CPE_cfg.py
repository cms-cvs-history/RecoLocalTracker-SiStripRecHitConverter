import FWCore.ParameterSet.Config as cms

process = cms.Process("Clusterizer")
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("RecoLocalTracker.SiStripRecHitConverter.test.RelValTTbar_cff")

process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.MagneticField_38T_cff")

process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load("Configuration.StandardSequences.Simulation_cff")

process.load("Configuration.StandardSequences.MixingNoPileUp_cff")

process.load("Configuration.StandardSequences.FakeConditions_cff")

process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff")

process.load("RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi")

process.load("Configuration.EventContent.EventContent_cff")

process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")

process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/relval/2008/6/22/RelVal-RelValQCD_Pt_50_80-1213987236-IDEAL_V2-2nd/0003/2613CED1-E340-DD11-98E7-000423D9880C.root')
)

process.CPEtest = cms.EDFilter("CPE",
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    associatePixel = cms.bool(False),
    ROUList = cms.vstring('g4SimHitsTrackerHitsTIBLowTof', 
        'g4SimHitsTrackerHitsTIBHighTof', 
        'g4SimHitsTrackerHitsTIDLowTof', 
        'g4SimHitsTrackerHitsTIDHighTof', 
        'g4SimHitsTrackerHitsTOBLowTof', 
        'g4SimHitsTrackerHitsTOBHighTof', 
        'g4SimHitsTrackerHitsTECLowTof', 
        'g4SimHitsTrackerHitsTECHighTof'),
    associateRecoTracks = cms.bool(False),
    associateStrip = cms.bool(True),
    VerbosityLevel = cms.untracked.int32(1),
    TkRecHits = cms.InputTag("TrackRefitter"),
    RecHitProducer = cms.string('siStripMatchedRecHits')
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('resolution_ttbar.root')
)

process.FEVT = cms.OutputModule("PoolOutputModule",
    process.FEVTSIMEventContent,
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RECO')
    ),
    fileName = cms.untracked.string('QCD.root')
)

process.p1 = cms.Path(process.mix*process.siStripMatchedRecHits*process.TrackRefitter*process.CPEtest)
process.TrackRefitter.TrajectoryInEvent = True


