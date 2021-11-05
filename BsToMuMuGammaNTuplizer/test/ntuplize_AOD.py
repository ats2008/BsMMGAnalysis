import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Configuration.Geometry.GeometryECALHCAL_cff")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('PhysicsTools.HepMCCandAlgos.genParticles_cfi')
process.load("Geometry.HcalEventSetup.CaloTowerTopology_cfi")
process.load("Configuration.Geometry.GeometryExtended2017_cff")
process.load("Configuration.Geometry.GeometryExtended2017Reco_cff")
process.load("RecoJets.Configuration.CaloTowersES_cfi")
process.load("Geometry.HcalEventSetup.hcalTopologyIdeal_cfi")

#process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
#process.GlobalTag.globaltag = cms.string('102X_upgrade2018_realistic_v15')

process.MessageLogger.cerr.FwkReport.reportEvery = 500

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.options = cms.untracked.PSet( numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True),

   )



process.source = cms.Source("PoolSource",
     duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
    fileNames = cms.untracked.vstring(
   'file:/eos/home-a/athachay/workarea/data/BsToMuMuGamma/eventModeling/MKT_ACCEP/bmmg_MKT_ACCEP_0.root'
   )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("bmmgAccepMKT.root")
)


process.decayfilter = cms.EDFilter("GenDecayKineFilter",
    SimGenParticle = cms.InputTag("genParticles"),
    DaughterIDs = cms.untracked.vint32(13, -13, 22),
    MaxEta = cms.untracked.vdouble(2.5, 2.5, 9999.0),
    MinEta = cms.untracked.vdouble(-2.5, -2.5, -9999.0),
    MinPt = cms.untracked.vdouble(3.5, 3.5, -99.0),
    NumberDaughters = cms.untracked.int32(3),
    ParticleID = cms.untracked.int32(531),
    verbose = cms.untracked.int32(20)
)

process.Ntuples = cms.EDAnalyzer("BsToMuMuGammaNTuplizer",
	muons   =cms.InputTag("muons"),
	beamSpot = cms.InputTag("offlineBeamSpot"),
	vertices = cms.InputTag("offlinePrimaryVertices"),
    	MustacheSCBarrelSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel"),
    	MustacheSCEndcapSrc= cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower"),
	GsfElectronSrc     = cms.InputTag("gedGsfElectrons"),
	muon_EtaMax      	= cms.untracked.double(1e3),	 
	muon_dcaMAX 		= cms.untracked.double(1e3),	
	muon_minPt  		= cms.untracked.double(1.0),	
	muon_zIPMax 		= cms.untracked.double(1e4),	
	muon_rIPMax 		= cms.untracked.double(1e4),	
	dimuon_minPt		= cms.untracked.double(0.0),
	dimuon_minInvMass 	= cms.untracked.double(-1e3),	
	dimuon_maxInvMass 	= cms.untracked.double(1e5),	
	dimuon_minVtxCL   	= cms.untracked.double(0.0),	
	dimuon_maxLStoBS  	= cms.untracked.double(1e5),	
	dimuon_maxDCAMuMu 	= cms.untracked.double(1e5),	
	dimuon_maxCosAlphaToBS 	= cms.untracked.double(1e5),	
        doHLT              = cms.bool(False),
    	doGenParticles     = cms.bool(True),
    	doFlatPt           = cms.bool(False),
   	doMuons            = cms.bool(False),
   	doPhotons          = cms.bool(False),
    	doPFPhotons        = cms.bool(False),
	PFPhoton_minPt     = cms.untracked.double(0.0),	
    	doSuperClusters    = cms.bool(False),
    	Run2_2018          = cms.bool(False),
    	genParticles       = cms.InputTag("genParticles"),
    	gedPhotonSrc       = cms.untracked.InputTag("gedPhotons"),
    	pfPhotonSrc        = cms.untracked.InputTag("particleFlow"),
	TriggerNames = cms.vstring("HLT_DoubleMu4_3_Bs_v14",
	                           "HLT_DoubleMu4_3_Jpsi_v2",
				   "HLT_DoubleMu4_JpsiTrk_Displaced_v15",
				   "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v15",
				   "HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v8"),
	HLTResult = cms.InputTag("TriggerResults","","HLT"),
	verbose  = cms.bool(False),
	doBsToMuMuGamma = cms.bool(False),
	isMC = cms.bool(True),
    	#caloParticleCollection          = cms.InputTag("mix","MergedCaloTruth"),
 	hbheRechitCollection            = cms.InputTag("reducedHcalRecHits","hbhereco","RECO"),
 	ebRechitCollection              = cms.InputTag("reducedEcalRecHitsEB","","RECO"),
   	eeRechitCollection              = cms.InputTag("reducedEcalRecHitsEE","","RECO"),
    	pfRechitCollection              = cms.InputTag("particleFlowRecHitECAL","","RECO"),
    	#pfClusterCollection            = cms.InputTag("particleFlowClusterECAL","","RECO"),
        doCompression                   = cms.bool(False),  #do the compression of floats
        nBits                           = cms.int32(23),   #nbits for float compression (<=23)
)


#process.p = cms.Path(process.decayfilter*process.Ntuples)
process.p = cms.Path(process.Ntuples)
