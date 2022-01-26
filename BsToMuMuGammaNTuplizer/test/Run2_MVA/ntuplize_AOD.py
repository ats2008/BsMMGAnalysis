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

process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v15_L1v1')
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
   'root://se01.indiacms.res.in//store/user/athachay/BsToMuMuGamma/Data/RunIISummer20UL18/106X_upgrade2018_realistic_v11_L1v1/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/crab_bs2mmg_RecoRun2UL18_v2/211227_041625/0000/BPH-RunIISummer20UL18RECO-00125_1.root',
   #'root://se01.indiacms.res.in//store/mc/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/30000/0D69C40E-AE25-9941-8F01-B43BBB86FF4D.root',
  #  'file:/eos/user/a/athachay/workarea/data/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
   #'file:/eos/user/a/athachay/workarea/data/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
 #  'file:/afs/cern.ch/work/a/athachay/public/BsToMuMuGamma/RunIIAutumn18DRPremix/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-evtgen-pythia8/AODSIM/102X_upgrade2018_realistic_v15-v1/606765BD-9BB4-9741-925C-A0C69B933039.root',      
   #'file:DoublePhotonGun/DoublePhoton0To40FlatPtAODSIM_HI_Reco_1.root',      
   #'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/Run2_analysis/CMSSW_10_6_20/src/BsMMGAnalysis/PhotonAnalyzer/test/AEFD418A-0A8D-414C-A8AC-86EE20287BDF.root'
   #'file:/eos/cms/store/group/phys_heavyions/rchudasa/lowPT_photonReco/SinglePhotonFlatPt1To20_GENSIM_Run2018/crab_SinglePhotonFlatPt1To20_default_RecoAOD_Run2018/210504_050451/0000/AODSIM_3.root'      
   #'root://se01.indiacms.res.in//store/user/rchudasa/BsMMG_2018UL/DoubleElectronFlatPt1To20_GENSIM_Run2018/crab_DoubleElectronFlatPt1To20_RecoAOD_Run2018/210816_180704/0000/AODSIM_1.root'
   #'file:EmptyBx_pp_2018_0B87F849-1E7B-A54C-AE8D-17C6D4549317.root'
   #'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/Run2_analysis/CMSSW_10_6_20/src/mcProduction/AODSIM_pion.root'
   #'/store/group/phys_heavyions/rchudasa/lowPT_photonReco/DoublePi0FlatPt1To20_GENSIM_Run2018_correctEta/crab_DoublePi0FlatPt1To20_RecoAOD_Run2018_correctEta/210908_093218/0000/AODSIM_pion_1.root'
   )
)
process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("muonNtuplizer.root")
    #fileName = cms.string("flatPtPhoton_ntuple_run2MVA.root")
    #fileName = cms.string("flatPtElectron_ntuple_pfIso.root")
    #fileName = cms.string("flatPtPi0_ntuple_pfIso.root")
    #fileName = cms.string("emptyBX_ntuple_pfIso.root")
    #fileName = cms.string("qcd_EmEnriched_ntuples.root")
    fileName = cms.string("bs2MuMuGamma_ntuples.root")
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
	#idName             = cms.string("cutBasedElectronID-Fall17-94X-V2-loose"),
	#idName             = cms.string("mvaPhoID-RunIIFall17-v2-wp90"),
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
    	doFlatPt           = cms.bool(True),
   	doMuons            = cms.bool(False),
   	doPhotons          = cms.bool(True),
    	doPFPhotons        = cms.bool(True),
	PFPhoton_minPt     = cms.untracked.double(0.0),	
    	doSuperClusters    = cms.bool(True),
    	Run2_2018          = cms.bool(True),
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
        doCompression                   = cms.bool(True),  #do the compression of floats
        nBits                           = cms.int32(23),   #nbits for float compression (<=23)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD

switchOnVIDPhotonIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.ntupler = cms.EDAnalyzer(
    'PhotonNtuplerVIDwithMVADemo',
    # The module automatically detects AOD vs miniAOD, so we configure both
    #
    # Common to all formats objects
    #
    # ... none ...
    #
    # Objects specific to AOD format
    #
    photons = cms.InputTag("gedPhotons"),
    genParticles = cms.InputTag("genParticles"),
    #
    # Objects specific to MiniAOD format
    #
    photonsMiniAOD = cms.InputTag("slimmedPhotons"),
    genParticlesMiniAOD = cms.InputTag("prunedGenParticles"),
    #
    # ID decisions (common to all formats)
    #
    # (the names of the ValueMaps for just decision and full info are the same,
    # they are distinguished by the type of the info)
    phoTightIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80"),
    phoTightIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80"),
    phoMediumIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
    phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
    # This is a fairly verbose mode if switched on, with full cut flow 
    # diagnostics for each candidate. Use it in a low event count test job.
    phoIdVerbose = cms.bool(False),
    #
    # ValueMaps with MVA results
    #
    mvaValuesMap     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1p1Values"),
    mvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v1p1Categories")
    )

#process.p = cms.Path(process.decayfilter*process.Ntuples)
process.p = cms.Path(process.decayfilter*process.Ntuples+process.egmPhotonIDSequence*process.ntupler)
