import FWCore.ParameterSet.Config as cms

def customizeBsToMMGTrigPathForTrigBasedFileGen(process):

    process.load('Configuration.EventContent.EventContent_cff')

    process.FEVTEventContent.outputCommands.append('drop *_*_*_SIM')
    process.FEVTEventContent.outputCommands.append('drop *_*_*_HLT')

    process.hltDoubleMu4BsToMMGL3Filtered.MaxDr = cms.vdouble( 20.0 )
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMin = cms.vdouble( 3.5 )
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMin = cms.vdouble( 3.5 )
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMax = cms.vdouble( 3.5 )
    
    process.hltDoubleMu4BsToMMGL3Filtered.MinInvMass = cms.vdouble( 1.0 )
    process.hltDoubleMu4BsToMMGL3Filtered.MaxInvMass = cms.vdouble( 10.0 )

    process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
    process.source.skipBadFiles = cms.untracked.bool(True)

    process.hltEgammaL1MatchFilterForBsToMMG = cms.EDFilter( "HLTEgammaL1TMatchFilterRegionalDummy",
        doIsolated = cms.bool( False ),
        endcap_end = cms.double( 2.65 ),
        region_phi_size = cms.double( 1.044 ),
        saveTags = cms.bool( True ),
        region_eta_size_ecap = cms.double( 1.0 ),
        barrel_end = cms.double( 1.4791 ),
        l1IsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        candIsolatedTag = cms.InputTag( "hltEgammaCandidatesForBsToMMG" ),
        l1CenJetsTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
        region_eta_size = cms.double( 0.522 ),
        L1SeedFilterTag = cms.InputTag( "hltL1sSingleEG15er2p5" ),
        candNonIsolatedTag = cms.InputTag( "" ),
        l1NonIsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        ncandcut = cms.int32( 1 ),
        l1TausTag = cms.InputTag( 'hltGtStage2Digis','Tau' )
       )
    
    # Output definition : L1 Pass/Fail Events
 
 
    process.Bs2MMG_L1Pass = cms.Path( process.HLTBeginSequence + process.hltL1sDoubleMu0er1p5SQOSdRMax1p4IorDoubleMu0er1p4SQOSdRMax1p4 + process.HLTEndSequence )

    process.RAWoutputL1Pass = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_L1Pass.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('Bs2MMG_L1Pass')
       ),
    maxSize = cms.untracked.int32(4000000),
    )

    process.RAWoutputL1Fail = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_L1Fail.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('!Bs2MMG_L1Pass')
       ),
    maxSize = cms.untracked.int32(4000000),
    ) 


    # Output definition : L3 Dimuon Pass/Fail Events

    process.Bs2MMG_L3DimuonPass= cms.Path( process.HLTBeginSequence + process.hltL1sDoubleMu0er1p5SQOSdRMax1p4IorDoubleMu0er1p4SQOSdRMax1p4 + process.hltPreDoubleMu43BsToMMG + process.hltL1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0 + process.HLTL2muonrecoSequence + cms.ignore(process.hltL2fL1sL1DoubleMu0er1p5SQOSdR1p4L1f0L2PreFiltered0) + process.HLTL3muonrecoSequence + cms.ignore(process.hltL1fForIterL3L1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0) + process.hltDoubleMu4BsToMMGL3Filtered + process.HLTEndSequence )

    process.RAWoutputL3DimuonPass = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_L3DimuonPass.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('Bs2MMG_L3DimuonPass')
       ),
    maxSize = cms.untracked.int32(4000000),
    )

    process.RAWoutputL3DimuonFail = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_L3DimuonFail.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('!Bs2MMG_L3DimuonPass')
       ),
    maxSize = cms.untracked.int32(4000000),
    ) 

    # Output definition : H/E Pass/Fail Events

    process.Bs2MMG_ECALHoverEPass = cms.Path( process.HLTBeginSequence + process.hltL1sDoubleMu0er1p5SQOSdRMax1p4IorDoubleMu0er1p4SQOSdRMax1p4 + process.hltPreDoubleMu43BsToMMG + process.hltL1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0 + process.HLTL2muonrecoSequence + cms.ignore(process.hltL2fL1sL1DoubleMu0er1p5SQOSdR1p4L1f0L2PreFiltered0) + process.HLTL3muonrecoSequence + cms.ignore(process.hltL1fForIterL3L1fL1sL1DoubleMu0er1p5SQOSdR1p4L1Filtered0) + process.hltDoubleMu4BsToMMGL3Filtered + process.hltDisplacedmumuVtxProducerDoubleMu4BsToMMG + process.hltDisplacedmumuFilterDoubleMu4BsToMMG + process.HLTDoFullUnpackingEgammaEcalMFSequence + process.HLTPFClusteringForEgammaForBsToMMG + process.hltEgammaCandidatesForBsToMMG + cms.ignore(process.hltEgammaL1MatchFilterForBsToMMG)+ process.HLTDoLocalHcalSequence + process.HLTFastJetForEgamma + process.hltEgammaHoverEForBsToMMG + process.hltPhotonHoverEFilterForBsToMMG  +process.HLTEndSequence )

    process.RAWoutputHoverEPass = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_ECALhOverEPass.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('Bs2MMG_ECALHoverEPass')
       ),
    maxSize = cms.untracked.int32(4000000),
    )

    process.RAWoutputHoverEFail = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),  
    fileName = cms.untracked.string('file:BsToMMG_GENSIM-DIGI-RAW_pu30To80CM_ECALhOverEFail.root'),
    outputCommands = process.FEVTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('!Bs2MMG_ECALHoverEPass')
       ),
    maxSize = cms.untracked.int32(4000000),
    ) 

    
    process.trigger_output_step = cms.EndPath(process.RAWoutputL1Pass+process.RAWoutputL1Fail+process.RAWoutputL1Fail+process.RAWoutputL3DimuonFail+process.RAWoutputL3DimuonPass+process.RAWoutputHoverEFail+process.RAWoutputHoverEPass)
    
    process.HLTSchedule+= cms.Schedule(*[process.trigger_output_step])
    

      
    
    return process


def customizeBsToMMGTrigPath(process):

    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMin = cms.vdouble( 3.5 )
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMax = cms.vdouble( 3.5 )
    
    process.hltDoubleMu4BsToMMGL3Filtered.MinInvMass = cms.vdouble( 1.0 )
    process.hltDoubleMu4BsToMMGL3Filtered.MaxInvMass = cms.vdouble( 6.0 )

    process.source.duplicateCheckMode = cms.untracked.string("noDuplicateCheck")

    process.hltEgammaL1MatchFilterForBsToMMG = cms.EDFilter( "HLTEgammaL1TMatchFilterRegionalDummy",
        doIsolated = cms.bool( False ),
        endcap_end = cms.double( 2.65 ),
        region_phi_size = cms.double( 1.044 ),
        saveTags = cms.bool( True ),
        region_eta_size_ecap = cms.double( 1.0 ),
        barrel_end = cms.double( 1.4791 ),
        l1IsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        candIsolatedTag = cms.InputTag( "hltEgammaCandidatesForBsToMMG" ),
        l1CenJetsTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
        region_eta_size = cms.double( 0.522 ),
        L1SeedFilterTag = cms.InputTag( "hltL1sSingleEG15er2p5" ),
        candNonIsolatedTag = cms.InputTag( "" ),
        l1NonIsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        ncandcut = cms.int32( 1 ),
        l1TausTag = cms.InputTag( 'hltGtStage2Digis','Tau' )
       )
    return process    


def customizeBsToMMGChageProducerModules(process):
    process.hltMuMuPhotonFilter = cms.EDFilter( "HLTPhotonMuonDR",
        saveTags = cms.bool( True ),
        electronTag = cms.InputTag( "hltEgammaGsfElectrons" ),
        originTag2 = cms.VInputTag( 'hltEgammaCandidatesForBsToMMG' ),
        MinPixHitsForDZ = cms.int32( 1 ),
        MinN = cms.int32( 2 ),
        originTag1 = cms.VInputTag( 'hltIterL3MuonCandidates' ),
        triggerType1 = cms.int32( 83 ),
        triggerType2 = cms.int32( 81 ),
        MinDR = cms.double( -1.0 ),
        MaxDZ = cms.double( 40.0 ),
        inputTag1 = cms.InputTag( "hltDisplacedmumuFilterDoubleMu4BsToMMG" ),
        checkSC = cms.bool( False ),
        inputTag2 = cms.InputTag( "hltPhotonHoverEFilterForBsToMMG" )
        )
    
    process.hltEgammaL1MatchFilterForBsToMMG = cms.EDFilter( "HLTEgammaL1TMatchFilterRegionalDummy",
        doIsolated = cms.bool( False ),
        endcap_end = cms.double( 2.65 ),
        region_phi_size = cms.double( 1.044 ),
        saveTags = cms.bool( True ),
        region_eta_size_ecap = cms.double( 1.0 ),
        barrel_end = cms.double( 1.4791 ),
        l1IsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        candIsolatedTag = cms.InputTag( "hltEgammaCandidatesForBsToMMG" ),
        l1CenJetsTag = cms.InputTag( 'hltGtStage2Digis','Jet' ),
        region_eta_size = cms.double( 0.522 ),
        L1SeedFilterTag = cms.InputTag( "hltL1sSingleEG15er2p5" ),
        candNonIsolatedTag = cms.InputTag( "" ),
        l1NonIsolatedTag = cms.InputTag( 'hltGtStage2Digis','EGamma' ),
        ncandcut = cms.int32( 1 ),
        l1TausTag = cms.InputTag( 'hltGtStage2Digis','Tau' )
       )
    return process

def addHLTEgammaExtraForBsToMMG(process):
    process.hltEgammaHLTExtraForBsToMMG = cms.EDProducer("EgammaHLTExtraProducer",
                                   egCands = cms.VPSet(
                                       cms.PSet(
                                           pixelSeeds = cms.InputTag("hltEgammaElectronPixelSeedsUnseeded"),
                                           ecalCands = cms.InputTag("hltEgammaCandidatesForBsToMMG"),
                                           gsfTracks = cms.InputTag("hltEgammaGsfTracksUnseeded"),
                                           label = cms.string('ForBsToMMG')
                                       ),
                                   ),                 
                                   ecal = cms.VPSet(
                                       cms.PSet(
                                           src= cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
                                           label = cms.string("EcalRecHitsEB")
                                       ),
                                       cms.PSet(
                                           src= cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
                                           label = cms.string("EcalRecHitsEE")
                                       )
                                   ),
                                   pfClusIso = cms.VPSet(
                                       cms.PSet(
                                           src = cms.InputTag("hltParticleFlowClusterECALBsToMMGSeeded"),
                                           label = cms.string("Ecal")
                                       ),
                                       cms.PSet(
                                           src = cms.InputTag("hltParticleFlowClusterHCAL"),
                                           label = cms.string("Hcal")
                                       ),
                                   ),
                                   hcal = cms.VPSet(cms.PSet(src=cms.InputTag("hltHbhereco"),label=cms.string(""))),
                                   trks = cms.VPSet(cms.PSet(src=cms.InputTag("hltMergedTracks"),label=cms.string(""))),                                  
                                   minPtToSaveHits = cms.double(8.),
                                   saveHitsPlusHalfPi = cms.bool(True),
                                   saveHitsPlusPi = cms.bool(False)
                                   
                        )

    for outmodname in process.outputModules_():
        outmod = process.outputModules_()[outmodname]
        if outmod.type_()=='PoolOutputModule':
            for endpathname in process.endpaths_():
                endpath = process.endpaths_()[endpathname]
                if outmodname in endpath.moduleNames():
                    endpath.insert(0,process.hltEgammaHLTExtraForBsToMMG)
    return process


def customizeBsToMMGTrigPathForEfficiencyFiles(process):
    
    process= customizeBsToMMGChageProducerModules(process)

    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMin  = cms.vdouble(1.0)       
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtMax  = cms.vdouble(1.0)
    process.hltDoubleMu4BsToMMGL3Filtered.MaxInvMass= cms.vdouble(100.0)
    process.hltDoubleMu4BsToMMGL3Filtered.MaxInvMass= cms.vdouble(0.0)
    process.hltDoubleMu4BsToMMGL3Filtered.MinPtPair = cms.vdouble(1.0)
    process.hltDoubleMu4BsToMMGL3Filtered.MaxDr     = 100.0
    process.hltDisplacedmumuVtxProducerDoubleMu4BsToMMG.matchToPrevious = False
    process.hltRecHitInRegionForBsToMMGMuonsECAL.regionEtaMargin =  1.4
    process.hltRecHitInRegionForBsToMMGMuonsECAL.regionPhiMargin =  1.4
    process.hltRecHitInRegionForBsToMMGMuonsES.regionEtaMargin =  1.4
    process.hltRecHitInRegionForBsToMMGMuonsES.regionPhiMargin =  1.4

    process.RAWoutput = cms.OutputModule("PoolOutputModule",
        dataset = cms.untracked.PSet(
            dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
            filterName = cms.untracked.string('')
        ),  
        fileName = cms.untracked.string('file:BsToMMG_hltDev.root'),
        outputCommands =cms.untracked.vstring(
                'drop *',
                'keep *_genParticles_*_*',
                'keep *_hltGtStage2Digis_*_*',
                'keep *_hltEgammaHLTExtraForBsToMMG_*_*',
                'keep *_hltIterL3MuonCandidates_*_*',
                'keep *_hltDoubleMu4BsToMMGL3Filtered_*_*',
                'keep *_hltDisplacedmumuVtxProducerDoubleMu4BsToMMG_*_*' 
                ),
        splitLevel = cms.untracked.int32(0),
        SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('HLT_DoubleMu4_3_BsToMMG_v0')
           ),
        maxSize = cms.untracked.int32(4000000)   )
    
    process.trigger_output_step = cms.EndPath(process.RAWoutput)
    process.HLTSchedule+= cms.Schedule(*[process.trigger_output_step])
    
    process = addHLTEgammaExtraForBsToMMG(process)    

    return process


def customizeBsToMMGTrigPathForLowPtReco(process):
    process.hltParticleFlowSuperClusterECALBsToMMGSeeded.thresh_SCEt = 1.0
    process.hltParticleFlowSuperClusterECALBsToMMGSeeded.thresh_PFClusterSeedBarrel = 0.50
    process.hltParticleFlowSuperClusterECALBsToMMGSeeded.thresh_PFClusterSeedEndcap = 0.50

    return process
