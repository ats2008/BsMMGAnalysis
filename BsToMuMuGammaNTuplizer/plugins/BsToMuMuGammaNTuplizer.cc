// -*- C++ -*-
//
// Package:    BsMMGAnalysis/BsToMuMuGammaNTuplizer
// Class:      BsToMuMuGammaNTuplizer
//
/**\class BsToMuMuGammaNTuplizer BsToMuMuGammaNTuplizer.cc BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author: 
//         Created: 
//
//


// system include files
#include <memory>
#include <map>
#include <string>
#include "TH1.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "EgammaAnalysis/ElectronTools/interface/SuperClusterHelper.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/Utils.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "Geometry/CaloTopology/interface/CaloTowerConstituentsMap.h"
#include <TTree.h>
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHadTower.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EGHcalRecHitSelector.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/pfIsoCalculator.h"

BsToMuMuGammaNTuplizer::BsToMuMuGammaNTuplizer(const edm::ParameterSet& iConfig) :

  doGenParticles_(iConfig.getParameter<bool>("doGenParticles")),
  doFlatPt_(iConfig.getParameter<bool>("doFlatPt")),
  doMuons_(iConfig.getParameter<bool>("doMuons")),
  doPhotons_(iConfig.getParameter<bool>("doPhotons")),
  doPFPhotons_(iConfig.getParameter<bool>("doPFPhotons")),
  doSuperClusters_(iConfig.getParameter<bool>("doSuperClusters")),
  doHLT(iConfig.getParameter<bool>("doHLT")),
  Run2_2018_(iConfig.getParameter<bool>("Run2_2018"))

{
  
  
  Utility= new Utils();
  
  if(doMuons_) muonToken_              = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  //caloPartToken_                 = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
  
  if(doPhotons_)    gedPhotonsCollection_       = consumes<std::vector<reco::Photon>>(iConfig.getUntrackedParameter<edm::InputTag>("gedPhotonSrc"));
  
    pfPhotonsCollection_        = consumes<edm::View<reco::PFCandidate>>(iConfig.getUntrackedParameter<edm::InputTag>("pfPhotonSrc"));
  
  if(doSuperClusters_){
    MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc"));
    MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"));
    gsfElectronToken_                       = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("GsfElectronSrc"));
    hbheRechitToken_               = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit>>>(iConfig.getParameter<edm::InputTag>("hbheRechitCollection"));
    ebRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
    eeRechitToken_                 = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
  }
  
  if(doHLT) {
    trigTable    =iConfig.getParameter<std::vector<std::string>>("TriggerNames");
    triggerBits_ =consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTResult"));
  }
  
  beamSpotToken_  	   =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  if(doGenParticles_)genParticlesCollection_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  
  etaMax_muon               =  iConfig.getUntrackedParameter<double>("muon_EtaMax")        ;
  dcaMax_muon_bs            =  iConfig.getUntrackedParameter<double>("muon_dcaMAX")        ;
  pTMinMuons 		      =  iConfig.getUntrackedParameter<double>("muon_minPt");
  pTMinPFPhotons               =  iConfig.getUntrackedParameter<double>("PFPhoton_minPt");
  trackIP_zMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_zIPMax")        ;
  trackIP_rMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_rIPMax")        ;
  
  minDimuon_pt              =  iConfig.getUntrackedParameter<double>("dimuon_minPt")      ;
  cl_dimuon_vtx             =  iConfig.getUntrackedParameter<double>("dimuon_minVtxCL")      ;
  ls_max_dimuonBS           =  iConfig.getUntrackedParameter<double>("dimuon_maxLStoBS")       ;
  dcaMax_dimuon_mumu        =  iConfig.getUntrackedParameter<double>("dimuon_maxDCAMuMu")        ;
  cosAlphaMax_dimuonBs      =  iConfig.getUntrackedParameter<double>("dimuon_maxCosAlphaToBS") ;
  minDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_minInvMass")    ;
  maxDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_maxInvMass")    ;
  
  printMsg=iConfig.getParameter<bool>("verbose");
  isMC=iConfig.getParameter<bool>("isMC");
  doBsToMuMuGamma=iConfig.getParameter<bool>("doBsToMuMuGamma");
  nBits_                         = iConfig.getParameter<int>("nBits"); 
  doCompression_                 = iConfig.getParameter<bool>("doCompression");  
  // initialize output TTree
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("EventTree", "Event data");
  
  theTree->Branch("run",    &run_);
  theTree->Branch("event",  &event_);
  theTree->Branch("lumis",  &lumis_);
  theTree->Branch("isData", &isData_);
 
  theTree->Branch("beamspot_x",         		 &beamspot_x_);
  theTree->Branch("beamspot_y",       	    	 &beamspot_y_);
  theTree->Branch("beamspot_z",       		     &beamspot_z_);
  theTree->Branch("beamspot_x_error",       	 &beamspot_x_error_);
  theTree->Branch("beamspot_y_error",       	 &beamspot_y_error_);
  theTree->Branch("beamspot_z_error",        	 &beamspot_z_error_);
  theTree->Branch("beamspot_covXX",         	 &beamspot_covXX);
  theTree->Branch("beamspot_covXY",         	 &beamspot_covXY);
  theTree->Branch("beamspot_covXZ",         	 &beamspot_covXZ);
  theTree->Branch("beamspot_covYY",         	 &beamspot_covYY);
  theTree->Branch("beamspot_covYZ",         	 &beamspot_covYZ);
  theTree->Branch("beamspot_covZZ",         	 &beamspot_covZZ);

  theTree->Branch("beamspot_dxdz",               &beamspot_dxdz_);
  theTree->Branch("beamspot_dydz",               &beamspot_dydz_);
  theTree->Branch("beamspot_sigmaZ",             &beamspot_sigmaZ_);
  theTree->Branch("beamspot_dxdz_error",         &beamspot_dxdz_error_);
  theTree->Branch("beamspot_dydz_error",         &beamspot_dydz_error_);
  theTree->Branch("beamspot_sigmaZError",        &beamspot_sigmaZError_);
  theTree->Branch("beamspot_beamWidthX",         &beamspot_beamWidthX_);
  theTree->Branch("beamspot_beamWidthY",         &beamspot_beamWidthY_);
  theTree->Branch("beamspot_beamWidthX_error",   &beamspot_beamWidthX_error_);
  theTree->Branch("beamspot_beamWidthY_error",   &beamspot_beamWidthY_error_);

  theTree->Branch("nPrimaryVertex",            &nPrimaryVertex_);
  theTree->Branch("primaryVertex_isFake",      &primaryVertex_isFake_);
  theTree->Branch("primaryVertex_x",           &primaryVertex_x_);
  theTree->Branch("primaryVertex_y",           &primaryVertex_y_);
  theTree->Branch("primaryVertex_z",           &primaryVertex_z_);
  theTree->Branch("primaryVertex_t",           &primaryVertex_t_);
  theTree->Branch("primaryVertex_covXX",       &primaryVertex_covXX);
  theTree->Branch("primaryVertex_covXY",       &primaryVertex_covXY);
  theTree->Branch("primaryVertex_covXZ",       &primaryVertex_covXZ);
  theTree->Branch("primaryVertex_covYY",       &primaryVertex_covYY);
  theTree->Branch("primaryVertex_covYZ",       &primaryVertex_covYZ);
  theTree->Branch("primaryVertex_covZZ",       &primaryVertex_covZZ);
  theTree->Branch("primaryVertex_x_error",     &primaryVertex_x_error_);
  theTree->Branch("primaryVertex_y_error",     &primaryVertex_y_error_);
  theTree->Branch("primaryVertex_z_error",     &primaryVertex_z_error_);
  theTree->Branch("primaryVertex_t_error",     &primaryVertex_t_error_);
  theTree->Branch("primaryVertex_ntracks",     &primaryVertex_ntracks_);
  theTree->Branch("primaryVertex_ndof",        &primaryVertex_ndof_);
  theTree->Branch("primaryVertex_chi2",        &primaryVertex_chi2_); 
  theTree->Branch("primaryVertex_normalizedChi2", &primaryVertex_normalizedChi2_);

  if (doHLT) {
    // ### Trigger ###
    //theTree->Branch("trigTable",     &TrigTable);
    TrigTable_store=nullptr;
    TrigResult_store=nullptr;
    TrigPrescales_store=nullptr;
    theTree->Branch("trigResult",    &trigResult);
    theTree->Branch("trigPrescales", &trigPrescales);
    theTree->Branch("l1Table",       &l1Table);
    theTree->Branch("l1Prescales",   &l1Prescales);
  
    SetupTriggerStorageVectors();
    SetupTriggerBranches();

  }


  if (doGenParticles_) {
    theTree->Branch("gen_nBs"			,&gen_nBs_);
    theTree->Branch("gen_Bs_pt"			,&gen_Bs_pt_);
    theTree->Branch("gen_Bs_energy"		,&gen_Bs_energy_);
    theTree->Branch("gen_Bs_eta"		,&gen_Bs_eta_);
    theTree->Branch("gen_Bs_phi"		,&gen_Bs_phi_);
    theTree->Branch("gen_Bs_pz"			,&gen_Bs_pz_);
    theTree->Branch("gen_Bs_pdgId"		,&gen_Bs_pdgId_);

    theTree->Branch("gen_nBsMuonM"		,&gen_nBsMuonM_);
    theTree->Branch("gen_BsMuonM_pt"		,&gen_BsMuonM_pt_);
    theTree->Branch("gen_BsMuonM_eta"		,&gen_BsMuonM_eta_);
    theTree->Branch("gen_BsMuonM_phi"		,&gen_BsMuonM_phi_);

    theTree->Branch("gen_nBsMuonP"		,&gen_nBsMuonP_);
    theTree->Branch("gen_BsMuonP_pt"		,&gen_BsMuonP_pt_);
    theTree->Branch("gen_BsMuonP_eta"		,&gen_BsMuonP_eta_);
    theTree->Branch("gen_BsMuonP_phi"		,&gen_BsMuonP_phi_);

    theTree->Branch("gen_nBsPhoton"		,&gen_nBsPhoton_);
    theTree->Branch("gen_BsPhoton_pt"		,&gen_BsPhoton_pt_);
    theTree->Branch("gen_BsPhoton_energy"	,&gen_BsPhoton_energy_);
    theTree->Branch("gen_BsPhoton_eta"		,&gen_BsPhoton_eta_);
    theTree->Branch("gen_BsPhoton_phi"		,&gen_BsPhoton_phi_);
 
    if(doFlatPt_){
    theTree->Branch("nMC",          &nMC_);
    theTree->Branch("mcPID",        &mcPID_);
    theTree->Branch("mcStatus",     &mcStatus_);
    theTree->Branch("mcVtx_x",      &mcVtx_x_);
    theTree->Branch("mcVtx_y",      &mcVtx_y_);
    theTree->Branch("mcVtx_z",      &mcVtx_z_);
    theTree->Branch("mcPt",         &mcPt_);
    theTree->Branch("mcEta",        &mcEta_);
    theTree->Branch("mcPhi",        &mcPhi_);
    theTree->Branch("mcE",          &mcE_);
    theTree->Branch("mcEt",         &mcEt_);
    theTree->Branch("mcMass",       &mcMass_);
   }
  }

  if(doMuons_){

    // ### mu- ###  
    theTree->Branch("nMuM",             &nMuM_);
    theTree->Branch("mumHighPurity",    &mumHighPurity_);
    theTree->Branch("mumPt",            &mumPt_);
    theTree->Branch("mumEta",           &mumEta_);
    theTree->Branch("mumPhi",           &mumPhi_);
    theTree->Branch("mumCL",            &mumCL_);
    theTree->Branch("mumNormChi2",      &mumNormChi2_);
    theTree->Branch("mumVx",            &mumVx_);
    theTree->Branch("mumVy",            &mumVy_);
    theTree->Branch("mumVz",            &mumVz_);
    theTree->Branch("mumDCABS",         &mumDCABS_);
    theTree->Branch("mumDCABSE",        &mumDCABSE_);
    theTree->Branch("mumFracHits",      &mumFracHits_);
    theTree->Branch("mumdxyBS",         &mumdxyBS_);
    theTree->Branch("mumdzBS",          &mumdzBS_);
    theTree->Branch("mumIdx",           &mumIdx_);
    theTree->Branch("mumNPixHits",      &mumNPixHits_);
    theTree->Branch("mumNPixLayers",    &mumNPixLayers_);
    theTree->Branch("mumNTrkHits",      &mumNTrkHits_);
    theTree->Branch("mumNTrkLayers",    &mumNTrkLayers_);
    theTree->Branch("mumNMuonHits",     &mumNMuonHits_);
    theTree->Branch("mumNMatchStation", &mumNMatchStation_);
    theTree->Branch("mum_isGlobalMuon", &mum_isGlobalMuon_);
    theTree->Branch("mum_isTrackerMuon",&mum_isTrackerMuon_);
    theTree->Branch("mum_StandAloneMuon",&mum_StandAloneMuon_);
    theTree->Branch("mum_isCaloMuon",    &mum_isCaloMuon_);
    theTree->Branch("mum_isPFMuon",      &mum_isPFMuon_);

    theTree->Branch("mum_selector"          	,  &mum_selector_); 
    theTree->Branch("mum_isIsolationValid"  	,  &mum_isIsolationValid_);
    theTree->Branch("mum_isPFIsolationValid"	,  &mum_isPFIsolationValid_);
   
    theTree->Branch("mum_isolationR03_trackSumPt"			,&mum_isolationR03_trackSumPt_);
    theTree->Branch("mum_isolationR03_trackEcalSumEt"			,&mum_isolationR03_trackEcalSumEt_);
    theTree->Branch("mum_isolationR03_trackHcalSumEt"			,&mum_isolationR03_trackHcalSumEt_);
    theTree->Branch("mum_isolationR03_trackHoSumEt"			,&mum_isolationR03_trackHoSumEt_);
    theTree->Branch("mum_isolationR03_trackNTracks"			,&mum_isolationR03_trackNTracks_);
    theTree->Branch("mum_isolationR03_trackNJets"			,&mum_isolationR03_trackNJets_);
    theTree->Branch("mum_isolationR03_trackerVetoSumPt"			,&mum_isolationR03_trackerVetoSumPt_);
    theTree->Branch("mum_isolationR03_emVetoSumEt"			,&mum_isolationR03_emVetoSumEt_);
    theTree->Branch("mum_isolationR03_hadVetoSumEt"			,&mum_isolationR03_hadVetoSumEt_);
    theTree->Branch("mum_isolationR03_hoVetoEt"				,&mum_isolationR03_hoVetoEt_);

    theTree->Branch("mum_isolationR05_trackSumPt"			,&mum_isolationR05_trackSumPt_);
    theTree->Branch("mum_isolationR05_trackEcalSumEt"			,&mum_isolationR05_trackEcalSumEt_);
    theTree->Branch("mum_isolationR05_trackHcalSumEt"			,&mum_isolationR05_trackHcalSumEt_);
    theTree->Branch("mum_isolationR05_trackHoSumEt"			,&mum_isolationR05_trackHoSumEt_);
    theTree->Branch("mum_isolationR05_trackNTracks"			,&mum_isolationR05_trackNTracks_);
    theTree->Branch("mum_isolationR05_trackNJets"			,&mum_isolationR05_trackNJets_);
    theTree->Branch("mum_isolationR05_trackerVetoSumPt"			,&mum_isolationR05_trackerVetoSumPt_);
    theTree->Branch("mum_isolationR05_emVetoSumEt"			,&mum_isolationR05_emVetoSumEt_);
    theTree->Branch("mum_isolationR05_hadVetoSumEt"			,&mum_isolationR05_hadVetoSumEt_);
    theTree->Branch("mum_isolationR05_hoVetoEt"				,&mum_isolationR05_hoVetoEt_);
                                                                                                                                                                                                                                                          
    theTree->Branch("mum_PFIsolationR03_sumChargedHadronPt"		,&mum_PFIsolationR03_sumChargedHadronPt_);
    theTree->Branch("mum_PFIsolationR03_sumChargedParticlePt"		,&mum_PFIsolationR03_sumChargedParticlePt_);
    theTree->Branch("mum_PFIsolationR03_sumNeutralHadronEt"		,&mum_PFIsolationR03_sumNeutralHadronEt_);
    theTree->Branch("mum_PFIsolationR03_sumPhotonEt"			,&mum_PFIsolationR03_sumPhotonEt_);
    theTree->Branch("mum_PFIsolationR03_sumNeutralHadronEtHighThreshold",&mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR03_sumPhotonEtHighThreshold"	,&mum_PFIsolationR03_sumPhotonEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR03_sumPUPt"			        ,&mum_PFIsolationR03_sumPUPt_);
                                                                                                                             
    theTree->Branch("mum_PFIsolationR04_sumChargedHadronPt"		,&mum_PFIsolationR04_sumChargedHadronPt_);
    theTree->Branch("mum_PFIsolationR04_sumChargedParticlePt"		,&mum_PFIsolationR04_sumChargedParticlePt_);
    theTree->Branch("mum_PFIsolationR04_sumNeutralHadronEt"		,&mum_PFIsolationR04_sumNeutralHadronEt_);
    theTree->Branch("mum_PFIsolationR04_sumPhotonEt"			,&mum_PFIsolationR04_sumPhotonEt_);
    theTree->Branch("mum_PFIsolationR04_sumNeutralHadronEtHighThreshold",&mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR04_sumPhotonEtHighThreshold"	,&mum_PFIsolationR04_sumPhotonEtHighThreshold_);
    theTree->Branch("mum_PFIsolationR04_sumPUPt"			        ,&mum_PFIsolationR04_sumPUPt_);

    // ### mu+ ###  
    theTree->Branch("nMuP",             &nMuP_);
    theTree->Branch("mupHighPurity",    &mupHighPurity_);
    theTree->Branch("mupPt",            &mupPt_);
    theTree->Branch("mupEta",           &mupEta_);
    theTree->Branch("mupPhi",           &mupPhi_);
    theTree->Branch("mupCL",            &mupCL_);
    theTree->Branch("mupNormChi2",      &mupNormChi2_);
    theTree->Branch("mupVx",            &mupVx_);
    theTree->Branch("mupVy",            &mupVy_);
    theTree->Branch("mupVz",            &mupVz_);
    theTree->Branch("mupDCABS",         &mupDCABS_);
    theTree->Branch("mupDCABSE",        &mupDCABSE_);
    theTree->Branch("mupFracHits",      &mupFracHits_);
    theTree->Branch("mupdxyBS",         &mupdxyBS_);
    theTree->Branch("mupdzBS",          &mupdzBS_);
    theTree->Branch("mupIdx_",          &mupIdx_);
    theTree->Branch("mupNPixHits",      &mupNPixHits_);
    theTree->Branch("mupNPixLayers",    &mupNPixLayers_);
    theTree->Branch("mupNTrkHits",      &mupNTrkHits_);
    theTree->Branch("mupNTrkLayers",    &mupNTrkLayers_);
    theTree->Branch("mupNMuonHits",     &mupNMuonHits_);
    theTree->Branch("mupNMatchStation", &mupNMatchStation_);
    theTree->Branch("mup_isGlobalMuon", &mup_isGlobalMuon_);
    theTree->Branch("mup_isTrackerMuon",&mup_isTrackerMuon_);
    theTree->Branch("mup_StandAloneMuon",&mup_StandAloneMuon_);
    theTree->Branch("mup_isCaloMuon",    &mup_isCaloMuon_);
    theTree->Branch("mup_isPFMuon",      &mup_isPFMuon_);

    theTree->Branch("mup_selector"          	,  &mup_selector_); 
    theTree->Branch("mup_isIsolationValid"  	,  &mup_isIsolationValid_);
    theTree->Branch("mup_isPFIsolationValid"	,  &mup_isPFIsolationValid_);
   
    theTree->Branch("mup_isolationR03_trackSumPt"			,&mup_isolationR03_trackSumPt_);
    theTree->Branch("mup_isolationR03_trackEcalSumEt"			,&mup_isolationR03_trackEcalSumEt_);
    theTree->Branch("mup_isolationR03_trackHcalSumEt"			,&mup_isolationR03_trackHcalSumEt_);
    theTree->Branch("mup_isolationR03_trackHoSumEt"			,&mup_isolationR03_trackHoSumEt_);
    theTree->Branch("mup_isolationR03_trackNTracks"			,&mup_isolationR03_trackNTracks_);
    theTree->Branch("mup_isolationR03_trackNJets"			,&mup_isolationR03_trackNJets_);
    theTree->Branch("mup_isolationR03_trackerVetoSumPt"			,&mup_isolationR03_trackerVetoSumPt_);
    theTree->Branch("mup_isolationR03_emVetoSumEt"			,&mup_isolationR03_emVetoSumEt_);
    theTree->Branch("mup_isolationR03_hadVetoSumEt"			,&mup_isolationR03_hadVetoSumEt_);
    theTree->Branch("mup_isolationR03_hoVetoEt"				,&mup_isolationR03_hoVetoEt_);

    theTree->Branch("mup_isolationR05_trackSumPt"			,&mup_isolationR05_trackSumPt_);
    theTree->Branch("mup_isolationR05_trackEcalSumEt"			,&mup_isolationR05_trackEcalSumEt_);
    theTree->Branch("mup_isolationR05_trackHcalSumEt"			,&mup_isolationR05_trackHcalSumEt_);
    theTree->Branch("mup_isolationR05_trackHoSumEt"			,&mup_isolationR05_trackHoSumEt_);
    theTree->Branch("mup_isolationR05_trackNTracks"			,&mup_isolationR05_trackNTracks_);
    theTree->Branch("mup_isolationR05_trackNJets"			,&mup_isolationR05_trackNJets_);
    theTree->Branch("mup_isolationR05_trackerVetoSumPt"			,&mup_isolationR05_trackerVetoSumPt_);
    theTree->Branch("mup_isolationR05_emVetoSumEt"			,&mup_isolationR05_emVetoSumEt_);
    theTree->Branch("mup_isolationR05_hadVetoSumEt"			,&mup_isolationR05_hadVetoSumEt_);
    theTree->Branch("mup_isolationR05_hoVetoEt"				,&mup_isolationR05_hoVetoEt_);
                                                                                                                                                                                                                                                          
    theTree->Branch("mup_PFIsolationR03_sumChargedHadronPt"		,&mup_PFIsolationR03_sumChargedHadronPt_);
    theTree->Branch("mup_PFIsolationR03_sumChargedParticlePt"		,&mup_PFIsolationR03_sumChargedParticlePt_);
    theTree->Branch("mup_PFIsolationR03_sumNeutralHadronEt"		,&mup_PFIsolationR03_sumNeutralHadronEt_);
    theTree->Branch("mup_PFIsolationR03_sumPhotonEt"			,&mup_PFIsolationR03_sumPhotonEt_);
    theTree->Branch("mup_PFIsolationR03_sumNeutralHadronEtHighThreshold",&mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR03_sumPhotonEtHighThreshold"	,&mup_PFIsolationR03_sumPhotonEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR03_sumPUPt"			        ,&mup_PFIsolationR03_sumPUPt_);
                                                                                                                             
    theTree->Branch("mup_PFIsolationR04_sumChargedHadronPt"		,&mup_PFIsolationR04_sumChargedHadronPt_);
    theTree->Branch("mup_PFIsolationR04_sumChargedParticlePt"		,&mup_PFIsolationR04_sumChargedParticlePt_);
    theTree->Branch("mup_PFIsolationR04_sumNeutralHadronEt"		,&mup_PFIsolationR04_sumNeutralHadronEt_);
    theTree->Branch("mup_PFIsolationR04_sumPhotonEt"			,&mup_PFIsolationR04_sumPhotonEt_);
    theTree->Branch("mup_PFIsolationR04_sumNeutralHadronEtHighThreshold",&mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR04_sumPhotonEtHighThreshold"	,&mup_PFIsolationR04_sumPhotonEtHighThreshold_);
    theTree->Branch("mup_PFIsolationR04_sumPUPt"			        ,&mup_PFIsolationR04_sumPUPt_);

    // ### mu+ mu- Mass ###
    theTree->Branch("mumuPt",    &mumuPt_);
    theTree->Branch("mumuEta",   &mumuEta_);
    theTree->Branch("mumuRapidity",&mumuRapidity_);
    theTree->Branch("mumuPhi",   &mumuPhi_);
    theTree->Branch("mumuMass",  &mumuMass_);
    theTree->Branch("mumuMassE", &mumuMassE_);
    theTree->Branch("mumuPx",    &mumuPx_);
    theTree->Branch("mumuPy",    &mumuPy_);
    theTree->Branch("mumuPz",    &mumuPz_);
    theTree->Branch("mumuCovPxPx_",    &mumuCovPxPx_);
    theTree->Branch("mumuCovPxPy_",    &mumuCovPxPy_);
    theTree->Branch("mumuCovPxPz_",    &mumuCovPxPz_);
    theTree->Branch("mumuCovPyPy_",    &mumuCovPyPy_);
    theTree->Branch("mumuCovPyPz_",    &mumuCovPyPz_);
    theTree->Branch("mumuCovPzPz_",    &mumuCovPzPz_);
    theTree->Branch("mumuDR",    &mumuDR_);
    theTree->Branch("mumuParentMuP",    &mumuParentMuP_);
    theTree->Branch("mumuParentMuM",    &mumuParentMuM_);

    // ### mu+ mu- Vtx ###
    theTree->Branch("mumuVtxCL",       &mumuVtxCL_);
    theTree->Branch("mumuVtxX",        &mumuVtxX_);
    theTree->Branch("mumuVtxY",        &mumuVtxY_);
    theTree->Branch("mumuVtxZ",        &mumuVtxZ_);
    theTree->Branch("mumuVtxCovXX_",   &mumuVtxCovXX_);
    theTree->Branch("mumuVtxCovXY_",   &mumuVtxCovXY_);
    theTree->Branch("mumuVtxCovXZ_",   &mumuVtxCovYY_);
    theTree->Branch("mumuVtxCovYY_",   &mumuVtxCovYY_);
    theTree->Branch("mumuVtxCovYZ_",   &mumuVtxCovYZ_);
    theTree->Branch("mumuVtxCovZZ_",   &mumuVtxCovZZ_);
    theTree->Branch("mumuVtxChi2",     &mumuVtxChi2_);
    theTree->Branch("mumuVtxNdof",     &mumuVtxNdof_);

    theTree->Branch("mumuCosAlphaBS",  &mumuCosAlphaBS_);
    theTree->Branch("mumuCosAlphaBSE", &mumuCosAlphaBSE_);
    theTree->Branch("mumuLBS",         &mumuLBS_);
    theTree->Branch("mumuLBSE",        &mumuLBSE_);
    theTree->Branch("mumuDCA",         &mumuDCA_);


  }
  
  if (doPhotons_) {
    theTree->Branch("nPho",                  &nPho_);
    theTree->Branch("phoE",                  &phoE_);
    theTree->Branch("phoEt",                 &phoEt_);
    theTree->Branch("phoEta",                &phoEta_);
    theTree->Branch("phoPhi",                &phoPhi_);
    theTree->Branch("phoSigmaE",               &phoSigmaE_);
    theTree->Branch("phoCalibE",               &phoCalibE_);
    theTree->Branch("phoCalibEt",              &phoCalibEt_);
    theTree->Branch("phoSCE",                  &phoSCE_);
    theTree->Branch("phoSCEt",                 &phoSCEt_);
    theTree->Branch("phoSCRawE",               &phoSCRawE_);
    theTree->Branch("phoESEnP1",               &phoESEnP1_);
    theTree->Branch("phoESEnP2",               &phoESEnP2_);
    theTree->Branch("phoSCEta",                &phoSCEta_);
    theTree->Branch("phoSCPhi",                &phoSCPhi_);
    theTree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
    theTree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
    theTree->Branch("phoSCBrem",               &phoSCBrem_);
    theTree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
    theTree->Branch("phoEleVeto",              &phoEleVeto_);
    theTree->Branch("phoR9",                   &phoR9_);
    theTree->Branch("phoHoverE",               &phoHoverE_);
    theTree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
    theTree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
    theTree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
    theTree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
    theTree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
    theTree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
    theTree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
    
    theTree->Branch("phoPFChIso",              &phoPFChIso_);
    theTree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
    theTree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
    theTree->Branch("phoEcalPFClusterIso",     &phoEcalPFClusterIso_);
    theTree->Branch("phoHcalPFClusterIso",     &phoHcalPFClusterIso_);
    theTree->Branch("phoIDMVA",                &phoIDMVA_);
   
    theTree->Branch("phoSeedTime",             &phoSeedTime_);
    theTree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
    theTree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);
    theTree->Branch("phoMIPChi2",                      &phoMIPChi2_);
    theTree->Branch("phoMIPSlope",                     &phoMIPSlope_);
    theTree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
    theTree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
    theTree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

  }

  if (doPFPhotons_) {
    theTree->Branch("nPFPho",                  &nPFPho_);
    theTree->Branch("phoPFE",                  &phoPFE_);
    theTree->Branch("phoPFEt",                 &phoPFEt_);
    theTree->Branch("phoPFEta",                &phoPFEta_);
    theTree->Branch("phoPFPhi",                &phoPFPhi_);
  }

  if (doSuperClusters_) {
    theTree->Branch("nSC",                  &nSC_);
    theTree->Branch("scE",                  &scE_);
    theTree->Branch("scEt",                 &scEt_);
    theTree->Branch("scRawE",               &scRawE_);
    theTree->Branch("scEta",                &scEta_);
    theTree->Branch("scPhi",                &scPhi_);
    theTree->Branch("scX",        &scX_);
    theTree->Branch("scY",        &scY_);
    theTree->Branch("scZ",        &scZ_);
    theTree->Branch("scEtaWidth", &scEtaWidth_);
    theTree->Branch("scPhiWidth", &scPhiWidth_);
    theTree->Branch("scRawEt",    &scRawEt_);
    theTree->Branch("scMinDrWithGsfElectornSC_",  &scMinDrWithGsfElectornSC_);
    theTree->Branch("scFoundGsfMatch_" ,        &scFoundGsfMatch_);

    theTree->Branch("scE5x5",   &scE5x5_);
    theTree->Branch("scE2x2Ratio",   &scE2x2Ratio_);
    theTree->Branch("scE3x3Ratio",   &scE3x3Ratio_);
    theTree->Branch("scEMaxRatio",   &scEMaxRatio_);
    theTree->Branch("scE2ndRatio",   &scE2ndRatio_);
    theTree->Branch("scETopRatio",   &scETopRatio_);
    theTree->Branch("scERightRatio",   &scERightRatio_);
    theTree->Branch("scEBottomRatio",   &scEBottomRatio_);
    theTree->Branch("scELeftRatio",   &scELeftRatio_);
    theTree->Branch("scE2x5MaxRatio",   &scE2x5MaxRatio_);
    theTree->Branch("scE2x5TopRatio",   &scE2x5TopRatio_);
    theTree->Branch("scE2x5RightRatio",   &scE2x5RightRatio_);
    theTree->Branch("scE2x5BottomRatio",   &scE2x5BottomRatio_); 
    theTree->Branch("scE2x5LeftRatio",   &scE2x5LeftRatio_); 
    theTree->Branch("scSwissCross",   &scSwissCross_); 
    theTree->Branch("scR9",   &scR9_);
    theTree->Branch("scSigmaIetaIeta",   &scSigmaIetaIeta_);
    theTree->Branch("scSigmaIetaIphi",   &scSigmaIetaIphi_);
    theTree->Branch("scSigmaIphiIphi",   &scSigmaIphiIphi_);
    theTree->Branch("scFull5x5_e5x5",   &scFull5x5_e5x5_);
    theTree->Branch("scFull5x5_e2x2Ratio",   &scFull5x5_e2x2Ratio_);
    theTree->Branch("scFull5x5_e3x3Ratio",   &scFull5x5_e3x3Ratio_);
    theTree->Branch("scFull5x5_eMaxRatio",   &scFull5x5_eMaxRatio_);
    theTree->Branch("scFull5x5_e2ndRatio",   &scFull5x5_e2ndRatio_);
    theTree->Branch("scFull5x5_eTopRatio",   &scFull5x5_eTopRatio_);
    theTree->Branch("scFull5x5_eRightRatio",   &scFull5x5_eRightRatio_);
    theTree->Branch("scFull5x5_eBottomRatio",   &scFull5x5_eBottomRatio_);
    theTree->Branch("scFull5x5_eLeftRatio",   &scFull5x5_eLeftRatio_);
    theTree->Branch("scFull5x5_e2x5MaxRatio",   &scFull5x5_e2x5MaxRatio_);
    theTree->Branch("scFull5x5_e2x5TopRatio",   &scFull5x5_e2x5TopRatio_);
    theTree->Branch("scFull5x5_e2x5RightRatio",   &scFull5x5_e2x5RightRatio_);
    theTree->Branch("scFull5x5_e2x5BottomRatio",   &scFull5x5_e2x5BottomRatio_); 
    theTree->Branch("scFull5x5_e2x5LeftRatio",   &scFull5x5_e2x5LeftRatio_); 
    theTree->Branch("scFull5x5_swissCross",   &scFull5x5_swissCross_); 
    theTree->Branch("scFull5x5_r9",   &scFull5x5_r9_);
    theTree->Branch("scFull5x5_sigmaIetaIeta",   &scFull5x5_sigmaIetaIeta_);
    theTree->Branch("scFull5x5_sigmaIetaIphi",   &scFull5x5_sigmaIetaIphi_);
    theTree->Branch("scFull5x5_sigmaIphiIphi",   &scFull5x5_sigmaIphiIphi_);
  
    theTree->Branch("nhcalRechit",      &nhcalRechit_);
    theTree->Branch("hcalRechitIEta",   &hcalRechitIEta_);
    theTree->Branch("hcalRechitIPhi",   &hcalRechitIPhi_);
    theTree->Branch("hcalRechitEnergy", &hcalRechitEnergy_);

    theTree->Branch("scPFChIso1",              &scPFChIso1_);
    theTree->Branch("scPFChIso2",              &scPFChIso2_);
    theTree->Branch("scPFChIso3",              &scPFChIso3_);
    theTree->Branch("scPFChIso4",              &scPFChIso4_);
    theTree->Branch("scPFChIso5",              &scPFChIso5_);
    
    theTree->Branch("scPFPhoIso1",             &scPFPhoIso1_);
    theTree->Branch("scPFPhoIso2",             &scPFPhoIso2_);
    theTree->Branch("scPFPhoIso3",             &scPFPhoIso3_);
    theTree->Branch("scPFPhoIso4",             &scPFPhoIso4_);
    theTree->Branch("scPFPhoIso5",             &scPFPhoIso5_);
    
    theTree->Branch("scPFNeuIso1",             &scPFNeuIso1_);
    theTree->Branch("scPFNeuIso2",             &scPFNeuIso2_);
    theTree->Branch("scPFNeuIso3",             &scPFNeuIso3_);
    theTree->Branch("scPFNeuIso4",             &scPFNeuIso4_);
    theTree->Branch("scPFNeuIso5",             &scPFNeuIso5_);

  }
}

//BsToMuMuGammaNTuplizer::~BsToMuMuGammaNTuplizer(){
//}


// ------------ method called for each event  ------------
void 
BsToMuMuGammaNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  // ## BEAMSOPT STUFF  ## //
  beamspot_x_  = 0.0   ;
  beamspot_y_  = 0.0   ;
  beamspot_z_  = 0.0   ;
  beamspot_x_error_  = 0.0   ;
  beamspot_y_error_  = 0.0   ;
  beamspot_z_error_  = 0.0   ;
  beamspot_dxdz_  = 0.0   ;
  beamspot_dydz_  = 0.0   ;
  beamspot_sigmaZ_  = 0.0   ;
  beamspot_dxdz_error_  = 0.0   ;
  beamspot_dydz_error_  = 0.0   ;
  beamspot_sigmaZError_  = 0.0   ;
  beamspot_beamWidthX_  = 0.0   ;
  beamspot_beamWidthY_  = 0.0   ;
  beamspot_beamWidthX_error_  = 0.0   ;
  beamspot_beamWidthY_error_  = 0.0   ;

  // # offlinePrimaryVertices # //
  nPrimaryVertex_ = 0;
  primaryVertex_isFake_.clear();
  primaryVertex_x_.clear();
  primaryVertex_y_.clear();
  primaryVertex_z_.clear();
  primaryVertex_t_.clear();
  primaryVertex_covXX.clear();
  primaryVertex_covXY.clear();
  primaryVertex_covXZ.clear();
  primaryVertex_covYY.clear();
  primaryVertex_covYZ.clear();
  primaryVertex_covZZ.clear();
  primaryVertex_x_error_.clear();
  primaryVertex_y_error_.clear();
  primaryVertex_z_error_.clear();
  primaryVertex_t_error_.clear();
  primaryVertex_ntracks_.clear();
  primaryVertex_ndof_.clear();
  primaryVertex_chi2_.clear();
  primaryVertex_normalizedChi2_.clear();

  if (doHLT) {

    // ### Trigger ###
    //TrigTable.clear();
    trigNames.clear();
    trigResult.clear();
    trigPrescales.clear();
    l1Table.clear();
    l1Prescales.clear();
    ClearTrggerStorages();
  
  }
  if (doGenParticles_) {
  
    gen_nBs_ = 0;
    gen_nBsMuonM_ = 0 ;
    gen_nBsMuonP_ = 0 ;
    gen_nBsPhoton_ = 0 ;

    gen_Bs_pt_.clear() ;
    gen_Bs_energy_.clear() ;
    gen_Bs_eta_.clear() ;
    gen_Bs_phi_.clear() ;
    gen_Bs_pz_.clear() ;
    gen_Bs_pdgId_.clear();
    gen_BsMuonM_pt_.clear() ;
    gen_BsMuonM_eta_.clear() ;
    gen_BsMuonM_phi_.clear();
    gen_BsMuonP_pt_.clear() ;
    gen_BsMuonP_eta_.clear() ;
    gen_BsMuonP_phi_.clear();
    gen_BsPhoton_pt_.clear() ;
    gen_BsPhoton_energy_.clear() ;
    gen_BsPhoton_eta_.clear() ;
    gen_BsPhoton_phi_.clear();
  
   if(doFlatPt_){
    nMC_ = 0;
    mcPID_                .clear();
    mcStatus_             .clear();
    mcVtx_x_              .clear();
    mcVtx_y_              .clear();
    mcVtx_z_              .clear();
    mcPt_                 .clear();
    mcEta_                .clear();
    mcPhi_                .clear();
    mcE_                  .clear();
    mcEt_                 .clear();
    mcMass_               .clear();
     }
   }

  if(doMuons_){
    nMuM_ = 0;
    nMuP_ = 0; 
    mumHighPurity_.clear();
    mumPt_.clear();
    mumEta_.clear();
    mumPhi_.clear();
    mumCL_.clear(); 
    mumNormChi2_.clear();
    mumVx_.clear();
    mumVy_.clear();
    mumVz_.clear();
 
    mumDCABS_.clear();
    mumDCABSE_.clear();
    mumFracHits_.clear();
    mumdxyBS_.clear();
    mumdzBS_.clear();
    mumIdx_.clear();
    mumCharge_.clear();
    mumNPixHits_.clear();
    mumNPixLayers_.clear();
    mumNTrkHits_.clear();
    mumNTrkLayers_.clear();
    mumNMuonHits_.clear();
    mumNMatchStation_.clear();
    mum_isGlobalMuon_.clear();
    mum_isTrackerMuon_.clear();
    mum_StandAloneMuon_.clear();
    mum_isCaloMuon_.clear();
    mum_isPFMuon_.clear();

    mum_selector_.clear();
 
    mum_isIsolationValid_.clear();
    mum_isPFIsolationValid_.clear();
   
    mum_isolationR03_trackSumPt_.clear();
    mum_isolationR03_trackEcalSumEt_.clear();
    mum_isolationR03_trackHcalSumEt_.clear();
    mum_isolationR03_trackHoSumEt_.clear();
    mum_isolationR03_trackNTracks_.clear();
    mum_isolationR03_trackNJets_.clear();
    mum_isolationR03_trackerVetoSumPt_.clear();
    mum_isolationR03_emVetoSumEt_.clear();
    mum_isolationR03_hadVetoSumEt_.clear();
    mum_isolationR03_hoVetoEt_.clear();
  
    mum_isolationR05_trackSumPt_.clear();
    mum_isolationR05_trackEcalSumEt_.clear();
    mum_isolationR05_trackHcalSumEt_.clear();
    mum_isolationR05_trackHoSumEt_.clear();
    mum_isolationR05_trackNTracks_.clear();
    mum_isolationR05_trackNJets_.clear();
    mum_isolationR05_trackerVetoSumPt_.clear();
    mum_isolationR05_emVetoSumEt_.clear();
    mum_isolationR05_hadVetoSumEt_.clear();
    mum_isolationR05_hoVetoEt_.clear();
  
    mum_PFIsolationR03_sumChargedHadronPt_.clear();
    mum_PFIsolationR03_sumChargedParticlePt_.clear();
    mum_PFIsolationR03_sumNeutralHadronEt_.clear();
    mum_PFIsolationR03_sumPhotonEt_.clear();
    mum_PFIsolationR03_sumNeutralHadronEtHighThreshold_.clear();
    mum_PFIsolationR03_sumPhotonEtHighThreshold_.clear();
    mum_PFIsolationR03_sumPUPt_.clear();
    
    mum_PFIsolationR04_sumChargedHadronPt_.clear();
    mum_PFIsolationR04_sumChargedParticlePt_.clear();
    mum_PFIsolationR04_sumNeutralHadronEt_.clear();
    mum_PFIsolationR04_sumPhotonEt_.clear();
    mum_PFIsolationR04_sumNeutralHadronEtHighThreshold_.clear();
    mum_PFIsolationR04_sumPhotonEtHighThreshold_.clear();
    mum_PFIsolationR04_sumPUPt_.clear();
    
    mupHighPurity_.clear();
    mupPt_.clear();
    mupEta_.clear();
    mupPhi_.clear();
    mupCL_.clear(); 
    mupNormChi2_.clear();
    mupVx_.clear();
    mupVy_.clear();
    mupVz_.clear();
    mupDCABS_.clear();
    mupDCABSE_.clear();
  
    mupFracHits_.clear();
    mupdxyBS_.clear();
    mupdzBS_.clear();
    mupIdx_.clear();
    mupCharge_.clear();
    mupNPixHits_.clear();
    mupNPixLayers_.clear();
    mupNTrkHits_.clear();
    mupNTrkLayers_.clear();
    mupNMuonHits_.clear();
    mupNMatchStation_.clear();
    mup_isGlobalMuon_.clear();
    mup_isTrackerMuon_.clear();
    mup_StandAloneMuon_.clear();
    mup_isCaloMuon_.clear();
    mup_isPFMuon_.clear();

    mup_selector_.clear();
 
    mup_isIsolationValid_.clear();
    mup_isPFIsolationValid_.clear();
   
    mup_isolationR03_trackSumPt_.clear();
    mup_isolationR03_trackEcalSumEt_.clear();
    mup_isolationR03_trackHcalSumEt_.clear();
    mup_isolationR03_trackHoSumEt_.clear();
    mup_isolationR03_trackNTracks_.clear();
    mup_isolationR03_trackNJets_.clear();
    mup_isolationR03_trackerVetoSumPt_.clear();
    mup_isolationR03_emVetoSumEt_.clear();
    mup_isolationR03_hadVetoSumEt_.clear();
    mup_isolationR03_hoVetoEt_.clear();
  
    mup_isolationR05_trackSumPt_.clear();
    mup_isolationR05_trackEcalSumEt_.clear();
    mup_isolationR05_trackHcalSumEt_.clear();
    mup_isolationR05_trackHoSumEt_.clear();
    mup_isolationR05_trackNTracks_.clear();
    mup_isolationR05_trackNJets_.clear();
    mup_isolationR05_trackerVetoSumPt_.clear();
    mup_isolationR05_emVetoSumEt_.clear();
    mup_isolationR05_hadVetoSumEt_.clear();
    mup_isolationR05_hoVetoEt_.clear();
  
    mup_PFIsolationR03_sumChargedHadronPt_.clear();
    mup_PFIsolationR03_sumChargedParticlePt_.clear();
    mup_PFIsolationR03_sumNeutralHadronEt_.clear();
    mup_PFIsolationR03_sumPhotonEt_.clear();
    mup_PFIsolationR03_sumNeutralHadronEtHighThreshold_.clear();
    mup_PFIsolationR03_sumPhotonEtHighThreshold_.clear();
    mup_PFIsolationR03_sumPUPt_.clear();
    
    mup_PFIsolationR04_sumChargedHadronPt_.clear();
    mup_PFIsolationR04_sumChargedParticlePt_.clear();
    mup_PFIsolationR04_sumNeutralHadronEt_.clear();
    mup_PFIsolationR04_sumPhotonEt_.clear();
    mup_PFIsolationR04_sumNeutralHadronEtHighThreshold_.clear();
    mup_PFIsolationR04_sumPhotonEtHighThreshold_.clear();
    mup_PFIsolationR04_sumPUPt_.clear();

    // ### mu+ mu- variables ###
    mumuPt_.clear();
    mumuEta_.clear();
    mumuRapidity_.clear();
    mumuPhi_.clear();
    mumuMass_.clear();
    mumuMassE_.clear();
    mumuPx_.clear();
    mumuPy_.clear();
    mumuPz_.clear();
    mumuCovPxPx_.clear();
    mumuCovPxPy_.clear();
    mumuCovPxPz_.clear();
    mumuCovPyPy_.clear();
    mumuCovPyPz_.clear();
    mumuCovPzPz_.clear();
    mumuDR_.clear();
    mumuParentMuP_.clear();
    mumuParentMuM_.clear();
  
    mumuVtxCL_.clear();
    mumuVtxX_.clear();
    mumuVtxY_.clear();
    mumuVtxZ_.clear();  
    mumuVtxCovXX_.clear();  
    mumuVtxCovXY_.clear();  
    mumuVtxCovXZ_.clear();  
    mumuVtxCovYY_.clear();  
    mumuVtxCovYZ_.clear();  
    mumuVtxCovZZ_.clear();  
    mumuVtxChi2_.clear();
    mumuVtxNdof_.clear();


    mumuCosAlphaBS_.clear();
    mumuCosAlphaBSE_.clear(); 
    mumuLBS_.clear();
    mumuLBSE_.clear();
    mumuDCA_.clear();
 
    
    nMuP_=0;
    nMuM_=0;
  }
  
  if (doPhotons_) {
    nPho_ = 0;
    phoE_                 .clear();
    phoEt_                .clear();
    phoEta_               .clear();
    phoPhi_               .clear();
    phoSigmaE_              .clear();
    phoCalibE_              .clear();
    phoCalibEt_             .clear();
    phoSCE_                 .clear();
    phoSCEt_                 .clear();
    phoSCRawE_              .clear();
    phoESEnP1_              .clear();
    phoESEnP2_              .clear();
    phoSCEta_               .clear();
    phoSCPhi_               .clear();
    phoSCEtaWidth_          .clear();
    phoSCPhiWidth_          .clear();
    phoSCBrem_              .clear();
    phohasPixelSeed_        .clear();
    phoEleVeto_             .clear();
    phoR9_                  .clear();
    phoHoverE_              .clear();
    phoESEffSigmaRR_        .clear();
    phoSigmaIEtaIEtaFull5x5_.clear();
    phoSigmaIEtaIPhiFull5x5_.clear();
    phoSigmaIPhiIPhiFull5x5_.clear();
    phoE2x2Full5x5_         .clear();
    phoE5x5Full5x5_         .clear();
    phoR9Full5x5_           .clear();
    phoPFChIso_             .clear();
    phoPFPhoIso_            .clear();
    phoPFNeuIso_            .clear();
    phoEcalPFClusterIso_    .clear();
    phoHcalPFClusterIso_    .clear();
    phoIDMVA_               .clear();
    phoSeedTime_            .clear();
    phoSeedEnergy_          .clear();
    phoMIPTotEnergy_        .clear();  
    phoMIPChi2_           .clear();
    phoMIPSlope_          .clear();
    phoMIPIntercept_      .clear();
    phoMIPNhitCone_       .clear();
    phoMIPIsHalo_         .clear();
  }
  
  if (doPFPhotons_) {
    nPFPho_ = 0;
    phoPFE_                 .clear();
    phoPFEt_                .clear();
    phoPFEta_               .clear();
    phoPFPhi_               .clear();
  }

  if (doSuperClusters_) {
    nSC_ = 0;
    scE_                  .clear();
    scEt_                 .clear();
    scRawE_               .clear();
    scEta_                .clear();
    scPhi_                .clear();
    scX_          .clear();      
    scY_          .clear();      
    scZ_          .clear();      
    scEtaWidth_   .clear();         
    scPhiWidth_   .clear();         
    scRawEt_      .clear();   
    scMinDrWithGsfElectornSC_.clear();
    scFoundGsfMatch_.clear();

    scE5x5_.clear();
    scE2x2Ratio_.clear();
    scE3x3Ratio_.clear();
    scEMaxRatio_.clear();
    scE2ndRatio_.clear();
    scETopRatio_.clear();
    scERightRatio_.clear();
    scEBottomRatio_.clear();
    scELeftRatio_.clear();
    scE2x5MaxRatio_.clear();
    scE2x5TopRatio_.clear();
    scE2x5RightRatio_.clear();
    scE2x5BottomRatio_.clear();
    scE2x5LeftRatio_.clear();
    scSwissCross_.clear();
    scR9_.clear();
    scSigmaIetaIeta_.clear(); 
    scSigmaIetaIphi_.clear(); 
    scSigmaIphiIphi_.clear(); 
    scFull5x5_e5x5_.clear();
    scFull5x5_e2x2Ratio_.clear();
    scFull5x5_e3x3Ratio_.clear();
    scFull5x5_eMaxRatio_.clear();
    scFull5x5_e2ndRatio_.clear();
    scFull5x5_eTopRatio_.clear();
    scFull5x5_eRightRatio_.clear();
    scFull5x5_eBottomRatio_.clear();
    scFull5x5_eLeftRatio_.clear();
    scFull5x5_e2x5MaxRatio_.clear();
    scFull5x5_e2x5TopRatio_.clear();
    scFull5x5_e2x5RightRatio_.clear();
    scFull5x5_e2x5BottomRatio_.clear();
    scFull5x5_e2x5LeftRatio_.clear();
    scFull5x5_swissCross_.clear();
    scFull5x5_r9_.clear();
    scFull5x5_sigmaIetaIeta_.clear(); 
    scFull5x5_sigmaIetaIphi_.clear(); 
    scFull5x5_sigmaIphiIphi_.clear();  
  
    nhcalRechit_ =0;
    hcalRechitIEta_.clear();
    hcalRechitIPhi_.clear();
    hcalRechitEnergy_.clear();

    scPFChIso1_             .clear();
    scPFChIso2_             .clear();
    scPFChIso3_             .clear();
    scPFChIso4_             .clear();
    scPFChIso5_             .clear();
    
    scPFPhoIso1_            .clear();
    scPFPhoIso2_            .clear();
    scPFPhoIso3_            .clear();
    scPFPhoIso4_            .clear();
    scPFPhoIso5_            .clear();
    
    scPFNeuIso1_            .clear();
    scPFNeuIso2_            .clear();
    scPFNeuIso3_            .clear();
    scPFNeuIso4_            .clear();
    scPFNeuIso5_            .clear();

  }

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();
  isData_ = iEvent.isRealData();

  /*edm::Handle<std::vector<reco::Vertex>> primaryVertexCollection;
  iEvent.getByToken(primaryVtxToken_, primaryVertexCollection);
 
  // adding  BEAMSOPT 
  beamspot_x_			= beamSpot.x0();  ;
  beamspot_y_			= beamSpot.y0();  ;
  beamspot_z_			= beamSpot.z0();  ;
  beamspot_x_error_		= beamSpot.x0Error();  ;
  beamspot_y_error_		= beamSpot.y0Error();  ;
  beamspot_z_error_		= beamSpot.z0Error();  ;
  beamspot_covXX        = beamSpot.covariance()(0,0);
  beamspot_covXY        = beamSpot.covariance()(0,1);
  beamspot_covXZ        = beamSpot.covariance()(0,2);
  beamspot_covYY        = beamSpot.covariance()(1,1);
  beamspot_covYZ        = beamSpot.covariance()(1,2);
  beamspot_covZZ        = beamSpot.covariance()(2,2);

  beamspot_dxdz_   		= beamSpot.dxdz();  ;
  beamspot_dydz_         	= beamSpot.dydz();  ;
  beamspot_sigmaZ_		= beamSpot.sigmaZ();  ;
  beamspot_dxdz_error_		= beamSpot.dxdzError();  ;
  beamspot_dydz_error_		= beamSpot.dydzError();  ;
  beamspot_sigmaZError_		= beamSpot.sigmaZ0Error();  ;
  beamspot_beamWidthX_		= beamSpot.BeamWidthX();  ;
  beamspot_beamWidthY_		= beamSpot.BeamWidthY();  ;
  beamspot_beamWidthX_error_	= beamSpot.BeamWidthXError();  ;
  beamspot_beamWidthY_error_	= beamSpot.BeamWidthXError();  ;

  reco::Vertex pv(math::XYZPoint(0, 0, -999), math::Error<3>::type()); 
  for(auto&  aVertex : *primaryVertexCollection){

    if( not aVertex.isValid() ) continue;
    
    // # offlinePrimaryVertices # //
    primaryVertex_isFake_ .push_back(   aVertex.isFake() );
    primaryVertex_x_ .push_back(   aVertex.x() );
    primaryVertex_y_ .push_back(   aVertex.y()  );
    primaryVertex_z_ .push_back(   aVertex.z()  );
    primaryVertex_t_ .push_back(   aVertex.t()  );
    primaryVertex_covXX.push_back( aVertex.covariance(0,0) );
    primaryVertex_covXY.push_back( aVertex.covariance(0,1) );
    primaryVertex_covXZ.push_back( aVertex.covariance(0,2) );
    primaryVertex_covYY.push_back( aVertex.covariance(1,1) );
    primaryVertex_covYZ.push_back( aVertex.covariance(1,2) );
    primaryVertex_covZZ.push_back( aVertex.covariance(2,2) );
    primaryVertex_x_error_ .push_back(   aVertex.xError()  );
    primaryVertex_y_error_ .push_back(   aVertex.yError() );
    primaryVertex_z_error_ .push_back(   aVertex.zError()  );
    primaryVertex_t_error_ .push_back(   aVertex.tError()  );
    primaryVertex_ntracks_ .push_back(   aVertex.nTracks() );
    primaryVertex_ndof_ .push_back(   aVertex.ndof() 	 	  );
    primaryVertex_chi2_ .push_back(   aVertex.chi2()  );
    primaryVertex_normalizedChi2_ .push_back(   aVertex.normalizedChi2()  );

    pv = aVertex;
    nPrimaryVertex_++;
  } // loop over primary vertex collection
*/
  // MC truth
  if (isMC) {
    fillGenParticles(iEvent);
  }

  if (doHLT) 		fillHLT(iEvent);
  if (doMuons_)     	fillMuons(iEvent, iSetup);
  if (doPhotons_)    	fillPhotons(iEvent, iSetup);
  if (doPFPhotons_) 	fillPFPhotons(iEvent, iSetup);
  //if (doSuperClusters_) fillSC(iEvent, iSetup, pv);
  theTree->Fill();
  
}



void BsToMuMuGammaNTuplizer::fillGenParticles(const edm::Event& iEvent)
{
  
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  iEvent.getByToken(genParticlesCollection_, genParticleCollection);
 
  /*if(doFlatPt_){
    for (auto p = genParticleCollection->begin(); p != genParticleCollection->end(); ++p) {
      mcPID_   .push_back(p->pdgId());
      mcStatus_.push_back(p->status());
      mcVtx_x_ .push_back(p->vx());
      mcVtx_y_ .push_back(p->vy());
      mcVtx_z_ .push_back(p->vz());
      mcPt_    .push_back(p->pt());
      mcEta_   .push_back(p->eta());
      mcPhi_   .push_back(p->phi());
      mcE_     .push_back(p->energy());
      mcEt_    .push_back(p->et());
      mcMass_  .push_back(p->mass());
      nMC_++;      
    } 
  } //if flat pT
  */
  int phoMul(-1),muMMul(-1),muPMul(-1);
  
  
  for(auto& aBsMeson : *genParticleCollection){
    
    if(abs(aBsMeson.pdgId())!=531) continue;
    
    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == -13) muPMul++;
      if(bsDaughter.pdgId() ==  13) muMMul++;
      if(bsDaughter.pdgId() ==  22) phoMul++;
    }
    if(muMMul <0 or muPMul <0 ) continue;
    if(phoMul <0 and doBsToMuMuGamma  ) continue;

    for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
      auto& bsDaughter = *(aBsMeson.daughter(j));
      if(bsDaughter.pdgId() == 13) {
	gen_BsMuonM_pt_.push_back(bsDaughter.pt());
	gen_BsMuonM_eta_.push_back(bsDaughter.eta());
	gen_BsMuonM_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonM_++;
      }
      if(bsDaughter.pdgId() == -13){
	gen_BsMuonP_pt_.push_back(bsDaughter.pt());
	gen_BsMuonP_eta_.push_back(bsDaughter.eta());
	gen_BsMuonP_phi_.push_back(bsDaughter.phi());
	gen_nBsMuonP_++;
      }
      if(bsDaughter.pdgId() ==  22){
	gen_BsPhoton_pt_.push_back(bsDaughter.pt());
	gen_BsPhoton_energy_.push_back(bsDaughter.energy());
	gen_BsPhoton_eta_.push_back(bsDaughter.eta());
	gen_BsPhoton_phi_.push_back(bsDaughter.phi());
	gen_nBsPhoton_++;
	
      }
    }  //number of daughters

    gen_Bs_pt_.push_back(aBsMeson.pt());
    gen_Bs_energy_.push_back(aBsMeson.energy());
    gen_Bs_eta_.push_back(aBsMeson.eta());
    gen_Bs_phi_.push_back(aBsMeson.phi());
    gen_Bs_pz_.push_back(aBsMeson.pz());
    gen_Bs_pdgId_.push_back(aBsMeson.pdgId());
    gen_nBs_++;
    
  } // genparticle collection
  
} // fill gen particles


void BsToMuMuGammaNTuplizer::fillMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

 } //fill muons

void BsToMuMuGammaNTuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es)
{

}

void BsToMuMuGammaNTuplizer::fillPFPhotons(const edm::Event& e, const edm::EventSetup& es)
{

}


void BsToMuMuGammaNTuplizer::fillSC(edm::Event const& e, const edm::EventSetup& es, reco::Vertex& pv) {

}



void BsToMuMuGammaNTuplizer::SetupTriggerStorageVectors()
{
  auto numTrigs = trigTable.size();
  TrigPrescales_store = new std::vector<int> [numTrigs];
  TrigResult_store    = new std::vector<bool>[numTrigs];
}

void BsToMuMuGammaNTuplizer::ClearTrggerStorages()
{
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      TrigResult_store[i].clear();
      TrigPrescales_store[i].clear();
    }
}

void BsToMuMuGammaNTuplizer::FillTrggerBranches()
{
	
}

void BsToMuMuGammaNTuplizer::SetupTriggerBranches()
{
  std::string branchName;
  for(uint32_t i=0;i<trigTable.size();i++)
    {
      branchName=trigTable[i]+"_result";
      theTree->Branch(branchName.c_str(),&(TrigResult_store[i]));
      branchName=trigTable[i]+"_prescale";
      theTree->Branch(branchName.c_str(),&(TrigPrescales_store[i]));
    }
}

void BsToMuMuGammaNTuplizer::fillHLT(edm::Event const& iEvent)
{
   
}

std::vector<float> BsToMuMuGammaNTuplizer::getShowerShapes(reco::CaloCluster* caloBC, const EcalRecHitCollection* recHits, const CaloTopology *topology)
{
  std::vector<float> shapes;
  return shapes; 
}

float BsToMuMuGammaNTuplizer::reduceFloat(float val, int bits)
{
  if(!doCompression_) return val;
  else return MiniFloatConverter::reduceMantissaToNbitsRounding(val,bits);
}

int BsToMuMuGammaNTuplizer::calDIPhi(int iPhi1, int iPhi2) {
  int dPhi = iPhi1 - iPhi2;
  if (dPhi > 72 / 2)
    dPhi -= 72;
  else if (dPhi < -72 / 2)
    dPhi += 72;
  return dPhi;
}

int BsToMuMuGammaNTuplizer::calDIEta(int iEta1, int iEta2) {
  int dEta = iEta1 - iEta2;
  if (iEta1 * iEta2 < 0) {  //-ve to +ve transistion and no crystal at zero
    if (dEta < 0)
      dEta++;
    else
      dEta--;
  }
  return dEta;
}

float BsToMuMuGammaNTuplizer::getMinEnergyHCAL_(HcalDetId id) const {
  if ( (id.subdetId() == HcalBarrel)  ) {
    if ( (Run2_2018_ == 1) )
      return 0.7;
    else if ( (Run2_2018_ == 0) ) { //means Run3
      if (id.depth() == 1)
	return 0.1;
      else if (id.depth() == 2)
	return 0.2;
      else
	return 0.3;
    }
    else //neither 2018 , nor Run3, not supported
      return 9999.0;
  } 

  else if (id.subdetId() == HcalEndcap) {
    if (id.depth() == 1)
      return 0.1;
    else
      return 0.2;
  } else
    return 9999.0;
}
//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
