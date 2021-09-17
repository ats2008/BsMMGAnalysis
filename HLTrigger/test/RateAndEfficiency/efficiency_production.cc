/***************************
 *
 * Usage  root -b -q efficiency_production.c
 *
 ***********************************/

#include "HLTNtupleTreeV3.h"
#include "efficiencyMeasurement.h"
R__LOAD_LIBRARY(HLTNtupleTreeV3_C.so)

#define MU_MASS 0.1056583755
#define PHO_MASS 0.0
#define L1_DR_MU_MAX 0.3
#define L1_DR_PHO_MAX 0.5
#define HLT_DR_MU_MAX 0.1
#define HLT_DR_PHO_MAX 0.1
#define ECAL_BARREL_ETA_MAX 1.48
#define ECAL_ECAP_ETA_MIN 1.48
#define ECAL_ECAP_ETA_MAX 3.0

const Double_t TWO_PI(3.141592653589*2) ;
const Double_t PI(3.141592653589)       ;

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) ;
Double_t getDPHI( Double_t phi1, Double_t phi2) ;
Double_t getDETA(Double_t eta1, Double_t eta2) ;

void produceTurnOns(string infile,string ofileName,string prefix="",Double_t DR_MAX_L1=1.4, Double_t ETA_RESTRIC_L1 = 1.54,Long64_t maxEvents=100000);

//void efficiency_production(string infile,string ofileName,string prefix="",Long64_t maxEvents=100000)
void efficiency_production(string infile="",string ofileName="",string prefix="",Long64_t maxEvents=100000)
{

    produceTurnOns( "/grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/Run314TevBsTommGHLTEfficiencyNtupleV3_0.root",
                    "HLT_Hists.root",
                    "run3MC_",
                   1.4,
                   1.54,
                  -100);
    return;   
    produceTurnOns( "/grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/Run314TevBsTommGHLTEfficiencyNtupleV3_lowPtCust_0.root",
                    "HLT_Hists_LowPtCust.root",
                    "run3MC_",
                   1.4,
                   1.54,
                  -100);
}


// Selection Functions


void produceTurnOns(string infile,string ofileName,string prefix="",Double_t DR_MAX_L1=1.4, Double_t ETA_RESTRIC_L1 = 154,Long64_t maxEvents=100000)
{

    std::map<string,efficiencyMeasurement * > efficiencyMap;
    std::map<string, TH1D * > histMap;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
   
    const int XBINS=1;
    Double_t xEdges[XBINS+1];
    for(int i=0;i<=XBINS;i++) xEdges[i]= 0.5 + i*1.0;
    const int XBINS_forPt=30;
    Double_t xEdges_forPt[XBINS_forPt+1];
    for(int i=0;i<=XBINS_forPt;i++) xEdges_forPt[i]= 0.0 + i*1.0;
    const int XBINS_forEta=60;
    Double_t xEdges_forEta[XBINS_forEta+1];
    for(int i=0;i<=XBINS_forEta;i++) xEdges_forEta[i]= -3.0 + i*0.1;
    const int XBINS_forDr=40;
    Double_t xEdges_forDr[XBINS_forDr+1];
    for(int i=0;i<=XBINS_forDr;i++) xEdges_forDr[i]= 0.0 + i*1.0/XBINS_forDr;
    const int XBINS_forPhi=60;
    Double_t xEdges_forPhi[XBINS_forPhi+1];
    for(int i=0;i<=XBINS_forPhi;i++) xEdges_forPhi[i]= -3.5 + i*0.1;
    const int XBINS_forMass=120;
    Double_t xEdges_forMass[XBINS_forMass+1];
    for(int i=0;i<=XBINS_forMass;i++) xEdges_forMass[i]= 0.0 + i*0.1;
    const int XBINS_forDiMuMass=60;
    Double_t xEdges_forDiMuMass[XBINS_forDiMuMass+1];
    for(int i=0;i<=XBINS_forDiMuMass;i++) xEdges_forDiMuMass[i]= 0.0 + i*0.1;
    const int XBINS_forR9=40;
    Double_t xEdges_forR9[XBINS_forR9+1];
    for(int i=0;i<=XBINS_forR9;i++) xEdges_forR9[i]= 0.0 + i*4.0/XBINS_forR9;
    const int XBINS_forSigmaIEtaIEta=80;
    Double_t xEdges_forSigmaIEtaIEta[XBINS_forSigmaIEtaIEta+1];
    for(int i=0;i<=XBINS_forSigmaIEtaIEta;i++) xEdges_forSigmaIEtaIEta[i]= 0.0 + i*0.08/XBINS_forSigmaIEtaIEta;
    const int XBINS_forHoverE=1000;
    Double_t xEdges_forHoverE[XBINS_forHoverE+1];
    for(int i=0;i<=XBINS_forHoverE;i++) xEdges_forHoverE[i]= 0.0 + i*50.0/XBINS_forHoverE;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
     
    //  Histogram bookings
    histMap["gen_BsMuM_pt"] = new TH1D("gen_BsMuM_pt","gen_BsMuM_pt",XBINS_forPt,xEdges_forPt);
    histMap["gen_BsMuM_eta"] = new TH1D("gen_BsMuM_eta","gen_BsMuM_eta",XBINS_forEta,xEdges_forEta);
    histMap["gen_BsMuM_phi"] = new TH1D("gen_BsMuM_phi","gen_BsMuM_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsMuM_pt"] = new TH1D("recoL1_BsMuM_pt","recoL1_BsMuM_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoL1_BsMuM_eta"] = new TH1D("recoL1_BsMuM_eta","recoL1_BsMuM_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoL1_BsMuM_phi"] = new TH1D("recoL1_BsMuM_phi","recoL1_BsMuM_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsMuM_dr"] = new TH1D("recoL1_BsMuM_dr","recoL1_BsMuM_dr",XBINS_forDr,xEdges_forDr);
    histMap["recoHLT_BsMuM_pt"] = new TH1D("recoHLT_BsMuM_pt","recoHLT_BsMuM_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsMuM_eta"] = new TH1D("recoHLT_BsMuM_eta","recoHLT_BsMuM_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoHLT_BsMuM_phi"] = new TH1D("recoHLT_BsMuM_phi","recoHLT_BsMuM_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoHLT_BsMuM_dr"] = new TH1D("recoHLT_BsMuM_dr","recoHLT_BsMuM_dr",XBINS_forDr,xEdges_forDr);
    
    histMap["gen_BsMuP_pt"] = new TH1D("gen_BsMuP_pt","gen_BsMuP_pt",XBINS_forPt,xEdges_forPt);
    histMap["gen_BsMuP_eta"] = new TH1D("gen_BsMuP_eta","gen_BsMuP_eta",XBINS_forEta,xEdges_forEta);
    histMap["gen_BsMuP_phi"] = new TH1D("gen_BsMuP_phi","gen_BsMuP_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsMuP_pt"] = new TH1D("recoL1_BsMuP_pt","recoL1_BsMuP_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoL1_BsMuP_eta"] = new TH1D("recoL1_BsMuP_eta","recoL1_BsMuP_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoL1_BsMuP_phi"] = new TH1D("recoL1_BsMuP_phi","recoL1_BsMuP_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsMuP_dr"] = new TH1D("recoL1_BsMuP_dr","recoL1_BsMuP_dr",XBINS_forDr,xEdges_forDr);
    histMap["recoHLT_BsMuP_pt"] = new TH1D("recoHLT_BsMuP_pt","recoHLT_BsMuP_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsMuP_eta"] = new TH1D("recoHLT_BsMuP_eta","recoHLT_BsMuP_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoHLT_BsMuP_phi"] = new TH1D("recoHLT_BsMuP_phi","recoHLT_BsMuP_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoHLT_BsMuP_dr"] = new TH1D("recoHLT_BsMuP_dr","recoHLT_BsMuP_dr",XBINS_forDr,xEdges_forDr);
    
    histMap["gen_BsPho_pt"] = new TH1D("gen_BsPho_pt","gen_BsPho_pt",XBINS_forPt,xEdges_forPt);
    histMap["gen_BsPho_eta"] = new TH1D("gen_BsPho_eta","gen_BsPho_eta",XBINS_forEta,xEdges_forEta);
    histMap["gen_BsPho_phi"] = new TH1D("gen_BsPho_phi","gen_BsPho_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsPho_pt"] = new TH1D("recoL1_BsPho_pt","recoL1_BsPho_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoL1_BsPho_eta"] = new TH1D("recoL1_BsPho_eta","recoL1_BsPho_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoL1_BsPho_phi"] = new TH1D("recoL1_BsPho_phi","recoL1_BsPho_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoL1_BsPho_dr"] = new TH1D("recoL1_BsPho_dr","recoL1_BsPho_dr",XBINS_forDr,xEdges_forDr);
    histMap["recoHLT_BsPho_pt"] = new TH1D("recoHLT_BsPho_pt","recoHLT_BsPho_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsPho_eta"] = new TH1D("recoHLT_BsPho_eta","recoHLT_BsPho_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoHLT_BsPho_phi"] = new TH1D("recoHLT_BsPho_phi","recoHLT_BsPho_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoHLT_BsPho_HForHoverE"]     = new TH1D("recoHLT_BsPho_HForHoverE","recoHLT_BsPho_HForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_HoverE"]     = new TH1D("recoHLT_BsPho_HoverE","recoHLT_BsPho_HoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_R9"]         = new TH1D("recoHLT_BsPho_R9","recoHLT_BsPho_R9",XBINS_forR9,xEdges_forR9);
    histMap["recoHLT_BsPho_NoiseCleanedSigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_NoiseCleanedSigmaIEtaIEta","recoHLT_BsPho_NoiseCleanedSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_SigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_SigmaIEtaIEta","recoHLT_BsPho_SigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_BarrelPt"] = new TH1D("recoHLT_BsPho_BarrelPt","recoHLT_BsPho_BarrelPt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsPho_ECapHForHoverE"]     = new TH1D("recoHLT_BsPho_ECapHForHoverE","recoHLT_BsPho_ECapHForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_ECapHoverE"]     = new TH1D("recoHLT_BsPho_ECapHoverE","recoHLT_BsPho_ECapHoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_ECapR9"]         = new TH1D("recoHLT_BsPho_ECapR9","recoHLT_BsPho_ECapR9",XBINS_forR9,xEdges_forR9);
    histMap["recoHLT_BsPho_ECapNoiseCleanedSigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_ECapNoiseCleanedSigmaIEtaIEta","recoHLT_BsPho_ECapNoiseCleanedSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_ECapSigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_ECapSigmaIEtaIEta","recoHLT_BsPho_ECapSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_ECapPt"] = new TH1D("recoHLT_BsPho_ECapPt","recoHLT_BsPho_ECapPt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsPho_BarrelHForHoverE"]     = new TH1D("recoHLT_BsPho_BarrelHForHoverE","recoHLT_BsPho_BarrelHForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_BarrelHoverE"]     = new TH1D("recoHLT_BsPho_BarrelHoverE","recoHLT_BsPho_BarrelHoverE",XBINS_forHoverE,xEdges_forHoverE);
    histMap["recoHLT_BsPho_BarrelR9"]         = new TH1D("recoHLT_BsPho_BarrelR9","recoHLT_BsPho_BarrelR9",XBINS_forR9,xEdges_forR9);
    histMap["recoHLT_BsPho_BarrelNoiseCleanedSigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_BarrelNoiseCleanedSigmaIEtaIEta","recoHLT_BsPho_BarrelNoiseCleanedSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_BarrelSigmaIEtaIEta"] = new TH1D("recoHLT_BsPho_BarrelSigmaIEtaIEta","recoHLT_BsPho_BarrelSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    histMap["recoHLT_BsPho_dr"] = new TH1D("recoHLT_BsPho_dr","recoHLT_BsPho_dr",XBINS_forDr,xEdges_forDr);
    
    histMap["gen_BsDiMu_mass"] = new TH1D("gen_BsDiMu_mass","gen_BsDiMu_mass",XBINS_forDiMuMass,xEdges_forDiMuMass);
    histMap["gen_BsDiMu_pt"] = new TH1D("gen_BsDiMu_pt","gen_BsDiMu_pt",XBINS_forPt,xEdges_forPt);
    histMap["gen_BsDiMu_eta"] = new TH1D("gen_BsDiMu_eta","gen_BsDiMu_eta",XBINS_forEta,xEdges_forEta);
    histMap["gen_BsDiMu_phi"] = new TH1D("gen_BsDiMu_phi","gen_BsDiMu_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoHLT_BsDiMu_mass"] = new TH1D("recoHLT_BsDiMu_mass","recoHLT_BsDiMu_mass",XBINS_forDiMuMass,xEdges_forDiMuMass);
    histMap["recoHLT_BsDiMu_pt"] = new TH1D("recoHLT_BsDiMu_pt","recoHLT_BsDiMu_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_BsDiMu_eta"] = new TH1D("recoHLT_BsDiMu_eta","recoHLT_BsDiMu_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoHLT_BsDiMu_phi"] = new TH1D("recoHLT_BsDiMu_phi","recoHLT_BsDiMu_phi",XBINS_forPhi,xEdges_forPhi);
    
    histMap["gen_Bs_mass"] = new TH1D("gen_Bs_mass","gen_Bs_mass",XBINS_forMass,xEdges_forMass);
    histMap["gen_Bs_pt"] = new TH1D("gen_Bs_pt","gen_Bs_pt",XBINS_forPt,xEdges_forPt);
    histMap["gen_Bs_eta"] = new TH1D("gen_Bs_eta","gen_Bs_eta",XBINS_forEta,xEdges_forEta);
    histMap["gen_Bs_phi"] = new TH1D("gen_Bs_phi","gen_Bs_phi",XBINS_forPhi,xEdges_forPhi);
    histMap["recoHLT_Bs_mass"] = new TH1D("recoHLT_Bs_mass","recoHLT_Bs_mass",XBINS_forMass,xEdges_forMass);
    histMap["recoHLT_Bs_pt"]  = new TH1D("recoHLT_Bs_pt" ,"recoHLT_Bs_pt",XBINS_forPt,xEdges_forPt);
    histMap["recoHLT_Bs_eta"] = new TH1D("recoHLT_Bs_eta","recoHLT_Bs_eta",XBINS_forEta,xEdges_forEta);
    histMap["recoHLT_Bs_phi"] = new TH1D("recoHLT_Bs_phi","recoHLT_Bs_phi",XBINS_forPhi,xEdges_forPhi);
    
   
    // Effciency bookings

    efficiencyMap["dimuonL1"] = new efficiencyMeasurement("dimuonL1",XBINS,xEdges);
    
    efficiencyMap["muonPRecoL1vsPt"]  = new efficiencyMeasurement("muonPRecoL1vsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["muonPRecoL1vsEta"] = new efficiencyMeasurement("muonPRecoL1vsEta",XBINS_forEta,xEdges_forEta);
    
    efficiencyMap["muonMRecoL1vsPt"]  = new efficiencyMeasurement("muonMRecoL1vsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["muonMRecoL1vsEta"] = new efficiencyMeasurement("muonMRecoL1vsEta",XBINS_forEta,xEdges_forEta);
    
    efficiencyMap["phoRecoL1vsPt"]  = new efficiencyMeasurement("phoRecoL1vsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["phoRecoL1vsEta"] = new efficiencyMeasurement("phoRecoL1vsEta",XBINS_forEta,xEdges_forEta);
    
    
    efficiencyMap["dimuonHLT"] = new efficiencyMeasurement("dimuonHLT",XBINS,xEdges);
    
    efficiencyMap["muonPRecoHLTvsPt"]  = new efficiencyMeasurement("muonPRecoHLTvsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["muonPRecoHLTvsEta"] = new efficiencyMeasurement("muonPRecoHLTvsEta",XBINS_forEta,xEdges_forEta);
    
    efficiencyMap["muonMRecoHLTvsPt"]  = new efficiencyMeasurement("muonMRecoHLTvsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["muonMRecoHLTvsEta"] = new efficiencyMeasurement("muonMRecoHLTvsEta",XBINS_forEta,xEdges_forEta);
    
    efficiencyMap["phoRecoHLTvsPt"]  = new efficiencyMeasurement("phoRecoHLTvsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["phoRecoHLTvsEta"] = new efficiencyMeasurement("phoRecoHLTvsEta",XBINS_forEta,xEdges_forEta);
    efficiencyMap["phoRecoHLTvsHForHoverE"] = new efficiencyMeasurement("phoRecoHForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsHoverE"] = new efficiencyMeasurement("phoRecoHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsR9"] = new efficiencyMeasurement("phoRecoHLTvsR9",XBINS_forR9,xEdges_forR9);
    efficiencyMap["phoRecoHLTvsNoiseCleanedSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    efficiencyMap["phoRecoHLTvsSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    
    efficiencyMap["phoRecoHLTvsBarrelPt"]  = new efficiencyMeasurement("phoRecoHLTvsBarrelPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["phoRecoHLTvsECapPt"]  = new efficiencyMeasurement("phoRecoHLTECapvsPt",XBINS_forPt,xEdges_forPt);
    efficiencyMap["phoRecoHLTvsBarrelHForHoverE"] = new efficiencyMeasurement("phoRecoHLTvsBarrelHForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsBarrelHoverE"] = new efficiencyMeasurement("phoRecoHLTvsBarrelHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsBarrelR9"] = new efficiencyMeasurement("phoRecoHLTvsBarrelR9",XBINS_forR9,xEdges_forR9);
    efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsBarrelSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    efficiencyMap["phoRecoHLTvsECapHForHoverE"] = new efficiencyMeasurement("phoRecoHLTvsECapHForHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsECapHoverE"] = new efficiencyMeasurement("phoRecoHLTvsECapHoverE",XBINS_forHoverE,xEdges_forHoverE);
    efficiencyMap["phoRecoHLTvsECapR9"] = new efficiencyMeasurement("phoRecoHLTvsECapR9",XBINS_forR9,xEdges_forR9);
    efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta"] = new efficiencyMeasurement("phoRecoHLTvsECapSigmaIEtaIEta",XBINS_forSigmaIEtaIEta,xEdges_forSigmaIEtaIEta);
    
    efficiencyMap["HLTDimuonMass"] = new efficiencyMeasurement("HLTDimuonMass",XBINS_forDiMuMass,xEdges_forDiMuMass);
    efficiencyMap["HLTmu_pTMin"] = new efficiencyMeasurement("HLTmu_pTMin",XBINS_forPt,xEdges_forPt);
    efficiencyMap["HLTmu_pTMax"] = new efficiencyMeasurement("HLTmu_pTMax",XBINS_forPt,xEdges_forPt);


    
    TChain *treeChain = new TChain("egHLTRun3Tree");
    treeChain->Add(infile.c_str());
    HLTNtupleTreeV3 ntupleRawTree(treeChain);

	Long64_t nentries = treeChain->GetEntries();
    cout<<" Available total "<<nentries<<" \n";
    if (maxEvents >0 ) nentries = nentries > maxEvents ? maxEvents : nentries;
    cout<<" Processing total "<<nentries<<" \n";
    Long64_t nb = 0,nbytes=0 ;

    UInt_t nMus;
    UInt_t GenCount(0);
    UInt_t L1PassCount(0);
    UInt_t HLTDiMuPassCount(0);
    UInt_t HLTPassCount(0);
    Double_t dr;
    Double_t   dra=0.2; Int_t ida=-1;
    Double_t   drb=0.2; Int_t idb=-1;
    Double_t   drc=0.2; Int_t idc=-1;
    TLorentzVector aLV,bLV,cLV;

    bool    muP = false ;
    bool    muM = false ;
    bool    pho = false ;
    bool    muL1Pass = false ;
    bool HLTDiMuPass = false ;
    bool HLTPass = false ;
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {
       if(jentry%10000 ==0 )
       {
            cout<<"Processing jentry : "<<jentry<<"\n";
       }
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
        dra=L1_DR_MU_MAX; ida=-1;
        drb=L1_DR_MU_MAX; idb=-1;
        drc=L1_DR_PHO_MAX; idc=-1;
        muP = false ;
        muM = false ;
        pho = false ;
        muL1Pass = false;
        HLTDiMuPass = false;
        HLTPass = false;

        // GEN SELECTION    
        if( ntupleRawTree.bsMuP_pt < 3.5 || ntupleRawTree.bsMuM_pt < 3.5 || ntupleRawTree.bsPho_pt < 4.0  ) continue;
        if( abs(ntupleRawTree.bsMuP_eta) > ETA_RESTRIC_L1 || abs(ntupleRawTree.bsMuM_eta) > ETA_RESTRIC_L1) continue;
        dr = getDR(ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuM_phi,ntupleRawTree.bsMuP_eta , ntupleRawTree.bsMuP_phi  ) ;
        if( dr > DR_MAX_L1) continue;
        GenCount++;
        
        aLV.SetPtEtaPhiM(ntupleRawTree.bsMuP_pt,ntupleRawTree.bsMuP_eta,ntupleRawTree.bsMuP_phi,MU_MASS);
        bLV.SetPtEtaPhiM(ntupleRawTree.bsMuM_pt,ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuM_phi,MU_MASS);
        cLV.SetPtEtaPhiM(ntupleRawTree.bsPho_pt,ntupleRawTree.bsPho_eta,ntupleRawTree.bsPho_phi,PHO_MASS);
        

        histMap["gen_Bs_mass" ]->Fill((aLV + bLV +cLV).M());     histMap["gen_BsDiMu_mass" ]->Fill((aLV + bLV ).M() );
        histMap["gen_Bs_pt" ]->Fill((aLV + bLV +cLV).Pt());     histMap["gen_BsDiMu_pt" ]->Fill((aLV + bLV ).Pt() );
        histMap["gen_Bs_eta"]->Fill((aLV + bLV +cLV).Eta());    histMap["gen_BsDiMu_eta"]->Fill((aLV + bLV ).Eta());
        histMap["gen_Bs_phi"]->Fill((aLV + bLV +cLV).Phi());    histMap["gen_BsDiMu_phi"]->Fill((aLV + bLV ).Phi());

        histMap["gen_BsMuP_pt"]->Fill(ntupleRawTree.bsMuP_pt);    histMap["gen_BsMuM_pt"]->Fill(ntupleRawTree.bsMuP_pt);
        histMap["gen_BsMuP_eta"]->Fill(ntupleRawTree.bsMuP_eta);  histMap["gen_BsMuM_eta"]->Fill(ntupleRawTree.bsMuP_eta);
        histMap["gen_BsMuP_phi"]->Fill(ntupleRawTree.bsMuP_phi);  histMap["gen_BsMuM_phi"]->Fill(ntupleRawTree.bsMuP_phi);
        
        histMap["gen_BsPho_pt"]->Fill(ntupleRawTree.bsPho_pt);   
        histMap["gen_BsPho_eta"]->Fill(ntupleRawTree.bsPho_eta); 
        histMap["gen_BsPho_phi"]->Fill(ntupleRawTree.bsPho_phi); 

        //std::cout<<"      ==   ==    \n";
        //std::cout<<"gen muP   : pt "<<ntupleRawTree.bsMuP_pt<<" , eta "<<ntupleRawTree.bsMuP_eta<<" , phi "<<ntupleRawTree.bsMuP_phi<<" : "<<muP<<"\n";
        //std::cout<<"gen muM   : pt "<<ntupleRawTree.bsMuM_pt<<" , eta "<<ntupleRawTree.bsMuM_eta<<" , phi "<<ntupleRawTree.bsMuM_phi<<" : "<<muM<<"\n";
        //std::cout<<"gen pho   : pt "<<ntupleRawTree.bsPho_pt<<" , eta "<<ntupleRawTree.bsPho_eta<<" , phi "<<ntupleRawTree.bsPho_phi<<" : "<<pho<<"\n";
       
       // L1 SELECTION
        
        nMus=ntupleRawTree.nrAllL1Muons;
      
       for(auto i=0;i<nMus;i++)
        {
            if( ntupleRawTree.allL1mu_charge[i]== 1) {
                dr = getDR(ntupleRawTree.allL1mu_eta[i],ntupleRawTree.allL1mu_phi[i],ntupleRawTree.bsMuP_eta , ntupleRawTree.bsMuP_phi  ) ;
                //std::cout<<i<<" L1 Mu: "<<ntupleRawTree.allL1mu_charge[i]<<" : "<<ntupleRawTree.allL1mu_pt[i]<<" , "<<ntupleRawTree.allL1mu_eta[i]<<" , "<<ntupleRawTree.allL1mu_phi[i]<<"  dr : "<<dr<<"\n";
                if( dra > dr )
                {
                    dra=dr;
                    ida=i;
                }
                }
            if( ntupleRawTree.allL1mu_charge[i]==-1)
            {
                dr = getDR(ntupleRawTree.allL1mu_eta[i],ntupleRawTree.allL1mu_phi[i],ntupleRawTree.bsMuM_eta , ntupleRawTree.bsMuM_phi  ) ;
                //std::cout<<i<<" L1 Mu: "<<ntupleRawTree.allL1mu_charge[i]<<" : "<<ntupleRawTree.allL1mu_pt[i]<<" , "<<ntupleRawTree.allL1mu_eta[i]<<" , "<<ntupleRawTree.allL1mu_phi[i]<<"  dr : "<<dr<<"\n";
                if( drb > dr )
                {
                    drb=dr;
                    idb=i;
                }
            }
        }
 
     for(auto i=0;i<ntupleRawTree.nrAllL1EGs;i++)
        {
            dr = getDR(ntupleRawTree.allL1eg_eta[i],ntupleRawTree.allL1eg_phi[i],ntupleRawTree.bsPho_eta , ntupleRawTree.bsPho_phi  ) ;
            //std::cout<<i<<"L1 Eg : "<<ntupleRawTree.allL1eg_pt[i]<<" , "<<ntupleRawTree.allL1eg_eta[i]<<" , "<<ntupleRawTree.allL1eg_phi[i]<<"  dr : "<<dr<<"\n";
            if( drc > dr )
            {
                drc=dr;
                idc=i;
            }

        }


        if(ida>-1) if( abs(ntupleRawTree.allL1mu_eta[ida]) < ETA_RESTRIC_L1 ) muP=true;
        if(idb>-1) if( abs(ntupleRawTree.allL1mu_eta[idb]) < ETA_RESTRIC_L1 ) muM=true;
        if(idc>-1) pho=true;
        

        dr=-1;

        if(muP)
        {
            
            histMap["recoL1_BsMuP_pt"]->Fill(ntupleRawTree.allL1mu_pt[ida]);
            histMap["recoL1_BsMuP_eta"]->Fill(ntupleRawTree.allL1mu_eta[ida]);
            histMap["recoL1_BsMuP_phi"]->Fill(ntupleRawTree.allL1mu_phi[ida]);
            histMap["recoL1_BsMuP_dr"]->Fill(dra);
        }
        if(muM)
        {
            histMap["recoL1_BsMuM_pt"]->Fill(ntupleRawTree.allL1mu_pt[idb]);
            histMap["recoL1_BsMuM_eta"]->Fill(ntupleRawTree.allL1mu_eta[idb]);
            histMap["recoL1_BsMuM_phi"]->Fill(ntupleRawTree.allL1mu_phi[idb]);
            histMap["recoL1_BsMuM_dr"]->Fill(drb);
        }
        if(pho)
        {
            histMap["recoL1_BsPho_pt" ]->Fill(ntupleRawTree.allL1eg_pt[idc]);
            histMap["recoL1_BsPho_eta"]->Fill(ntupleRawTree.allL1eg_eta[idc]);
            histMap["recoL1_BsPho_phi"]->Fill(ntupleRawTree.allL1eg_phi[idc]);
            histMap["recoL1_BsPho_dr" ]->Fill(drc);
        }



        if( muP and muM )
        {
            dr = getDR(ntupleRawTree.allL1mu_eta[ida],ntupleRawTree.allL1mu_phi[ida],ntupleRawTree.allL1mu_eta[idb] , ntupleRawTree.allL1mu_phi[idb]  ) ;
            if( 
                dr < DR_MAX_L1 &&
                abs(ntupleRawTree.allL1mu_eta[ida]) < ETA_RESTRIC_L1 &&
                abs(ntupleRawTree.allL1mu_eta[idb]) < ETA_RESTRIC_L1 
              )
                {
                    muL1Pass = true;
                    L1PassCount+=1;
                }
        }
        
       efficiencyMap["muonPRecoL1vsPt"]->fill(  muP , 1.0  );
       efficiencyMap["muonPRecoL1vsEta"]->fill( muP , 1.0  );
       efficiencyMap["muonMRecoL1vsPt"]->fill(  muM , 1.0  );
       efficiencyMap["muonMRecoL1vsEta"]->fill( muM , 1.0  );
       efficiencyMap["phoRecoL1vsPt"]->fill ( pho , 1.0  );
       efficiencyMap["phoRecoL1vsEta"]->fill( pho , 1.0  );
       
       efficiencyMap["dimuonL1"]->fill( muL1Pass , 1.0  );
       
       //std::cout<<"L1 dr MuP= "<<dra<<"  MuM : "<<drb<<" : \t\tL1 RSLT : "<<muL1Pass<<"\n";
        
        if( not muL1Pass ) continue ;

        nMus=ntupleRawTree.nrL3Muons;
       
       // HLT Efficiency
        dra=HLT_DR_MU_MAX; ida=-1;
        drb=HLT_DR_MU_MAX; idb=-1;
        drc=HLT_DR_PHO_MAX; idc=-1;
        muP = false ;
        muM = false ;
        pho = false ;
       for(auto i=0;i<nMus;i++)
        {
            dr=-1.0;
            //std::cout<<"HLT L3Mu : "<<i<<" : "<<ntupleRawTree.l3Mu_charge[i]<<" : "<<ntupleRawTree.l3Mu_pt[i]<<" , "<<ntupleRawTree.l3Mu_eta[i]<<" , "<<ntupleRawTree.l3Mu_phi[i]<<" dr :";
            if(ntupleRawTree.l3Mu_charge[i]==-1) {
            dr = getDR(ntupleRawTree.l3Mu_eta[i],ntupleRawTree.l3Mu_phi[i],ntupleRawTree.bsMuM_eta , ntupleRawTree.bsMuM_phi  ) ;
            if( dra > dr )
            {
                dra=dr;
                ida=i;
            }
            }
                
            if(ntupleRawTree.l3Mu_charge[i]==1) {
            dr = getDR(ntupleRawTree.l3Mu_eta[i],ntupleRawTree.l3Mu_phi[i],ntupleRawTree.bsMuP_eta , ntupleRawTree.bsMuP_phi  ) ;
            if( drb > dr )
            {
                drb=dr;
                idb=i;
            }
            }
            //std::cout<<dr<< "\n";
        }
 
     for(auto i=0;i<ntupleRawTree.nrEgs;i++)
        {
            //std::cout<<i<<" HLT Eg : "<<ntupleRawTree.eg_et[i]<<" , "<<ntupleRawTree.eg_eta[i]<<" , "<<ntupleRawTree.eg_phi[i]<<"\n";
            dr = getDR(ntupleRawTree.eg_eta[i],ntupleRawTree.eg_phi[i],ntupleRawTree.bsPho_eta , ntupleRawTree.bsPho_phi  ) ;
            if( drc > dr )
            {
                drc=dr;
                idc=i;
            }
        }

        if(ida>-1)  muP=true;
        if(idb>-1)  muM=true;
        if(idc>-1)  pho=true;

        if(muP)
        {
            histMap["recoHLT_BsMuP_pt"]->Fill(ntupleRawTree.l3Mu_pt[ida]);
            histMap["recoHLT_BsMuP_eta"]->Fill(ntupleRawTree.l3Mu_eta[ida]);
            histMap["recoHLT_BsMuP_phi"]->Fill(ntupleRawTree.l3Mu_phi[ida]);
            histMap["recoHLT_BsMuP_dr"]->Fill(dra);
        }
        if(muM)
        {
            histMap["recoHLT_BsMuM_pt"]->Fill(ntupleRawTree.l3Mu_pt[idb]);
            histMap["recoHLT_BsMuM_eta"]->Fill(ntupleRawTree.l3Mu_eta[idb]);
            histMap["recoHLT_BsMuM_phi"]->Fill(ntupleRawTree.l3Mu_phi[idb]);
            histMap["recoHLT_BsMuM_dr"]->Fill(drb);
        }
        if(pho)
        {
            histMap["recoHLT_BsPho_pt" ]->Fill(ntupleRawTree.eg_et[idc]);
            histMap["recoHLT_BsPho_eta"]->Fill(ntupleRawTree.eg_eta[idc]);
            histMap["recoHLT_BsPho_phi"]->Fill(ntupleRawTree.eg_phi[idc]);
            histMap["recoHLT_BsPho_HForHoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]);
            histMap["recoHLT_BsPho_HoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]/ntupleRawTree.eg_energy[idc]);
            histMap["recoHLT_BsPho_NoiseCleanedSigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEtaNoise[idc]);
            histMap["recoHLT_BsPho_SigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEta[idc]);
            histMap["recoHLT_BsPho_R9" ]->Fill(ntupleRawTree.eg_r9Val[idc]);
            histMap["recoHLT_BsPho_dr" ]->Fill(drc);
            
            if(abs(ntupleRawTree.bsPho_eta) < ECAL_BARREL_ETA_MAX)
            {
                histMap["recoHLT_BsPho_BarrelPt" ]->Fill(ntupleRawTree.eg_et[idc]);
                histMap["recoHLT_BsPho_BarrelHForHoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]);
                histMap["recoHLT_BsPho_BarrelHoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]/ntupleRawTree.eg_energy[idc]);
                histMap["recoHLT_BsPho_BarrelNoiseCleanedSigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEtaNoise[idc]);
                histMap["recoHLT_BsPho_BarrelSigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEta[idc]);
                histMap["recoHLT_BsPho_BarrelR9" ]->Fill(ntupleRawTree.eg_r9Val[idc]);
                efficiencyMap["phoRecoHLTvsBarrelHForHoverE"]->fill(  pho , ntupleRawTree.eg_HoverE[idc]  );
                efficiencyMap["phoRecoHLTvsBarrelHoverE"]->fill(  pho , ntupleRawTree.eg_HoverE[idc]/ntupleRawTree.eg_energy[idc]  );
                efficiencyMap["phoRecoHLTvsBarrelR9"]->fill(      pho , ntupleRawTree.eg_r9Val[idc]  );
                efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta"]->fill( pho , ntupleRawTree.eg_sigmaIEtaIEtaNoise[idc]  );
                efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta"]->fill( pho , ntupleRawTree.eg_sigmaIEtaIEta[idc]  );
            }

            if( (abs(ntupleRawTree.bsPho_eta) >ECAL_ECAP_ETA_MIN  ) and  ( abs(ntupleRawTree.bsPho_eta) < ECAL_ECAP_ETA_MAX) )
            {
                histMap["recoHLT_BsPho_ECapPt" ]->Fill(ntupleRawTree.eg_et[idc]);
                histMap["recoHLT_BsPho_ECapHForHoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]);
                histMap["recoHLT_BsPho_ECapHoverE" ]->Fill(ntupleRawTree.eg_HoverE[idc]/ntupleRawTree.eg_energy[idc]);
                histMap["recoHLT_BsPho_ECapNoiseCleanedSigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEtaNoise[idc]);
                histMap["recoHLT_BsPho_ECapSigmaIEtaIEta" ]->Fill(ntupleRawTree.eg_sigmaIEtaIEta[idc]);
                histMap["recoHLT_BsPho_ECapR9" ]->Fill(ntupleRawTree.eg_r9Val[idc]);
                efficiencyMap["phoRecoHLTvsECapHForHoverE"]->fill(  pho , ntupleRawTree.eg_HoverE[idc]  );
                efficiencyMap["phoRecoHLTvsECapHoverE"]->fill(  pho , ntupleRawTree.eg_HoverE[idc] / ntupleRawTree.eg_energy[idc]  );
                efficiencyMap["phoRecoHLTvsECapR9"]->fill(      pho , ntupleRawTree.eg_r9Val[idc]  );
                efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta"]->fill( pho , ntupleRawTree.eg_sigmaIEtaIEtaNoise[idc]  );
                efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta"]->fill( pho , ntupleRawTree.eg_sigmaIEtaIEta[idc]  );
            }
        }

        dr=-1;
        if( muP and muM )
         {
            
            aLV.SetPtEtaPhiM(ntupleRawTree.l3Mu_pt[ida],ntupleRawTree.l3Mu_eta[ida],ntupleRawTree.l3Mu_phi[ida],MU_MASS);
            bLV.SetPtEtaPhiM(ntupleRawTree.l3Mu_pt[idb],ntupleRawTree.l3Mu_eta[idb],ntupleRawTree.l3Mu_phi[idb],MU_MASS);
            
            histMap["recoHLT_BsDiMu_mass" ]->Fill( ( aLV + bLV ).M() );  
            histMap["recoHLT_BsDiMu_pt" ]->Fill((aLV + bLV ).Pt() );  
            histMap["recoHLT_BsDiMu_eta"]->Fill((aLV + bLV ).Eta());
            histMap["recoHLT_BsDiMu_phi"]->Fill((aLV + bLV ).Phi());
            HLTDiMuPass = true;   
            HLTDiMuPassCount++;
            
            if(pho)
            {
                cLV.SetPtEtaPhiM(ntupleRawTree.eg_et[idb],ntupleRawTree.eg_eta[idb],ntupleRawTree.eg_phi[idb],PHO_MASS);
            
                histMap["recoHLT_Bs_mass" ]->Fill( ( aLV + bLV + cLV ).M() );  
                histMap["recoHLT_Bs_pt" ]->Fill((aLV + bLV + cLV ).Pt() );  
                histMap["recoHLT_Bs_eta"]->Fill((aLV + bLV + cLV ).Eta());
                histMap["recoHLT_Bs_phi"]->Fill((aLV + bLV + cLV ).Phi());
                HLTPass = true;   
                HLTPassCount++;
            

            }
         
         }

        //std::cout<<"HLTDiMu dr MuP= "<<dra<<"  MuM : "<<drb<<" : \t\tHLT DiMuRSLT : "<<HLTDiMuPass<<"\n";
        //std::cout<<"HLT dr Pho : "<<drc<<" : \t\t HLT RSLT : "<<HLTPass<<"\n";
         
        efficiencyMap["muonPRecoHLTvsPt"]->fill( muP , ntupleRawTree.bsMuP_pt  );
        efficiencyMap["muonPRecoHLTvsEta"]->fill( muP , ntupleRawTree.bsMuP_eta  );
        
        efficiencyMap["muonMRecoHLTvsPt"]->fill( muM , ntupleRawTree.bsMuM_pt  );
        efficiencyMap["muonMRecoHLTvsEta"]->fill( muM , ntupleRawTree.bsMuM_eta  );
        
        efficiencyMap["phoRecoHLTvsPt"]->fill( pho , ntupleRawTree.bsPho_pt  );
        efficiencyMap["phoRecoHLTvsEta"]->fill( pho , ntupleRawTree.bsPho_eta  );
        
        
        if(abs(ntupleRawTree.bsPho_eta) < ECAL_BARREL_ETA_MAX)
        {
         efficiencyMap["phoRecoHLTvsBarrelPt"]->fill( pho , ntupleRawTree.bsPho_pt  );
        }
        if( abs(ntupleRawTree.bsPho_eta) > ECAL_ECAP_ETA_MIN  and abs(ntupleRawTree.bsPho_eta) < ECAL_BARREL_ETA_MAX)
        {
         efficiencyMap["phoRecoHLTvsECapPt"]->fill( pho , ntupleRawTree.bsPho_pt  );
        }

        if(pho)
        {
        } 

        aLV.SetPtEtaPhiM(ntupleRawTree.bsMuP_pt,ntupleRawTree.bsMuP_eta,ntupleRawTree.bsMuP_phi,MU_MASS);
        bLV.SetPtEtaPhiM(ntupleRawTree.bsMuM_pt,ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuM_phi,MU_MASS);
         
         efficiencyMap["HLTDimuonMass"]->fill( HLTDiMuPass , (aLV + bLV).M() );
         auto ptmax = aLV.Pt() > bLV.Pt() ? aLV.Pt() : bLV.Pt() ;
         auto ptmin = aLV.Pt() < bLV.Pt() ? aLV.Pt() : bLV.Pt() ;
         efficiencyMap["HLTmu_pTMin"]->fill( HLTDiMuPass , ptmax );
         efficiencyMap["HLTmu_pTMax"]->fill( HLTDiMuPass , ptmin );
         efficiencyMap["dimuonHLT"]->fill( HLTDiMuPass , 1.0  );
        
	}
    

    // Making the Min/Max Threshold efficiencies
    efficiencyMap["phoRecoL1vsPt_minThr"]  	                =  efficiencyMap["phoRecoL1vsPt"]  	         ->Clone("phoRecoL1vsPt_minThr"); 
    efficiencyMap["muonPRecoHLTvsPt_minThr"]  	            =  efficiencyMap["muonPRecoHLTvsPt"]  	     ->Clone("muonPRecoHLTvsPt_minThr");
    efficiencyMap["muonMRecoHLTvsPt_minThr"]  	            =  efficiencyMap["muonMRecoHLTvsPt"]  	     ->Clone("muonMRecoHLTvsPt_minThr");
    efficiencyMap["phoRecoHLTvsPt_minThr"]  	            =  efficiencyMap["phoRecoHLTvsPt"]  	     ->Clone("phoRecoHLTvsPt_minThr");
    efficiencyMap["phoRecoHLTvsBarrelPt_minThr"]  	            =  efficiencyMap["phoRecoHLTvsBarrelPt"]  	     ->Clone("phoRecoHLTvsBarrelPt_minThr");
    efficiencyMap["phoRecoHLTvsECapPt_minThr"]  	            =  efficiencyMap["phoRecoHLTvsECapPt"]  	     ->Clone("phoRecoHLTvsECapPt_minThr");
    efficiencyMap["phoRecoHLTvsBarrelHForHoverE_maxThr"]        =  efficiencyMap["phoRecoHLTvsBarrelHForHoverE"]       ->Clone("phoRecoHLTvsBarrelHForHoverE_maxThr");
    efficiencyMap["phoRecoHLTvsECapHForHoverE_maxThr"]          =  efficiencyMap["phoRecoHLTvsECapHForHoverE"]       ->Clone("phoRecoHLTvsECapHForHoverE_maxThr");
    efficiencyMap["phoRecoHLTvsBarrelHoverE_maxThr"]        =  efficiencyMap["phoRecoHLTvsBarrelHoverE"]       ->Clone("phoRecoHLTvsBarrelHoverE_maxThr");
    efficiencyMap["phoRecoHLTvsECapHoverE_maxThr"]          =  efficiencyMap["phoRecoHLTvsECapHoverE"]       ->Clone("phoRecoHLTvsECapHoverE_maxThr");
    efficiencyMap["phoRecoHLTvsBarrelR9_minThr"] 	        =  efficiencyMap["phoRecoHLTvsBarrelR9"] 	         ->Clone("phoRecoHLTvsBarrelR9_minThr");
    efficiencyMap["phoRecoHLTvsECapR9_minThr"] 	            =  efficiencyMap["phoRecoHLTvsECapR9"] 	         ->Clone("phoRecoHLTvsECapR9_minThr");
    efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta_maxThr"] =  efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta"]->Clone("phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta_maxThr");
    efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta_maxThr"] =  efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta"]->Clone("phoRecoHLTvsBarrelSigmaIEtaIEta_maxThr");
    efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta_maxThr"] 	=  efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta"]->Clone("phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta_maxThr");
    efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta_maxThr"] 	=  efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta"]->Clone("phoRecoHLTvsECapSigmaIEtaIEta_maxThr");
    efficiencyMap["HLTDimuonMass_minThr"] 	                        =  efficiencyMap["HLTDimuonMass"] 	         ->Clone("HLTDimuonMass_minThr");
    efficiencyMap["HLTmu_pTMin_minThr"] 	                        =  efficiencyMap["HLTmu_pTMin"] 	         ->Clone("HLTmu_pTMin_minThr");
    efficiencyMap["HLTmu_pTMax_minThr"] 	                        =  efficiencyMap["HLTmu_pTMax"] 	         ->Clone("HLTmu_pTMax_minThr");

    efficiencyMap["phoRecoL1vsPt_minThr"]  	               ->sumForMinCut(); 
    efficiencyMap["muonPRecoHLTvsPt_minThr"]  	           ->sumForMinCut(); 
    efficiencyMap["muonMRecoHLTvsPt_minThr"]  	           ->sumForMinCut(); 
    efficiencyMap["phoRecoHLTvsPt_minThr"]  	           ->sumForMinCut(); 
    efficiencyMap["phoRecoHLTvsBarrelHForHoverE_maxThr"]       ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsBarrelHoverE_maxThr"]       ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta_maxThr"]->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta_maxThr"]->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsBarrelR9_minThr"] 	       ->sumForMinCut(); 
    efficiencyMap["phoRecoHLTvsECapHForHoverE_maxThr"]         ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsECapHoverE_maxThr"]         ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta_maxThr"]  ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta_maxThr"]  ->sumForMaxCut(); 
    efficiencyMap["phoRecoHLTvsECapR9_minThr"] 	           ->sumForMinCut(); 
    efficiencyMap["HLTDimuonMass_minThr"] 	                       ->sumForMinCut(); 
    efficiencyMap["HLTmu_pTMin_minThr"] 	                       ->sumForMinCut(); 
    efficiencyMap["HLTmu_pTMax_minThr"] 	                       ->sumForMinCut(); 

    efficiencyMap["phoRecoHLTvsBarrelHForHoverE_maxThr"]       ->divByTotalForEfficiency = true;  
    efficiencyMap["phoRecoHLTvsBarrelHoverE_maxThr"]       ->divByTotalForEfficiency = true;  
    efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta_maxThr"]->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta_maxThr"]->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsBarrelR9_minThr"] 	       ->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsECapHForHoverE_maxThr"]         ->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsECapHoverE_maxThr"]         ->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsECapNoiseCleanedSigmaIEtaIEta_maxThr"]  ->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsECapSigmaIEtaIEta_maxThr"]  ->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsECapR9_minThr"] 	           ->divByTotalForEfficiency = true;

    efficiencyMap["phoRecoHLTvsBarrelNoiseCleanedSigmaIEtaIEta_maxThr"]->divByTotalForEfficiency = true;
    efficiencyMap["phoRecoHLTvsBarrelSigmaIEtaIEta_maxThr"]->divByTotalForEfficiency = true;

    auto ofile=new TFile((prefix+ofileName).c_str(),"RECREATE");
    std::cout<<" Writing to file "<<ofile->GetName()<<"\n";
    ofile->cd();

    for (std::map<string,efficiencyMeasurement * >::iterator it=efficiencyMap.begin() ; it!=efficiencyMap.end(); ++it)
    {
        cout<<"Processing Efficiency "<<it->first<<"  \n";
        auto &curTrigger = *(it->second); 
        curTrigger.calculateEfficiency();
        curTrigger.Write();
    }
    
    for (std::map<string,TH1D * >::iterator it=histMap.begin() ; it!=histMap.end(); ++it)
    {
        cout<<"Processing Hist "<<it->first<<"  \n";
        auto &curTrigger = *(it->second); 
        curTrigger.Write();
    }
    
    cout<<"\n\n";
    cout<<" L1  stats : Total      :"<<L1PassCount<<" / "<<GenCount<<"  : "<<1.0*L1PassCount/GenCount<<"\n";
    
    cout<<"\n\n";
    cout<<" HLT DiMu stats : Total      :"<<HLTDiMuPassCount<<" / "<<GenCount<<"  : "<<1.0*HLTDiMuPassCount/GenCount<<"\n";
    cout<<"                  Individual :"<<HLTDiMuPassCount<<" / "<<L1PassCount<<"  : "<<1.0*HLTDiMuPassCount/L1PassCount<<"\n";

    cout<<"\n\n";
    cout<<" Photon Reco Stats after DiMuSelection  Individual :"<<HLTPassCount<<" / "<<HLTDiMuPassCount<<"  : "<<1.0*HLTPassCount/HLTDiMuPassCount<<"\n";
    
    cout<<"\n\n";
    cout<<" HLT stats : Total      :"<<HLTPassCount<<" / "<<GenCount<<"  : "<<1.0*HLTPassCount/GenCount<<"\n";
    cout<<"             Individual :"<<HLTPassCount<<" / "<<L1PassCount<<"  : "<<1.0*HLTPassCount/L1PassCount<<"\n";

    cout<<"\n\n";
   
   ofile->Purge();
   ofile->Close();
}


Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) {


    Double_t de = getDETA(eta1,eta2);
    Double_t dp = getDPHI(phi1,phi2); 
    return sqrt(de*de + dp*dp);

}


Double_t getDPHI( Double_t phi1, Double_t phi2) {

  Double_t dphi = phi1 - phi2;
  
  while( dphi > PI)
        dphi-= TWO_PI; 
  while( dphi < -PI)
        dphi += TWO_PI; 

  //std::cout<<"phi1  : "<<phi1<<" , phi2 : "<<phi2<<" dphi : "<<dphi<<"\n";
  
  if ( TMath::Abs(dphi) > 3.141592653589 ) {
    cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 "<< dphi << endl;
  }
  
  return TMath::Abs(dphi);
  //return dphi;
}



Double_t getDETA(Double_t eta1, Double_t eta2){
  return TMath::Abs(eta1 - eta2);
}

