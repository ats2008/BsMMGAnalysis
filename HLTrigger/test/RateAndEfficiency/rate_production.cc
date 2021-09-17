/*
 *
 *          USAGE
 *
 * root -b -q 'rate_production.cc("YOURCONFIG.CFG")'
 *
 *
 * */

#include "HLTNtupleTreeV3.h"
R__LOAD_LIBRARY(HLTNtupleTreeV3_C.so)

#include "efficiencyMeasurement.h"

#define HLT_PRESCALE 1100.0
#define LUMIS_LENGTH 23.31

#define MU_MASS 0.1056583755
#define PHO_MASS 0.0

#define L1_MATCHING_DR 0.3
#define ECAL_BARREL_ETA_MAX 1.48
#define ECAL_ECAP_ETA_MIN 1.48
#define ECAL_ECAP_ETA_MAX 3.0

const Double_t TWO_PI(3.141592653589*2) ;
const Double_t PI(3.141592653589)       ;

Double_t getDR( Double_t eta1, Double_t phi1,Double_t eta2 ,Double_t phi2) ;
Double_t getDPHI( Double_t phi1, Double_t phi2) ;
Double_t getDETA(Double_t eta1, Double_t eta2) ;

void produceTurnOns(string infile,string ofileName,string prefix="",Double_t DR_MAX_L1=1.4, Double_t ETA_RESTRIC_L1 = 1.54,Long64_t maxEvents=100000);

void rate_production(string fname)
{
	fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;
	
	std::vector<string> InFileList;
	std::vector<string> fSonList;
	
    Double_t aDouble;

    Long64_t maxEvents(10);
    string ofileName("output.root");
    string prefix("");

    // GenSelection
    bool doGenSelection;
    //L1 : Dimuon Selection
    vector<Double_t> min_AbsEtaL1Mu		{ 0.0 } ;
    vector<Double_t> max_AbsEtaL1Mu		{ 2.4 } ;
    vector<Double_t> min_PtL1Mu		{ 3.5 } ;
    vector<Double_t> max_DrL1MuMu	{ 1.2 } ;
    vector<Double_t> max_DEtaL1MuMu	{ 1e9 } ;

    //HLT :  Muon and Dimuon selectiuon
    vector<Double_t> min_PtMaxMu		{ 4.0 } ;
    vector<Double_t> min_PtMinMu		{ 3.0 } ;
    vector<Double_t> min_AbsEtaMaxMu		{ 0.0 } ;
    vector<Double_t> max_AbsEtaMaxMu		{ 2.4 } ;
    vector<Double_t> min_AbsEtaMinMu		{ 0.0 } ;
    vector<Double_t> max_AbsEtaMinMu		{ 2.4 } ;
    vector<Double_t> min_PtDimu		    { 3.0 } ;
    vector<Double_t> max_DeltaRMuMu	    { 3.0 } ;
    vector<Double_t> min_InvMassDimuon	{ 4.5 } ;
    vector<Double_t> max_InvMassDimuon	{ 6.0 } ;
    Double_t min_FitProbablityDiMuVertex (0.005 ) ;
    
    //HLT :  Photon Selection
    Double_t max_DeltaRPhoDimu	 { 1.4  } ;
    Double_t min_PtPho		     { 4.0  } ;
    Double_t max_HoverEBarrelPho { 0.15 } ;
    Double_t max_HoverEECapPho	 { 0.10 } ;
    Double_t max_SigmaIEtaIEta5x5BarrelPho { 0.14 } ;
    Double_t max_SigmaIEtaIEta5x5ECapPho { 0.035 } ;
    Double_t min_R9BarrelPho { 0.5 } ;
    Double_t min_R9ECapPho	 { 0.8 } ;
    
    //HLT : BsMMG candidate selection
    Double_t min_MMGInvMass	{ 4.5 } ;
    Double_t max_MMGInvMass	{ 6.0 } ;
    
    /**********************************************/
    
                    #include "readConfig.h"
    
    /**********************************************/
    
    
    // Trigger Logics 
    std::vector<TLorentzVector>  dimuonVectors;
    std::set<UInt_t> LumiSections;
    std::set<UInt_t>::iterator itrLumiSec;
    std::pair<std::set<UInt_t>::iterator,bool> lumi_it_pair;//=LumiSections.begin();
    TChain *treeChain = new TChain("egHLTRun3Tree");
    for(auto i=0;i<InFileList.size();i++)
    {
        treeChain->Add(InFileList[i].c_str());
    }
    

    HLTNtupleTreeV3 ntupleRawTree(treeChain);

	Long64_t nentries = treeChain->GetEntries();
    cout<<" Available total "<<nentries<<" \n";
    if (maxEvents >0 ) nentries = nentries > maxEvents ? maxEvents : nentries;
    cout<<" Processing total "<<nentries<<" \n";
    Long64_t nb = 0,nbytes=0 ;

    UInt_t nMus;
    UInt_t EventCount(0);
    UInt_t GenEventCount(0);
    UInt_t L1PassCount(0);
    UInt_t HLTDiMuPassCount(0);
    UInt_t HLTDiMuVtxPassCount(0);
    UInt_t HLTPassCount(0);
    Double_t dr;
    Double_t   dra=0.2; Int_t ida=-1;
    Double_t   drb=0.2; Int_t idb=-1;
    Double_t   drc=0.2; Int_t idc=-1;
    TLorentzVector aLV,bLV,cLV,dLV;

    bool    muP = false ;
    bool    muM = false ;
    bool    pho = false ;
    bool    muL1Pass = false ;
    bool HLTDiMuPass = false ;
    bool HLTDiMuVtxPass = false ;
    bool HLTPass = false ;
    bool tmpStat=false;
    
    std::set<UInt_t> l1SeedingMuonIdx;
    std::set<UInt_t>::iterator l1SeedingMuonItr;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++)
    {
       if(jentry%10000 ==0 )
       {
            cout<<"Processing jentry : "<<jentry<<"\n";
       }
	
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       EventCount++;
       
       if(doGenSelection)
       {
            tmpStat=true;
            if(ntupleRawTree.bsPho_pt <min_PtPho) tmpStat=false;

            if(not tmpStat) continue;
            
            dr = getDR(ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuP_phi,ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuM_phi);
            for(auto k=0;k < min_PtL1Mu.size() ;k++)
            {
              tmpStat=true;

              if( ntupleRawTree.bsMuP_pt < min_PtL1Mu[k] ) tmpStat =false;
              if( ntupleRawTree.bsMuM_pt < min_PtL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.bsMuP_eta) < min_AbsEtaL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.bsMuM_eta) < min_AbsEtaL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.bsMuP_eta) > max_AbsEtaL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.bsMuM_eta) > max_AbsEtaL1Mu[k] ) tmpStat =false;
     //         std::cout<<"abs(ntupleRawTree.bsMuP_eta - ntupleRawTree.bsMuM_eta) > max_DEtaL1MuMu[k]"
     //                   <<abs(ntupleRawTree.bsMuP_eta - ntupleRawTree.bsMuM_eta)<<" < "<<max_DEtaL1MuMu[k]<<"\n";
              
              if( abs(ntupleRawTree.bsMuP_eta - ntupleRawTree.bsMuM_eta) > max_DEtaL1MuMu[k] ) tmpStat =false;
              if( dr > max_DrL1MuMu[k] ) tmpStat =false;
              if(tmpStat) break;
            }
            
            if(not tmpStat) continue;
            
            aLV.SetPtEtaPhiM(ntupleRawTree.bsMuP_pt,ntupleRawTree.bsMuP_eta,ntupleRawTree.bsMuP_phi,MU_MASS);
            bLV.SetPtEtaPhiM(ntupleRawTree.bsMuM_pt,ntupleRawTree.bsMuM_eta,ntupleRawTree.bsMuM_phi,MU_MASS);
            
            for(auto k=0; k<lenDimuParVec; k++)
            {
                    tmpStat=true;
                    if ( ntupleRawTree.bsMuM_pt < ntupleRawTree.bsMuP_pt )
                    {
                        if (     ntupleRawTree.bsMuM_pt < min_PtMinMu[k] ) tmpStat = false;
                        if (     ntupleRawTree.bsMuP_pt < min_PtMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuM_eta) < min_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuM_eta) > max_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuP_eta) < min_AbsEtaMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuP_eta) > max_AbsEtaMaxMu[k] ) tmpStat = false;
                    }
                    else
                    {
                        if (     ntupleRawTree.bsMuP_pt < min_PtMinMu[k] ) tmpStat = false;
                        if (     ntupleRawTree.bsMuM_pt < min_PtMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuP_eta) < min_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuP_eta) > max_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuM_eta) < min_AbsEtaMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.bsMuM_eta) > max_AbsEtaMaxMu[k] ) tmpStat = false;
                    }
                    
                    if( dr > max_DeltaRMuMu[k] )   tmpStat = false ; 
                    if( not tmpStat ) continue;
                    
                    if( (aLV + bLV).Pt() < min_PtDimu[k] ) tmpStat = false ; 
                    
                    if( (aLV + bLV).M() < min_InvMassDimuon[k] ) tmpStat = false ; 
                    if( (aLV + bLV).M() > max_InvMassDimuon[k] ) tmpStat = false ; 
                    
                    if(tmpStat)
                    {
                        dimuonVectors.push_back(cLV);
                        break;
                    }
             }

            if(not tmpStat) continue;
           
       }
       GenEventCount++;
    
       lumi_it_pair=LumiSections.insert(ntupleRawTree.runnr*3000+ntupleRawTree.lumiSec); 
       
       muL1Pass = false;
       HLTDiMuPass = false;
       HLTDiMuVtxPass = false;
       HLTPass = false;
       
       // L1 SELECTION
       l1SeedingMuonIdx.clear();
       nMus=ntupleRawTree.nrAllL1Muons;
       for(auto i=0;i<nMus;i++)
       {
            //std::cout<<i<<"/"<<nMus<<" L1 Mu: "<<ntupleRawTree.allL1mu_charge[i]<<" : "<<ntupleRawTree.allL1mu_pt[i]<<" , "<<ntupleRawTree.allL1mu_eta[i]<<" , "<<ntupleRawTree.allL1mu_phi[i]<<"\n";
            for(auto k=0;k < min_PtL1Mu.size() ;k++)
            {
              tmpStat=true;
              if( ntupleRawTree.allL1mu_pt[i] < min_PtL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.allL1mu_eta[i]) < min_AbsEtaL1Mu[k] ) tmpStat =false;
              if( abs(ntupleRawTree.allL1mu_eta[i]) > max_AbsEtaL1Mu[k] ) tmpStat =false;
              if(tmpStat) break;
            }

            if(not tmpStat) continue;

            for(auto j=i+1;j<nMus;j++)
            {
                  if( ntupleRawTree.allL1mu_charge[i]*ntupleRawTree.allL1mu_charge[j] !=-1 ) continue;
                 dr = getDR(ntupleRawTree.allL1mu_eta[i],ntupleRawTree.allL1mu_phi[i],ntupleRawTree.allL1mu_eta[j],ntupleRawTree.allL1mu_phi[j]);
                //std::cout<<"\t\t -> Mu2 :"<<j<<" "<<ntupleRawTree.allL1mu_charge[j]<<","<<ntupleRawTree.allL1mu_pt[j]<<","<<ntupleRawTree.allL1mu_eta[i]<<","<<ntupleRawTree.allL1mu_phi[i]<<"  dr : "<<dr<<"\n";
                for(auto k=0;k < min_PtL1Mu.size() ;k++)
                {
                  tmpStat=true;
                  if( ntupleRawTree.allL1mu_pt[j] < min_PtL1Mu[k] ) tmpStat =false;
                  if( abs(ntupleRawTree.allL1mu_eta[j]) < min_AbsEtaL1Mu[k] ) tmpStat =false;
                  if( abs(ntupleRawTree.allL1mu_eta[j]) > max_AbsEtaL1Mu[k] ) tmpStat =false;
                  if( abs(ntupleRawTree.allL1mu_eta[j] - ntupleRawTree.allL1mu_eta[i]) > max_DEtaL1MuMu[k] ) tmpStat =false;
                  if( dr > max_DrL1MuMu[k] ) tmpStat =false;
                  if(tmpStat) break;
                }
              
              if(not tmpStat) continue;

             muL1Pass = true;
             l1SeedingMuonIdx.insert(i);
             l1SeedingMuonIdx.insert(j);
            }
            
      }

      //std::cout<<" muL1Pass : "<<muL1Pass<<"\n";
      
      if( not muL1Pass ) continue ;
      L1PassCount++; 
      
      // HLT Dimuon Selection

        //std::cout<<"      ==   ==    \n";
        //std::cout<<"gen muP   : pt "<<ntupleRawTree.bsMuP_pt<<" , eta "<<ntupleRawTree.bsMuP_eta<<" , phi "<<ntupleRawTree.bsMuP_phi<<" : "<<muP<<"\n";
        //std::cout<<"gen muM   : pt "<<ntupleRawTree.bsMuM_pt<<" , eta "<<ntupleRawTree.bsMuM_eta<<" , phi "<<ntupleRawTree.bsMuM_phi<<" : "<<muM<<"\n";
        //std::cout<<"gen pho   : pt "<<ntupleRawTree.bsPho_pt<<" , eta "<<ntupleRawTree.bsPho_eta<<" , phi "<<ntupleRawTree.bsPho_phi<<" : "<<pho<<"\n";
        for( l1SeedingMuonItr =l1SeedingMuonIdx.begin() ; l1SeedingMuonItr != l1SeedingMuonIdx.end() ; l1SeedingMuonItr++)
        {
             //std::cout<<"L1 Seed : "<<ntupleRawTree.allL1mu_charge[*l1SeedingMuonItr]<<" , "<<ntupleRawTree.allL1mu_pt[*l1SeedingMuonItr]<<" , "\
                    <<ntupleRawTree.allL1mu_eta[*l1SeedingMuonItr]<<","<<ntupleRawTree.allL1mu_phi[*l1SeedingMuonItr]<<"\n";
        }
        nMus=ntupleRawTree.nrL3Muons;
        dimuonVectors.clear();
        for(auto i=0;i<nMus;i++)
        {
               tmpStat=false;
               for( l1SeedingMuonItr =l1SeedingMuonIdx.begin() ; l1SeedingMuonItr != l1SeedingMuonIdx.end() ; l1SeedingMuonItr++)
               {
                    dr=getDR(ntupleRawTree.allL1mu_eta[*l1SeedingMuonItr],ntupleRawTree.allL1mu_phi[*l1SeedingMuonItr],
                             ntupleRawTree.l3Mu_eta[i],ntupleRawTree.l3Mu_phi[i]);
                    if(dr < L1_MATCHING_DR ) { tmpStat=true; break ;}

               }
            //std::cout<<"HLT L3Mu : "<<i<<"/"<<nMus<<" : "<<ntupleRawTree.l3Mu_charge[i]<<" : "<<ntupleRawTree.l3Mu_pt[i]<<" , "<<ntupleRawTree.l3Mu_eta[i]<<" , "<<ntupleRawTree.l3Mu_phi[i]<<" , L1 dr : "<<dr<<"\n";
               if( not tmpStat ) continue;
            
            for(auto j=i+1;j<nMus;j++)
            {
               tmpStat=false;
               for( l1SeedingMuonItr =l1SeedingMuonIdx.begin() ; l1SeedingMuonItr != l1SeedingMuonIdx.end() ; l1SeedingMuonItr++)
               {
                    dr=getDR(ntupleRawTree.allL1mu_eta[*l1SeedingMuonItr],ntupleRawTree.allL1mu_phi[*l1SeedingMuonItr],
                             ntupleRawTree.l3Mu_eta[j],ntupleRawTree.l3Mu_phi[j]);
                    if(dr < L1_MATCHING_DR ) { tmpStat=true; break ;}

               }
                if( not tmpStat ) continue;
                
                if(ntupleRawTree.l3Mu_charge[i]*ntupleRawTree.l3Mu_charge[j]!=-1) continue ;
                dr = getDR(ntupleRawTree.l3Mu_eta[i],ntupleRawTree.l3Mu_phi[i],ntupleRawTree.l3Mu_eta[j],ntupleRawTree.l3Mu_phi[j]) ;
                aLV.SetPtEtaPhiM(ntupleRawTree.l3Mu_pt[i],ntupleRawTree.l3Mu_eta[i],ntupleRawTree.l3Mu_phi[i],MU_MASS);
                bLV.SetPtEtaPhiM(ntupleRawTree.l3Mu_pt[j],ntupleRawTree.l3Mu_eta[j],ntupleRawTree.l3Mu_phi[j],MU_MASS);
                cLV = aLV + bLV;    
                //std::cout<<"   Mu2   : "<<j<<"/"<<nMus<<" : "<<ntupleRawTree.l3Mu_charge[j]<<" : "<<ntupleRawTree.l3Mu_pt[j]<<" , "<<ntupleRawTree.l3Mu_eta[j]<<" , "<<ntupleRawTree.l3Mu_phi[j]<<" di mu : pt "<<cLV.Pt()<<" , mass : "<<cLV.M()<<", dr : "<<dr;
               
                bool tmpStat = false     ;

                for(auto k=0; k<lenDimuParVec; k++)
                {
                    tmpStat=true;

                    
                    if ( ntupleRawTree.l3Mu_pt[i] < ntupleRawTree.l3Mu_pt[j] )
                    {
                        if ( ntupleRawTree.l3Mu_pt[i] < min_PtMinMu[k] ) tmpStat = false;
                        if ( ntupleRawTree.l3Mu_pt[j] < min_PtMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[i]) < min_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[i]) > max_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[j]) < min_AbsEtaMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[j]) > max_AbsEtaMaxMu[k] ) tmpStat = false;
                    }
                    if ( ntupleRawTree.l3Mu_pt[i] >= ntupleRawTree.l3Mu_pt[j] )
                    {
                        if ( ntupleRawTree.l3Mu_pt[j] < min_PtMinMu[k] ) tmpStat = false;
                        if ( ntupleRawTree.l3Mu_pt[i] < min_PtMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[j]) < min_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[j]) > max_AbsEtaMinMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[i]) < min_AbsEtaMaxMu[k] ) tmpStat = false;
                        if ( abs(ntupleRawTree.l3Mu_eta[i]) > max_AbsEtaMaxMu[k] ) tmpStat = false;
                    }
                   
                    
                   if( dr > max_DeltaRMuMu[k] )   tmpStat = false ; 
                   
                  if( not tmpStat ) continue;
                    
                    if( (cLV.Pt() < min_PtDimu[k]) ) tmpStat = false ; 
                    

                    if( (cLV.M() < min_InvMassDimuon[k]) ) tmpStat = false ; 
                    if( (cLV.M() > max_InvMassDimuon[k]) ) tmpStat = false ; 
                    
                    
                    if(tmpStat)
                    {
                        dimuonVectors.push_back(cLV);
                        break;
                    }
                }
                //std::cout<<dr<< "\n";
            }
        }
    
        
    if(dimuonVectors.size() > 0 )
    {
        HLTDiMuPass=true;
        HLTDiMuPassCount++;
    }

    //std::cout<<" HLTDiMuPass : "<<HLTDiMuPass<<"\n";

    if( not HLTDiMuPass ) continue ;

    // HLT DisplacedVertex Selection
    tmpStat=false;
    for(auto i=0 ; i<ntupleRawTree.nrDiMuVertex ; i++ )
    {
       if( TMath::Prob(ntupleRawTree.DiMuVertex_chi2[i],ntupleRawTree.DiMuVertex_ndof[i]) > min_FitProbablityDiMuVertex ) tmpStat=true;
    }

    if(not tmpStat ) continue;
    HLTDiMuVtxPass = true;
    HLTDiMuVtxPassCount++;
    
   // HLT Photon Selection

   for( auto i= 0 ; i<dimuonVectors.size() ; i++)
    {
       for(auto j=0;j<ntupleRawTree.nrEgs;j++)
       {
           if( ntupleRawTree.eg_et[j] < min_PtPho ) continue ;
           dr = getDR(  ntupleRawTree.eg_eta[j],ntupleRawTree.eg_phi[j], dimuonVectors[i].Eta(), dimuonVectors[i].Phi()  ) ;
           
           if( dr > max_DeltaRPhoDimu ) continue;
            
           if( abs(ntupleRawTree.eg_eta[j]) < ECAL_BARREL_ETA_MAX ) 
           {
                if( ntupleRawTree.eg_HoverE[j]/ntupleRawTree.eg_energy[j] > max_HoverEBarrelPho ) continue;
                //if( ntupleRawTree.eg_sigmaIEtaIEta[j] > max_SigmaIEtaIEta5x5BarrelPho   ) continue;
                if( ntupleRawTree.eg_sigmaIEtaIEtaNoise[j] > max_SigmaIEtaIEta5x5BarrelPho   ) continue;
                if( ntupleRawTree.eg_r9Val[j] < min_R9BarrelPho ) continue;
           }
           if( abs(ntupleRawTree.eg_eta[j]) > ECAL_ECAP_ETA_MIN and abs(ntupleRawTree.eg_eta[j]) < ECAL_ECAP_ETA_MAX ) 
           {
                if( ntupleRawTree.eg_HoverE[j]/ntupleRawTree.eg_energy[j] > max_HoverEBarrelPho ) continue;
                //if( ntupleRawTree.eg_sigmaIEtaIEta[j] > max_SigmaIEtaIEta5x5BarrelPho   ) continue;
                if( ntupleRawTree.eg_sigmaIEtaIEtaNoise[j] > max_SigmaIEtaIEta5x5BarrelPho   ) continue;
                if( ntupleRawTree.eg_r9Val[j] < min_R9ECapPho ) continue;
           }
           
           aLV.SetPtEtaPhiM(ntupleRawTree.eg_et[i],ntupleRawTree.eg_eta[i],ntupleRawTree.eg_phi[i],PHO_MASS);
           cLV = dimuonVectors[i] + aLV;
           
           if( cLV.M() < min_MMGInvMass or cLV.M() > max_MMGInvMass ) continue;

           HLTPass=true;
           HLTPassCount++;
           if(HLTPass) break; 
       }
       if(HLTPass) break; 
     }
     //std::cout<<" HLTPass : "<<HLTPass<<"\n";
  }

    auto nLumi= LumiSections.size();
    double factor = HLT_PRESCALE/(nLumi*LUMIS_LENGTH);
    double selection;
    cout<<"\n\n";
    cout<<" Runs Processed     : "<<" ** \n";
    cout<<" Lumisections Processed     : "<<nLumi<<"\n";
    cout<<" Events Processed     : "<<EventCount <<"\n";
    if(doGenSelection)
    {  
      cout<<" Events Passing Gen Selection   : "<<GenEventCount <<"\n";
      selection=GenEventCount;
    }
    else
    {
      selection=EventCount;
    }
    cout<<"\n";
    cout<<" Selection         : "<<selection          <<" \t\t [         L1 / Selection    = "<<selection*1.0/EventCount  <<" ] "<<"\n";
    cout<<" L1  Pass          : "<<L1PassCount        <<" \t\t [         L1 / Selection    = "<<L1PassCount*1.0/selection  <<" ] "<<"\n";
    cout<<" HLT DiMu Pass     : "<<HLTDiMuPassCount   <<" \t\t [    HLTDiMu / Selection    = "<<HLTDiMuPassCount*1.0/selection<<" ] "<<"\n";
    cout<<" HLT DiMuVtx Pass  : "<<HLTDiMuVtxPassCount<<" \t\t [ HLTDiMuVtx / Selection    = "<<HLTDiMuVtxPassCount*1.0/selection<<" ] "<<"\n";
    cout<<" HLT Pass          : "<<HLTPassCount       <<" \t\t [        HLT / Selection    = "<<HLTPassCount*1.0/selection<<" ] "<<"\n";
    cout<<"\n";
    cout<<" HLT Prescale         :  "<<HLT_PRESCALE<<"   "<<"\n";
    cout<<" Lumisection Length   :  "<<LUMIS_LENGTH<<" s "<<"\n";
    cout<<" L1 Rate              :  "<<L1PassCount*factor<<"   "<<"\n";
    cout<<" HLT DiMu Rate        :  "<<HLTDiMuPassCount*factor<<"   "<<"\n";
    cout<<" HLT MMG  Rate        :  "<<HLTPassCount*factor<<"   "<<"\n";
    cout<<"\n\n";
   

    
    fstream ofile((prefix+ofileName).c_str(),ios::out);
    ofile<<"\n\n";
    ofile<<" Runs Processed             : "<<" ** "<<"\n";
    ofile<<" Lumisections Processed     : "<<nLumi<<"\n";
    ofile<<" Events Processed           : "<<EventCount<<"\n";
    if(doGenSelection)
    {  
      ofile<<" Events Passing Gen Selection   : "<<GenEventCount <<"\n";
      selection=GenEventCount;
    }
    else
    {
      selection=EventCount;
    }
    ofile<<"\n";
    ofile<<" Selection         : "<<selection          <<"\t\t [         L1 / Selection    = "<<selection*1.0/EventCount  <<" ] "<<"\n";
    ofile<<" L1  Pass          : "<<L1PassCount        <<" \t\t [         L1 / Selection    = "<<L1PassCount*1.0/selection  <<" ] "<<"\n";
    ofile<<" HLT DiMu Pass     : "<<HLTDiMuPassCount   <<" \t\t [    HLTDiMu / Selection    = "<<HLTDiMuPassCount*1.0/selection<<" ] "<<"\n";
    ofile<<" HLT DiMuVtx Pass  : "<<HLTDiMuVtxPassCount<<" \t\t [ HLTDiMuVtx / Selection    = "<<HLTDiMuVtxPassCount*1.0/selection<<" ] "<<"\n";
    ofile<<" HLT Pass          : "<<HLTPassCount       <<" \t\t [        HLT / Selection    = "<<HLTPassCount*1.0/selection<<" ] "<<"\n";
    ofile<<"\n";
    ofile<<" HLT Prescale         :  "<<HLT_PRESCALE<<"   "<<"\n";
    ofile<<" Lumisection Length   :  "<<LUMIS_LENGTH<<" s "<<"\n";
    ofile<<" L1 Rate              :  "<<L1PassCount*factor<<"   "<<"\n";
    ofile<<" HLT DiMu Rate        :  "<<HLTDiMuPassCount*factor<<"   "<<"\n";
    ofile<<" HLT MMG  Rate        :  "<<HLTPassCount*factor<<"   "<<"\n";
    ofile<<"\n\n";
    
    ofile.close();
   

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

