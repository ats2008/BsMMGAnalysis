void TreeMaker::Pi0ParticleSCMaker()
{
    AddSCTree("mergedPi0_SCTree");
    AddSCTree("leadPi0_SCTree");
    AddSCTree("subLeadPi0Gamma_SCTree");
    AddSCHistos("mergedPi0_");
    AddSCHistos("leadPi0GammaSC_");
    AddSCHistos("subLeadPi0GammaSC_");
    Double_t dr;
    
    std::cout<<"\nBegining Pi0 Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands(0),nMergedCands(0),nLeadGammaCands(0),nSubLeadGammaCands(0);
    Long64_t nCandsFromInclusion(0), nCandsFormSameSC(0);
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    bool goodRunLumi = false;

    Double_t drMin;
    Int_t scMatchIdx,scG1MatchIdx,scG2MatchIdx,tempI;
    Int_t g1MCIdx,g2MCIdx;
    Bool_t foundmatch=false;
    TLorentzVector g1,g2,pi0;
    Double_t g1g2DR,drMinG1,drMinG2,tempD;
    Bool_t isMerged;

    th1fStore["gen_p0daugterGammaGammaDR"] = new TH1F("gen_p0daugterGammaGammaDR","gen_p0daugterGammaGammaDR",100,0.0,2.0);
    th1fStore["gen_Pi0Gamma1DR"] = new TH1F("gen_Pi0Gamma1DR","gen_Pi0Gamma1DR",100,0.0,2.0);
    th1fStore["gen_Pi0Gamma2DR"] = new TH1F("gen_Pi0Gamma2DR","gen_Pi0Gamma2DR",100,0.0,2.0);
    

    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       eventGenMultiplicity=0;
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%10000 == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       EventCount++;
       foundmatch=false;
        
       if(isMC)
       for(int i=0;i < ntupleRawTree.nMC ; i++)
       {
           // std::cout<<"\tpdg id : "<<(ntupleRawTree.mcPID)->at(i)<<" pt : "<<(ntupleRawTree.mcPt)->at(i)<<" eta : "<<(ntupleRawTree.mcEta)->at(i)<<" phi : "<<ntupleRawTree.mcPhi->at(i)<<"\n";
            if(ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
            if(genParticleIsStable) if( ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                
            pi0.SetPtEtaPhiM(ntupleRawTree.mcPt->at(i) , ntupleRawTree.mcEta->at(i) , ntupleRawTree.mcPhi->at(i) , ntupleRawTree.mcMass->at(i) );
            drMin=1e9;
            for(int idx=0;idx < ntupleRawTree.nMC ; idx++)
            {
                if(ntupleRawTree.mcPID->at(idx) != 22 ) continue;
                g1.SetPtEtaPhiM(ntupleRawTree.mcPt->at(idx) , ntupleRawTree.mcEta->at(idx) , ntupleRawTree.mcPhi->at(idx) , ntupleRawTree.mcMass->at(idx) );
                // std::cout<<"\t\tpdg id : "<<(ntupleRawTree.mcPID)->at(idx)<<" pt : "<<(ntupleRawTree.mcPt)->at(idx)<<" eta : "<<(ntupleRawTree.mcEta)->at(idx)<<" phi : "<<ntupleRawTree.mcPhi->at(idx)<<"\n";
                for(int jdx=idx+1;jdx < ntupleRawTree.nMC ; jdx++)
                {
                    if(ntupleRawTree.mcPID->at(jdx) != 22 ) continue;
                    g2.SetPtEtaPhiM(ntupleRawTree.mcPt->at(jdx) , ntupleRawTree.mcEta->at(jdx) , ntupleRawTree.mcPhi->at(jdx) , ntupleRawTree.mcMass->at(jdx) );
                    dr=pi0.DeltaR(g1+g2);
                   // std::cout<<"\t\t\t\tpdg id : "<<(ntupleRawTree.mcPID)->at(jdx)<<" pt : "<<(ntupleRawTree.mcPt)->at(jdx)<<" eta : "<<(ntupleRawTree.mcEta)->at(jdx)<<" phi : "<<ntupleRawTree.mcPhi->at(jdx)<<" dr : "<<dr<<"\n";
                    if(dr<drMin)
                    {
                       drMin=dr;
                       g1MCIdx = idx; 
                       g2MCIdx = jdx; 
                    }
                    
                }
            }


            
            if( drMin < 0.001 ) {

              g1.SetPtEtaPhiM(ntupleRawTree.mcPt->at(g1MCIdx) , ntupleRawTree.mcEta->at(g1MCIdx) , ntupleRawTree.mcPhi->at(g1MCIdx) , ntupleRawTree.mcMass->at(g1MCIdx) );
              g2.SetPtEtaPhiM(ntupleRawTree.mcPt->at(g2MCIdx) , ntupleRawTree.mcEta->at(g2MCIdx) , ntupleRawTree.mcPhi->at(g2MCIdx) , ntupleRawTree.mcMass->at(g2MCIdx) );
              g1g2DR=g1.DeltaR(g2);
              th1fStore["gen_p0daugterGammaGammaDR"]->Fill(g1g2DR);
              th1fStore["gen_Pi0Gamma1DR"]->Fill(pi0.DeltaR(g1));
              th1fStore["gen_Pi0Gamma2DR"]->Fill(pi0.DeltaR(g2));
            }
            else
            {
                continue;
            }

 //           continue;
            
            drMinG1=drGenMatchMin;
            drMinG2=drGenMatchMin;
            scG1MatchIdx=-1;
            scG2MatchIdx=-1;
            
            isMerged=false;

            for( Int_t j =0 ;j< ntupleRawTree.nSC ;j++)
            {
                if(abs(ntupleRawTree.scEta->at(j)) < scAbsEtaMin ) continue;
                if(abs(ntupleRawTree.scEta->at(j)) > scAbsEtaMax ) continue;

                    dr=getDR(ntupleRawTree.mcEta->at(g1MCIdx),ntupleRawTree.mcPhi->at(g2MCIdx), ntupleRawTree.scEta->at(j),ntupleRawTree.scPhi->at(j));
                    if(dr<drMinG1)
                    {
                        drMinG1=dr;
                        scG1MatchIdx=j;
                    }
                    dr=getDR(ntupleRawTree.mcEta->at(g2MCIdx),ntupleRawTree.mcPhi->at(g2MCIdx), ntupleRawTree.scEta->at(j),ntupleRawTree.scPhi->at(j));
                    if(dr<drMinG2)
                    {
                        drMinG2=dr;
                        scG2MatchIdx=j;
                    }
            }
            if( (scG1MatchIdx==scG2MatchIdx) and (scG2MatchIdx !=-1) )
            {
                isMerged = true;
                nCandsFormSameSC++;
            }
            else
            {
                if( scG1MatchIdx > -1 )
                {
                    if(scG2MatchIdx > -1)
                    {
                        if( ntupleRawTree.scEt->at(scG1MatchIdx) < ntupleRawTree.scEt->at(scG2MatchIdx))
                        {
                            tempI = scG1MatchIdx;
                            scG1MatchIdx=scG2MatchIdx;
                            scG2MatchIdx=tempI;

                            tempD = drMinG1 ;
                            drMinG1 = drMinG2;
                            drMinG2= tempD;

                            tempI=g1MCIdx;
                            g1MCIdx=g2MCIdx;
                            g2MCIdx=tempI;
                        }
                        isMerged=false;
                    }
                }
                else
                {
                    if(scG2MatchIdx > -1)
                    {
                            tempI = scG1MatchIdx;
                            scG1MatchIdx=scG2MatchIdx;
                            scG2MatchIdx=tempI;

                            tempD = drMinG1 ;
                            drMinG1 = drMinG2;
                            drMinG2= tempD;

                            tempI=g1MCIdx;
                            g1MCIdx=g2MCIdx;
                            g2MCIdx=tempI;
                    }
                }

                if(scG1MatchIdx > -1 and scG2MatchIdx < 0 )
                {
                    if( 0.5*( sqrt( 
                                ntupleRawTree.scEtaWidth->at(scG1MatchIdx)*ntupleRawTree.scEtaWidth->at(scG1MatchIdx) 
                                + ntupleRawTree.scPhiWidth->at(scG1MatchIdx)*ntupleRawTree.scPhiWidth->at(scG1MatchIdx) 
                                ))  
                          > 
                           getDR(ntupleRawTree.scEta->at(scG1MatchIdx) , ntupleRawTree.scPhi->at(scG1MatchIdx) , ntupleRawTree.mcEta->at(g2MCIdx) , ntupleRawTree.mcPhi->at(g2MCIdx) ) 
                      ) 
                    {
                            nCandsFromInclusion++;
                           isMerged=true; 
                    }
                    else
                    {
                           isMerged=false;
                    }
                }
                else
                {
                      isMerged=false;
                }
                
            }

            if(isMerged) 
            { 
                fill_scHists(scG1MatchIdx,"mergedPi0_",drMinG1);
                fillSCVariablesToOutTree(scG1MatchIdx,"mergedPi0_SCTree");
                nMergedCands++;
                foundmatch=true;
            }
            else
            {
                if(scG1MatchIdx > -1 ) 
                {
                    fill_scHists(scG1MatchIdx,"leadPi0GammaSC_",drMinG1);
                    fillSCVariablesToOutTree(scG1MatchIdx,"leadPi0_SCTree");
                    nLeadGammaCands++;
                    foundmatch=true;
                }
                if(scG2MatchIdx > -1 ) 
                {
                    fill_scHists(scG2MatchIdx,"subLeadPi0GammaSC_",drMinG2);
                    fillSCVariablesToOutTree(scG2MatchIdx,"subLeadPi0Gamma_SCTree");
                    nSubLeadGammaCands++;
                    foundmatch=true;
                }
                if(scG1MatchIdx > -1 and scG2MatchIdx > -1 ) 
                {
                    nCands++;
                }
            }
        
            fill_genHists(i);
       }
       if(foundmatch==true)     EventCountWithCand++;
       
       for( Int_t j =0 ;j< ntupleRawTree.nSC ;j++)
       {
            if(abs(ntupleRawTree.scEta->at(j)) < scAbsEtaMin ) continue;
            if(abs(ntupleRawTree.scEta->at(j)) > scAbsEtaMax ) continue;

            drMin=2.99;
            if(isMC)
            for(int i=0;i < ntupleRawTree.nMC ; i++)
            {
                 if( ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
                 if(genParticleIsStable) if( ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                 dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.scEta->at(j),ntupleRawTree.scPhi->at(j));
                 if(dr<drMin)
                 {
                        drMin=dr;
                 }
            }
            fill_scHists(j,"allSC_",dr);
       }
      fill_eventHists();
    }

    std::cout<<" Number of Evnets processed                      : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates                : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of Events with merged candidates         : "<<nMergedCands<<"\n";
    std::cout<<"            single SC Matched candidates         : "<<nCandsFormSameSC<<"\n";
    std::cout<<"          Incusive SC Matched candidates         : "<<nCandsFromInclusion<<"\n";
    std::cout<<" Number of Events with lead candidates     [NM]  : "<<nLeadGammaCands<<"\n";
    std::cout<<" Number of Events with sub lead candidates [NM]  : "<<nSubLeadGammaCands<<"\n";
    std::cout<<" Number of Events with both candidates     [NM]  : "<<nCands<<"\n";
}

void TreeMaker::genParticleSCMaker()
{
    AddSCHistos("genMatchedSC_");
    AddSCTree("genMatchedSCTree");
    
    Double_t dr;
    
    std::cout<<"\nBegining Analysis Script !";
    if (maxEvents >0 ) maxEvents = nentries > maxEvents ? maxEvents : nentries;
    cout<<"\nProcessing total "<<maxEvents<<" events \n\n";
   
    Long64_t EventCount=0;
    Long64_t EventCountWithCand=0;
    Long64_t nCands=0;
    
    Long64_t nb = 0,nbytes=0 ;

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    bool goodRunLumi = false;

    Double_t drMin;
    Int_t scMatchIdx;
    Bool_t foundmatch=false;
    for (Long64_t jentry=0; jentry<maxEvents; jentry++)
    {  
       eventGenMultiplicity=0;
       Long64_t ientry_evt = ntupleRawTree.LoadTree(jentry);
       if (ientry_evt < 0) break;
       nb = ntupleRawTree.fChain->GetEntry(jentry);   nbytes += nb;
       
       if(jentry%10000 == 0 )
       {
             t_end = std::chrono::high_resolution_clock::now();
             std::cout<<"Processing Entry in event loop : "<<jentry<<" / "<<maxEvents<<"  [ "<<100.0*jentry/maxEvents<<"  % ]  "
                      << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                      <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents - jentry)/(1e-9 + jentry)* 0.001
                      <<std::endl;
       
       }
       
       EventCount++;
       foundmatch=false;
       if(isMC)
       for(int i=0;i < ntupleRawTree.nMC ; i++)
       {
         //   std::cout<<"\tpdg id : "<<(ntupleRawTree.mcPID)->at(i)<<" pt : "  \
                                    <<(ntupleRawTree.mcPt)->at(i)<<" eta : "<<(ntupleRawTree.mcEta)->at(i)<<" phi : "<<ntupleRawTree.mcPhi->at(i) \
                                    <<" status : "<<ntupleRawTree.mcStatus->at(i) \
                                    <<"\n";
            if(ntupleRawTree.mcPID->at(i) != genParticlePDGID )      continue;
            if(genParticleIsStable)   if( ntupleRawTree.mcStatus->at(i) != 1 )          continue;
            
            drMin=drGenMatchMin;
            scMatchIdx=-1;
            for( Int_t j =0 ;j< ntupleRawTree.nSC ;j++)
            {
                if(abs(ntupleRawTree.scEta->at(j)) < scAbsEtaMin ) continue;
                if(abs(ntupleRawTree.scEta->at(j)) > scAbsEtaMax ) continue;

                    dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.scEta->at(j),ntupleRawTree.scPhi->at(j));
                    if(dr<drMin)
                    {
                        drMin=dr;
                        scMatchIdx=j;
                    }
            }
            
            fill_genHists(i);
            if( scMatchIdx > -1)
            {
                fill_scHists(scMatchIdx,"genMatchedSC_",drMin);
                fillSCVariablesToOutTree(scMatchIdx,"genMatchedSCTree");
                nCands++;
                foundmatch=true;
            }
       }
       if(foundmatch==true)     EventCountWithCand++;
       
       for( Int_t j =0 ;j< ntupleRawTree.nSC ;j++)
       {
            if(abs(ntupleRawTree.scEta->at(j)) < scAbsEtaMin ) continue;
            if(abs(ntupleRawTree.scEta->at(j)) > scAbsEtaMax ) continue;

            drMin=2.99;
            if(isMC)
            for(int i=0;i < ntupleRawTree.nMC ; i++)
            {
                 if( ntupleRawTree.mcPID->at(i) != genParticlePDGID ) continue;
                 if(genParticleIsStable) if(not ntupleRawTree.mcStatus->at(i) != 1 ) continue;
                 dr=getDR(ntupleRawTree.mcEta->at(i),ntupleRawTree.mcPhi->at(i), ntupleRawTree.scEta->at(j),ntupleRawTree.scPhi->at(j));
                 if(dr<drMin)
                 {
                        drMin=dr;
                 }
            }
            fill_scHists(j,"allSC_",dr);
            

       }
      fill_eventHists();
    }

    std::cout<<" Number of Evnets processed        : "<<EventCount<<"\n";
    std::cout<<" Number of Evnets with candidates  : "<<EventCountWithCand<<"\n";
    std::cout<<" Number of candidates              : "<<nCands<<"\n";
}

void TreeMaker::setupOutputSCTree()
{
}


void TreeMaker::AddSCTree(TString SCTreeName)
{

    
    auto outSC_Tree = new TTree(SCTreeName,"variables for the PhotonID from SC");

    Int_t idx(0), offset(candidateMapDouble["SCTreeStorage"]);
	outSC_Tree->Branch("scE"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scEt"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scRawE"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scEta"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPhi"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scX"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scY"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scZ"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scEtaWidth"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPhiWidth"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scRawEt"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scMinDrWithGsfElectornSC_",&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFoundGsfMatch_"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE5x5"		            ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2x2Ratio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE3x3Ratio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scEMaxRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2ndRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scETopRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scERightRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scEBottomRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scELeftRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("sce2x5_MaxRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2x5TopRatio"		    ,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2x5RightRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2x5BottomRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scE2x5LeftRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scSwissCross"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scR9"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scsigmaIetaIeta"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scsigmaIetaIphi"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scsigmaIphiIphi"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e5x5"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x2Ratio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e3x3Ratio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_eMaxRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2ndRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_eTopRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_eRightRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_eBottomRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_eLeftRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x5MaxRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x5TopRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x5RightRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x5BottomRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_e2x5LeftRatio"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_swissCross"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_r9"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_sigmaIetaIeta"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_sigmaIetaIphi"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scFull5x5_sigmaIphiIphi"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFChIso1"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFChIso2"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFChIso"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFChIso4"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFChIso5"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFPhoIso1"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFPhoIso2"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFPhoIso"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFPhoIso4"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFPhoIso5"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFNeuIso1"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFNeuIso2"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFNeuIso"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFNeuIso4"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;
	outSC_Tree->Branch("scPFNeuIso5"		,&storageArrayDouble[ idx +  offset ]); idx+=1 ;

    treeStore[SCTreeName]=outSC_Tree;
}

void TreeMaker::fillSCVariablesToOutTree(Int_t scIDX,TString SCTreeName)
{
    Int_t idx(0) , offset(candidateMapDouble["SCTreeStorage"]);

	storageArrayDouble[idx + offset ] =ntupleRawTree.scE->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scEt->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scRawE->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scEta->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPhi->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scX->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scY->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scZ->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scEtaWidth->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPhiWidth->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scRawEt->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scMinDrWithGsfElectornSC_->at(scIDX)  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFoundGsfMatch_->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE5x5->at(scIDX)		              ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x2Ratio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE3x3Ratio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scEMaxRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2ndRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scETopRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scERightRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scEBottomRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scELeftRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x5MaxRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x5TopRatio->at(scIDX)		      ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x5RightRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x5BottomRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scE2x5LeftRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scSwissCross->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scR9->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scSigmaIetaIeta->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scSigmaIetaIphi->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scSigmaIphiIphi->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e5x5->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x2Ratio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e3x3Ratio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_eMaxRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2ndRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_eTopRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_eRightRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_eBottomRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_eLeftRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x5MaxRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x5TopRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x5RightRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x5BottomRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_e2x5LeftRatio->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_swissCross->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_r9->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_sigmaIetaIeta->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_sigmaIetaIphi->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scFull5x5_sigmaIphiIphi->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFChIso1->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFChIso2->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFChIso3->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFChIso4->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFChIso5->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFPhoIso1->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFPhoIso2->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFPhoIso3->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFPhoIso4->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFPhoIso5->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFNeuIso1->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFNeuIso2->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFNeuIso3->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFNeuIso4->at(scIDX)		  ;idx+=1 ;
	storageArrayDouble[idx + offset ] =ntupleRawTree.scPFNeuIso5->at(scIDX)		  ;idx+=1 ;
    
    treeStore[SCTreeName]->Fill();
}

void TreeMaker::bookHistograms()
{
     const int   pTBins=50;
     const float pTmin(0.0);
     const float pTmax(50.0);

     const int   etaBins=70;
     const float etamin(-3.5);
     const float etamax(3.5);

     const int   deltaRNBins=30;
     const float deltaRMin(0.0);
     const float deltaRMax(3.0);
     // All sc 
     for( TString tag : {"allSC_"} )
     {
        AddSCHistos(tag);
     }

     for( TString tag : {"gen_" })
     {
        th1fStore[tag+"Multiplicity"            ]  = new TH1F(tag + "Multiplicity" ,"Multiplicity", 200 , -0.5  , 199.5  );
        th1fStore[tag+"Pt"                      ]  = new TH1F(tag + "Pt" ,"Pt of GEN", pTBins , pTmin  , pTmax  );
        th1fStore[tag+"Eta"                     ]  = new TH1F(tag + "Eta","Eta of GEN", etaBins , etamin  , etamax  );
        th1fStore[tag+"phi"                     ]  = new TH1F(tag + "phi","Phi of GEN", 64 , -3.20  , 3.20  );
        th1fStore[tag+"E"		                ]  = new TH1F(tag + "E"  ,"E", 600 , 0.0 ,150.0 ) ;                  
     }
     
     th1fStore["sc_Multiplicity"            ]  = new TH1F("scMultiplicity" ,"Multiplicity", 200 , -0.5  , 199.5  );
}

void TreeMaker::AddSCHistos(TString tag)
{ 
     const int   pTBins=50;
     const float pTmin(0.0);
     const float pTmax(50.0);

     const int   etaBins=70;
     const float etamin(-3.5);
     const float etamax(3.5);

     const int   deltaRNBins=30;
     const float deltaRMin(0.0);
     const float deltaRMax(3.0);

        th1fStore[tag+"deltaR"                  ]  = new TH1F(tag + "deltaR","deltaR of sc Candidate", deltaRNBins , deltaRMin  , deltaRMax  );
        th1fStore[tag+"Eta"                     ]  = new TH1F(tag + "Eta","Eta of sc Candidate", etaBins , etamin  , etamax  );
        th1fStore[tag+"phi"                     ]  = new TH1F(tag + "phi","Phi of sc Candidate", 64 , -3.20  , 3.20  );
        th1fStore[tag+"E"		                ] = new TH1F(tag + "E"		                ,   "E"		         , 600 , 0.0 ,150.0 ) ;                  
	    th1fStore[tag+"Et"		                ] = new TH1F(tag + "Et"		                ,   "Et"		     , 250 , 0.0 ,50.0) ;
	    th1fStore[tag+"RawE"		            ] = new TH1F(tag + "RawE"		                ,   "RawE"		     , 600 , 0.0 ,150.0) ; 
	    th1fStore[tag+"X"		                ] = new TH1F(tag + "X"		                ,   "X"		    , 640 , -160.0,160.0) ;
	    th1fStore[tag+"Y"		                ] = new TH1F(tag + "Y"		                ,   "Y"		    , 640 , -160.0,160.0);   
	    th1fStore[tag+"Z"		                ] = new TH1F(tag + "Z"		                ,   "Z"		    , 800 , -400.0 ,400.0);
	    th1fStore[tag+"EtaWidth"		        ] = new TH1F(tag + "EtaWidth"		            ,   "EtaWidth"	, 180 , 0.0 ,0.045);	
	    th1fStore[tag+"PhiWidth"		        ] = new TH1F(tag + "PhiWidth"		            ,   "PhiWidth"	, 100 , 0.0 ,0.25);	
	    th1fStore[tag+"RawEt"		            ] = new TH1F(tag + "RawEt"		            ,   "RawEt"		, 250.0 , 0.0, 50.0);
	    th1fStore[tag+"MinDrWithGsfElectornSC_" ] = new TH1F(tag + "MinDrWithGsfElectornSC_"  ,   "MinDrWithGsfElectornSC_",2,0.0,1.0);
	    th1fStore[tag+"FoundGsfMatch_"		    ] = new TH1F(tag + "FoundGsfMatch_"		    ,   "FoundGsfMatch_",2,-0.5,1.5);		
	    th1fStore[tag+"E5x5"		            ] = new TH1F(tag + "E5x5"		                ,   "E5x5"		    ,480,0.0,120.0);       
	    th1fStore[tag+"E2x2Ratio"		        ] = new TH1F(tag + "E2x2Ratio"		        ,   "E2x2Ratio"		,120 ,-0.1,1.1) ;    
	    th1fStore[tag+"E3x3Ratio"		        ] = new TH1F(tag + "E3x3Ratio"		        ,   "E3x3Ratio"		,120 ,-0.1,1.1) ;      
	    th1fStore[tag+"EMaxRatio"		        ] = new TH1F(tag + "EMaxRatio"		        ,   "EMaxRatio"		,120 ,-0.1,1.1) ;    
	    th1fStore[tag+"E2ndRatio"		        ] = new TH1F(tag + "E2ndRatio"		        ,   "E2ndRatio"		,120 ,-0.1,1.1) ;    
	    th1fStore[tag+"ETopRatio"		        ] = new TH1F(tag + "ETopRatio"		        ,   "ETopRatio"		,120 ,-0.1,1.1) ;    
	    th1fStore[tag+"ERightRatio"		        ] = new TH1F(tag + "ERightRatio"		        ,   "ERightRatio"	,120 ,-0.1,1.1) ;	    
	    th1fStore[tag+"EBottomRatio"		    ] = new TH1F(tag + "EBottomRatio"		        ,   "EBottomRatio"	,120 ,-0.1,1.1) ;	    
	    th1fStore[tag+"ELeftRatio"		        ] = new TH1F(tag + "ELeftRatio"		        ,   "ELeftRatio"	,120 ,-0.1,1.1) ;	    
	    th1fStore[tag+"e2x5_MaxRatio"		    ] = new TH1F(tag + "e2x5_MaxRatio"		    ,   "e2x5_MaxRatio"	,120 ,-0.1,1.1) ;	    
	    th1fStore[tag+"E2x5TopRatio"		    ] = new TH1F(tag + "E2x5TopRatio"		        ,   "E2x5TopRatio"	,120 ,-0.1,1.1) ;	    
	    th1fStore[tag+"E2x5RightRatio"		    ] = new TH1F(tag + "E2x5RightRatio"		    ,   "E2x5RightRatio",120 ,-0.1,1.1) ;		
	    th1fStore[tag+"E2x5BottomRatio"		    ] = new TH1F(tag + "E2x5BottomRatio"		    ,   "E2x5BottomRatio",120 ,-0.1,1.1) ;		
	    th1fStore[tag+"E2x5LeftRatio"		    ] = new TH1F(tag + "E2x5LeftRatio"		    ,   "E2x5LeftRatio"	,120 ,-0.1,1.1) ;	
	    th1fStore[tag+"SwissCross"		        ] = new TH1F(tag + "SwissCross"		        ,   "SwissCross"	,120 ,-0.1,1.1) ;	
	    th1fStore[tag+"R9"		                ] = new TH1F(tag + "R9"		                ,   "R9"		, 420 ,-0.2,4.0) ;
	    th1fStore[tag+"sigmaIetaIeta"		    ] = new TH1F(tag + "sigmaIetaIeta"		    ,   "sigmaIetaIeta", 184 ,-0.001,0.045);  		
	    th1fStore[tag+"sigmaIetaIphi"		    ] = new TH1F(tag + "sigmaIetaIphi"		    ,   "sigmaIetaIphi", 800 , -0.0008, 0.0008);		
	    th1fStore[tag+"sigmaIphiIphi"		    ] = new TH1F(tag + "sigmaIphiIphi"		    ,   "sigmaIphiIphi", 142 ,-0.01 , 0.07 );		
	    th1fStore[tag+"Full5x5_e5x5"		    ] = new TH1F(tag + "Full5x5_e5x5"		        ,   "Full5x5_e5x5" , 600 , 0.0 , 150.0);		
	    th1fStore[tag+"Full5x5_e2x2Ratio"		] = new TH1F(tag + "Full5x5_e2x2Ratio"		,   "Full5x5_e2x2Ratio"		,120 ,-0.1,1.1) ;    
	    th1fStore[tag+"Full5x5_e3x3Ratio"		] = new TH1F(tag + "Full5x5_e3x3Ratio"		,   "Full5x5_e3x3Ratio"		,120 ,-0.1,1.1) ;
	    th1fStore[tag+"Full5x5_eMaxRatio"		] = new TH1F(tag + "Full5x5_eMaxRatio"		,   "Full5x5_eMaxRatio"		,120 ,-0.1,1.1) ;
	    th1fStore[tag+"Full5x5_e2ndRatio"		] = new TH1F(tag + "Full5x5_e2ndRatio"		,   "Full5x5_e2ndRatio"		,120 ,-0.1,1.1) ;
	    th1fStore[tag+"Full5x5_eTopRatio"		] = new TH1F(tag + "Full5x5_eTopRatio"		,   "Full5x5_eTopRatio"		,120 ,-0.1,1.1) ;
	    th1fStore[tag+"Full5x5_eRightRatio"	    ] = new TH1F(tag + "Full5x5_eRightRatio"		,   "Full5x5_eRightRatio"	,120 ,-0.1,1.1) ;	     	
	    th1fStore[tag+"Full5x5_eBottomRatio"	] = new TH1F(tag + "Full5x5_eBottomRatio"		,   "Full5x5_eBottomRatio"	,120 ,-0.1,1.1) ;	   	
	    th1fStore[tag+"Full5x5_eLeftRatio"		] = new TH1F(tag + "Full5x5_eLeftRatio"		,   "Full5x5_eLeftRatio"	,120 ,-0.1,1.1) ;	
	    th1fStore[tag+"Full5x5_e2x5MaxRatio"	] = new TH1F(tag + "Full5x5_e2x5MaxRatio"		,   "Full5x5_e2x5MaxRatio"	,120 ,-0.1,1.1) ;	   	
	    th1fStore[tag+"Full5x5_e2x5TopRatio"	] = new TH1F(tag + "Full5x5_e2x5TopRatio"		,   "Full5x5_e2x5TopRatio"	,120 ,-0.1,1.1) ;	   	
	    th1fStore[tag+"Full5x5_e2x5RightRatio"	] = new TH1F(tag + "Full5x5_e2x5RightRatio"	,   "Full5x5_e2x5RightRatio",120 ,-0.1,1.1) ;		   	
	    th1fStore[tag+"Full5x5_e2x5BottomRatio" ] = new TH1F(tag + "Full5x5_e2x5BottomRatio"	,   "Full5x5_e2x5BottomRatio",120 ,-0.1,1.1) ;		 		
	    th1fStore[tag+"Full5x5_e2x5LeftRatio"	] = new TH1F(tag + "Full5x5_e2x5LeftRatio"	,   "Full5x5_e2x5LeftRatio"	,120 ,-0.1,1.1) ;	     	
	    th1fStore[tag+"Full5x5_swissCross"		] = new TH1F(tag + "Full5x5_swissCross"		,   "Full5x5_swissCross"	,120 ,-0.1,1.1) ;	
	    th1fStore[tag+"Full5x5_r9"		        ] = new TH1F(tag + "Full5x5_r9"		        ,   "Full5x5_r9"		    , 420 , -0.2 , 4.0);
	    th1fStore[tag+"Full5x5_sigmaIetaIeta"	] = new TH1F(tag + "Full5x5_sigmaIetaIeta"	,   "Full5x5_sigmaIetaIeta" , 184 ,-0.001,0.045);  		  		     	
	    th1fStore[tag+"Full5x5_sigmaIetaIphi"	] = new TH1F(tag + "Full5x5_sigmaIetaIphi"	,   "Full5x5_sigmaIetaIphi"	, 800 , -0.0008, 0.0008);			     	
	    th1fStore[tag+"Full5x5_sigmaIphiIphi"	] = new TH1F(tag + "Full5x5_sigmaIphiIphi"	,   "Full5x5_sigmaIphiIphi"	, 142 ,-0.01 , 0.07 );		   	     	
	    th1fStore[tag+"PFChIso1"		        ] = new TH1F(tag + "PFChIso1"		            ,   "PFChIso1"	, 1000 , 0.0 , 0.20 ) ; 
	    th1fStore[tag+"PFChIso2"		        ] = new TH1F(tag + "PFChIso2"		            ,   "PFChIso2"	, 1000 , 0.0 , 0.20 ) ;	
	    th1fStore[tag+"PFChIso"		            ] = new TH1F(tag + "PFChIso3"		            ,   "PFChIso"	, 1000 , 0.0 , 0.20 ) ;	
	    th1fStore[tag+"PFChIso4"		        ] = new TH1F(tag + "PFChIso4"		            ,   "PFChIso4"	, 1000 , 0.0 , 0.20 ) ;	
	    th1fStore[tag+"PFChIso5"		        ] = new TH1F(tag + "PFChIso5"		            ,   "PFChIso5"	, 1000 , 0.0 , 0.20 ) ;	
	    th1fStore[tag+"PFPhoIso1"		        ] = new TH1F(tag + "PFPhoIso1"		        ,   "PFPhoIso1"	, 1000 , 0.0 , 20.0 ) ;	
	    th1fStore[tag+"PFPhoIso2"		        ] = new TH1F(tag + "PFPhoIso2"		        ,   "PFPhoIso2" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFPhoIso"		        ] = new TH1F(tag + "PFPhoIso3"	            ,   "PFPhoIso"  , 1000 , 0.0 , 20.0 ) ; 	
	    th1fStore[tag+"PFPhoIso4"		        ] = new TH1F(tag + "PFPhoIso4"		        ,   "PFPhoIso4" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFPhoIso5"		        ] = new TH1F(tag + "PFPhoIso5"		        ,   "PFPhoIso5" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFNeuIso1"		        ] = new TH1F(tag + "PFNeuIso1"		        ,   "PFNeuIso1" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFNeuIso2"		        ] = new TH1F(tag + "PFNeuIso2"		        ,   "PFNeuIso2" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFNeuIso"		        ] = new TH1F(tag + "PFNeuIso3"	            ,   "PFNeuIso3" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFNeuIso4"		        ] = new TH1F(tag + "PFNeuIso4"		        ,   "PFNeuIso4" , 1000 , 0.0 , 20.0 ) ;		
	    th1fStore[tag+"PFNeuIso5"		        ] = new TH1F(tag + "PFNeuIso5"		        ,   "PFNeuIso5" , 1000 , 0.0 , 20.0 ) ;		


}
void TreeMaker::fill_genHists(Int_t idx)
{
        TString tag="gen_";
        th1fStore[tag+"Pt"    ]->Fill(ntupleRawTree.mcPt->at(idx));   
        th1fStore[tag+"Eta"   ]->Fill(ntupleRawTree.mcEta->at(idx));  
        th1fStore[tag+"phi"   ]->Fill(ntupleRawTree.mcPhi->at(idx));  
        th1fStore[tag+"E"	  ]->Fill(ntupleRawTree.mcPhi->at(idx));  
}

void TreeMaker::fill_eventHists()
{
        th1fStore["gen_Multiplicity"    ]->Fill(eventGenMultiplicity);   
        th1fStore["sc_Multiplicity"    ]->Fill(ntupleRawTree.nSC);   
}

void TreeMaker::fill_scHists(Int_t scIDX,TString tag,Double_t dr)
{
        th1fStore[tag+"deltaR"                  ]->Fill( dr );
        th1fStore[tag+"Eta"                     ]->Fill(ntupleRawTree.scEta->at(scIDX)		      );
        th1fStore[tag+"phi"                     ]->Fill(ntupleRawTree.scPhi->at(scIDX)		      );
        th1fStore[tag+"E"		                ]->Fill(ntupleRawTree.scE->at(scIDX)		      );
	    th1fStore[tag+"Et"		                ]->Fill(ntupleRawTree.scEt->at(scIDX)		      );
	    th1fStore[tag+"RawE"		            ]->Fill(ntupleRawTree.scRawE->at(scIDX)		  );
	    th1fStore[tag+"X"		                ]->Fill(ntupleRawTree.scX->at(scIDX)		      );
	    th1fStore[tag+"Y"		                ]->Fill(ntupleRawTree.scY->at(scIDX)		      );
	    th1fStore[tag+"Z"		                ]->Fill(ntupleRawTree.scZ->at(scIDX)		      );
	    th1fStore[tag+"EtaWidth"		        ]->Fill(ntupleRawTree.scEtaWidth->at(scIDX)		  );
	    th1fStore[tag+"PhiWidth"		        ]->Fill(ntupleRawTree.scPhiWidth->at(scIDX)		  );
	    th1fStore[tag+"RawEt"		            ]->Fill(ntupleRawTree.scRawEt->at(scIDX)		  );
	    th1fStore[tag+"MinDrWithGsfElectornSC_" ]->Fill(ntupleRawTree.scMinDrWithGsfElectornSC_->at(scIDX)  );
	    th1fStore[tag+"FoundGsfMatch_"		    ]->Fill(ntupleRawTree.scFoundGsfMatch_->at(scIDX)		  );
	    th1fStore[tag+"E5x5"		            ]->Fill(ntupleRawTree.scE5x5->at(scIDX)		              );
	    th1fStore[tag+"E2x2Ratio"		        ]->Fill(ntupleRawTree.scE2x2Ratio->at(scIDX)		      );
	    th1fStore[tag+"E3x3Ratio"		        ]->Fill(ntupleRawTree.scE3x3Ratio->at(scIDX)		      );
	    th1fStore[tag+"EMaxRatio"		        ]->Fill(ntupleRawTree.scEMaxRatio->at(scIDX)		      );
	    th1fStore[tag+"E2ndRatio"		        ]->Fill(ntupleRawTree.scE2ndRatio->at(scIDX)		      );
	    th1fStore[tag+"ETopRatio"		        ]->Fill(ntupleRawTree.scETopRatio->at(scIDX)		      );
	    th1fStore[tag+"ERightRatio"		        ]->Fill(ntupleRawTree.scERightRatio->at(scIDX)		      );
	    th1fStore[tag+"EBottomRatio"		    ]->Fill(ntupleRawTree.scEBottomRatio->at(scIDX)		      );
	    th1fStore[tag+"ELeftRatio"		        ]->Fill(ntupleRawTree.scELeftRatio->at(scIDX)		      );
	    th1fStore[tag+"e2x5_MaxRatio"		    ]->Fill(ntupleRawTree.scE2x5MaxRatio->at(scIDX)		      );
	    th1fStore[tag+"E2x5TopRatio"		    ]->Fill(ntupleRawTree.scE2x5TopRatio->at(scIDX)		      );
	    th1fStore[tag+"E2x5RightRatio"		    ]->Fill(ntupleRawTree.scE2x5RightRatio->at(scIDX)		  );
	    th1fStore[tag+"E2x5BottomRatio"		    ]->Fill(ntupleRawTree.scE2x5BottomRatio->at(scIDX)		  );
	    th1fStore[tag+"E2x5LeftRatio"		    ]->Fill(ntupleRawTree.scE2x5LeftRatio->at(scIDX)		  );
	    th1fStore[tag+"SwissCross"		        ]->Fill(ntupleRawTree.scSwissCross->at(scIDX)		  );
	    th1fStore[tag+"R9"		                ]->Fill(ntupleRawTree.scR9->at(scIDX)		  );
	    th1fStore[tag+"sigmaIetaIeta"		    ]->Fill(ntupleRawTree.scSigmaIetaIeta->at(scIDX)		  );
	    th1fStore[tag+"sigmaIetaIphi"		    ]->Fill(ntupleRawTree.scSigmaIetaIphi->at(scIDX)		  );
	    th1fStore[tag+"sigmaIphiIphi"		    ]->Fill(ntupleRawTree.scSigmaIphiIphi->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e5x5"		    ]->Fill(ntupleRawTree.scFull5x5_e5x5->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x2Ratio"		]->Fill(ntupleRawTree.scFull5x5_e2x2Ratio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e3x3Ratio"		]->Fill(ntupleRawTree.scFull5x5_e3x3Ratio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_eMaxRatio"		]->Fill(ntupleRawTree.scFull5x5_eMaxRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2ndRatio"		]->Fill(ntupleRawTree.scFull5x5_e2ndRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_eTopRatio"		]->Fill(ntupleRawTree.scFull5x5_eTopRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_eRightRatio"	    ]->Fill(ntupleRawTree.scFull5x5_eRightRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_eBottomRatio"	]->Fill(ntupleRawTree.scFull5x5_eBottomRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_eLeftRatio"		]->Fill(ntupleRawTree.scFull5x5_eLeftRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x5MaxRatio"	]->Fill(ntupleRawTree.scFull5x5_e2x5MaxRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x5TopRatio"	]->Fill(ntupleRawTree.scFull5x5_e2x5TopRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x5RightRatio"	]->Fill(ntupleRawTree.scFull5x5_e2x5RightRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x5BottomRatio" ]->Fill(ntupleRawTree.scFull5x5_e2x5BottomRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_e2x5LeftRatio"	]->Fill(ntupleRawTree.scFull5x5_e2x5LeftRatio->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_swissCross"		]->Fill(ntupleRawTree.scFull5x5_swissCross->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_r9"		        ]->Fill(ntupleRawTree.scFull5x5_r9->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_sigmaIetaIeta"	]->Fill(ntupleRawTree.scFull5x5_sigmaIetaIeta->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_sigmaIetaIphi"	]->Fill(ntupleRawTree.scFull5x5_sigmaIetaIphi->at(scIDX)		  );
	    th1fStore[tag+"Full5x5_sigmaIphiIphi"	]->Fill(ntupleRawTree.scFull5x5_sigmaIphiIphi->at(scIDX)		  );
	    th1fStore[tag+"PFChIso1"		        ]->Fill(ntupleRawTree.scPFChIso1->at(scIDX)		  );
	    th1fStore[tag+"PFChIso2"		        ]->Fill(ntupleRawTree.scPFChIso2->at(scIDX)		  );
	    th1fStore[tag+"PFChIso"		            ]->Fill(ntupleRawTree.scPFChIso3->at(scIDX)		  );
	    th1fStore[tag+"PFChIso4"		        ]->Fill(ntupleRawTree.scPFChIso4->at(scIDX)		  );
	    th1fStore[tag+"PFChIso5"		        ]->Fill(ntupleRawTree.scPFChIso5->at(scIDX)		  );
	    th1fStore[tag+"PFPhoIso1"		        ]->Fill(ntupleRawTree.scPFPhoIso1->at(scIDX)		  );
	    th1fStore[tag+"PFPhoIso2"		        ]->Fill(ntupleRawTree.scPFPhoIso2->at(scIDX)		  );
	    th1fStore[tag+"PFPhoIso"		        ]->Fill(ntupleRawTree.scPFPhoIso3->at(scIDX)		  );
	    th1fStore[tag+"PFPhoIso4"		        ]->Fill(ntupleRawTree.scPFPhoIso4->at(scIDX)		  );
	    th1fStore[tag+"PFPhoIso5"		        ]->Fill(ntupleRawTree.scPFPhoIso5->at(scIDX)		  );
	    th1fStore[tag+"PFNeuIso1"		        ]->Fill(ntupleRawTree.scPFNeuIso1->at(scIDX)		  );
	    th1fStore[tag+"PFNeuIso2"		        ]->Fill(ntupleRawTree.scPFNeuIso2->at(scIDX)		  );             
	    th1fStore[tag+"PFNeuIso"		        ]->Fill(ntupleRawTree.scPFNeuIso3->at(scIDX)		  );             
	    th1fStore[tag+"PFNeuIso4"		        ]->Fill(ntupleRawTree.scPFNeuIso4->at(scIDX)		  );
	    th1fStore[tag+"PFNeuIso5"		        ]->Fill(ntupleRawTree.scPFNeuIso5->at(scIDX)		  );
}                                                       
                                                        
