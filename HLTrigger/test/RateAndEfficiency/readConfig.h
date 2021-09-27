
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
    //std::cout<<"Reading Params "<<"\n";
    
    int gg=0;

	while(std::getline(cfgFile,line))
	{
	   if(line=="#PARAMS_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#PARAMS_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
       strStream.clear();
       strStream.str(line);
       while (getline(strStream, field,'='))
       {
            if(field.compare("OutputFile")==0){
                 getline(strStream, field);
                 ofileName=field;
                 //std::cout<<" setting ofileName = "<<ofileName<<"\n";
            }
            if(field.compare("OutputPrefix")==0){
                 getline(strStream, field);
                 prefix=field;
                 cout<<" setting prefix = "<<prefix<<"\n";
            }
            if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("DoGenSelection")==0){
                 getline(strStream, field);
                 if(std::atoi(field.c_str())==1)  doGenSelection=1;
                 else doGenSelection = 0;
                 cout<<" setting doGenSelection  = "<<doGenSelection<<"\n";
            }
            if(field.compare("MinPtL1Mu")==0){
                 min_PtL1Mu.clear();
                 cout<<" setting min_PtL1Mu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_PtL1Mu.push_back(aDouble);
                 }
                 cout<<" }\n";
            } 
            if(field.compare("MinAbsEtaL1Mu")==0){
                 min_AbsEtaL1Mu.clear();
                 cout<<" setting min_AbsEtaL1Mu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_AbsEtaL1Mu.push_back(aDouble);
                 }
                 cout<<" }\n";
            } 
            if(field.compare("MaxAbsEtaL1Mu")==0){
                 max_AbsEtaL1Mu.clear();
                 cout<<" setting max_AbsEtaL1Mu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_AbsEtaL1Mu.push_back(aDouble);
                 }
                 cout<<" }\n";
            } 
            if(field.compare("MaxDEtaL1MuMu")==0){
                 max_DEtaL1MuMu.clear();
                 cout<<" setting max_DEtaL1MuMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_DEtaL1MuMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }  
            if(field.compare("MaxDrL1MuMu")==0){
                 max_DrL1MuMu.clear();
                 cout<<" setting max_DrL1MuMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_DrL1MuMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinPtMaxMu")==0){
                 min_PtMaxMu.clear();
                 cout<<" setting min_PtMaxMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_PtMaxMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinPtMinMu")==0){
                 min_PtMinMu.clear();
                 cout<<" setting min_PtMinMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_PtMinMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MaxAbsEtaMaxMu")==0){
                 min_PtMinMu.clear();
                 cout<<" setting max_AbsEtaMaxMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_AbsEtaMaxMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MaxAbsEtaMinMu")==0){
                 min_PtMinMu.clear();
                 cout<<" setting max_AbsEtaMinMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_AbsEtaMinMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            
            if(field.compare("MinAbsEtaMaxMu")==0){
                 min_PtMinMu.clear();
                 cout<<" setting min_AbsEtaMaxMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_AbsEtaMaxMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }}
            if(field.compare("MinAbsEtaMinMu")==0){
                 min_PtMinMu.clear();
                 cout<<" setting min_AbsEtaMinMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_AbsEtaMinMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinPtDimu")==0){
                 min_PtDimu.clear();
                 cout<<" setting min_PtDimu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_PtDimu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinDeltaRMuMu")==0){
                 max_DeltaRMuMu.clear();
                 cout<<" setting max_DeltaRMuMu  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_DeltaRMuMu.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinInvMassDimuon")==0){
                 min_InvMassDimuon.clear();
                 cout<<" setting min_InvMassDimuon  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    min_InvMassDimuon.push_back(aDouble);
                 }
                 cout<<" }\n";
            } 
            if(field.compare("MaxInvMassDimuon")==0){
                 max_InvMassDimuon.clear();
                 cout<<" setting max_InvMassDimuon  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_InvMassDimuon.push_back(aDouble);
                 }
                 cout<<" }\n";
            }
            if(field.compare("MinFitProbablityDiMuVertex")==0){
                 cout<<" setting min_FitProbablityDiMuVertex  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<"";
                    min_FitProbablityDiMuVertex=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MaxDeltaRPhoDimu")==0){
                 cout<<" setting max_DeltaRPhoDimu  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_DeltaRPhoDimu = aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MinPtPho")==0){
                 cout<<" setting min_PtPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    min_PtPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MaxHoverEBarrelPho")==0){
                 cout<<" setting min_HoverEBarrelPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_HoverEBarrelPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MaxHoverEECapPho")==0){
                 cout<<" setting max_HoverEECapPho  = { ";
                 while( strStream >> aDouble )
                 {
                    cout<<aDouble<<", ";
                    max_HoverEECapPho=aDouble;
                 }
                 cout<<" }\n";
            }

            if(field.compare("MaxSigmaIEtaIEta5x5BarrelPho")==0){
                 cout<<" setting max_SigmaIEtaIEta5x5BarrelPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    max_SigmaIEtaIEta5x5BarrelPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MaxSigmaIEtaIEta5x5ECapPho")==0){
                 cout<<" setting max_SigmaIEtaIEta5x5ECapPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    max_SigmaIEtaIEta5x5ECapPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MinR9BarrelPho")==0){
                 cout<<" setting min_R9BarrelPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    min_R9BarrelPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MinR9ECapPho")==0){
                 cout<<" setting min_R9ECapPho  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    min_R9ECapPho=aDouble;
                 }
                 cout<<"\n";
            }

            if(field.compare("MinMMGInvMass")==0){
                 cout<<" setting min_MMGInvMass  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    min_MMGInvMass=aDouble;
                 }
                 cout<<"\n";
            } 

            if(field.compare("MaxMMGInvMass")==0){
                 cout<<" setting max_MMGInvMass  =  ";
                 while( strStream >> aDouble )
                 {
                    cout<<" "<<aDouble;
                    max_MMGInvMass=aDouble;
                 }
                 cout<<"\n";
            }
       }
	
    }

    // Sanity check for configuration
    auto lenDimuParVec= min_PtMaxMu.size() ;
    auto l = lenDimuParVec;
    //std::cout<<"lenDimuParVec = "<<lenDimuParVec<<"\n";
    for(auto k=0; k<lenDimuParVec; k++)
         {
            std::cout<<min_PtMaxMu       .size() <<" / k = "<<k<<"   min_PtMaxMu[k]       : "<<" "<<min_PtMaxMu[k]      <<"\n";      
            std::cout<<min_PtMinMu       .size() <<" / k = "<<k<<"   min_PtMinMu[k]       : "<<" "<<min_PtMinMu[k]      <<"\n";      
            std::cout<<min_PtDimu        .size() <<" / k = "<<k<<"   min_PtDimu[k]        : "<<" "<<min_PtDimu[k]       <<"\n";       
            std::cout<<max_DeltaRMuMu    .size() <<" / k = "<<k<<"   max_DeltaRMuMu[k]    : "<<" "<<max_DeltaRMuMu[k]   <<"\n";   
            std::cout<<min_InvMassDimuon .size() <<" / k = "<<k<<"   min_InvMassDimuon[k] : "<<" "<<min_InvMassDimuon[k]<<"\n";
            std::cout<<max_InvMassDimuon .size() <<" / k = "<<k<<"   max_InvMassDimuon[k] : "<<" "<<max_InvMassDimuon[k]<<"\n";
         }
    if( l != min_PtMinMu.size()     )    {  cout<<"l = "<<l<<" min_PtMinMu.size() = "<<min_PtMinMu.size()<<"\n";             }
    if( l != min_PtDimu.size()      )    {  cout<<"l = "<<l<<" min_PtDimu.size() = "<<min_PtDimu.size()<<"\n";               }
    if( l != max_DeltaRMuMu.size()  )    {  cout<<"l = "<<l<<" max_DeltaRMuMu.size() = "<<max_DeltaRMuMu.size()<<"\n";       }
    if( l != min_InvMassDimuon.size())   {  cout<<"l = "<<l<<" min_InvMassDimuon.size() = "<<min_InvMassDimuon.size()<<"\n"; }
    if( l != max_InvMassDimuon.size())   {  cout<<"l = "<<l<<" max_InvMassDimuon.size() = "<<max_InvMassDimuon.size()<<"\n"; }


	// getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#FILELIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#FILELIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   //std::cout<<line<<"|\n";
	   InFileList.push_back(line);
	}

	//std::cout<<"Read File List with  "<<InFileList.size()<<" entries \n";
    
