import itertools

cfgTxt="\n\
#PARAMS_BEG\n\
InputFilePrefix=\n\
OutputPrefix=\n\
OutputFile=hltResult_@@IDX.txt\n\
MaxEvents=-1000\n\
DoGenSelection=@@ISMC\n\
MinAbsEtaL1Mu= @@MinAbsEtaL1Mu \n\
MaxAbsEtaL1Mu= @@MaxAbsEtaL1Mu \n\
MinPtL1Mu= @@MinPtL1Mu \n\
MaxDrL1MuMu= @@MaxDrL1MuMu  \n\
MaxDEtaL1MuMu= @@MaxDEtaL1MuMu \n\
MinDeltaRMuMu=     @@MinDeltaRMuMu      \n\
MinAbsEtaMaxMu=    @@MinAbsEtaMaxMu     \n\
MaxAbsEtaMaxMu=    @@MaxAbsEtaMaxMu     \n\
MinAbsEtaMinMu=    @@MinAbsEtaMinMu     \n\
MaxAbsEtaMinMu=    @@MaxAbsEtaMinMu     \n\
MinPtMinMu=        @@MinPtMinMu         \n\
MinPtMaxMu=        @@MinPtMaxMu         \n\
MinPtDimu=         @@MinPtDimu          \n\
MinInvMassDimuon=  @@MinInvMassDimuon   \n\
MaxInvMassDimuon=  @@MaxInvMassDimuon   \n\
MinFitProbablityDiMuVertex= 0.005\n\
MaxDeltaRPhoDimu= @@MaxDeltaRPhoDimu\n\
MinPtPho= @@MinPtPho\n\
MaxHoverEBarrelPho= @@MaxHoverEBarrelPho\n\
MaxHoverEECapPho= @@MaxHoverEECapPho\n\
MaxSigmaIEtaIEta5x5BarrelPho= @@MaxSigmaIEtaIEta5x5BarrelPho\n\
MaxSigmaIEtaIEta5x5ECapPho= @@MaxSigmaIEtaIEta5x5ECapPho\n\
MinR9BarrelPho=@@MinR9BarrelPho\n\
MinR9ECapPho=@@MinR9ECapPho\n\
MinMMGInvMass=@@MinMMGInvMass\n\
MaxMMGInvMass=@@MaxMMGInvMass\n\
#PARAMS_END\n\
#FILELIST_BEG\n\
@@FNAMES\n\
#FILELIST_END\n\
#EOF\n\
"

def getConfigTemplateForData():
    return cfgTxt.replace("@@ISMC",'0')

def getConfigTemplateForMC():
    return cfgTxt.replace("@@ISMC",'1')

def getPSetForPhotonID() :
    L1Sel=['L1A','L1B']
    MinDeltaRMuMu=  [1.4]  
    MinAbsEtaMaxMu= [0.0] 
    MaxAbsEtaMaxMu= [2.4] 
    MinAbsEtaMinMu= [0.0] 
    MaxAbsEtaMinMu= [2.4] 
    MinPtXMu=    [[4.0,3.5]]  
    MinPtDimu=     [4.9]  
    MinInvMassDimuon=[1.0]
    MaxInvMassDimuon=[6.0]
    MaxDeltaRPhoDimu=[0.7,0.9,1.1,1.3]
    MinPtPho=[4.0]
    
    MaxHoverEBarrelPho=[0.15,0.30 ,0.60,4.0]
    MaxHoverEECapPho = [0.10,0.24 ,0.50,3.2]
    MaxSigmaIEtaIEta5x5BarrelPho=[0.01  , 0.014 ,0.020,0.025]
    MaxSigmaIEtaIEta5x5ECapPho  =[0.030 , 0.035 ,0.040,0.05]
    MinR9BarrelPho=[ 0.4 ,0.4 ,0.4,0.5,0.6]
    MinR9ECapPho = [ 0.4 ,0.5 ,0.6,0.8,0.9]
    
    MaxHoverE=[[i,j] for i, j in zip(MaxHoverEBarrelPho,MaxHoverEECapPho)]
    MaxSigmaIEtaIEta5x5=[[i,j] for i, j in zip(MaxSigmaIEtaIEta5x5BarrelPho,MaxSigmaIEtaIEta5x5ECapPho)]
    MinR9=[[i,j] for i, j in zip(MinR9BarrelPho,MinR9ECapPho)]
    MinMMGInvMass=[4.0]
    MaxMMGInvMass=[6.5]
    MMGInv=[[i,j] for i, j in zip(MinMMGInvMass,MaxMMGInvMass)]
    
    dimuList=[MinDeltaRMuMu,MinPtXMu,MinAbsEtaMaxMu,MaxAbsEtaMaxMu,MinAbsEtaMinMu,MaxAbsEtaMinMu,MinPtDimu,MinInvMassDimuon,MaxInvMassDimuon]
    diMuPSet = list(itertools.product(*dimuList))
    
    photonleg=[L1Sel,diMuPSet,MaxDeltaRPhoDimu,MinPtPho,MMGInv,MaxHoverE,MaxSigmaIEtaIEta5x5,MinR9]
    
    total=1
    for l in photonleg:
        print(len(l)," * ",end="")
        total*=len(l)
    print( "   = ",total)
    parset = list(itertools.product(*photonleg))
    print("Parameterset made for : ",len(parset)," configs")

    return parset

def getPSetForL1AcceptanceStudy() :
    L1Sel=['L1A','L1B','L1C','L1D','L1E']
    #L1Sel=['L1E']
    MinDeltaRMuMu=  [1.4]  
    MinAbsEtaMaxMu= [0.0] 
    MaxAbsEtaMaxMu= [2.4] 
    MinAbsEtaMinMu= [0.0] 
    MaxAbsEtaMinMu= [2.4] 
    
    MinPtXMu=    [[3.5,3.5],[4.0,3.5],[4.0,4.0]]  
    MinPtDimu=     [4.9]  
    MinInvMassDimuon=[0.0]
    MaxInvMassDimuon=[6.0]
    MaxDeltaRPhoDimu=[1.4]
    MinPtPho=[4.0]
    
    MaxHoverEBarrelPho=[0.15 ,0.50]
    MaxHoverEECapPho = [0.10 ,0.40]
    MaxSigmaIEtaIEta5x5BarrelPho=[0.014 ,0.25]
    MaxSigmaIEtaIEta5x5ECapPho  =[0.035 ,0.05]
    MinR9BarrelPho=[0.5,0.5]
    MinR9ECapPho  =[0.8,0.5]
    
    MaxHoverE=[[i,j] for i, j in zip(MaxHoverEBarrelPho,MaxHoverEECapPho)]
    MaxSigmaIEtaIEta5x5=[[i,j] for i, j in zip(MaxSigmaIEtaIEta5x5BarrelPho,MaxSigmaIEtaIEta5x5ECapPho)]
    MinR9=[[i,j] for i, j in zip(MinR9BarrelPho,MinR9ECapPho)]
    MinMMGInvMass=[4.0,4.5]
    MaxMMGInvMass=[6.5,6.0]
    MMGInv=[[i,j] for i, j in zip(MinMMGInvMass,MaxMMGInvMass)]
    
    dimuList=[MinDeltaRMuMu,MinPtXMu,MinAbsEtaMaxMu,MaxAbsEtaMaxMu,MinAbsEtaMinMu,MaxAbsEtaMinMu,MinPtDimu,MinInvMassDimuon,MaxInvMassDimuon]
    diMuPSet = list(itertools.product(*dimuList))
    
    photonleg=[L1Sel,diMuPSet,MaxDeltaRPhoDimu,MinPtPho,MMGInv,MaxHoverE,MaxSigmaIEtaIEta5x5,MinR9]
    
    total=1
    for l in photonleg:
        print(len(l)," * ",end="")
        total*=len(l)
    print( "   = ",total)
    parset = list(itertools.product(*photonleg))
    print("Parameterset made for : ",len(parset)," configs")
    return parset




def getPSetForRateAndEfficJPsiWindow() :
    L1Sel=['L1A','L1B','L1E']
    MinDeltaRMuMu=  [1.4]  
    MinAbsEtaMaxMu= [0.0] 
    MaxAbsEtaMaxMu= [2.4] 
    MinAbsEtaMinMu= [0.0] 
    MaxAbsEtaMinMu= [2.4] 
    
    MinPtXMu=    [[3.5,3.5],[4.0,3.5],[4.0,4.0]]  
    MinPtDimu=     [4.9]  
    MinInvMassDimuon=[2.9]
    MaxInvMassDimuon=[3.3]
    MaxDeltaRPhoDimu=[1.4]
    MinPtPho=[4.0]
    
    MaxHoverEBarrelPho=[0.15 ,0.50]
    MaxHoverEECapPho = [0.10 ,0.40]
    MaxSigmaIEtaIEta5x5BarrelPho=[0.014 ,0.25]
    MaxSigmaIEtaIEta5x5ECapPho  =[0.035 ,0.05]
    MinR9BarrelPho=[0.5,0.5]
    MinR9ECapPho  =[0.8,0.5]
    
    
    MaxHoverE=[[i,j] for i, j in zip(MaxHoverEBarrelPho,MaxHoverEECapPho)]
    MaxSigmaIEtaIEta5x5=[[i,j] for i, j in zip(MaxSigmaIEtaIEta5x5BarrelPho,MaxSigmaIEtaIEta5x5ECapPho)]
    MinR9=[[i,j] for i, j in zip(MinR9BarrelPho,MinR9ECapPho)]
    MinMMGInvMass=[4.0,4.5]
    MaxMMGInvMass=[6.5,6.0]
    MMGInv=[[i,j] for i, j in zip(MinMMGInvMass,MaxMMGInvMass)]
    
    dimuList=[MinDeltaRMuMu,MinPtXMu,MinAbsEtaMaxMu,MaxAbsEtaMaxMu,MinAbsEtaMinMu,MaxAbsEtaMinMu,MinPtDimu,MinInvMassDimuon,MaxInvMassDimuon]
    diMuPSet = list(itertools.product(*dimuList))
    
    photonleg=[L1Sel,diMuPSet,MaxDeltaRPhoDimu,MinPtPho,MMGInv,MaxHoverE,MaxSigmaIEtaIEta5x5,MinR9]

    total=1
    for l in photonleg:
        print(len(l)," * ",end="")
        total*=len(l)
    print( "   = ",total)
    parset = list(itertools.product(*photonleg))
    print("Parameterset made for : ",len(parset)," configs")
    return parset


def getPSetForRateAndEffic() :
    L1Sel=['L1A','L1B','L1E']
    MinDeltaRMuMu=  [1.4]  
    MinAbsEtaMaxMu= [0.0] 
    MaxAbsEtaMaxMu= [2.4] 
    MinAbsEtaMinMu= [0.0] 
    MaxAbsEtaMinMu= [2.4] 
    
    MinPtXMu=    [[3.5,3.5],[4.0,3.5],[4.0,4.0]]  
    MinPtDimu=     [4.9]  
    MinInvMassDimuon=[0.0,1.0,1.5,2.0,2.4,2.6,2.9,3.3,4.5]
    MaxInvMassDimuon=[6.0]
    MaxDeltaRPhoDimu=[1.4]
    MinPtPho=[4.0]
    
    MaxHoverEBarrelPho=[0.15 ,0.50]
    MaxHoverEECapPho = [0.10 ,0.40]
    MaxSigmaIEtaIEta5x5BarrelPho=[0.014 ,0.25]
    MaxSigmaIEtaIEta5x5ECapPho  =[0.035 ,0.05]
    MinR9BarrelPho=[0.5,0.5]
    MinR9ECapPho  =[0.8,0.5]
    
    
    MaxHoverE=[[i,j] for i, j in zip(MaxHoverEBarrelPho,MaxHoverEECapPho)]
    MaxSigmaIEtaIEta5x5=[[i,j] for i, j in zip(MaxSigmaIEtaIEta5x5BarrelPho,MaxSigmaIEtaIEta5x5ECapPho)]
    MinR9=[[i,j] for i, j in zip(MinR9BarrelPho,MinR9ECapPho)]
    MinMMGInvMass=[4.0,4.25]
    MaxMMGInvMass=[6.5,6.25]
    MMGInv=[[i,j] for i, j in zip(MinMMGInvMass,MaxMMGInvMass)]
    
    dimuList=[MinDeltaRMuMu,MinPtXMu,MinAbsEtaMaxMu,MaxAbsEtaMaxMu,MinAbsEtaMinMu,MaxAbsEtaMinMu,MinPtDimu,MinInvMassDimuon,MaxInvMassDimuon]
    diMuPSet = list(itertools.product(*dimuList))
    
    photonleg=[L1Sel,diMuPSet,MaxDeltaRPhoDimu,MinPtPho,MMGInv,MaxHoverE,MaxSigmaIEtaIEta5x5,MinR9]

    total=1
    for l in photonleg:
        print(len(l)," * ",end="")
        total*=len(l)
    print( "   = ",total)
    parset = list(itertools.product(*photonleg))
    print("Parameterset made for : ",len(parset)," configs")
    return parset

def customizeConfigFile(custTxt,pset,maskJPsi):
    l1Code=pset[0]
    diMusel=pset[1]
    
    custTxt= 'JPsiMask='+str(maskJPsi)+'\n'+custTxt
    
    if l1Code=='L1A':
       custTxt= 'L1MODE=L1A\n'+custTxt
       custTxt=custTxt.replace('@@MinAbsEtaL1Mu','0.0')
       custTxt=custTxt.replace('@@MaxAbsEtaL1Mu','1.4')
       custTxt=custTxt.replace('@@MinPtL1Mu','0.0')
       custTxt=custTxt.replace('@@MaxDrL1MuMu','1.4')
       custTxt=custTxt.replace('@@MaxDEtaL1MuMu','1e9')

    if l1Code=='L1B':
       custTxt= 'L1MODE=L1B\n'+custTxt
       custTxt=custTxt.replace('@@MinAbsEtaL1Mu','0.0 0.0')
       custTxt=custTxt.replace('@@MaxAbsEtaL1Mu','1.4 3.0')
       custTxt=custTxt.replace('@@MinPtL1Mu','0.0 4.5')
       custTxt=custTxt.replace('@@MaxDrL1MuMu','1.4 1.2')
       custTxt=custTxt.replace('@@MaxDEtaL1MuMu','1e9 1e9')
    
    # L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6
    if l1Code=='L1C':
       custTxt= 'L1MODE=L1C\n'+custTxt
       custTxt=custTxt.replace('@@MinAbsEtaL1Mu','0.0')
       custTxt=custTxt.replace('@@MaxAbsEtaL1Mu','2.0')
       custTxt=custTxt.replace('@@MinPtL1Mu','0.0')
       custTxt=custTxt.replace('@@MaxDrL1MuMu','1e9')
       custTxt=custTxt.replace('@@MaxDEtaL1MuMu','1.6')
    
   # L1_DoubleMuOpen_er1p4_OS_dEta_Max1p6 OR L1_DoubleMu3_er2p0_SQ_OS_dR_Max1p4 OR  L1_DoubleMu4_SQ_OS_dR_Max1p2
    if l1Code=='L1D':
       custTxt= 'L1MODE=L1D\n'+custTxt
       custTxt=custTxt.replace('@@MinAbsEtaL1Mu','0.0 0.0 0.0')
       custTxt=custTxt.replace('@@MaxAbsEtaL1Mu','1.4 2.0 3.0')
       custTxt=custTxt.replace('@@MinPtL1Mu','0.0 3.0 4.0')
       custTxt=custTxt.replace('@@MaxDrL1MuMu','1e9 1.4 1.2')
       custTxt=custTxt.replace('@@MaxDEtaL1MuMu','1.6 1e9 1e9')

    if l1Code=='L1E':
       custTxt= 'L1MODE=L1E\n'+custTxt
       custTxt=custTxt.replace('@@MinAbsEtaL1Mu','0.0')
       custTxt=custTxt.replace('@@MaxAbsEtaL1Mu','1e9')
       custTxt=custTxt.replace('@@MinPtL1Mu','4.0')
       custTxt=custTxt.replace('@@MaxDrL1MuMu','1.2')
       custTxt=custTxt.replace('@@MaxDEtaL1MuMu','1e9')

    mul=1

    if maskJPsi and diMusel[7] <3.3:
        custTxt=custTxt.replace('@@MinInvMassDimuon',str(diMusel[7])+' 3.3')
        custTxt=custTxt.replace('@@MaxInvMassDimuon','2.9 '+str(diMusel[8]))
        mul=2
    else:
        custTxt=custTxt.replace('@@MinInvMassDimuon',str(diMusel[7]))
        custTxt=custTxt.replace('@@MaxInvMassDimuon',str(diMusel[8]))
        mul=1

    custTxt=custTxt.replace('@@MinDeltaRMuMu',( str(diMusel[0])+' ')*mul)
    custTxt=custTxt.replace('@@MinAbsEtaMaxMu',( str(diMusel[2])+' ')*mul)
    custTxt=custTxt.replace('@@MaxAbsEtaMaxMu',( str(diMusel[3])+' ')*mul)
    custTxt=custTxt.replace('@@MinAbsEtaMinMu',( str(diMusel[4])+' ')*mul)
    custTxt=custTxt.replace('@@MaxAbsEtaMinMu',( str(diMusel[5])+' ')*mul)
    custTxt=custTxt.replace('@@MinPtDimu',( str(diMusel[6])+' ')*mul)
    custTxt=custTxt.replace('@@MinPtMaxMu',( str(diMusel[1][0])+' ')*mul)
    custTxt=custTxt.replace('@@MinPtMinMu',( str(diMusel[1][1])+' ')*mul)
    
    custTxt=custTxt.replace('@@MaxDeltaRPhoDimu',str(pset[2]))
    custTxt=custTxt.replace('@@MinPtPho',str(pset[3]))
    custTxt=custTxt.replace('@@MinMMGInvMass',str(pset[4][0]))
    custTxt=custTxt.replace('@@MaxMMGInvMass',str(pset[4][1]))
    
    custTxt=custTxt.replace('@@MaxHoverEBarrelPho', str(pset[5][0]))
    custTxt=custTxt.replace('@@MaxHoverEECapPho', str(pset[5][1]))
    custTxt=custTxt.replace('@@MaxSigmaIEtaIEta5x5BarrelPho', str(pset[6][0]))
    custTxt=custTxt.replace('@@MaxSigmaIEtaIEta5x5ECapPho', str(pset[6][1]))
    custTxt=custTxt.replace('@@MinR9BarrelPho', str(pset[7][0]))
    custTxt=custTxt.replace('@@MinR9ECapPho', str(pset[7][1]))
    
    return custTxt
