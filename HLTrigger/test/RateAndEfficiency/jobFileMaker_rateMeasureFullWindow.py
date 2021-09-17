#!/usr/bin/env python3
import os
import itertools
import sys
import parameterSetup as pSetup

print("sys.argv"," : ",sys.argv)
prefix=''
if(len(sys.argv) > 1):
    prefix=sys.argv[1]
prefix='run2FullWindow'
NJOBS=3000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=1
PU_NEVENTS=-1
pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']


maskJPsi=False
destination='/grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/efficencyANDrates/Run2FullWindow'
FNAME = "/grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/Run2HLT2018DFullILumiG1p6RateStudyNtuple_0.root"
cfgTxt=pSetup.getConfigTemplateForData()
parset=pSetup.getPSetForRateAndEffic()

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)cdr.stdout\n\
error = $Fp(filename)cdr.stderr\n\
log = $Fp(filename)cdr.log\n\
"
condorScript=open(prefix+'_job.sub','w')
condorScript.write(condorScriptString)

runScriptTxt="\
#!/bin/bash\n\
set -x\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME=/home/athachay\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd /home/athachay/t3store3/bs2mumug/run2studies/CMSSW_11_2_2_patch1/src\n\
cd "+pwd+"/@@DIRNAME \n\
cp /grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/*.h . \n\
cp /grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/rate_production.cc . \n\
cp /grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/rate_production.cc . \n\
cp /grid_mnt/t3storage3/athachay/bs2mumug/hltDev/CMSSW_11_3_0/src/Analysis/HLTAnalyserPy/Run3HLTNtuplizer/HLTNtupleTreeV3* . \n\
root -b -q 'rate_production.cc(\"@@CFGFILENAME\")'\n\
if [ $? -eq 0 ]; then \n\
    echo OK\n\
    cp hlt* "+destination+"\n\
    mv @@RUNSCRIPTNAME @@RUNSCRIPTNAME.sucess\n\
else\n\
    echo FAIL\n\
fi\n\
"

baseDir=prefix+'Jobs/'
if not os.path.exists(baseDir):
    os.system('mkdir '+baseDir)
print("Making ",NJOBS," Jobs ")

NJOBS = min(NJOBS,len(parset))
nid=int((NJOBS/10) + 1)

for ii in range(NJOBS):
    i=ii+ZERO_OFFSET
    
    dirName= baseDir+'Job_'+str(i)
    if(i%nid == 0):
        print(i," Job Made")
    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir '+dirName)
    cfgFileName='hltConfig_'+str(i)+".cfg"
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    cfgTxttmp=cfgTxt.replace("@@IDX",str(i))
    cfgTxttmp=cfgTxttmp.replace("@@FNAMES",FNAME)
    cfgTxttmp=pSetup.customizeConfigFile(cfgTxttmp,parset[ii],maskJPsi)
    cfgFile.write(cfgTxttmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/run'+str(i)+'.sh'
    runScriptNameFull=pwd+'/'+dirName+'/run'+str(i)+'.sh'
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@IDX",str(i))
    tmp=tmp.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    tmp=tmp.replace("@@RUNSCRIPTNAME",runScriptNameFull)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")

condorScript.close()
