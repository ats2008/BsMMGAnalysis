#!/usr/bin/env python3
import os
import sys,time
version='v1'
"""
    Usage
    ./makeCondorJobForAnalysis.py <EXECUTABLE> <InputFileListFname> <CFG_TEMPLATE> <analysisOption> <destination> <NJOBS> <FILES_PER_JOB> <jobPrefix>"

"""

def loadConfigTemplate(fname):
    f=open(fname,'r')
    txt=f.readlines()
    f.close()
    cfgString="".join(txt)
    return cfgString    


MVAsToTrain=['BDT','BDTG','KNN','MLP','SVM','BDTD','BDTB','BDTF','RuleFit','DNN_CPU']
executable='analysisNtupleMaker.exe'


NJOBS=20000
NEVENTS_PER_JOB = -1
ZERO_OFFSET=0
FILES_PER_JOB=40


pwd=os.environ['PWD']
proxy_path=os.environ['X509_USER_PROXY']
HOME=os.environ['HOME']
xrdRedirector="root://cms-xrd-global.cern.ch/"

FileSource ="bmm5FileList.txt"
destination='/grid_mnt/t3storage3/athachay/bs2mumug/run2studies/CMSSW_10_6_19_patch2/src/BsMMGAnalysis/MergeWithBMMNtuples/RunLumiEventFileMaker/runLumiList/'
cfgTemplateFname="analysis2018Temp.cfg"
tag=""
argC=1
hText="""
 ./mkaeJob.py <exe> <cfg template>  <destination> <tag>
"""
if len(sys.argv) > 1:
    executable=sys.argv[argC]  
else:
    print(hText)
    sys.exit(1)
argC+=1
if len(sys.argv) > argC:
    print(argC," : ",sys.argv[argC])
    cfgTemplateFname=sys.argv[argC]  
argC+=1
if len(sys.argv) > argC:
    print(argC," : ",sys.argv[argC])
    destination=sys.argv[argC]  
argC+=1
if len(sys.argv) > argC:
    print(argC," : ",sys.argv[argC])
    tag=sys.argv[argC]  

if(not os.path.exists(destination)):
    os.system("mkdir -p "+destination)
destination=os.path.abspath(destination)

print("Executable ",executable)
print("CFG template  file ",cfgTemplateFname)
print("destination : ",destination)
print("tag : ",tag)

configurationTxt= loadConfigTemplate(cfgTemplateFname)
configurationTxt= configurationTxt.replace("@@PWD",pwd)
configurationTxt= configurationTxt.replace("@@VERSION",version)
configurationTxt= configurationTxt.replace("@@TAG",tag)

condorScriptString="\
executable = $(filename)\n\
output = $Fp(filename)run.$(Cluster).stdout\n\
error = $Fp(filename)run.$(Cluster).stderr\n\
log = $Fp(filename)run.$(Cluster).log\n\
+JobFlavour = \"espresso\"\n\
"

runScriptTxt="\
#!/bin/bash\n\
set -x\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
export HOME="+HOME+"\n\
export X509_USER_PROXY="+proxy_path+"\n\
cd @@DIRNAME \n\
eval `scramv1 runtime -sh`\n\
TMPDIR=`mktemp -d`\n\
cd $TMPDIR\n\
cp  "+pwd+"/"+executable+" .\n\
cp @@DIRNAME/@@CFGFILENAME .\n\
mv @@RUNSCRIPT @@RUNSCRIPT.busy \n\
./"+executable+" @@CFGFILENAME \n\
if [ $? -eq 0 ]; then \n\
    rm *.exe\n\
    mkdir -p  @@DESTN\n\
    mv * @@DESTN\n\
    if [ $? -eq 0 ] ; then \n\
        mv @@CFGFILENAME @@DESTN\n\
        mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess \n\
        echo OK\n\
    else\n\
        mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
        echo FAIL\n\
    fi\n\
else\n\
    mv @@RUNSCRIPT.busy @@RUNSCRIPT \n\
    echo FAIL\n\
fi\n\
"

head='Condor/Jobs'+tag
if not os.path.exists(head ):
    os.system('mkdir -p '+head)

condorScriptName=head+'/job'+tag+'.sub'
condorScript=open(condorScriptName,'w')
condorScript.write(condorScriptString)


njobs=0
tm=int(time.time()*100)
for mva in MVAsToTrain:
    tm+=1
    dirName= pwd+'/'+head+'/Job_'+str(mva)
    
    if os.path.exists(dirName):
        k=True
    else:
        os.system('mkdir -p '+dirName)
    
    cfgFileName='MVATrainer_'+tag+'_'+mva+'.cfg'
    cfgFile=open(dirName+'/'+cfgFileName,'w')
    tmp=configurationTxt.replace("@@MVA",mva)
    cfgFile.write(tmp)
    cfgFile.close()   
    
    runScriptName=dirName+'/'+tag+'run'+mva+'.sh'
    if os.path.exists(runScriptName+'.sucess'):
       os.system('rm '+runScriptName+'.sucess')
    runScript=open(runScriptName,'w')
    tmp=runScriptTxt.replace("@@DIRNAME",dirName)
    tmp=tmp.replace("@@CFGFILENAME",cfgFileName)
    tmp=tmp.replace("@@RUNSCRIPT",runScriptName)
    dtn=destination+'/Results_'+tag+'_'+mva+'_'+str(tm)
    tmp=tmp.replace("@@DESTN",dtn)
    runScript.write(tmp)
    runScript.close()
    os.system('chmod +x '+runScriptName)
    condorScript.write("queue filename matching ("+runScriptName+")\n")
    njobs+=1
print()
print(" Number of jobs made : ", njobs)
condorScript.close()

print()
print("Job Files are in  : ",head)
print("Destination of the results are  : ",destination)
print("condor submit file is : ",condorScriptName)
