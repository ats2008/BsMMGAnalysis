#!/bin/bash
set -x
#for Full menu 
#hltGetConfiguration /users/athachay/HLT_Bs2MMG_v0/V3 \

#for slimed menu
hltGetConfiguration /users/athachay/HLT_11_3_0_bmmgDev_v0/V26 \
   --setup /dev/CMSSW_11_3_0/GRun \
   --globaltag 113X_dataRun3_HLT_v1 \
   --input root://cms-xrd-global.cern.ch//eos/cms/store/data/Run2018A/HLTPhysics/RAW/v1/000/316/944/00000/E400BE7E-7E61-E811-BAA5-FA163E4E5866.root \
   --mc \
   --process MYHLT \
   --full \
   --offline \
   --path HLTriggerFirstPath,HLTriggerFinalPath,HLT_DoubleMu4_3_BsToMMG_v0 \
   --max-events 10 \
   --type GRun \
   --open \
   > hltMenuBs2MMGOnly_egamma.py 
   
   #--path HLTriggerFirstPath,HLTriggerFinalPath,HLT_MyTriggerPath_v1 \

hltGetConfiguration /users/athachay/HLT_11_3_0_bmmgDev_v0/V26 \
   --setup /dev/CMSSW_11_3_0/GRun \
   --globaltag 113X_dataRun3_HLT_v1 \
   --input root://cms-xrd-global.cern.ch//eos/cms/store/data/Run2018A/HLTPhysics/RAW/v1/000/316/944/00000/E400BE7E-7E61-E811-BAA5-FA163E4E5866.root \
   --mc \
   --process MYHLT \
   --full \
   --offline \
   --max-events 10 \
   --type GRun\
   > hltMenuBs2MMG_egammaAndOthers.py 

set +x  
