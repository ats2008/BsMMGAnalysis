#!/bin/bash
set -x
#for Full menu 
#hltGetConfiguration /users/athachay/HLT_Bs2MMG_v0/V3 \

#for slimed menu
hltGetConfiguration /users/athachay/HLT_11_3_0_bmmgDev_v0/V26 \
   --setup /dev/CMSSW_11_3_0/GRun \
   --globaltag 113X_dataRun3_HLT_v1 \
   --input root://cms-xrd-global.cern.ch//eos/cms/store/data/Run2018A/HLTPhysics/RAW/v1/000/316/944/00000/E400BE7E-7E61-E811-BAA5-FA163E4E5866.root \
   --data \
   --process MYHLT \
   --full \
   --offline \
   --l1-emulator uGT \
   --l1 L1Menu_Collisions2018_v2_1_0-d1_xml \
   --path HLTriggerFirstPath,HLTriggerFinalPath,HLT_DoubleMu4_3_BsToMMG_v0,HLT_Ele32_WPTight_Gsf_v15 \
   --max-events 10 \
   --type GRun \
   --open \
   > hltMenuBs2MMGOnly_egamma.py 
set +x  
