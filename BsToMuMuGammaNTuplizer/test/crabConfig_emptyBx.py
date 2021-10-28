from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_EmptyBX_onlySCInfo'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True


config.Data.inputDataset ='/EmptyBX/Run2018D-PromptReco-v2/AOD'
#config.Data.lumiMask = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt"
#config.Data.lumiMask = "/afs/cern.ch/work/r/rchudasa/private/hiforest_1034/CMSSW_10_3_4/src/HeavyIonsAnalysis/PhotonAnalysis/test/crab_projects/crab_ntuples_data/results/notFinishedLumis.json" #with seed info run on 13 Apr 2021 and unfinished jobs on 15th Apr 2021 
#config.Data.lumiMask = "/afs/cern.ch/work/j/jniedzie/private/lbl_ntuplizer/CMSSW_10_3_4/src/HeavyIonsAnalysis/PhotonAnalysis/test/crab_projects/crab_ntuples_data_lbl/results/notFinishedLumis.json"
config.Data.splitting = 'FileBased'
#config.Data.runRange = '326381-327564'
config.Data.unitsPerJob = 1

#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco/'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_EmptyBX_onlySCInfo'
config.Site.storageSite = 'T2_CH_CERN'
