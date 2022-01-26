from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_BsMMG_2018UL'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.inputDataset ='/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/athachay-crab_bs2mmg_RecoRun2UL18_v2-30c2008a2fc4ea0206d76ce8a5759bfe/USER'

config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 100
config.Data.outLFNDirBase = '/store/group/phys_bphys/rchudasa/'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_BsMMG_2018UL'
config.Site.storageSite = 'T2_CH_CERN'
