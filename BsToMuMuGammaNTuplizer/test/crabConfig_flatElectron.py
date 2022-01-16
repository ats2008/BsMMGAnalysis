from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'ntuples_flatElectron_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ntuplize_AOD.py'
config.JobType.maxMemoryMB = 4000
config.Data.inputDBS = 'phys03'
config.Data.inputDataset ='/DoubleElectronFlatPt1To20_GENSIM_Run2018/rchudasa-crab_DoubleElectronFlatPt1To20_RecoAOD_Run2018-648c701074ba72c4e16eb15dae8f7a8f/USER'

config.JobType.allowUndistributedCMSSW = True

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1

config.Data.outLFNDirBase = '/store/user/athachay/BsToMuMuGamma/Run2Studies/PhotonIDNtuples'
config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.outputDatasetTag = 'ntuples_defaultReco_flatElectron_pfIso'
config.Site.storageSite = 'T2_IN_TIFR'
