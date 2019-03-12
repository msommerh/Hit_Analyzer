from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'btagHits'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'config.py'

config.Data.inputDataset = '/eos/cms/store/express/Run2018E/ExpressPhysics/FEVT/Express-v1/000/325/308/00000' #'/JetHT/Run2016G-PromptReco-v1/RECO'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/msommerh/'# % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/t3home/msommerh/samples/'
config.Data.publication = False
config.Data.outputDatasetTag = 'btagHits_Ana_v3'

config.Site.storageSite = 'T3_CH_PSI'
