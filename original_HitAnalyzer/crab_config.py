from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

datasets_bg = ["/MC_QCD_170to300/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_300to470/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_470to600/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_600to800/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_800to1000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_1000to1400/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_1400to1800/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_1800to2400/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_2400to3200/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_QCD_3200toInf/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER"]
datasets_2017 = ["/MC_ZPrime_to_BBar_2017_M1000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M1200/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER", 
		"/MC_ZPrime_to_BBar_2017_M1400/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M1600/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M1800/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M2000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M2500/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M3000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M3500/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M4000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M4500/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M5000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M5500/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER",
		"/MC_ZPrime_to_BBar_2017_M6000/msommerh-AOD-86331fa3e5403a58492e6b8f4cc1a457/USER"]
datasets_2018 = ["/MC_ZPrime_to_BBar_2018_M1000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M1200/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER", 
		"/MC_ZPrime_to_BBar_2018_M1400/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M1600/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M1800/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M2000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M2500/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M3000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M3500/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M4000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M4500/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M5000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M5500/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER",
		"/MC_ZPrime_to_BBar_2018_M6000/msommerh-AOD-2fd59cbde119ecab78af65e08efe8aae/USER"]
bg_bins = ['170to300', '300to470', '470to600', '600to800', '800to1000', '1000to1400', '1400to1800', '1800to2400', '2400to3200', '3200toInf']
signal_bins = ['1000', '1200', '1400', '1600', '1800', '2000', '2500', '3000', '3500', '4000', '4500', '5000', '5500', '6000']

data_nr = 9

config.General.requestName = 'updated_matching_nocuts_QCD_'+bg_bins[data_nr]+"_0"
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
#config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'config.py'

#config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 2
config.Data.ignoreLocality = True
config.Site.whitelist = ["T2_CH*", "T2_FR*", "T2_IT*", "T2_DE*", "T2_AT*", "T2_BE*", "T2_ES*"]

#config.Data.inputDataset = '/eos/cms/store/express/Run2018E/ExpressPhysics/FEVT/Express-v1/000/325/308/00000' #'/JetHT/Run2016G-PromptReco-v1/RECO'
config.Data.inputDataset = datasets_bg[data_nr]
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
#totalUnits = 3
config.Data.outLFNDirBase = '/store/user/msommerh/HitAnalyzer_updated_nocuts/QCD'
#config.Data.outLFNDirBase = '/afs/cern.ch/work/m/msommerh/public/MC_samples/ZPrime_to_BBar_M4000/flat_tuples/'
#config.Data.outLFNDirBase = '/t3home/msommerh/samples/'
config.Data.publication = False
config.Data.outputDatasetTag = 'updated_matching_nocuts_QCD_'+bg_bins[data_nr]

config.Site.storageSite = 'T3_CH_PSI'
