import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("Demo")

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.Services_cff")


process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('/afs/cern.ch/work/m/msommerh/public/test_samples/flatTuple.root')
                                   )


process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )
                     
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v10', '') #For upgrade2017
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v11', '')   #For ZprimeBB

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#"/store/user/msommerh/MC_ZPrime_to_BBar_2018_M4000/AOD/190415_071252/0000/AOD_ZPrime_to_BBar_2018_M1000_7.root" #signal test sample (2018, M=4000)
	"/store/user/msommerh/MC_QCD_2400to3200/AOD/190415_075003/0000/QCD_2400to3200_19.root"	#background test sample (QCD, pt: 2400 to 3200)
	)
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.demo = cms.EDAnalyzer('HitAnalyzer',
    Verbosity = cms.untracked.bool(False),
    phase1 = cms.untracked.bool(True),
    isMC = cms.untracked.bool(True),
    pT_cut = cms.untracked.double(200),
    src = cms.InputTag("siPixelClusters"),
)

process.p = cms.Path(process.demo)
