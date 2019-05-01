import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.patTemplate_cfg import *

process = cms.Process("Demo")

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.Services_cff")


#Filters:

## HBHE Noise Filter
process.load("CommonTools.RecoAlgos.HBHENoiseFilter_cfi")
process.load("CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi")

## The good primary vertex filter ____________________________________________||
process.primaryVertexFilter = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
    filter = cms.bool(True)
    )

## The beam scraping filter __________________________________________________||
process.noscraping = cms.EDFilter(
    "FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
    )

## The CSC beam halo tight filter ____________________________________________||
#process.load('RecoMET.METFilters.CSCTightHaloFilter_cfi')

## The HCAL laser filter _____________________________________________________||
#process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
#process.hcalLaserEventFilter.vetoByRunEventNumber=cms.bool(False)
#process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.bool(True)

### The ECAL dead cell trigger primitive filter _______________________________||
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

### The EE bad SuperCrystal filter ____________________________________________||
process.load('RecoMET.METFilters.eeBadScFilter_cfi')


process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('flatTuple.root')
                                   )



process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )
                     
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '') #Doesnt work with pixel geom
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')  #For upgrade 2023 samples
# process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v10', '') #For upgrade2017
# process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_v3', '')   #For tt
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

filepath = "file:/eos/user/m/msommerh/MC_samples/QCD/"
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/RunIISpring16reHLT80/ZprimeToTTJet_M-4000_TuneCUETP8M1_13TeV-amcatnlo-pythia8/AODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/90000/00E593DA-7139-E611-86E4-0CC47A4D99A4.root')
    fileNames = cms.untracked.vstring(
	#'file:/eos/cms/store/express/Run2018E/ExpressPhysics/FEVT/Express-v1/000/325/308/00000/023F1B34-4E2E-A343-8E2C-09C411E86530.root'  #this data sample works
	#'file:/afs/cern.ch/work/t/thaarres/public/bTag_ntracks/CMSSW_9_4_0_patch1/src/bTag_nHits/HitAnalyzer/ZprimeBBbar_M4000_GENSIMDIGIRECO_1.root' #this MC file works
	#file_path+'AOD_ZPrime_to_BBar_M400_2018_0.root'
	#'file:/afs/cern.ch/user/m/msommerh/CMSSW_10_2_7/src/MCprodForBc2JPsilnu/QCD/AOD_QCD_GRID_test3200toInf_3200toInf_2.root'
	#'file:/afs/cern.ch/user/m/msommerh/CMSSW_10_2_7/src/MC_production/crab_projects/crab_MC_QCD_DR2_test_1/results/DR_step2_default_5to10_0_1.root'
	#'root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/msommerh/MC_samples/ZPrime_to_BBar_2018/M4000/AOD_ZPrime_to_BBar_2018_M4000_17.root'
	#"file:/eos/user/m/msommerh/MC_samples/ZPrime_to_BBar_2018/M4000/AOD_ZPrime_to_BBar_2018_M4000_17.root"
	"file:AOD_ZPrime_to_BBar_2018_M4000_1.root"
	#"file:AOD_ZPrime_to_BBar_2017_M4000_1.root"
	#"file:QCD_1800to2400_1.root"
	)
)

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.demo = cms.EDAnalyzer('HitAnalyzer',
    Verbosity = cms.untracked.bool(False),
    phase1 = cms.untracked.bool(True),
    isMC = cms.untracked.bool(True),
    pT_cut = cms.untracked.double(200),
    nJets_cut = cms.untracked.int32(2),
    leading_jet_eta = cms.untracked.double(2.5),
    loose_jets_cut = cms.untracked.bool(True),
    tight_jets_cut = cms.untracked.bool(True),
    MET_over_sumEt_cut = cms.untracked.double(0.5),
    src = cms.InputTag("siPixelClusters"),
    HLTtriggers = cms.InputTag("TriggerResults", "", "HLT")
)

process.p = cms.Path(
	process.noscraping *
	process.HBHENoiseFilterResultProducer * 
	process.HBHENoiseFilter * 
	#process.CSCTightHaloFilter *
	#process.hcalLaserEventFilter *
	process.EcalDeadCellTriggerPrimitiveFilter *
	process.eeBadScFilter *
	process.primaryVertexFilter *
	process.demo)
#process.p = cms.Path(process.demo)
