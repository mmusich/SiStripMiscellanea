import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptySource")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Express_v7', '')

process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.demo = cms.EDAnalyzer('SiStripApvGainInspector',
                              inputFile = cms.untracked.string("/tmp/musich/CMSSW_10_1_5/src/SiStripMiscellanea/SiStripApvGainInspector/test/DQM_V0001_R000999999__StreamExpress__Run2018B-PromptCalibProdSiStripGainsAAG-Express-v1-317279-317340__ALCAPROMPT.root")
)


process.p = cms.Path(process.demo)
