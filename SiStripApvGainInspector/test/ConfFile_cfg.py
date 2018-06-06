import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.source = cms.Source("EmptyIOVSource",
                            firstValue = cms.uint64(317340),
                            lastValue  = cms.uint64(317340),
                            timetype  = cms.string('runnumber'),
                            interval = cms.uint64(1)
                            )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Express_v7', '')

process.load("Configuration.Geometry.GeometryRecoDB_cff")

####################################################################
# Output file
####################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string("APVGainsTree.root")
                                   )                                    

process.demo = cms.EDAnalyzer('SiStripApvGainInspector',
                              inputFile = cms.untracked.string("DQM_V0001_R000999999__StreamExpress__Run2018B-PromptCalibProdSiStripGainsAAG-Express-v1-317279-317340__ALCAPROMPT.root")
)


process.p = cms.Path(process.demo)
