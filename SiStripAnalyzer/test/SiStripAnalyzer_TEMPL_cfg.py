import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
                            fileNames = readFiles
                            )

readFiles.extend(XXX_FILES_XXX)

###################################################################
# The track refitter
###################################################################
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'92X_dataRun2_Express_v4','')

###################################################################
# The BeamSpot
###################################################################
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

###################################################################
# The track refitter
###################################################################
from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone(src = cms.InputTag("ALCARECOSiStripCalMinBiasAAG"))

###################################################################
# The module
###################################################################
process.demo = cms.EDAnalyzer('SiStripAnalyzer',
                              lumiScalers = cms.InputTag("scalersRawToDigi"),
                              Tracks      = cms.InputTag("tracksRefit")
                              )

###################################################################
# Output name
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("SiStripInfo.root")
                                   )


process.p = cms.Path(process.offlineBeamSpot*process.MeasurementTrackerEvent*process.tracksRefit*process.demo)
