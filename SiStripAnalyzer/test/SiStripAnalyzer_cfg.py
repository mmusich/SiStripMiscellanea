import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:SiStripCalMinBias_oldNoise.root'
        #'file:SiStripCalMinBias.root'
        )
                            )

###################################################################
# The track refitter
###################################################################
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'102X_dataRun2_Prompt_v1','')
# process.GlobalTag.toGet = cms.VPSet(
#     cms.PSet(record = cms.string("SiStripNoisesRcd"),
#              tag = cms.string("myTag"),
#              connect = cms.string('sqlite_file:../../../myNoise.db')
#              )
#     )

###################################################################
# The BeamSpot
###################################################################
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

###################################################################
# The track refitter
###################################################################
from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone(src = cms.InputTag("ALCARECOSiStripCalMinBias"))

###################################################################
# The module
###################################################################
process.demo = cms.EDAnalyzer('SiStripClusterAnalyzer',
                              lumiScalers = cms.InputTag("scalersRawToDigi"),
                              Tracks      = cms.InputTag("tracksRefit"),
                              verbose_fit = cms.bool(True)
                              )

###################################################################
# Output name
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_oldNoise.root")
                                   #fileName = cms.string("test_newNoise.root")
                                   )


process.p = cms.Path(process.offlineBeamSpot*process.MeasurementTrackerEvent*process.tracksRefit*process.demo)
