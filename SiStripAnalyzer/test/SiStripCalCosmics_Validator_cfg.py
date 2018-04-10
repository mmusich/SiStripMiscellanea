import FWCore.ParameterSet.Config as cms

process = cms.Process("SiStripGainValidator")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring('/store/group/alca_trackeralign/musich/test_out/test_Commissioning2018_DATA_ReReco_SiStripCalCosmics/myTest_ReRecoCosmics18_SiStripCalCosmics_313052_1.root') 
                            )


from HLTrigger.HLTfilters.triggerResultsFilter_cfi import *
process.AAGFilter = triggerResultsFilter.clone(triggerConditions = cms.vstring("HLT_ZeroBias_FirstCollisionAfterAbortGap_*"),
                                               hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                               l1tResults = cms.InputTag( "" ),
                                               throw = cms.bool(False)
                                               )

###################################################################
# The track refitter
###################################################################
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,'101X_dataRun2_Express_v5','')
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("SiStripApvGain3Rcd"),
             tag = cms.string("SiStripGainFromParticles"),
             connect = cms.string("sqlite_file:/afs/cern.ch/user/m/musich/public/forStripDBcontacts/G2_initial2017/Gains_G2_299061_Sqlite.db")
             )
    )

###################################################################
# The BeamSpot
###################################################################
process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

###################################################################
# The track refitter
###################################################################
from RecoTracker.TrackProducer.TrackRefitter_cfi import TrackRefitter
process.load('RecoTracker.TrackProducer.TrackRefitters_cff')
process.tracksRefit = TrackRefitter.clone(src = cms.InputTag("ALCARECOSiStripCalCosmics"))

###################################################################
# The module
###################################################################
process.SiStripGainsValidator = cms.EDAnalyzer('SiStripGainsValidator',
                                               Tracks      = cms.InputTag("tracksRefit")
                                               )

###################################################################
# Output name
###################################################################
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("SiStripGainValidation.root")
                                   )

process.p = cms.Path(process.offlineBeamSpot*process.MeasurementTrackerEvent*process.tracksRefit*process.SiStripGainsValidator)
