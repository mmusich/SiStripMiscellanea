import FWCore.ParameterSet.Config as cms

process = cms.Process("SiStripGainValidator")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/34E2149B-3485-E711-AF48-02163E01A332.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D6564C87-3B85-E711-BB6B-02163E019C3B.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6421B77D-3B85-E711-9A58-02163E01A5BD.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/B4A0BF91-3B85-E711-B01F-02163E011E55.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/88309B7A-3B85-E711-9097-02163E019DCF.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/323534F0-3C85-E711-B196-02163E01A3CE.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/DEF43A19-3D85-E711-86B7-02163E014138.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6CB08228-3D85-E711-8C76-02163E011A01.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/AAE49908-3D85-E711-BEF4-02163E012A7E.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D414AF2F-4285-E711-9A79-02163E019CCB.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/C448402F-4285-E711-BD41-02163E01432C.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/828BFA2D-4285-E711-B18D-02163E01440E.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7E85E420-4285-E711-8885-02163E01A4F5.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/B0FF2920-4285-E711-86EA-02163E019CD0.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/C481C723-4285-E711-88B9-02163E01A3D2.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/D6CF070F-4385-E711-8CCD-02163E0141EA.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/14B185E5-4285-E711-99BE-02163E01A1E4.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/0A126D10-4385-E711-A9FB-02163E0142B1.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/96081603-4385-E711-85D2-02163E01A1D3.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/F888CC3D-4285-E711-A81F-02163E019DDD.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/4C04EC35-4285-E711-98EA-02163E01A674.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7C235249-4285-E711-A75F-02163E01262C.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/58BC536E-4285-E711-85AE-02163E01A20D.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/20036652-4285-E711-86BC-02163E01A2C3.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/22602AA2-4285-E711-A472-02163E011E2B.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/1CA79F25-4585-E711-942A-02163E0136F1.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/6A67EE36-4585-E711-80F6-02163E01A5CE.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/DEA89555-4785-E711-8014-02163E019CD2.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7E940E56-4785-E711-AAEC-02163E012510.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/5E048D5A-4785-E711-9D46-02163E01A5D1.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/98CB945B-4785-E711-80E8-02163E019DF3.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/52981554-4785-E711-8F2C-02163E01200E.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/940AB340-4A85-E711-8EBD-02163E0142EC.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7A3A5042-4A85-E711-BD60-02163E01473A.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/428E1863-4A85-E711-A5D9-02163E01A7A4.root',
'/store/express/Run2017C/StreamExpress/ALCARECO/SiStripCalMinBias-Express-v3/000/301/461/00000/7A44C144-4A85-E711-BC19-02163E019C13.root',
 )                            
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
process.GlobalTag = GlobalTag(process.GlobalTag,'92X_dataRun2_Express_v4','')
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
process.tracksRefit = TrackRefitter.clone(src = cms.InputTag("ALCARECOSiStripCalMinBias"))

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
