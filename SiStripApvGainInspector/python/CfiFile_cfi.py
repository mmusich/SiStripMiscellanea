import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('SiStripApvGainInspector'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
