import glob
import FWCore.ParameterSet.Config as cms
process = cms.Process("SiStripDigisAnalysis")

###################################################################
# Messages
###################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

###################################################################
# Source
###################################################################
readFiles = cms.untracked.vstring()
#readFiles.extend(['/store/relval/CMSSW_11_2_0_pre8/RelValPREMIXUP18_PU25/PREMIX/PU_112X_upgrade2018_realistic_v4-v1/00000/ffce1120-e325-409b-8382-0747045edf08.root'])

the_files=[]
file_list = glob.glob("/eos/cms/store/relval/CMSSW_11_2_0_pre8/RelValPREMIXUP18_PU25/PREMIX/PU_112X_upgrade2018_realistic_v4-v1/00000/*")
for f in file_list:
    the_files.append(f.replace("/eos/cms",""))
readFiles.extend(the_files)
print(the_files)
process.source = cms.Source("PoolSource",fileNames = readFiles)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

###################################################################
# Analyzer
###################################################################
process.myanalysis = cms.EDAnalyzer("SiStripDigisAnalyzer",
                                    src = cms.InputTag("simSiStripDigis","ZeroSuppressed")
                                    )

###################################################################
# Path
###################################################################
process.p1 = cms.Path(process.myanalysis)
