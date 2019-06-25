# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections

import FWCore.ParameterSet.Config as cms

# CMS process initialization
process = cms.Process('jecprocess')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# connect to global tag
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')
process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'

# setup JetCorrectorDBReader 
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))
process.source = cms.Source('EmptySource')
process.ak5 = cms.EDAnalyzer('JetCorrectorDBReader', 
                             payloadName=cms.untracked.string('AK5PF'),
                             globalTag=cms.untracked.string('FT_53_LV5_AN1'),  
                             printScreen=cms.untracked.bool(False),
                             createTextFile=cms.untracked.bool(True))
process.p = cms.Path(process.ak5)
