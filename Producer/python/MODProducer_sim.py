import json
import os
import sys
import time

import FWCore.ParameterSet.Config as cms
from FWCore.MessageService.MessageLogger_cfi import *

# grab options from command line
argc = len(sys.argv)
dataset_ptmin = sys.argv[2]
aod_filename  = sys.argv[3]
data_type     = 'sim'
global_tag    = 'START53_LV6A1'
data_dir      = '/media/sf_disk2/CMS2011A_QCDSimDatasets'

# check that data path exists
assert os.path.exists(data_dir), 'Path ' + data_dir + ' does not exist'

# make paths to subfolders
paths = {}
for subfolder in ['aod', 'mod', 'jec', 'log']:
    paths[subfolder] = (os.path.join(data_dir, subfolder) if subfolder == 'jec' else 
                        os.path.join(data_dir, subfolder, dataset_ptmin))
    assert os.path.exists(paths[subfolder]), 'Path ' + paths[subfolder] + ' does not exist'

# check that aod file exists
aod_filepath = os.path.join(paths['aod'], aod_filename)
assert os.path.exists(aod_filepath), 'File ' + aod_filepath + ' does not exist'

# get basename, which is aod_filename without the suffix
base_filename = aod_filename.split('.')[0]

# construct output filepath
mod_filepath = os.path.join(paths['mod'], base_filename + '.mod')

# construct log filepath
log_filepath = os.path.join(paths['log'], base_filename + '.log')

# begin log info
print '\nLog file:', log_filepath
with open(log_filepath, 'w') as f:
    f.write(' '.join(sys.argv) + '\n' +
            time.asctime(time.gmtime()) + '\n\n'
            'Data Type:    ' + data_type + '\n'
            'Log Filepath: ' + log_filepath + '\n'
            'AOD Filepath: ' + aod_filepath + '\n'
            'MOD Filepath: ' + mod_filepath + '\n'
            'JEC Path:     ' + paths['jec'] + '\n\n')

# make cms process
process = cms.Process('MODProduction')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# handle global tag 
process.load('Configuration.StandardSequences.Services_cff') # NOTE: this was not done in MODProducer
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/%s.db' % global_tag)
process.GlobalTag.globaltag = '%s::All' % global_tag

# setup source
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring('file:' + aod_filepath))

# -1 here means all events
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# MODProducer
process.MODProducer = cms.EDProducer('MODProducer',
                                      data_type=cms.string(data_type),
                                      jec_path=cms.string(os.path.join(paths['jec'], global_tag + '_')),
                                      mod_filepath=cms.string(mod_filepath),
                                      cert_str=cms.vstring([]),
                                      log_filepath=cms.string(log_filepath),
                                      print_every=cms.uint32(1000),
                                      min_jetpt_jet=cms.double(10.),
                                      min_jetpt_particles=cms.double(300.))

process.producer = cms.Path(process.MODProducer)
process.schedule = cms.Schedule(process.producer)
