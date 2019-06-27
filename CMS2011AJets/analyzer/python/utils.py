from collections import OrderedDict
import os
import time

import energyflow as ef
import numpy as np

with open(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'DATAPATH'), 'r') as f:
    DATAPATH = f.read().strip()

DATASETS = {'cms': 'JetPrimaryDataset', 'sim': 'QCDSimDatasets'}
DATASET_PATHS = {k: os.path.join(DATAPATH, v) for k,v in DATASETS.items()}
FILEFORMATS = {'mod', 'jets'}

def get_filenames(path=None, dataset='cms', subdir='mod', ptmin=None, remove_ending=False, 
                  include_path=False, must_include=None):

    path = DATAPATH if path is None else path
    path = os.path.join(path, DATASETS[dataset], '' if subdir is None else subdir)

    if ptmin is not None:
            path = os.path.join(path, str(ptmin))

    raw_filenames = os.listdir(path)
    if must_include is None:
        filenames = sorted(raw_filenames)
    else:
        filenames = sorted(filter(lambda x: must_include in x, raw_filenames))

    if include_path:
        filenames = [os.path.join(path, filename) for filename in filenames]

    if remove_ending:
        filenames = [filename.split('.')[0] for filename in filenames]

    return filenames

def get_sim_filenames_dict(path=None, subdir='mod', remove_ending=False, include_path=False, must_include=None):

    # sort these according to numerical value
    sim_ptmins = sorted(get_filenames(path=path, dataset='sim', subdir=subdir), key=int)

    # get ordered dictionary of sim filenames
    sim_filenames = OrderedDict([(ptmin, get_filenames(path=path, dataset='sim', subdir=subdir, ptmin=ptmin, 
                                                       remove_ending=remove_ending, include_path=include_path,
                                                       must_include=must_include))
                                for ptmin in sim_ptmins])
    return sim_filenames

def separate_particle_arrays(particles, particles_index):
    
    # array to hold particles
    particles_array = np.zeros(len(particles_index) - 1, dtype='O')
    
    # iterate over indices
    for j,pi in enumerate(particles_index[:-1]):
        particles_array[j] = np.asarray(particles[pi:particles_index[j+1]])
        
    return particles_array

def concat_npz_files(filepaths, dataset, store_pfcs=True, store_gens=True):
    
    # check for validity
    assert dataset in {'cms', 'sim', 'gen'}, "dataset should be one of 'cms', 'sim', 'gen'"
    
    start = time.time()
    
    # arrays
    jets_float, jets_int = [], []
    pfcs, pfcs_index = [], []
    gens, gens_index = [], []
    
    nprev_pfcs = nprev_gens = 0
    
    # iterate over filepaths
    for i,filepath in enumerate(filepaths):
        
        # load
        f_npz = np.load(filepath)
        
        # check that we have any jets
        jets_int_i = f_npz['jets_int']
        if not len(jets_int_i):
            continue
        
        # append jets
        jets_float.append(f_npz['jets_float'])
        jets_int.append(jets_int_i)
        
        # pfcs
        if dataset != 'gen' and store_pfcs:
            pfcs.append(f_npz['pfcs'])
            pfcs_index.append(f_npz['pfcs_index'][:-1] + nprev_pfcs)
            nprev_pfcs += len(pfcs[-1])

            if len(pfcs[-1]) != f_npz['pfcs_index'][-1]:
                print('pfcs index error', filepath)
            
        # gens
        if dataset != 'cms' and store_gens:
            gens.append(f_npz['gens'])
            gens_index.append(f_npz['gens_index'][:-1] + nprev_gens)
            nprev_gens += len(gens[-1])

            if len(gens[-1]) != f_npz['gens_index'][-1]:
                print('gens index error', filepath)
            
        if (i+1) % 100 == 0:
            print('  [{}/{}] files processed in {:.3f}s'.format(i+1, len(filepaths), time.time() - start))
            
    # concatenate
    jets_float = np.concatenate(jets_float, axis=0)
    jets_int = np.concatenate(jets_int, axis=0)
    
    if dataset != 'gen' and store_pfcs:
        pfcs = np.concatenate(pfcs, axis=0)
        pfcs_index.append([nprev_pfcs])
        pfcs_index = np.concatenate(pfcs_index)
        
    if dataset != 'cms' and store_gens:
        gens = np.concatenate(gens, axis=0)
        gens_index.append([nprev_gens])
        gens_index = np.concatenate(gens_index)
        
    print('Finished and concatenated in {:.3f}s'.format(time.time() - start))
    
    return jets_float, jets_int, pfcs, pfcs_index, gens, gens_index

def filter_particles(particles, which='all', pt_cut=None, chs=False, pt_i=0, pid_i=4, vertex_i=5):
    
    mask = np.ones(len(particles), dtype=bool)
    
    # pt cut
    if pt_cut is not None:
        mask &= (particles[:,pt_i] >= pt_cut)
        
    # select specified particles
    if which != 'all':
        chrg_mask = ef.ischrgd(particles[:,pid_i])
        
        if which == 'charged':
            mask &= chrg_mask
        else:
            mask &= ~chrg_mask
            
    # apply chs
    if chs:
        mask &= (particles[:,vertex_i] <= 0)
        
    return mask