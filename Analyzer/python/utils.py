from collections import OrderedDict
import os
import time

import energyflow as ef
import numpy as np

PATHS = {'cms': '/mnt/disk0/CMS2011A/JetPrimaryDataset', 
         'sim': '/mnt/disk0/CMS2011A/QCDSimDatasets'}

DATASETS = {'cms': 0, 'sim': 1}
FILEFORMATS = {'umod', 'jets'}

def get_filenames(path=None, dataset='cms', subdir='umod', ptmin=None, remove_ending=False, 
                  include_path=False, must_include=None, **kwargs):
    
    if path is None:
        path = os.path.join(PATHS[dataset], '' if subdir is None else subdir)

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
        return remove_endings(filenames)

    return filenames

def remove_endings(filename_array):
    return [filename.split('.')[0] for filename in filename_array]

def get_sim_filenames_dict(subdir='umod', remove_ending=False, include_path=False, must_include=None):
    sim_ptmins = sorted(get_filenames(dataset='sim', subdir=subdir), key=int)
    sim_filenames = OrderedDict([(ptmin, get_filenames(dataset='sim', subdir=subdir, ptmin=ptmin, 
                                                       remove_ending=remove_ending, include_path=include_path,
                                                       must_include=must_include))
                                for ptmin in sim_ptmins])
    return sim_filenames

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