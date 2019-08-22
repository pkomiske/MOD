from collections import OrderedDict
import os
import time

import energyflow as ef
import numpy as np

with open(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, 'DATAPATH'), 'r') as f:
    DATAPATH = f.read().strip()

DATASETS = {'cms': 'JetPrimaryDataset', 'sim': 'QCDSimDatasets'}
DATASET_PATHS = {k: os.path.join(DATAPATH, v) for k,v in DATASETS.items()}
FILEFORMATS = {'mod', 'jet'}

def path(dataset, *args):
    return os.path.join(DATASET_PATHS[dataset], *args)

def get_filenames(dataset='cms', subdir='mod', ptmin=None, remove_ending=False, 
                  include_path=False, must_include=None, must_ignore=None):

    local_path = path(dataset, '' if subdir is None else subdir)

    if ptmin is not None:
        local_path = os.path.join(local_path, str(ptmin))

    filenames = os.listdir(local_path)
    if must_include is None:
        filenames = sorted(filenames)
    else:
        filenames = sorted(filter(lambda x: must_include in x, filenames))

    if must_ignore is not None:
        filenames = sorted(filter(lambda x: must_ignore not in x, filenames))

    if include_path:
        filenames = [os.path.join(local_path, filename) for filename in filenames]

    if remove_ending:
        filenames = [filename.split('.')[0] for filename in filenames]

    return filenames

def get_sim_filenames_dict(subdir='mod', remove_ending=False, include_path=False, must_include=None, must_ignore=None):

    # sort these according to numerical value
    sim_ptmins = sorted(get_filenames(dataset='sim', subdir=subdir, must_ignore='.tar.gz'), key=int)

    # get ordered dictionary of sim filenames
    sim_filenames = OrderedDict([(ptmin, get_filenames(dataset='sim', subdir=subdir, ptmin=ptmin, 
                                                       remove_ending=remove_ending, include_path=include_path,
                                                       must_include=must_include, must_ignore=must_ignore))
                                for ptmin in sim_ptmins])
    return sim_filenames

def concat_npz_files(filepaths, dataset):
    
    # check for validity
    assert dataset in {'cms', 'sim', 'gen'}, "dataset should be one of 'cms', 'sim', 'gen'"
    
    start = time.time()
    
    # arrays
    jets_float, jets_int = [], []
    pfcs, pfcs_index = [], []
    gens, gens_index = [], []
    nprev_pfcs = nprev_gens = 0

    # counters
    nevs = nevs_valid = nevs_trig_fired = 0
    njets_gen = njets_gen_pass_ptmin = 0
    njets = njets_pass_ptmin = njets_ak5_match = 0

    prefix = 'gen_' if dataset == 'gen' else ''
    
    # iterate over filepaths
    for i,filepath in enumerate(filepaths):
        
        # load
        f_npz = np.load(filepath)

        # update counters
        if dataset == 'gen':
            njets_gen += f_npz['njets_gen']
            njets_gen_pass_ptmin += f_npz['njets_gen_pass_ptmin']
        else:
            nevs += f_npz['nevs']
            nevs_valid += f_npz['nevs_valid']
            nevs_trig_fired += f_npz['nevs_trig_fired']
            njets += f_npz['njets']
            njets_pass_ptmin += f_npz['njets_pass_ptmin']
            njets_ak5_match += f_npz['njets_ak5_match']
        
        # check that we have any jets
        jets_int_i = f_npz['{}jets_int'.format(prefix)]
        if not len(jets_int_i):
            continue
        
        # append jets
        jets_float.append(f_npz['{}jets_float'.format(prefix)])
        jets_int.append(jets_int_i)
        
        # pfcs
        if dataset != 'gen':
            pfcs.append(separate_particle_arrays(f_npz['pfcs'], f_npz['pfcs_index']))
            
        # gens
        if dataset != 'cms':
            p = 'sim_' if dataset == 'sim' else 'gen_'
            gens.append(separate_particle_arrays(f_npz[p + 'gens'], f_npz[p + 'gens_index']))
            
        if (i+1) % 100 == 0:
            print('  [{}/{}] files processed in {:.3f}s'.format(i+1, len(filepaths), time.time() - start))
            
    # concatenate
    jets_float = np.concatenate(jets_float, axis=0)
    jets_int = np.concatenate(jets_int, axis=0)
    
    if dataset != 'gen':
        pfcs = np.concatenate(pfcs, axis=0)
        
    if dataset != 'cms':
        gens = np.concatenate(gens, axis=0)
        
    print('Finished and concatenated in {:.3f}s'.format(time.time() - start))

    if dataset != 'gen':
        print('Num. Events:', nevs)
        print('Num. Events Valid:', nevs_valid)
        print('Num. Events Trig Fired:', nevs_trig_fired)
        print('Num. Hardest 2 Jets:', njets)
        print('Num. Jets > 375 GeV:', njets_pass_ptmin)
        print('Num. Jets AK5 Matched:', njets_ak5_match)
    else:
        print('Num. Gen Jets:', njets_gen)
        print('Num. Gen Jets > 375 GeV:', njets_gen_pass_ptmin)
    
    return jets_float, jets_int, pfcs, gens

def separate_particle_arrays(particles, particles_index):
    
    # array to hold particles
    particles_array = np.zeros(len(particles_index) - 1, dtype='O')
    
    # iterate over indices
    for j,pi in enumerate(particles_index[:-1]):
        particles_array[j] = np.asarray(particles[pi:particles_index[j+1]])
        
    return particles_array

def make_particles_index(particle_arrays):
    
    # list of indices
    index = [0]

    # iterate over all particles
    for particles in particle_arrays:
        index.append(index[-1] + len(particles))

    # convert to numpy array with proper dtype
    return np.asarray(index, dtype=np.uint32)

def write_large_object_array_to_h5(hf, name, arr, dtype, ncols=None, chunksize=10**5, **compression):

    nrows = sum([len(x) for x in arr])
    ncols = arr[0].shape[1] if ncols is None else ncols

    dataset = hf.create_dataset(name, (nrows, ncols), dtype=dtype, **compression)

    begin = end = ind = 0
    while end < len(arr):
        end = min(len(arr), end + chunksize)

        arr_chunk = np.concatenate(arr[begin:end], axis=0)
        dataset[ind:ind+len(arr_chunk)] = arr_chunk

        begin = end
        ind += len(arr_chunk)
        del arr_chunk

    return dataset

try:
    rwf = np.load(path('sim', 'ReweightingFactors.npz'))
    kf_x, kf_y = rwf['kfactor_x'], rwf['kfactor_y']
    npv_hist_ratios = rwf['npv_hist_ratios']
    residual_factor = rwf['residual_factor']
    rwf.close()

except:
    pass

def k_factors(pts):
    return np.interp(pts, kf_x, kf_y)

def npv_factors(npvs):
    return npv_hist_ratios[np.asarray(npvs, dtype=int)]

def sim_factors(pts, npvs):
    return k_factors(pts) * npv_factors(npvs) * residual_factor

def gen_factors(pts):
    return k_factors(pts) * residual_factor

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