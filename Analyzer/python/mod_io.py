import gc
import h5py
import multiprocessing
import os
import pickle
import time
import warnings

import numpy as np

from utils import *

def separate_particle_arrays(particles, particles_index):
    
    # array to hold particles
    particles_array = np.zeros(len(particles_index) - 1, dtype='O')
    
    # iterate over indices
    for j,pi in enumerate(particles_index[:-1]):
        particles_array[j] = np.asarray(particles[pi:particles_index[j+1]])
        
    return particles_array

class MODDataset(object):

    def __init__(self, filename, dataset, njets=-1, store_jets=True, store_pfcs=True, store_gens=True):

        # check for valid dataset argument
        assert dataset in {'cms', 'sim', 'gen'}, "dataset must be one of 'cms', 'sim', 'gen'"

        # store parameters
        self.cms = (dataset == 'cms')
        self.sim = (dataset == 'sim')
        self.gen = (dataset == 'gen')
        self.store_jets = store_jets
        self.store_pfcs = store_pfcs and not self.gen
        self.store_gens = store_gens and not self.cms

        # get filepath
        if not filename.endswith('.h5'):
            filename += '.h5'
        if self.cms:
            self.filepath = os.path.join(PATHS['cms'], 'cuts', filename)
        else:
            self.filepath = os.path.join(PATHS['sim'], 'cuts', filename)

        # open h5 file
        self.hf = h5py.File(self.filepath, 'r')

        if self.store_jets:

            # load selected jets
            jets_slice = (slice(None) if njets == -1 else slice(njets))
            self.jets_i = self.hf['jets_i'][jets_slice]
            self.jets_f = self.hf['jets_f'][jets_slice]

            # store jets cols
            self.store_cols('jets_i')
            self.store_cols('jets_f')

        if self.store_pfcs:

            # determine pfcs_index
            pfcs_index_slice = slice(None) if njets == -1 else slice(njets + 1)
            self.pfcs_index = self.hf['pfcs_index'][pfcs_index_slice]

            # store pfcs as separate arrays
            self.all_pfcs = self.hf['pfcs'][:self.pfcs_index[-1]]
            self.pfcs = separate_particle_arrays(self.all_pfcs, self.pfcs_index)

            # store pfcs cols
            self.store_cols('pfcs')

            # set pfcs as particls
            if not self.gen:
                self.particles = self.pfcs
                self.particles_cols = self.pfcs_cols

        if self.store_gens:

            # determine gens_index
            gens_index_slice = slice(None) if njets == -1 else slice(njets + 1)
            self.gens_index = self.hf['gens_index'][gens_index_slice]

            # store gens as separate arrays
            self.all_gens = self.hf['gens'][:self.gens_index[-1]]
            self.gens = separate_particle_arrays(self.all_gens, self.gens_index)

            # store gens cols
            if self.gen:
                self.store_cols('gens')

                self.particles = self.gens
                self.particle_cols = self.gens_cols

        # store filenames
        self.filenames = self.hf['filenames'][:].astype('U')

    def store_cols(self, dset):

        try:

            # get cols from file
            cols = self.hf[dset].attrs['cols'].astype('U')
            setattr(self, dset + '_cols', cols)

            for i,col in enumerate(cols):

                # ensure cols are unique
                if hasattr(self, col):
                    raise RuntimeError("Repeat instances of col '{}', check file validity".format(col))

                # store column index
                else:
                    setattr(self, col, i)
        except:
            print(cols)
            raise

    def close(self):
        self.hf.close()

    def __del__(self):
        self.close()
        gc.collect()

class uMODReader(object):

    def __init__(self, filename, path=None, ptmin=None, dataset='cms', fileformat='umod',
                                 store_ak5s=True, store_pfcs=True,
                                 store_hards=True, store_gens=True):

        # ensure valid dataset
        assert dataset in DATASETS, 'dataset {} not in {}'.format(dataset, DATASETS)
        self.dataset = dataset
        self.sim = (self.dataset == 'sim')
        if not self.sim:
            store_gens = store_hards = False
        
        # ensure valid file format
        assert fileformat in FILEFORMATS, 'fileformat {} not in {}'.format(fileformat, FILEFORMATS)
        self.fileformat = fileformat
        self.umod = (self.fileformat == 'umod')
        
        # store filepath information
        ending = '.' + fileformat
        self.filename = filename if filename.endswith(ending) else filename + ending
        self.base_name = self.filename.split('.')[0]
        self.ptmin = str(ptmin) if (self.sim and ptmin is not None) else ''
        self.path = os.path.join(PATHS[self.dataset], self.fileformat, self.ptmin) if path is None else path
        
        # store options
        self.store_ak5s = store_ak5s
        self.store_pfcs = store_pfcs
        self.store_hards = store_hards
        self.store_gens = store_gens
        
        # dictionary to hold everything
        self.data = {'base_name': self.base_name, 'ptmin': self.ptmin, 'dataset': self.dataset, 'fileformat': self.fileformat}
        
        # time the reading in of the uMOD file
        self.start_time = time.time()
        self.read()
        self.duration = time.time() - self.start_time

        # get lumi block event totals if in jets 
        if not self.umod:
            self.get_lb_nevs()

        # check event total
        lb_tot = sum((lb['nev'] for lb in self.lbs))
        if lb_tot != self.nev_good:
            warnings.warn('Number of events mismatch: {} in file vs. {} expected'.format(lb_tot, self.nev_good))

    def __getattr__(self, attr):
        return self.data[attr]

    def save(self, path=None, endstr=''):
        
        # check that the file was actually done
        assert 'umodproducer_duration' in self.data
        
        # handle paths
        if path is None:
            path = os.path.join(PATHS[self.dataset], 'pickle', self.ptmin if self.sim else '')
        os.makedirs(path, exist_ok=True)
        
        # save data
        filename = '{}_{}{}.pickle'.format(self.base_name, self.fileformat, endstr)
        with open(os.path.join(path, filename), 'wb') as f:
            pickle.dump(self.data, f)

    def get_lb_nevs(self):
        for lb_i,count in zip(*np.unique(self.lbs_i, return_counts=True)):
            self.lbs[lb_i]['nev'] = count

    def read(self):

        # event info containers
        evs_i, evns, npvs, times, rhos, ftrigs, lbs_i = [], [], [], [], [], [], []

        # event data containers
        ak5s, pfcs, hards, gens = [], [], [], []

        # lumi blocks
        lbs = []

        # matching
        if not self.umod:
            ak52pfc, ak52gen, hard2gen = [], [], []

        with open(os.path.join(self.path, self.filename), 'r') as f:

            # containers for this event
            ev_ak5s, ev_pfcs, ev_hards, ev_gens = [], [], [], []

            in_lumi_block, in_event = False, False
            for line in f:
                parts = line.split()

                # handle empty line
                if len(parts) == 0:

                    # final actions for lumi block
                    if in_lumi_block:
                        in_lumi_block = False
                        lb['prescales_1'] = np.asarray(lb['prescales_1'], dtype=float)
                        lb['prescales_2'] = np.asarray(lb['prescales_2'], dtype=float)
                        lb['prescales'] = lb['prescales_1']*lb['prescales_2']
                        if self.umod:
                            lb['counts'] = np.zeros(len(lb['triggers']))
                        lbs.append(lb)
                        
                    # final actions for event
                    elif in_event:
                        in_event = False

                        if self.umod:
                            lbs_i.append(len(lbs) - 1)
                            lb['nev'] += 1

                        if self.store_ak5s:
                            ak5s.append(np.asarray(ev_ak5s, dtype=float))
                            ev_ak5s = []

                        if self.store_pfcs:
                            pfcs.append(np.asarray(ev_pfcs, dtype=float))
                            ev_pfcs = []

                        if self.store_hards:
                            hards.append(np.asarray(ev_hards, dtype=float))
                            ev_hards = []

                        if self.store_gens:
                            gens.append(np.asarray(ev_gens, dtype=float))
                            ev_gens = []

                        if not self.umod and not has_matches:
                            ak52pfc.append(np.array([]))

                            if self.sim:
                                ak52gen.append(np.array([]))
                                hard2gen.append(np.array([]))

                    continue

                # beginning of line tells you what it is
                key = parts[0]

                # check keys in decreasing order of frequency
                # having the continue statements seems to be a bit faster

                if in_event:

                    # particle flow candidate
                    if key == 'P':
                        if self.store_pfcs:
                            ev_pfcs.append(parts[1:])
                        continue

                    # sim specific keys
                    if self.sim:

                        # gen particle
                        if key == 'G':
                            if self.store_gens:
                                ev_gens.append(parts[1:])
                            continue

                        # hard process particle
                        if key == 'H':
                            if self.store_hards:
                                ev_hards.append(parts[1:])
                            continue

                        # match CMS jets to gen jets
                        if key == 'MatchAK52GEN':
                            ak52gen.append(np.asarray(parts[1:], dtype=int))
                            continue
                            
                        # match hard partons to get jets
                        if key == 'MatchHARD2GEN':
                            hard2gen.append(np.asarray(parts[1:], dtype=int))
                            continue

                    # antikt (R=0.5) jet
                    if key == 'AK5':
                        if self.store_ak5s:
                            ev_ak5s.append(parts[1:])
                        continue

                    # match CMS jets to PFC jets
                    if key == 'MatchAK52PFC':
                        has_matches = True
                        ak52pfc.append(np.asarray(parts[1:], dtype=int))
                        continue
                    
                    # event number
                    if key == 'EV':
                        evns.append(parts[1])
                        continue
                        
                    # number of primary vertices
                    if key == 'NPV':
                        npvs.append(parts[1])
                        continue
                        
                    # time of event in seconds, umod only
                    if key == 't(s)':
                        times.append(float(parts[1]))
                        continue
                        
                    # time of event offset in microseconds, umod only
                    if key == 't(us)':
                        times[-1] += float(parts[1])/10**6
                        continue

                    # time of event in seconds, jets only
                    if key == 't':
                        times.append(parts[1])
                        continue
                        
                    # event-wide rho
                    if key == 'RHO':
                        rhos.append(parts[1])
                        continue

                    # lumi block index, jets only
                    if key == 'iLB':
                        lbs_i.append(parts[1])
                        continue
                        
                    # triggers for event
                    if key == 'TR':
                        fired_triggers = np.asarray(parts[1:], dtype=int)
                        ftrigs.append(fired_triggers)
                        if self.umod:
                            lb['counts'][fired_triggers] += 1
                        continue

                # signifier of start of event
                if key == 'i':
                    in_event = True
                    has_matches = False
                    evs_i.append(parts[1])
                    continue

                if in_lumi_block:

                    # trigger names and prescales
                    if key.startswith('HLT'):
                        lb['triggers'].append(key.split('_v')[0] if self.umod else key)
                        lb['prescales_1'].append(parts[1])
                        lb['prescales_2'].append(parts[2])
                        continue
                    
                    # run number
                    if key == 'RN':
                        lb['rn'] = int(parts[1])
                        continue

                    # cross section information
                    if key == 'XS':
                        lb['cross_section'] = float(parts[1])
                        continue

                    # fired counts, jets only
                    if key == 'FCs':
                        lb['counts'] = np.asarray(parts[1:], dtype=int)
                        continue

                    # error on unknown lumi block key
                    raise ValueError('Unknown key ' + key + ' in lumi block from ' + self.filename)

                # signifier of start of lumi block
                if key == 'LB':
                    in_lumi_block = True
                    lb = {'triggers': [], 'prescales_1': [], 'prescales_2': [], 
                          'nev': 0, 'j': len(lbs), 'lb': int(parts[1])}
                    continue

                # summary info
                if key == 'NumBadEvents':
                    self.data['nev_bad'] = int(parts[1])
                elif key == 'NumGoodEvents':
                    self.data['nev_good'] = int(parts[1])
                elif key == 'NumTotalEvents':
                    self.data['nev_total'] = int(parts[1])
                elif key == 'Duration(s)':
                    self.data['umodproducer_duration'] = float(parts[1])
                elif key.startswith('#'):
                    pass
                        
                # error on unknown key
                else:
                    raise ValueError('Unknown key ' + key + ' from ' + self.filename)

        # store in dictionary
        self.data['evs_i'] = np.asarray(evs_i, dtype=int)
        self.data['evns'] = np.asarray(evns, dtype=int)
        self.data['npvs'] = np.asarray(npvs, dtype=int)
        self.data['times'] = np.asarray(times, dtype=float)
        self.data['rhos'] = np.asarray(rhos, dtype=float)
        self.data['ftrigs'] = np.asarray(ftrigs)
        self.data['lbs_i'] = np.asarray(lbs_i, dtype=int)
        self.data['ak5s'] = np.asarray(ak5s)
        self.data['pfcs'] = np.asarray(pfcs)
        self.data['lbs'] = np.asarray(lbs)

        if not self.umod:
            self.data['ak52pfc'] = np.asarray(ak52pfc)

        if self.sim:
            self.data['gens'] = np.asarray(gens)
            self.data['hards'] = np.asarray(hards)
            
            if not self.umod:
                self.data['ak52gen'] = np.asarray(ak52gen)
                self.data['hard2gen'] = np.asarray(hard2gen)
