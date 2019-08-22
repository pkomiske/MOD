import gc
import h5py
import multiprocessing
import os
import pickle
import time
import warnings

import numpy as np
import six

# for import as part of package
try:
    from . import utils

# for standalone import
except ImportError:
    import utils

class MODReader(object):

    def __init__(self, filename, path=None, ptmin=None, dataset='cms', fileformat='mod', store_particles=True):

        # ensure valid dataset
        assert dataset in utils.DATASETS, 'dataset {} not in {}'.format(dataset, utils.DATASETS)
        self.dataset = dataset
        self.sim = (self.dataset == 'sim')
        if not self.sim:
            store_gens = store_hards = False
        
        # ensure valid file format
        assert fileformat in utils.FILEFORMATS, 'fileformat {} not in {}'.format(fileformat, utils.FILEFORMATS)
        self.fileformat = fileformat
        self.mod = (self.fileformat == 'mod')
        
        # store filepath information
        ending = '.' + fileformat
        self.filename = filename if filename.endswith(ending) else filename + ending
        self.base_name = self.filename.split('.')[0]
        self.ptmin = str(ptmin) if (self.sim and ptmin is not None) else ''
        self.path = utils.path(self.dataset, self.fileformat, self.ptmin) if path is None else path
        
        # store options
        self.store_particles = store_particles
        
        # dictionary to hold everything
        self.data = {'base_name': self.base_name, 'ptmin': self.ptmin, 'dataset': self.dataset, 'fileformat': self.fileformat}
        
        # time the reading in of the MOD file
        self.start_time = time.time()
        self.read()
        self.duration = time.time() - self.start_time

        # get lumi block event totals if in jets 
        if not self.mod:
            self.get_lb_nevs()

        # check event total
        lb_tot = sum((lb['nev'] for lb in self.lbs))
        if lb_tot != self.nev_good:
            warnings.warn('Number of events mismatch: {} in file vs. {} expected'.format(lb_tot, self.nev_good))

    def __getattr__(self, attr):
        return self.data[attr]

    def save(self, path=None, endstr='', save_particles=True):
        
        # check that the file was actually done
        assert 'modproducer_duration' in self.data
        
        # handle paths
        if path is None:
            path = os.path.join(self.path, os.pardir, 'pickle', self.ptmin if self.sim else '')
        os.makedirs(path, exist_ok=True)
        
        # save data
        filename = '{}_{}{}.pickle'.format(self.base_name, self.fileformat, endstr)
        with open(os.path.join(path, filename), 'wb') as f:
            data = self.data.copy()

            # remove all particle related arrays
            if not save_particles:
                del data['pfcs']

                if not self.mod:
                    del data['ak52pfc']

                if self.sim:
                    del data['gens'], data['hards']

                    if not self.mod:
                        del data['ak52gen'], data['hard2gen']

            pickle.dump(data, f)

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
        if not self.mod:
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
                        if self.mod:
                            lb['counts'] = np.zeros(len(lb['triggers']))
                        lbs.append(lb)
                        
                    # final actions for event
                    elif in_event:
                        in_event = False

                        if self.mod:
                            lbs_i.append(len(lbs) - 1)
                            lb['nev'] += 1

                        ak5s.append(np.asarray(ev_ak5s, dtype=float))
                        ev_ak5s = []

                        if self.store_particles:
                            pfcs.append(np.asarray(ev_pfcs, dtype=float))
                            ev_pfcs = []

                            hards.append(np.asarray(ev_hards, dtype=float))
                            ev_hards = []

                            gens.append(np.asarray(ev_gens, dtype=float))
                            ev_gens = []

                        if not self.mod and not has_matches:
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
                        if self.store_particles:
                            ev_pfcs.append(parts[1:])
                        continue

                    # sim specific keys
                    if self.sim:

                        # gen particle
                        if key == 'G':
                            if self.store_particles:
                                ev_gens.append(parts[1:])
                            continue

                        # hard process particle
                        if key == 'H':
                            if self.store_particles:
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
                        
                    # time of event in seconds, mod only
                    if key == 't(s)':
                        times.append(float(parts[1]))
                        continue
                        
                    # time of event offset in microseconds, mod only
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
                        if self.mod:
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
                        lb['triggers'].append(key.split('_v')[0] if self.mod else key)
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
                    self.data['modproducer_duration'] = float(parts[1])
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

        if not self.mod:
            self.data['ak52pfc'] = np.asarray(ak52pfc)

        if self.sim:
            self.data['gens'] = np.asarray(gens)
            self.data['hards'] = np.asarray(hards)
            
            if not self.mod:
                self.data['ak52gen'] = np.asarray(ak52gen)
                self.data['hard2gen'] = np.asarray(hard2gen)
