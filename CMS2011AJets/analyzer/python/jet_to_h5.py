import itertools
import multiprocessing
import os
import pickle
import time

import energyflow as ef
import energyflow.datasets.mod
import h5py
import numpy as np

__all__ = ['make_filename_arrays', 
           'cms_jets_to_npzs',
           'sim_jets_to_npzs',
           'cms_npzs_to_h5',
           'sim_npzs_to_h5',
           'gen_npzs_to_h5']

PTMIN = 375
TRIG_NAME = 'HLT_Jet300'
COMPRESSION = 6
NJETS_PER_FILE = 100000

def make_filename_arrays():

    # cms
    cms_filenames = utils.get_filenames(subdir='jet', remove_ending=True)

    print('Saving CMS FilenameArray, {} files'.format(len(cms_filenames)))
    np.save(utils.path('cms', 'FilenameArray.npy'), np.asarray(cms_filenames))

    print('Saving CMS FilenameMap')
    filename_map = {name: i for i,name in enumerate(cms_filenames)}
    with open(utils.path('cms', 'FilenameMap.pickle'), 'wb') as f:
        pickle.dump(filename_map, f)

    # sim
    filename_map = {}
    filename_array = []
    tot = 0

    # get sim filenames and iterate
    sim_names_dict = utils.get_sim_filenames_dict(subdir='jet', remove_ending=True)
    sim_ptmins = list(sim_names_dict.keys())
    for ptmin in sim_ptmins:
        filename_map[ptmin] = {}

        # iterate over individual names
        i = 0
        for name in sim_names_dict[ptmin]:
            filename_map[ptmin][name] = i + tot
            filename_array.append([ptmin, name])
            i += 1

        tot += i

    print('Saving SIM FilenameArray, {} files'.format(tot))
    np.save(utils.path('sim', 'FilenameArray.npy'), np.asarray(filename_array))

    print('Saving SIM FilenameMap')
    with open(utils.path('sim', 'FilenameMap.pickle'), 'wb') as f:
        pickle.dump(filename_map, f)

def process_gen_jet(ev_gens, gen_i, gen2hard, hards, phi=None):

    # get particles for this gen jet and append to overall list
    ev_gens_i = ev_gens[ev_gens[:,0] == gen_i, 1:]

    # calculate kinematic quantities for this gen jet
    gen_pt, gen_y, gen_phi, gen_m = ef.sum_ptyphims(ev_gens_i)
    gen_eta = np.arcsinh(np.sinh(gen_y)*np.sqrt(1. + (gen_m/gen_pt)**2))

    # phi fix
    if phi is not None:
        gen_phi = ef.phi_fix(gen_phi, phi, copy=False)
        ev_gens_i[:,2] = ef.phi_fix(ev_gens_i[:,2], phi, copy=False)

    # determine if gen jet matched to a hard parton
    hard_i = gen2hard.get(gen_i, -1)
    if hard_i == -1:
        hard_pt = hard_y = hard_phi = -1.
        hard_pid = 0  
    else:
        hard_pt, hard_y, hard_phi = hards[hard_i,1:4]
        hard_pid = hards[hard_i,5]

        # phi fix
        hard_phi = ef.phi_fix(hard_phi, gen_phi, copy=False)

    return (gen_pt, gen_y, gen_phi, gen_m, gen_eta, ev_gens_i, hard_pt, hard_y, hard_phi, hard_pid)

def jet_to_npz(arg):

    # unpack arguments
    dataset, filepath = arg
    sim, cms = (dataset != 'cms'), (dataset == 'cms')

    # grab globals and store locally
    local_weight = weight
    ptmin = PTMIN

    # load file
    data = np.load(filepath, allow_pickle=True)

    # arrays
    jets_float, jets_int = [], []
    pfcs, pfcs_index = [], [0]

    # arrays for gen
    if sim:
        sim_gens, sim_gens_index = [], [0]
        gen_gens, gen_gens_index = [], [0]
        gen_jets_float, gen_jets_int = [], []

    # counters
    nevs = nevs_valid = nevs_trig_fired = 0
    njets_gen = njets_gen_pass_ptmin = 0
    njets = njets_pass_ptmin = njets_ak5_match = 0

    # grab attributes of data that will be used repeatedly
    filename_index, f_lbs = filename_map[dataset][data['base_name']], data['lbs']

    # process lumiblocks for this file
    for lb in f_lbs:

        # note: make prescales negative in order to indicate "off" by default
        lb_triggers, lb_prescales = lb['triggers'], -lb['prescales']

        # the indices where the jet triggers are for this lumi block, -1 if not present
        lb['jtn_ind'] = lb_triggers.index(TRIG_NAME) if TRIG_NAME in lb_triggers else -1

    # sim only args
    if sim:
        more_args = (data['gens'], data['hards'], data['ak52gen'], data['hard2gen'])
    else:
        more_args = 4*[itertools.repeat(None)]

    # zipped iterable
    z = zip(data['ak5s'], data['pfcs'], data['ak52pfc'], data['ftrigs'], 
            data['lbs_i'], data['evns'], data['npvs'], *more_args)

    # iterate over events
    for ak5s, ev_pfcs, ak52pfc, ftrigs, lb_i, evn, npv, ev_gens, hards, ak52gen, hard2gen in z:

        # count total events
        nevs += 1

        # get lumiblock for this event
        lb = f_lbs[lb_i]
        rn, lbn = lb['rn'], lb['lb']

        # skip event if cms lumiblock invalid
        if cms:
            rnlb = (rn, lbn)
            if rnlb not in lbs or lbs[rnlb]['rec_lumi'] == 0.:
                continue

        # count valid events
        nevs_valid += 1

        # skip if trigger didn't fire
        if (lb['jtn_ind'] not in ftrigs):
            continue

        nevs_trig_fired += 1

        # skip if no jets
        if not len(ak5s):
            continue

        # count total number of jets
        njets += min(2, len(ak5s))

        # determine if we have particles
        have_particles = False
        if ev_pfcs.size > 0:
            have_particles = True
            
            # make mapping from gen index to hard index
            if sim:
                gen2hard = {gen: hard for hard,gen in enumerate(hard2gen)}

        # handle gen jets
        if sim and have_particles:

            # hardest two gen jets
            for gen_i in range(2):

                # process gen jet
                (gen_pt, gen_y, gen_phi, gen_m, gen_eta, ev_gens_i, 
                 hard_pt, hard_y, hard_phi, hard_pid) = process_gen_jet(ev_gens, gen_i, gen2hard, hards)

                # count gen_jets above 10 GeV
                if gen_pt > 10.:
                    njets_gen += 1

                # impose cuts on gen
                if gen_pt < ptmin:
                    continue

                njets_gen_pass_ptmin += 1

                # append information
                gen_jets_float.append((gen_pt, gen_y, gen_phi, gen_m, gen_eta,
                                       hard_pt, hard_y, hard_phi, local_weight))
                gen_jets_int.append((filename_index, rn, lbn, evn, hard_pid))
                gen_gens.append(ev_gens_i)
                gen_gens_index.append(gen_gens_index[-1] + len(ev_gens_i))

        # apply JECs
        pts_corrected = ak5s[:,1] * ak5s[:,5]
        argsort = np.argsort(pts_corrected)[::-1]

        # iterate over hardest two jets
        for i in argsort[:2]:

            # impose cuts on jet
            if pts_corrected[i] < ptmin:
                continue

            # record that jet passed cuts
            njets_pass_ptmin += 1

            # get matching index
            pfc_i = ak52pfc[i]
            if pfc_i == -1:
                continue

            # record that jet matched
            njets_ak5_match += 1

            # get kinematic quantities for the jet
            ak5 = ak5s[i]
            pt, eta, phi, m, jec, area = ak5[1:7]
            neutral_frac, quality = max(ak5[9:11]), ak5[13]

            # calculate rapidity of jet
            y = np.arcsinh(np.sinh(eta)/np.sqrt(1. + (m/pt)**2))

            # if we don't have particles, -1 (or other default) for everything
            if not have_particles:

                # set particles indexes to last position (to make empty array)
                pfcs_index_i = pfcs_index[-1]

                if sim:
                    gen_pt = gen_y = gen_phi = gen_m = gen_eta = hard_pt = hard_y = hard_phi = -1.
                    hard_pid = 0
                    sim_gens_index_i = sim_gens_index[-1]

            # we've got particles!
            else:
                
                # get pfcs for this jet and append to overall list
                ev_pfcs_i = ev_pfcs[ev_pfcs[:,0] == pfc_i, 1:]
                pfcs.append(ev_pfcs_i)
                pfcs_index_i = pfcs_index[-1] + len(ev_pfcs_i)

                # sim only
                if sim:
                
                    # determine matching gen jet and add to gen_inds
                    gen_i = ak52gen[i]

                    # if not matched, -1 for everything gen related
                    if gen_i == -1:
                        gen_pt = gen_y = gen_phi = gen_m = gen_eta = hard_pt = hard_y = hard_phi = -1.
                        hard_pid =  0
                        sim_gens_index_i = sim_gens_index[-1]

                    # we have a gen jet associated to this ak5
                    else:

                        # process gen jet
                        (gen_pt, gen_y, gen_phi, gen_m, gen_eta, ev_gens_i, 
                         hard_pt, hard_y, hard_phi, hard_pid) = process_gen_jet(ev_gens, gen_i, gen2hard, hards, phi=phi)

                        # store gen information for sim
                        sim_gens.append(ev_gens_i)
                        sim_gens_index_i = sim_gens_index[-1] + len(ev_gens_i)

            # append pfc index
            pfcs_index.append(pfcs_index_i)

            # jets info cms
            if cms:
                jets_float.append((pt, y, phi, m, eta, jec, area, neutral_frac, local_weight))
                jets_int.append((filename_index, rn, lbn, evn, npv, quality))

            # jets info sim
            else:
                jets_float.append((pt, y, phi, m, eta, jec, area, neutral_frac,
                                   gen_pt, gen_y, gen_phi, gen_m, gen_eta,
                                   hard_pt, hard_y, hard_phi, local_weight))
                jets_int.append((filename_index, rn, lbn, evn, npv, quality, hard_pid))
                sim_gens_index.append(sim_gens_index_i)

    # convert to numpy
    jets_float = np.asarray(jets_float, dtype=np.float64)
    jets_int = np.asarray(jets_int, dtype=np.int64)
    pfcs_index = np.asarray(pfcs_index, dtype=np.uint32)
    pfcs = (np.concatenate(pfcs, axis=0) if len(pfcs) else np.zeros((0, 6))).astype(np.float64)

    if sim:
        sim_gens_index = np.asarray(sim_gens_index, dtype=np.uint32)
        sim_gens = (np.concatenate(sim_gens, axis=0) if len(sim_gens) else np.zeros((0, 5))).astype(np.float64)

        # handle gen arrays
        gen_jets_float = np.asarray(gen_jets_float, dtype=np.float64)
        gen_jets_int = np.asarray(gen_jets_int, dtype=np.int64)
        gen_gens_index = np.asarray(gen_gens_index, dtype=np.uint32)
        gen_gens = (np.concatenate(gen_gens, axis=0) if len(gen_gens) else np.zeros((0, 5))).astype(np.float64)

    # arrays to be saved for cms/sim
    arrs = {
        'jets_float': jets_float, 'jets_int': jets_int,
        'pfcs_index': pfcs_index, 'pfcs': pfcs,
        'nevs': nevs,
        'nevs_valid': nevs_valid,
        'nevs_trig_fired': nevs_trig_fired,
        'njets_gen': njets_gen,
        'njets_gen_pass_ptmin': njets_gen_pass_ptmin,
        'njets': njets,
        'njets_pass_ptmin': njets_pass_ptmin,
        'njets_ak5_match': njets_ak5_match,
    }

    # save cms
    if cms:
        np.savez(utils.path('cms', 'npz', data['base_name'] + '.npz'), **arrs)

    else:

        # save sim
        arrs.update({'sim_gens_index': sim_gens_index, 'sim_gens': sim_gens})
        np.savez(utils.path('sim', 'sim_npz', dataset, data['base_name'] + '.npz'), **arrs)

        # save gen
        gen_arrs = {'gen_jets_float': gen_jets_float, 'gen_jets_int': gen_jets_int,
                    'gen_gens_index': gen_gens_index, 'gen_gens': gen_gens,
                    'njets_gen': njets_gen, 'njets_gen_pass_ptmin': njets_gen_pass_ptmin}
        np.savez(utils.path('sim', 'gen_npz', dataset, data['base_name'] + '.npz'), **gen_arrs)

def parallelize_jets_to_npz(args):

    start = time.time()
    with multiprocessing.Pool(processes=None) as pool:

        # iterate over files as they're done
        for i,_ in enumerate(pool.imap_unordered(jet_to_npz, args, chunksize=1)):

            # print progress
            if (i+1) % 100 == 0:
                print('  [{}/{}] files processed in {:.3f}s'.format(i+1, len(args), time.time() - start))

        print('  [{}/{}] files processed in {:.3f}s'.format(i+1, len(args), time.time() - start))

def cms_jets_to_npzs():

    # declare global variables
    global lbs, weight, filename_map

    # luminosity blocks
    lbs = np.load(utils.path('cms', 'lumis', 'ValidatedLumiblocks.pickle'), allow_pickle=True)

    # load effective luminosities for weight
    trig_leffs = np.load(utils.path('cms', 'lumis', 'EffectiveLuminositiesByTrigger.pickle'), allow_pickle=True)
    weight = 1.0/trig_leffs[TRIG_NAME]

    # load filename_map
    filename_map = {'cms': np.load(utils.path('cms', 'FilenameMap.pickle'), allow_pickle=True)}

    # ensure output directory exists
    os.makedirs(utils.path('cms', 'npz'), exist_ok=True)

    print('CMS Jet Primary Dataset: jet ptmin {} GeV'.format(PTMIN))

    # parallelize
    args = [('cms', f) for f in utils.get_filenames(subdir='pickle', must_ignore='noparticles', include_path=True)]
    parallelize_jets_to_npz(args)

def sim_jets_to_npzs():

    # declare global variables
    global weight, filename_map

    # load cross sections for weight
    cross_sections = np.load(utils.path('sim', 'CrossSections.pickle'), allow_pickle=True)

    # load filename map
    filename_map = np.load(utils.path('sim', 'FilenameMap.pickle'), allow_pickle=True)

    sim_filepaths = utils.get_sim_filenames_dict(subdir='pickle', must_ignore='noparticles', include_path=True)
    for ptmin,filepaths in sim_filepaths.items():

        # ensure output directories exists
        os.makedirs(utils.path('sim', 'sim_npz', ptmin), exist_ok=True)
        os.makedirs(utils.path('sim', 'gen_npz', ptmin), exist_ok=True)
        
        print('SIM {} Dataset: jet ptmin {} GeV'.format(ptmin, PTMIN))

        # convert from pb to nb
        weight = cross_sections[ptmin]['xs']/cross_sections[ptmin]['nev']/10**3

        args = [(ptmin, filepath) for filepath in filepaths]
        parallelize_jets_to_npz(args)

def cms_npzs_to_h5():

    # cms files
    filepaths = utils.get_filenames(include_path=True, subdir='npz')

    # cms filename array
    filename_array = np.load(utils.path('cms', 'FilenameArray.npy'))

    # load and concatenate all arrays
    jets_f, jets_i, pfcs, _ = utils.concat_npz_files(filepaths, 'cms')

    # cms cols
    jets_f_cols = ['jet_pt', 'jet_y', 'jet_phi', 'jet_m', 'jet_eta', 'jec', 'jet_area', 'jet_max_nef', 'weight']
    jets_i_cols = ['fn', 'rn', 'lbn', 'evn', 'npv', 'quality']
    pfcs_cols = ['pt', 'y', 'phi', 'm', 'pid', 'vertex']

    # make moddataset
    arrays = {'jets_f': jets_f, 'jets_i': jets_i, 'pfcs': pfcs,
              'jets_f_cols': jets_f_cols, 'jets_i_cols': jets_i_cols, 'pfcs_cols': pfcs_cols,
              'filenames': filename_array}

    name = 'CMS_Jet300_pT375-infGeV'
    filepath = utils.path('cms', 'h5', *((1 if NJETS_PER_FILE == -1 else 2)*[name]))

    cms = ef.datasets.mod.MODDataset(_dataset='cms', _arrays=arrays)
    cms.save(filepath, npf=NJETS_PER_FILE, compression=COMPRESSION, n_jobs=-1)
    del cms, jets_f, jets_i, pfcs

def sim_npzs_to_h5():

    # sim files
    sim_filepaths = utils.get_sim_filenames_dict(include_path=True, subdir='sim_npz')

    # sim filename array
    filename_array = np.load(utils.path('sim', 'FilenameArray.npy'))

    # sim cols
    jets_f_cols = ['jet_pt', 'jet_y', 'jet_phi', 'jet_m', 'jet_eta', 'jec', 'jet_area', 'jet_max_nef',
                   'gen_jet_pt', 'gen_jet_y', 'gen_jet_phi', 'gen_jet_m', 'gen_jet_eta', 
                   'hard_pt', 'hard_y', 'hard_phi', 'weight']
    jets_i_cols = ['fn', 'rn', 'lbn', 'evn', 'npv', 'quality', 'hard_pid']
    pfcs_cols = ['pt', 'y', 'phi', 'm', 'pid', 'vertex']

    # iterate over datasets
    for ptmin,filepaths in sim_filepaths.items():

        if ptmin not in ['120', '170']: 
            continue

        print('SIM {} Dataset'.format(ptmin))

        # load and concatenate all arrays
        jets_f, jets_i, pfcs, gens = utils.concat_npz_files(filepaths, 'sim')

        # add vertex column to gens
        gens = np.asarray([np.hstack((g, np.zeros((len(g), 1)))) for g in gens])

        # make moddataset
        arrays = {'jets_f': jets_f, 'jets_i': jets_i, 'pfcs': pfcs, 'gens': gens,
                  'jets_f_cols': jets_f_cols, 'jets_i_cols': jets_i_cols, 
                  'pfcs_cols': pfcs_cols, 'gens_cols': pfcs_cols,
                  'filenames': filename_array}

        name = 'SIM{}_Jet300_pT375-infGeV'.format(ptmin)
        filepath = utils.path('sim', 'h5', *((1 if NJETS_PER_FILE == -1 else 2)*[name]))

        sim = ef.datasets.mod.MODDataset(_dataset='sim', _arrays=arrays)
        sim.save(filepath, npf=NJETS_PER_FILE, compression=COMPRESSION, n_jobs=-1)
        del sim, jets_f, jets_i, pfcs, gens

        print()

def gen_npzs_to_h5():
    
    # gen files
    gen_filepaths = utils.get_sim_filenames_dict(include_path=True, subdir='gen_npz')

    # sim filename array
    filename_array = np.load(utils.path('sim', 'FilenameArray.npy'))

    # gen cols
    jets_f_cols = ['jet_pt', 'jet_y', 'jet_phi', 'jet_m', 'jet_eta', 
                   'hard_pt', 'hard_y', 'hard_phi', 'weight']
    jets_i_cols = ['fn', 'rn', 'lbn', 'evn', 'hard_pid']
    gens_cols = ['pt', 'y', 'phi', 'm', 'pid', 'vertex']

    # iterate over datasets
    for ptmin,filepaths in gen_filepaths.items():

        print('GEN {} Dataset'.format(ptmin))

        # load and concatenate all arrays
        jets_f, jets_i, _, gens = utils.concat_npz_files(filepaths, 'gen')

        # add vertex column to gens
        gens = np.asarray([np.hstack((g, np.zeros((len(g), 1)))) for g in gens])

        # make moddataset
        arrays = {'jets_f': jets_f, 'jets_i': jets_i, 'gens': gens,
                  'jets_f_cols': jets_f_cols, 'jets_i_cols': jets_i_cols, 'gens_cols': gens_cols,
                  'filenames': filename_array}

        name = 'GEN{}_pT375-infGeV'.format(ptmin)
        filepath = utils.path('sim', 'h5', *((1 if NJETS_PER_FILE == -1 else 2)*[name]))

        gen = ef.datasets.mod.MODDataset(_dataset='gen', _arrays=arrays)
        gen.save(filepath, npf=NJETS_PER_FILE, compression=COMPRESSION, n_jobs=-1)
        del gen, jets_f, jets_i, gens

        print()

if __name__ == '__main__':

    import mod_io
    import utils

    make_filename_arrays()
    cms_jets_to_npzs()
    sim_jets_to_npzs()
    cms_npzs_to_h5()
    sim_npzs_to_h5()
    gen_npzs_to_h5()

else:

    from . import mod_io
    from . import utils