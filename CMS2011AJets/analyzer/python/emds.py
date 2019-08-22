import time
import warnings

import energyflow as ef
import energyflow.emd
import energyflow.datasets.mod
import numpy as np

__all__ = ['calc_emds', 'calc_corrdims', 'calc_qg_corrdims']

num4emd_str = '40k'
ptranges = [(399, 401), (540, 560)]

dsets = ['cms', 'sim', 'gen']
names = [
    'CMS_Jet300_pT{}-{}GeV_EtaMax19_Quality2_{}',
    'SIM_Jet300_pT{}-{}GeV_EtaMax19_Quality2_{}',
    'GEN_pT{}-{}GeV_EtaMax19_{}',
]

rotations = [True, False]
filters = {
    'All': {'which': 'all', 'pt_cut': None, 'chs': True},
    'PFCpTCut': {'which': 'all', 'pt_cut': 1.0, 'chs': True},
    'Tracks': {'which': 'charged', 'pt_cut': None, 'chs': True},
    'TracksPFCpTCut': {'which': 'charged', 'pt_cut': 1.0, 'chs': True},
}

def calc_emds():

    for ptmin,ptmax in ptranges:
        for rotate in rotations:
            rot_str = '_Rotated' if rotate else ''

            for filtname,filt in filters.items():
                for dset,name in zip(dsets, names):
                    start = time.time()

                    dset_path = (dset if dset != 'gen' else 'sim')
                    name = name.format(ptmin, ptmax, num4emd_str)

                    print('Loading', name)
                    modds = ef.datasets.mod.MODDataset(name, path=utils.path(dset_path, 'h5'))

                    proc_ps = np.asarray([ps[utils.filter_particles(ps, **filt)] for ps in modds.particles])
                    proc_ps = np.asarray([ef.center_ptyphims(ps, copy=False, center='escheme') for ps in proc_ps])
                    if rotate:
                        proc_ps = np.asarray([ef.rotate_ptyphims(ps, copy=False, rotate='ptscheme') 
                                              for ps in proc_ps])
                    print('Processed particles in {:.3f}s'.format(time.time() - start))
                    
                    ptm = np.mean([ptmin, ptmax])
                    emds = ptm * ef.emd.emds(proc_ps, R=0.5, norm=True, gdim=2, mask=True, verbose=1, 
                                             print_every=5*10**7, n_jobs=82, empty_policy=-1)

                    emd_name = '{}_{}{}'.format(name, filtname, rot_str)
                    np.save(utils.path(dset_path, 'emds', emd_name), emds)
                    print('Saved {} EMDs, done in {:.3f}s overall'.format(emd_name, time.time() - start))
                    print()
                    
                    del emds, proc_ps, modds

def calc_corrdims():

    bins = 10**np.linspace(0, np.log10(500), 100)
    reg = 10**-30
    midbins = (bins[1:] + bins[:-1])/2
    dmidbins = np.log(midbins[1:]) - np.log(midbins[:-1])
    midbins2 = (midbins[:-1] + midbins[1:])/2

    for ptmin,ptmax in ptranges[1:]:
        for rotate in rotations[1:]:
            rot_str = '_Rotated' if rotate else ''

            for filtname,filt in filters.items():
                if filtname == 'All' or filtname == 'PFCpTCut':
                    print('skipping', filtname)
                    continue
                for dset,name in zip(dsets, names):
                    start = time.time()

                    name = name.format(ptmin, ptmax, num4emd_str)
                    dset_path = (dset if dset != 'gen' else 'sim')
                    modds = ef.datasets.mod.MODDataset(name, path=utils.path(dset_path, 'h5'), 
                                                       store_pfcs=False, store_gens=False)

                    emd_name = '{}_{}{}'.format(name, filtname, rot_str)
                    emds = np.load(utils.path(dset_path, 'emds', emd_name + '.npy'))

                    print('Loaded MODDataset and EMDs for {} in {:.3f}s'.format(emd_name, time.time() - start))
                    start = time.time()

                    emds = np.triu(emds)
                    mask = (emds > 0)
                    emds = emds[mask]
                    weights = (modds.weights[:,np.newaxis] * modds.weights)[mask]
                    del mask

                    hist, errs, _ = modplot.calc_hist(emds, bins=bins, weights=weights)
                    del emds, weights
                    print('Computed histogram in {:.3f}s'.format(time.time() - start))

                    counts = np.cumsum(hist) + reg
                    counts_errs = np.sqrt(np.cumsum(errs**2))

                    dims = (np.log(counts[1:]) - np.log(counts[:-1]))/dmidbins
                    dims_errs = np.sqrt((counts_errs[1:]/counts[1:])**2 + (counts_errs[:-1]/counts[:-1])**2)/dmidbins

                    fp = utils.path(dset_path, 'plotdata', 'EMDCorrDims_{}'.format(emd_name))
                    np.savez(fp, dims=dims, dims_errs=dims_errs, midbins2=midbins2)

def calc_qg_corrdims():

    bins = 10**np.linspace(0, np.log10(500), 100)
    reg = 10**-30
    midbins = (bins[1:] + bins[:-1])/2
    dmidbins = np.log(midbins[1:]) - np.log(midbins[:-1])
    midbins2 = (midbins[:-1] + midbins[1:])/2

    for ptmin,ptmax in ptranges[1:]:
        for rotate in rotations[:2]:
            rot_str = '_Rotated' if rotate else ''

            for filtname,filt in filters.items():
                for dset,name in zip(dsets, names):
                    if dset == 'cms': 
                        continue

                    start = time.time()

                    name = name.format(ptmin, ptmax, num4emd_str)
                    dset_path = 'sim'
                    modds = ef.datasets.mod.MODDataset(name, path=utils.path(dset_path, 'h5'), 
                                                       store_pfcs=False, store_gens=False)

                    emd_name = '{}_{}{}'.format(name, filtname, rot_str)
                    emds = np.load(utils.path(dset_path, 'emds', emd_name + '.npy'))

                    print('Loaded MODDataset and EMDs for {} in {:.3f}s'.format(emd_name, time.time() - start))
                    start = time.time()
                    
                    q_pids = np.array([1,2,3,4,5,6])
                    masks = {'q': np.isin(np.abs(modds.hard_pids), q_pids),
                             'g': (modds.hard_pids == 21)}
                    
                    emds = np.triu(emds)
                    dims_dict, dims_errs_dict = {}, {}
                    for p,mask in masks.items():

                        mask = (emds > 0) & (mask[:,np.newaxis] & mask)
                        uemds = emds[mask]
                        weights = (modds.weights[:,np.newaxis] * modds.weights)[mask]
                        del mask

                        hist, errs, _ = modplot.calc_hist(uemds, bins=bins, weights=weights)
                        del uemds, weights
                        print('Computed histogram in {:.3f}s'.format(time.time() - start))

                        counts = np.cumsum(hist) + reg
                        counts_errs = np.sqrt(np.cumsum(errs**2))

                        dims_dict[p+'dims'] = (np.log(counts[1:]) - np.log(counts[:-1]))/dmidbins
                        dims_errs_dict[p+'errs'] = np.sqrt((counts_errs[1:]/counts[1:])**2 + 
                                                           (counts_errs[:-1]/counts[:-1])**2)/dmidbins

                    fp = utils.path(dset_path, 'plotdata', 'EMDCorrDimsQG_{}'.format(emd_name))
                    np.savez(fp, midbins2=midbins2, **dims_dict, **dims_errs_dict)


if __name__ == '__main__':

    import modplot
    import utils

else:

    from . import modplot
    from . import utils
    