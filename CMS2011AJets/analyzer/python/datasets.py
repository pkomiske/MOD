import os
import time

from energyflow.datasets import mod
import numpy as np

__all__ = ['cms_375up',
           'sim_375to425',
           'sim_425to700',
           'sim_700up',
           'gen_375to425',
           'gen_425to700',
           'gen_700up',
           'sim_noparticles',
           'gen_noparticles',
           'cms_subsample_datasets_for_emds',
           'sim_subsample_datasets_for_emds',
           'gen_subsample_datasets_for_emds']

COMPRESSION = None
NJETS_PER_FILE = -1

def dataset_from_pt_range(dataset, ptmin, ptmax=None):

    s = {'dataset': dataset, 'ptmin': ptmin, 'ptmax': ptmax}

    def tmpfunc():

        dataset, ptmin, ptmax = s['dataset'], s['ptmin'], s['ptmax']

        specs = [('corr_jet_pt>=', ptmin)]
        if ptmax is not None:
            specs.append(('corr_jet_pt<', ptmax))
        else:
            ptmax = 'inf'

        dset = mod.load(*specs, dataset=dataset, amount=-1, verbose=1)

        if dataset == 'gen':
            name = 'GEN_pT{}-{}GeV'.format(ptmin, ptmax)
        elif dataset == 'cms':
            name = 'CMS_Jet300_pT{}-{}GeV'.format(ptmin, ptmax)
        elif dataset == 'sim':
            name = 'SIM_Jet300_pT{}-{}GeV'.format(ptmin, ptmax)
        else:
            raise ValueError('bad dataset')

        filepath = utils.path('sim' if dataset == 'gen' else dataset, 'h5', name)
        dset.save(filepath, npf=NJETS_PER_FILE, compression=COMPRESSION)

        del dset

    return tmpfunc

cms_375up = dataset_from_pt_range('cms', 375)

sim_375to425 = dataset_from_pt_range('sim', 375, 425)
sim_425to700 = dataset_from_pt_range('sim', 425, 700)
sim_700up = dataset_from_pt_range('sim', 700)

gen_375to425 = dataset_from_pt_range('gen', 375, 425)
gen_425to700 = dataset_from_pt_range('gen', 425, 700)
gen_700up = dataset_from_pt_range('gen', 700)

def sim_noparticles():
    
    dset = mod.load(dataset='sim', amount=-1, store_gens=False, store_pfcs=False, verbose=1)
    dset.save(utils.path('sim', 'h5', 'SIM_Jet300_pT375-infGeV_noparticles'))

def gen_noparticles():

    dset = mod.load(dataset='gen', amount=-1, store_gens=False, store_pfcs=False, verbose=1)
    dset.save(utils.path('sim', 'h5', 'GEN_pT375-infGeV_noparticles'))

nums4emd = [5000, 40000]

def cms_subsample_datasets_for_emds():

    for num4emd in nums4emd:
        num4emd_str = str(num4emd//1000) + 'k'
        subs = np.load(utils.path('sim', 'emds', 'SubsampleWeightsInds_{}.pickle'.format(num4emd_str)), allow_pickle=True)

        for k,v in subs.items():
            specs = v['sim_specs']

            # load CMS dataset
            cms = mod.MODDataset(utils.path('cms', 'h5', 'CMS_Jet300_pT375-infGeV'), *specs,
                                 num=num4emd, shuffle=True)

            print('CMS', k)

            # save new dataset
            filepath = utils.path('cms', 'h5', 'CMS_Jet300_pT{}-{}GeV_EtaMax{}_Quality{}_{}')

            ptmin, ptmax = specs[0][0], specs[0][2]
            absetamax = specs[1][1]
            quality = specs[2][1]
            cms.save(filepath.format(ptmin, ptmax, int(absetamax*10), quality, num4emd_str))
            del cms

def sim_subsample_datasets_for_emds():

    for num4emd in nums4emd:
        num4emd_str = str(num4emd//1000) + 'k'
        subs = np.load(utils.path('sim', 'emds', 'SubsampleWeightsInds_{}.pickle'.format(num4emd_str)), allow_pickle=True)

        for ptrange,v in subs.items():
            start = time.time()
            specs = v['sim_specs']

            sim = mod.MODDataset(utils.path('sim', 'h5', 'SIM_Jet300_pT{}-{}GeV'.format(*ptrange)), *specs)

            mask = np.zeros(len(sim), dtype=bool)
            mask[v['sim_inds']] = True

            sim.apply_mask(mask)
            print('SIM {} processed in {:.3f}s'.format(ptrange, time.time() - start))

            sim.jets_f[:,sim.weight] = v['sim_sub_weights']

            filepath = utils.path('sim', 'h5', 'SIM_Jet300_pT{}-{}GeV_EtaMax{}_Quality{}_{}')
            ptmin, ptmax = specs[0][0], specs[0][2]
            absetamax = specs[1][1]
            quality = specs[2][1]
            sim.save(filepath.format(ptmin, ptmax, int(absetamax*10), quality, num4emd_str))
            del sim

def gen_subsample_datasets_for_emds():

    for num4emd in nums4emd:
        num4emd_str = str(num4emd//1000) + 'k'
        subs = np.load(utils.path('sim', 'emds', 'SubsampleWeightsInds_{}.pickle'.format(num4emd_str)), allow_pickle=True)

        for ptrange,v in subs.items():
            start = time.time()
            specs = v['gen_specs']

            gen = mod.MODDataset(utils.path('sim', 'h5', 'GEN_pT{}-{}GeV'.format(*ptrange)), *specs)

            mask = np.zeros(len(gen), dtype=bool)
            mask[v['gen_inds']] = True

            gen.apply_mask(mask)
            print('GEN {} processed in {:.3f}s'.format(ptrange, time.time() - start))

            gen.jets_f[:,gen.weight] = v['gen_sub_weights']

            filepath = utils.path('sim', 'h5', 'GEN_pT{}-{}GeV_EtaMax{}_{}')

            ptmin, ptmax = specs[0][0], specs[0][2]
            absetamax = specs[1][1]
            gen.save(filepath.format(ptmin, ptmax, int(absetamax*10), num4emd_str))
            del gen

if __name__ == '__main__':

    import mod_io
    import utils

    cms_375to425()
    cms_425to700()
    cms_700up()
    sim_375to425()
    sim_425to700()
    sim_700up()
    gen_375to425()
    gen_425to700()
    gen_700up()
    sim_noparticles()
    gen_noparticles()
    subsample_datasets_for_emds()

else:

    from . import mod_io
    from . import utils