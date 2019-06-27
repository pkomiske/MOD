from collections import Counter
import json
import os
import pickle
import time

import numpy as np

__all__ = ['process_lumisbyls', 
           'process_jet_primary_dataset_lbs', 
           'count_jet_primary_dataset_lbs',
           'extract_sim_cross_sections']

def process_lumisbyls():

    lumipath = os.path.join(utils.DATASET_PATHS['cms'], 'lumis')

    # names of files
    cert_file = 'Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt'
    lumi_infile = '2011lumibyls.csv'
    lumi_outfile = '2011Lumis.pickle'

    # read cert file
    cert_path = os.path.join(lumipath, cert_file)
    print('Reading cert file from {}'.format(cert_path))

    with open(cert_path, 'r') as f:
        cert_lumis = {int(k): v for k,v in json.load(f).items()}

    # read lumis file
    lumi_path = os.path.join(lumipath, lumi_infile)
    print('Reading lumis from {}'.format(lumi_path))

    lumis = {}
    with open(lumi_path, 'r') as f:
        for line in f:
            if line.startswith ('#'):
                print(line[:-1])
                continue
                
            parts = line.split(',')
            rn = int(parts[0].split(':')[0])
            lb = int(parts[1].split(':')[0])
            del_lumi = float(parts[5])
            rec_lumi = float(parts[6])
            t = time.mktime(time.strptime(parts[2][:6] + '20' + parts[2][6:], '%m/%d/%Y %H:%M:%S'))
            
            lumis[(rn, lb)] = [del_lumi, rec_lumi, t]
            
            # check cert file
            gl_ranges = cert_lumis[rn]
            good = False
            for gl_range in gl_ranges:
                if gl_range[0] <= lb <= gl_range[1]:
                    good = True
                    break
                    
            if not good:
                print(line, 'in lumisbyls but not good according to Cert')

    # pickle lumis dictionary
    lumi_path = os.path.join(lumipath, lumi_outfile)
    print('Pickling lumis to {}'.format(lumi_path))
    with open(lumi_path, 'wb') as f:
        pickle.dump(lumis, f)

    print()

def process_jet_primary_dataset_lbs():

    lumifilepath = os.path.join(utils.DATASET_PATHS['cms'], 'lumis', '2011Lumis.pickle')

    # load 2011lumisbyls
    lumis = np.load(lumifilepath, allow_pickle=True)

    # dict for lbs
    lbs = {}

    # iterate over files
    start = time.time()
    filepaths = utils.get_filenames(subdir='pickle', include_path=True, must_include='noparticles')
    for i,filepath in enumerate(filepaths):

        # load file
        data = np.load(filepath, allow_pickle=True)

        # iterate over lumiblocks in this file
        for lb in data['lbs']:
            rnlb = (lb['rn'], lb['lb'])
            
            # skip lbs not in lumis
            if rnlb not in lumis:
                continue
                
            # add lumi info
            lb['del_lumi'], lb['rec_lumi'], lb['time'] = lumis[rnlb]
            
            # add file info
            file = os.path.basename(filepath)
            lb['file'] = file.split('.')[0]
            
            # check for existing lb
            if rnlb in lbs:
                print('{} with {} events in files {} and {}'.format(rnlb, lb['nev'], file, lbs[rnlb]['file']))
                print('    stored as {}'.format(rnlb + (1,)))
                lbs[rnlb + (1,)] = lb
            else:
                lbs[rnlb] = lb
            
        # check that each event had a trigger fire
        ntrs = [len(trs) for trs in data['ftrigs']]
        if len(ntrs) and min(ntrs) == 0:
            print('file {} has an event where no triggers fired'.format(file))

        if (i + 1) % 100 == 0:
            print('  [{}/{}] files processed in {:.3f}s'.format(i + 1, len(filepaths), time.time() - start))


    # handle the one split lumi block (verified all these operations make sense previously)
    assert lbs[(163255, 603)]['triggers'] == lbs[(163255, 603, 1)]['triggers']
    lbs[(163255, 603)]['counts'] += lbs[(163255, 603, 1)]['counts']
    lbs[(163255, 603)]['nev'] += lbs[(163255, 603, 1)]['nev']
    lbs[(163255, 603)]['additional_files'] = [lbs[(163255, 603)]['file']]
    lbs[(163255, 603)]['file'] = lbs[(163255, 603, 1)]['file']
    del lbs[(163255, 603, 1)]

    lbvalfp = os.path.join(utils.DATASET_PATHS['cms'], 'lumis', 'ValidatedLumiblocks.pickle')
    print('Saving CMS Jet Primary Dataset lumiblock information to {}'.format(lbvalfp))
    with open(lbvalfp, 'wb') as f:
        pickle.dump(lbs, f)

def count_jet_primary_dataset_lbs():

    lumipath = os.path.join(utils.DATASET_PATHS['cms'], 'lumis')

    lumis = np.load(os.path.join(lumipath, '2011Lumis.pickle'), allow_pickle=True)
    lbs = np.load(os.path.join(lumipath, 'ValidatedLumiblocks.pickle'), allow_pickle=True)

    min_2011a_rn, max_2011a_rn = 160404, 173692

    rns = set([k[0] for k in lbs.keys()])
    print('Number of unique run numbers:', len(rns))

    # missing lumi blocks
    missing_lbs = {}
    for rnlb,lumi in lumis.items():

        # it's in run A
        if min_2011a_rn <= rnlb[0] and rnlb[0] <= max_2011a_rn:
            
            # but not in the jet primary dataset
            if rnlb not in lbs:
                missing_lbs[rnlb] = lumi
    print('Number of missing lumi blocks:', len(missing_lbs))

    missing_rec_lumis = [lb[1] for lb in missing_lbs.values()]
    print('Missing lumiblocks sorted by recorded luminosity:')
    print(sorted(missing_rec_lumis))
    print('\nNumber of missing blocks that have zero recorded luminosity:', Counter(missing_rec_lumis)[0])

    mlbs = np.asarray(list(missing_lbs.values()))
    nz_mlbs = mlbs[mlbs[:,1] > 0]
    x = np.sum(nz_mlbs[:,0]/nz_mlbs[:,1] > 10) + np.sum((mlbs[:,0] > 0) & (mlbs[:,1] == 0))
    print('Number of missing lumi blocks where recorded was order of magnitude less than delivered:', x)

    y = np.sum((mlbs[:,0] == 0) & (mlbs[:,1] == 0))
    print('Number of missing lumi blocks where both delivered and reocrded is zero:', y)

    zero_rl_lbs = {}
    for rnlb,lumi in lumis.items():

        # it's in run A
        if min_2011a_rn <= rnlb[0] and rnlb[0] <= max_2011a_rn:
            
            # and has zero recorded luminosity
            if lumi[1] == 0:
                zero_rl_lbs[rnlb] = lumi

    print('Number of lumiblocks in run 2011A with zero recorded lumi:', len(zero_rl_lbs))

    zero_rl_present_lbs = [rnlb for rnlb in lbs if rnlb in zero_rl_lbs]
    print('Number of lumiblocks present in run 2011A with zero recorded lumi:', len(zero_rl_present_lbs))
    print()

def extract_sim_cross_sections():

    sim_filepaths = utils.get_sim_filenames_dict(subdir='pickle', include_path=True, must_include='noparticles')

    lbs = {}
    for ptmin,filepaths in sim_filepaths.items(): 
        lbs[ptmin] = {}
        
        for filepath in filepaths:
            lbs[ptmin][os.path.basename(filepath)] = np.load(filepath, allow_pickle=True)['lbs']

    cross_sections = {}
    for ptmin,lbs_pt in lbs.items():
        cross_sections[ptmin] = {'xs': set(), 'nev': 0}

        for filename,lbs_f in lbs_pt.items():
            for lb in lbs_f:

                cross_sections[ptmin]['xs'].add(lb['cross_section'])
                cross_sections[ptmin]['nev'] += lb['nev']
                
        xs = cross_sections[ptmin]['xs'].pop()
        if len(cross_sections[ptmin]['xs']):
            raise ValueError('cross sections differ across lumi block in the same sim dataset!')
        cross_sections[ptmin]['xs'] = xs

    fp = os.path.join(utils.DATASET_PATHS['sim'], 'CrossSections.pickle')
    print('SIM Cross Sections by ptmin:', cross_sections)
    print('Saved to', fp)
    with open(fp, 'wb') as f:
        pickle.dump(cross_sections, f)

if __name__ == '__main__':

    import utils

    #process_lumisbyls()
    #process_jet_primary_dataset_lbs()
    #count_jet_primary_dataset_lbs()
    extract_sim_cross_sections()

else:

    from . import utils
