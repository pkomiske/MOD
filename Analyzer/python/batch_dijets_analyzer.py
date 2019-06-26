import multiprocessing
import os
import subprocess
import time

import utils
import umod_io

cms_filepaths = utils.get_filenames(dataset='cms', include_path=True)
sim_filepaths = utils.get_sim_filenames_dict(include_path=True)

cms_args = [(filepath, 'cms', None) for filepath in cms_filepaths]
sim_args = {ptmin: [(filepath, 'sim', ptmin) for filepath in filepaths] for ptmin,filepaths in sim_filepaths.items()}

def run_dijets_analyzer(args):
    filepath, dataset, ptmin = args
    outfilepath = filepath.replace('umod', 'jets')
    
    command = '../bin/dijets_analyzer {} {} {}'.format(filepath, outfilepath, dataset)
    #print(command)
    #subprocess.run(command.split())

    umod_io.uMODReader(os.path.basename(outfilepath), dataset=dataset, ptmin=ptmin, fileformat='jets').save()

def analyze_dataset(args, name=''):
    print('Analyzing dataset', name)
    with multiprocessing.Pool(processes=None) as pool:
        start = time.time()
        for i,_ in enumerate(pool.imap_unordered(run_dijets_analyzer, args, chunksize=1)):
            if (i+1) % 100 == 0:
                print('  [{}/{}] done in {:.3f}s'.format(i+1, len(args), time.time() - start))

        print('Done in {:.3f}s'.format(time.time() - start))
        print()

# CMS dataset
analyze_dataset(cms_args, name='cms')

# SIM datasets
for ptmin,args in sim_args.items():
    analyze_dataset(args, name='sim {}'.format(ptmin))
