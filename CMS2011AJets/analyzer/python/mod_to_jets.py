import multiprocessing
import os
import subprocess
import time

# options controlling what this module does
do_cms = True
do_sim = True
do_dijets_analyzer = False
save_everything = True
save_noparticles = True

__all__ = ['batch_mod_to_jets']

dijets_analyzer = os.path.relpath(os.path.join(os.path.dirname(__file__), os.pardir, 'bin', 'dijets_analyzer'))

def run_mod_to_jets(args):

    filename, dataset, ptmin = args
    path = os.path.join(utils.DATASET_PATHS[dataset])
    
    indir = os.path.join(path, 'mod', '' if ptmin is None else ptmin)
    outdir = os.path.join(path, 'jets', '' if ptmin is None else ptmin)
    os.makedirs(outdir, exist_ok=True)
    
    if do_dijets_analyzer:
        command = '{} {}.mod {}.jets {}'.format(dijets_analyzer, os.path.join(indir, filename), 
                                                os.path.join(outdir, filename), dataset)
        print(command)
        subprocess.run(command.split())

    r = mod_io.MODReader(filename, path=outdir, dataset=dataset, ptmin=ptmin, fileformat='jets', store_particles=save_everything)
    nevs = (r.nev_good, r.nev_bad, r.nev_total)

    pickle_path = os.path.join(path, 'pickle', '' if ptmin is None else ptmin)
    if save_everything:
        r.save(path=pickle_path)

    if save_noparticles:
        r.save(path=pickle_path, endstr='_noparticles', save_particles=False)

    return nevs

def analyze_dataset(args, name=''):

    print('Analyzing dataset', name)
    with multiprocessing.Pool(processes=None) as pool:
        
        start = time.time()
        good_events = bad_events = total_events = 0
        for i,nevs in enumerate(pool.imap_unordered(run_mod_to_jets, args, chunksize=1)):
            if (i+1) % 100 == 0:
                print('{0}  [{1}/{2}] done in {3:.3f}s{0}'.format('\n' if do_dijets_analyzer else '', i+1, len(args), time.time() - start))

            good_events += nevs[0]
            bad_events += nevs[1]
            total_events += nevs[2]

        print('Good Events:', good_events)
        print('Bad Events:', bad_events)
        print('Total Events:', total_events)
        print('Done in {:.3f}s'.format(time.time() - start))
        print()

def batch_mod_to_jets():

    if do_cms:
        cms_filenames = utils.get_filenames(remove_ending=True)
        cms_args = [(filename, 'cms', None) for filename in cms_filenames]
        analyze_dataset(cms_args, name='CMS')

    if do_sim:
        sim_filenames = utils.get_sim_filenames_dict(remove_ending=True)
        sim_args = {ptmin: [(filename, 'sim', ptmin) for filename in filenames] for ptmin,filenames in sim_filenames.items()}
        for ptmin,args in sim_args.items():
            analyze_dataset(args, name='SIM {}'.format(ptmin))

if __name__ == '__main__':

    import utils
    import mod_io

    batch_dijets_analyzer()

else:

    from . import utils
    from . import mod_io
