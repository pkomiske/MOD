import multiprocessing
import pickle
import time

import energyflow as ef
import numpy as np

__all__ = ['make_filename_arrays']

def make_filename_arrays():

    # cms
    cms_filenames = utils.get_filenames(subdir='jet')

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
    sim_names_dict = utils.get_sim_filenames_dict(subdir='jet')
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



if __name__ == '__main__':

    import mod_io
    import utils

    make_filename_arrays()

else:

    from . import mod_io
    from . import utils