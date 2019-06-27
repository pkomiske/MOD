import importlib
import os
import shutil

# path where data will be stored
DATAPATH = '/disk0'

# list of data collections
COLLECTIONS = {
    'CMS2011AJets': {
        'datadirs': ['JetPrimaryDataset', 'QCDSimDatasets']
    }
}

# get path to git repo
GITREPO = os.path.dirname(__file__)

# function to copy each data collection to location specified by DATAPATH
def copy_collections():

    # iterate over data collections
    for collection,opts in COLLECTIONS.items():

        # make file recording datapath for this collection
        with open(os.path.join(GITREPO, collection, 'DATAPATH'), 'w') as f:
            f.write(os.path.join(DATAPATH, collection))

        # determine paths
        src_path = os.path.join(GITREPO, collection)
        dst_path = os.path.join(DATAPATH, collection)

        # handle if dst_path exists
        if os.path.exists(dst_path):

            # if pointing to same location, skip
            if os.path.samefile(src_path, dst_path):
                print('No need to copy {}, already in place'.format(collection))
                print()
                continue

            response = input('{} already exists at {}, do you want to remove it and reinitialize? [y/n]: '.format(collection, dst_path))
            print()

            if response == 'y':

                # SAFEGUARD
                print('You actually do not want to remove this since it has important files!')
                continue

                #print('Removing {}'.format(dst_path))
                #shutil.rmtree(dst_path, ignore_errors=True)

            else:
                print('Leaving {} as is'.format(collection))
                print()
                continue

        # make collection directory
        os.mkdir(dst_path)

        for datadir in opts['datadirs']:

            data_src_path = os.path.join(src_path, datadir)
            data_dst_path = os.path.join(dst_path, datadir)

            # copy data collection to data path
            print('Copying {} to {}'.format(data_src_path, data_dst_path))
            shutil.copytree(data_src_path, data_dst_path)

            # make plotdata dir
            os.mkdir(os.path.join(data_dst_path, 'plotdata'))

        print()

def init_collections():

     # iterate over data collections
    for collection,opts in COLLECTIONS.items():

        # call init functions for this collection
        module = importlib.import_module(collection)

        for func in module.init_funcs:
            print('Running {}.{}'.format(collection, func))
            getattr(module, func)()

if __name__ == '__main__':

    # copy collection
    copy_collections()

    # initialize each collection
    init_collections()
