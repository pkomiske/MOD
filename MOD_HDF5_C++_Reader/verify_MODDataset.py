import os
import subprocess

import energyflow as ef
import numpy as np

def read_jets_file(filepath, dtype):
    with open(filepath, 'r') as f:
        return np.asarray([row.split() for row in f], dtype=dtype)
    
def read_particles_file(filepath):
    with open(filepath, 'r') as f:
        events, event = [], []
        for row in f:
            parts = row.split()
            if len(parts):
                event.append(parts)
            else:
                events.append(np.asarray(event, dtype=float))
                event = []
        return np.asarray(events)
    
def compare_float_arrays(arr0, arr1, tol=10**-11, reg=10**-50):
    return np.all(2*np.abs(arr0 - arr1)/np.abs(arr0 + arr1 + reg) <= tol)

def compare_particle_arrays(arr0, arr1, tol=10**-11):
    return np.all([compare_float_arrays(x, y, tol=tol) for x,y in zip(arr0, arr1)])

def match_files(modds, outpath, index, verbose=0):

    jets_i = read_jets_file(os.path.join(outpath, 'jets_i_{}.txt'.format(index)), int)
    jets_f = read_jets_file(os.path.join(outpath, 'jets_f_{}.txt'.format(index)), float)

    # adjust weights (since energyflow adjusts this automatically to get the cross section right)
    jets_f[:,-1] *= modds.jets_f[0,-1]/jets_f[0,-1]

    jets_i_match = np.all(jets_i == modds.jets_i)
    jets_f_match = compare_float_arrays(jets_f, modds.jets_f)

    if modds.store_pfcs:
        pfcs = read_particles_file(os.path.join(outpath, 'pfcs_{}.txt'.format(index)))
        pfcs_match  = compare_particle_arrays(pfcs, modds.pfcs)
    else:
        pfcs_match = True

    if modds.store_gens:
        gens = read_particles_file(os.path.join(outpath, 'gens_{}.txt'.format(index)))
        gens_match  = compare_particle_arrays(gens, modds.gens)
    else:
        gens_match = True
    
    all_match = jets_i_match and jets_f_match and pfcs_match and gens_match

    if verbose > 0:
        print('jets_i match:', jets_i_match)
        print('jets_f match:', jets_f_match)
        if modds.store_pfcs:
            print('pfcs match:', pfcs_match)
        if modds.store_gens:
            print('gens match:', gens_match)
        print('all match:', all_match)

    return all_match

if __name__ == '__main__':

    # get mod file
    modds = ef.mod.load(amount=1, dataset='sim', subdatasets=['SIM300_Jet300_pT375-infGeV'], verbose=1)

    # compile example
    command = 'g++ MODDatasetTest.cc -o MODDatasetTest -std=c++11 -O3 -Wall -g -lhdf5_cpp -lhdf5'
    print(command)
    subprocess.run(command.split())

    # run example
    outpath = os.path.expanduser('~/.energyflow/datasets')
    h5fp = os.path.join(outpath, 'CMS2011AJets', 'SIM300_Jet300_pT375-infGeV', 'SIM300_Jet300_pT375-infGeV_0_compressed.h5')
    command = './MODDatasetTest {} {} {}'.format(h5fp, outpath, 0)
    print(command)
    subprocess.run(command.split())

    # verify file
    match_files(modds, outpath, 0, verbose=1)
