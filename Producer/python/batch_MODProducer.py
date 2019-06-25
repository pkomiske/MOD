import os
import sys
import multiprocessing
import subprocess

def cmsRun_MODProducer(filename):
    commands = 'cmsRun python/MODProducer.py ' + filename
    print commands
    subprocess.call(commands.split())

def cmsRun_MODProducer_sim(filename):
    commands = 'cmsRun python/MODProducer_sim.py ' + filename
    print commands
    subprocess.call(commands.split())

if __name__ == '__main__':

    i = sys.argv[1]
    nprocs = int(sys.argv[2]) if len(sys.argv) > 2 else 3

    with open('/media/sf_disk2/CMS2011A_QCDSimDatasets/aux/CMS2011A_QCDSim_' + i + '.txt', 'r') as f:
        files = [line.strip() for line in f]

    pool = multiprocessing.Pool(processes=nprocs)
    pool.map(cmsRun_MODProducer_sim, files, chunksize=1)
