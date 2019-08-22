import os

import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

__all__ = ['evaluate_triggers_nb',
           'evaluate_kfactor_nb',
           'evaluate_emd_subsampling_nb']

def run_and_save_notebook(filename):

    path = os.path.dirname(__file__)
    filepath = os.path.join(path, filename)

    # read notebook
    with open(filepath) as f:
        nb = nbformat.read(f, as_version=4)

    # run
    ep = ExecutePreprocessor(timeout=None, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': path}})

    # save
    with open(filepath, 'w') as f:
        nbformat.write(nb, f)

def evaluate_triggers_nb():
    run_and_save_notebook('Triggers.ipynb')

def evaluate_fileinfo_nb():
    run_and_save_notebook('FileInfo.ipynb')

def evaluate_kfactor_nb():
    run_and_save_notebook('K-Factor.ipynb')

def evaluate_emd_subsampling_nb():
    run_and_save_notebook('EMD-Subsampling.ipynb')