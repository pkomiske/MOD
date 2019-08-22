from . import datasets
from . import emds
from . import jet_to_h5
from . import lumis
from . import mod_to_jet
from . import notebooks

# since files and functions have the same names, need to put this here
__all__ = (datasets.__all__ + 
           emds.__all__ +
           jet_to_h5.__all__ + 
           lumis.__all__ + 
           mod_to_jet.__all__ + 
           notebooks.__all__)

from .datasets import *
from .emds import *
from .jet_to_h5 import *
from .lumis import *
from .mod_to_jet import *
from .notebooks import *

# this list controls which functions are run at initialization time
init_funcs = [

    #'process_lumisbyls', 
    #'batch_mod_to_jet',
    #'process_jet_primary_dataset_lbs',
    #'count_jet_primary_dataset_lbs',
    #'extract_sim_cross_sections',
    #'evaluate_triggers_nb',
    #'make_filename_arrays',
    #'cms_jets_to_npzs',
    #'sim_jets_to_npzs',
    #'cms_npzs_to_h5',
    #'sim_npzs_to_h5',
    #'gen_npzs_to_h5',
    #'sim_noparticles',
    #'gen_noparticles',
    #'cms_375up',
    #'sim_375to425',
    #'sim_425to700',
    #'sim_700up',
    #'gen_375to425',
    #'gen_425to700',
    #'gen_700up',
    #'evaluate_kfactor_nb',
    #'evaluate_emd_subsampling_nb',
    #'cms_subsample_datasets_for_emds',
    #'sim_subsample_datasets_for_emds',
    #'gen_subsample_datasets_for_emds',
    #'calc_emds',
    'calc_corrdims',
    #'calc_qg_corrdims'
]

__all__ += ['init_funcs']
