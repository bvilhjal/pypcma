"""
Re-implementation of LD score regression (Bulik-Sullivan et al., Nat Genet 2015).

1. Implement LD calculation in the 1K genomes
    - With and without adjustment for population Structure.
    - Store local LD structure.
2. Implement LD score regression for a given set of SNPs.
3. Implement heritability estimation.
4. Implement genetic correlation estimation.
"""
from itertools import izip
import cPickle
import os
import pdb
import random
import sys
import time

from numpy import linalg 
from plinkio import plinkfile
import h5py
import plinkio
import kgenome


try: 
    import scipy as sp
except Exception:
    import numpy as sp
    print 'Using Numpy instead of Scipy.'
    




def generate_1k_LD_scores(input_genotype_file, output_file, ld_radius=200,):
    """
    Generates 1k genomes LD scores and stores in the given file
    """
    pass


def pre_calculate_everything(input_genotype_file, output_file, kinship_pca_file, ld_radius=200, ancestry='EUR'):
    """
    Generates population structure adjusted 1k genomes LD scores and stores in the given file.
    """
    
    
    
    # 6. a) Calculate LD score.
    # 6. b) Calculate population structure adjusted LD score.
    # 7. Store everything.  EVERYTHING!
    
    

    
def estimate_heritability():
    pass

def estimate_genetic_correlation():   
    pass 
    
    
def generate_genetic_corr_figure():
    pass

def generate_herit_figure():
    pass



def get_1K_genomes_LD_scores(snps, ld_radius=200):
    """
    Calculates LD tables, and the LD scores in one go...
    """
    
    # First estimate the genotype covariance matrix with and
    
    pass
#     return ret_dict


def get_genotype_cov_mat():
    """
    Estimate genotype covariance matrix from the provided SNPs.
    """
    

def get_ld_tables(snps, ld_radius=500, ld_window_size=0, gm=None, gm_ld_radius=None):
    """
    Calculates LD tables, and the LD scores in one go, with or without a genetic map.
    """
    
    ld_dict = {}
    m, n = snps.size
    print m, n
    ld_scores = sp.ones(m)
    ret_dict = {}
    if gm_ld_radius is None:
        for snp_i, snp in enumerate(snps):
            # Calculate D
            start_i = max(0, snp_i - ld_radius)
            stop_i = min(m, snp_i + ld_radius + 1)
            X = snps[start_i: stop_i]
            D_i = sp.dot(snp, X.T) / n
            r2s = D_i ** 2
            ld_dict[snp_i] = D_i
            lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
            # lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
            ld_scores[snp_i] = lds_i
    else:
        assert gm is not None, 'Genetic map is missing.'
        window_sizes = []
        ld_boundaries = []
        for snp_i, snp in enumerate(snps):
            curr_cm = gm[snp_i] 
            
            # Now find lower boundary
            start_i = snp_i
            min_cm = gm[snp_i]
            while start_i > 0 and min_cm > curr_cm - gm_ld_radius:
                start_i = start_i - 1
                min_cm = gm[start_i]
            
            # Now find the upper boundary
            stop_i = snp_i
            max_cm = gm[snp_i]
            while stop_i > 0 and max_cm < curr_cm + gm_ld_radius:
                stop_i = stop_i + 1
                max_cm = gm[stop_i]
            
            ld_boundaries.append([start_i, stop_i])    
            curr_ws = stop_i - start_i
            window_sizes.append(curr_ws)
            assert curr_ws > 0, 'Some issues with the genetic map'

            X = snps[start_i: stop_i]
            D_i = sp.dot(snp, X.T) / n
            r2s = D_i ** 2
            ld_dict[snp_i] = D_i
            lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
#             lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
            ld_scores[snp_i] = lds_i
        
        avg_window_size = sp.mean(window_sizes)
        print 'Average # of SNPs in LD window was %0.2f' % avg_window_size
        if ld_window_size == 0:
            ld_window_size = avg_window_size * 2
        ret_dict['ld_boundaries'] = ld_boundaries
    ret_dict['ld_dict'] = ld_dict
    ret_dict['ld_scores'] = ld_scores

    return ret_dict

