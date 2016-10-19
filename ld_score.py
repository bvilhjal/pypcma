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
import ld
import time

from numpy import linalg 
from plinkio import plinkfile
import h5py
import plinkio
import kgenome
import h5py_util as hu


try: 
    import scipy as sp
except Exception:
    import numpy as sp
    print 'Using Numpy instead of Scipy.'
    
ok_nts = sp.array(['A', 'T', 'G', 'C'])



def generate_1k_LD_scores(input_genotype_file, chrom_snp_trans_mats,
                          maf_thres=0.01, ld_radius=200, debug_filter_frac=0.01,
                          indiv_filter=None , snp_filter=None):
    """
    Generates 1k genomes LD scores and stores in the given file
    """
    
    chrom_ld_scores_dict = {}
    ld_score_sum = 0
    struct_adj_ld_score_sum = 0
    num_snps = 0
    print 'Calculating LD information w. radius %d' % ld_radius
    in_h5f = h5py.File(input_genotype_file)
        
    print 'Calculating local LD'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        
        print 'Loading SNPs'
        g_dict = kgenome.get_genotype_data(in_h5f, chrom, maf_thres, indiv_filter=indiv_filter,
                                   snp_filter=snp_filter, randomize_sign=False, snps_signs=None,
                                   return_snps_info=True, debug_filter_frac=debug_filter_frac)
        
        norm_snps = g_dict['norm_snps']
        
        ret_dict = ld.get_ld_scores(norm_snps, ld_radius=ld_radius)
        avg_ld_score = sp.mean(ret_dict['ld_scores'])
        g_dict['ld_scores'] = ret_dict['ld_scores']
        g_dict['avg_ld_score'] = avg_ld_score
        ld_score_sum += sp.sum(ret_dict['ld_scores'])
        
        
        print 'Un-adjusted average LD score was: %0.3f' % avg_ld_score

        if chrom_snp_trans_mats is not None:
            snp_trans_mat = chrom_snp_trans_mats[chrom_str]
            norm_snps = sp.dot(norm_snps, snp_trans_mat.T)
            
            # Need to re-normalize?
            snp_means = sp.mean(norm_snps, 1)
            snp_means.shape = (len(snp_means), 1)
            snp_stds = sp.std(norm_snps, 1)
            snp_stds.shape = (len(snp_stds), 1)
            norm_snps = sp.array((norm_snps - snp_means) / snp_stds)

            ret_dict = ld.get_ld_scores(norm_snps, ld_radius=ld_radius)
            
            avg_ld_score = sp.mean(ret_dict['ld_scores'])
            print 'Pop-structure adjusted average LD score was: %0.3f' % avg_ld_score

            g_dict['struct_adj_ld_scores'] = ret_dict['ld_scores']
            g_dict['avg_struct_adj_ld_score'] = avg_ld_score
            struct_adj_ld_score_sum += sp.sum(ret_dict['ld_scores'])
        
        del g_dict['norm_snps']
        del g_dict['snp_means']
        del g_dict['snp_stds']
        chrom_ld_scores_dict[chrom_str] = g_dict
        num_snps += len(norm_snps)
    
    avg_gw_ld_score = ld_score_sum / float(num_snps)
    avg_gw_struct_adj_ld_score = ld_score_sum / float(num_snps)
    ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score,
                      'avg_gw_struct_adj_ld_score':avg_gw_struct_adj_ld_score,
                      'chrom_dict':chrom_ld_scores_dict}    
    
    print 'Done calculating the LD table and LD scores.'
    return ld_scores_dict
    
    

def calculate(input_genotype_file, input_ld_pruned_genotype_file,
              ld_score_file, kinship_pca_file, ld_radius=200, maf_thres=0.01, snp_filter_frac=0.05, debug_filter_frac=1):
    """
    Generates population structure adjusted 1k genomes LD scores and stores in the given file.
    """
    
    # Kinship
    kinship_pca_dict = kgenome.get_kinship_pca_dict(input_ld_pruned_genotype_file, kinship_pca_file,
                                                    maf_thres=maf_thres, snp_filter_frac=snp_filter_frac)
    chrom_snp_trans_mats = {}
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        chrom_snp_trans_mats[chrom_str] = kinship_pca_dict[chrom_str]['cholesky_decomp_inv_snp_cov']    
    
    # bla
    lds_dict = generate_1k_LD_scores(input_genotype_file, chrom_snp_trans_mats,
                                    maf_thres=maf_thres, ld_radius=ld_radius,
                                    debug_filter_frac=debug_filter_frac)
    
    
    # Store LD scores
    lds_h5f = h5py.File(ld_score_file)
    hu.dict_to_hdf5(lds_dict, lds_h5f)
    lds_h5f.close()
    

    
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
    





if __name__ == '__main__':
    calculate('/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated_ld_pruned.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/ld_scores.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1kgenomes_kinship_pca_f0.95.hdf5',
              ld_radius=1000, maf_thres=0.01, snp_filter_frac=0.95, debug_filter_frac=0.02)
