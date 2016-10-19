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


try: 
    import scipy as sp
except Exception:
    import numpy as sp
    print 'Using Numpy instead of Scipy.'
    
ok_nts = sp.array(['A', 'T', 'G', 'C'])



def generate_1k_LD_scores(input_genotype_file, chrom_snp_trans_mats,
                          gm_ld_radius=None, maf_thres=0.01, ld_radius=200, debug_filter=0.01):
    """
    Generates 1k genomes LD scores and stores in the given file
    """
    
    chrom_ld_scores_dict = {}
    if gm_ld_radius is not None:
        chrom_ld_boundaries = {}
    ld_score_sum = 0
    struct_adj_ld_score_sum = 0
    num_snps = 0
    print 'Calculating LD information w. radius %d' % ld_radius

    in_h5f = h5py.File(input_genotype_file)
    indiv_ids = in_h5f['indiv_ids'][...] 
    
    std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))
    
    print 'Calculating local LD'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        
        print 'Loading SNPs'
        snps = in_h5f[chrom_str]['snps'][...]
        snp_stds = in_h5f[chrom_str]['snp_stds'][...]
        snp_means = in_h5f[chrom_str]['snp_means'][...]
        nts = in_h5f[chrom_str]['nts'][...]
        print 'snps.shape: %s, snp_stds.shape: %s, snp_means.shape: %s' % (str(snps.shape), str(snp_stds.shape), str(snp_means.shape))
        
        if debug_filter < 1:
            debug_snp_filter = sp.random.random(len(snps)) < debug_filter
            snps = snps[debug_snp_filter]
            snp_stds = snp_stds[debug_snp_filter]
            snp_means = snp_means[debug_snp_filter]
            nts = nts[debug_snp_filter]

        
        print 'Filtering SNPs with MAF <', maf_thres
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        snp_means = snp_means[maf_filter]
        nts = nts[maf_filter]
        
        nt_filter = sp.in1d(nts, ok_nts)
        if not sp.all(nt_filter):
            print 'Removing SNPs with missing NT information'
            snps = snps[nt_filter]
            snp_stds = snp_stds[nt_filter]
            snp_means = snp_means[nt_filter]
        
        print '%d SNPs remaining' % len(snps)
        
        print 'Normalizing SNPs'
        norm_snps = (snps - snp_means) / snp_stds
    
        
        ret_dict = ld.get_ld_scores(norm_snps, ld_radius=ld_radius)
        avg_ld_score = sp.mean(ret_dict['ld_scores'])
        chrom_ld_scores_dict[chrom_str] = {'ld_scores':ret_dict['ld_scores'], 'avg_ld_score':avg_ld_score}
        ld_score_sum += sp.sum(ret_dict['ld_scores'])
        
        print 'Un-adjusted average LD score was: %0.3f' % avg_ld_score

        if chrom_snp_trans_mats is not None:
            snp_trans_mat = chrom_snp_trans_mats[chrom_str]
            norm_snps = sp.dot(norm_snps, snp_trans_mat.T)
    
            # Need to re-normalize?
            snp_means = sp.mean(snps, 1)
            snp_means.shape = (len(snp_means), 1)
            snp_stds = sp.std(snps, 1)
            snp_stds.shape = (len(snp_stds), 1)
            norm_snps = sp.array((norm_snps - snp_means) / snp_stds)

            ret_dict = ld.get_ld_scores(norm_snps, ld_radius=ld_radius)
            
            avg_ld_score = sp.mean(ret_dict['ld_scores'])
            print 'Pop-structure adjusted average LD score was: %0.3f' % avg_ld_score

            chrom_ld_scores_dict[chrom_str]['struct_adj_ld_scores'] = ret_dict['ld_scores']
            chrom_ld_scores_dict[chrom_str]['avg_struct_adj_ld_score'] = avg_ld_score
            struct_adj_ld_score_sum += sp.sum(ret_dict['ld_scores'])
            
        num_snps += len(norm_snps)
    
    avg_gw_ld_score = ld_score_sum / float(num_snps)
    ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score, 'chrom_dict':chrom_ld_scores_dict}    
    
    print 'Done calculating the LD table and LD scores.'
    if gm_ld_radius is not None:
        ld_scores_dict['chrom_ld_boundaries'] = chrom_ld_boundaries 
        
    return ld_scores_dict
    
    
def get_popadj_snp_trans_mat(kinship):
    return 

def calculate(input_genotype_file, input_ld_pruned_genotype_file,
              pca_adj_ld_score_file, ld_score_file, kinship_pca_file,
              ld_radius=200, maf_thres=0.01, snp_filter_frac=0.05):
    """
    Generates population structure adjusted 1k genomes LD scores and stores in the given file.
    """
    
    kinship_pca_dict = kgenome.get_kinship_pca_dict(input_ld_pruned_genotype_file, kinship_pca_file,
                                                    maf_thres=maf_thres, snp_filter_frac=snp_filter_frac)
    chrom_snp_trans_mats = {}
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        chrom_snp_trans_mats[chrom_str] = kinship_pca_dict[chrom_str]['cholesky_decomp_inv_snp_cov']    
    
    ld_dict = generate_1k_LD_scores(input_genotype_file, chrom_snp_trans_mats,
                                    maf_thres=maf_thres, ld_radius=ld_radius, debug_filter=1)
    
    
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
    





if __name__ == '__main__':
    calculate('/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated_ld_pruned.hdf5'
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/pca_adj_ld_scores.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/ld_scores.hdf5',
              '/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1kgenomes_kinship_pca_f0.95.hdf5',
              ld_radius=1000, maf_thres=0.01, snp_filter_frac=0.95)
