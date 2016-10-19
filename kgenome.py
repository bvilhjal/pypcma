"""
Methods for analysing 1000 genomes data.
"""

from itertools import izip
import cPickle
import gzip
import os
import ld

from scipy import linalg
import h5py
import h5py_util as hu

import scipy as sp
import time

__updated__ = '2016-10-19'

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }
ok_nts = sp.array(['A', 'T', 'G', 'C'])
 

    
def gen_unrelated_eur_1k_data(input_file='/home/bjarni/TheHonestGene/faststorage/1Kgenomes/phase3/1k_genomes_hg.hdf5' ,
                              out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
                              maf_thres=0.01, max_relatedness=0.05, K_thinning_frac=0.1, debug=False):
    h5f = h5py.File(input_file)
    num_indivs = len(h5f['indivs']['continent'])
    eur_filter = h5f['indivs']['continent'][...] == 'EUR'
    num_eur_indivs = sp.sum(eur_filter)
    print 'Number of European individuals: %d', num_eur_indivs
    K = sp.zeros((num_eur_indivs, num_eur_indivs), dtype='float64')
    num_snps = 0
    std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))

    print 'Calculating kinship'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        print 'Loading SNPs and data'
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')

        print 'Loading NTs'
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]
        
        print 'Filtering multi-allelic SNPs'
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]


        if K_thinning_frac < 1:
            print 'Thinning SNPs for kinship calculation'
            thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
            snps = snps[thinning_filter]
            alt_nts = alt_nts[thinning_filter]
            ref_nts = ref_nts[thinning_filter]

        print 'Filter SNPs with missing NT information'
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]

        print 'Filtering non-European individuals'
        snps = snps[:, eur_filter]

        print 'Filtering SNPs with MAF <', maf_thres
        snp_stds = sp.std(snps, 1)
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        
        print '%d SNPs remaining after all filtering steps.' % len(snps)

        print 'Normalizing SNPs'
        snp_means = sp.mean(snps, 1)
        norm_snps = (snps - snp_means[sp.newaxis].T) / snp_stds[sp.newaxis].T
        
        print 'Updating kinship'        
        K += sp.dot(norm_snps.T, norm_snps)
        num_snps += len(norm_snps)
        assert sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0)

    K = K / float(num_snps)
    print 'Kinship calculation done using %d SNPs\n' % num_snps
    
    # Filter individuals
    print 'Filtering individuals'
    keep_indiv_set = set(range(num_eur_indivs))
    for i in range(num_eur_indivs):
        if i in keep_indiv_set:
            for j in range(i + 1, num_eur_indivs):
                if K[i, j] > max_relatedness:
                    if j in keep_indiv_set:
                        keep_indiv_set.remove(j)
    keep_indivs = list(keep_indiv_set)
    keep_indivs.sort()
    print 'Retained %d individuals\n' % len(keep_indivs)
    
    # Checking that everything is ok!
    K_ok = K[keep_indivs]
    K_ok = K_ok[:, keep_indivs]
    assert (K_ok - sp.tril(K_ok)).max() < max_relatedness

    indiv_filter = sp.zeros(num_indivs, dtype='bool8')
    indiv_filter[(sp.arange(num_indivs)[eur_filter])[keep_indivs]] = 1
    
    assert sp.sum(indiv_filter) == len(keep_indivs)
    
    # Store in new file
    print 'Now storing data.'
    oh5f = h5py.File(out_file, 'w')
    indiv_ids = h5f['indivs']['indiv_ids'][indiv_filter]
    oh5f.create_dataset('indiv_ids', data=indiv_ids)    
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        print 'Loading SNPs and data'
        snps = sp.array(h5f[chrom_str]['calldata']['snps'][...], dtype='int8')
        snp_ids = h5f[chrom_str]['variants']['ID'][...]
        positions = h5f[chrom_str]['variants']['POS'][...]

        print 'Loading NTs'
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]
        
        print 'Filtering multi-allelic SNPs'
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]
        positions = positions[multi_allelic_filter]
        snp_ids = snp_ids[multi_allelic_filter]

        print 'Filter individuals'
        snps = snps[:, indiv_filter]
        
        print 'Filter SNPs with missing NT information'
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]
            ref_nts = ref_nts[nt_filter]
            alt_nts = alt_nts[nt_filter]
            positions = positions[nt_filter]
            snp_ids = snp_ids[nt_filter]
        
        print 'filter monomorphic SNPs'
        snp_stds = sp.std(snps, 1)
        mono_morph_filter = snp_stds > 0
        snps = snps[mono_morph_filter]
        ref_nts = ref_nts[mono_morph_filter]
        alt_nts = alt_nts[mono_morph_filter]
        positions = positions[mono_morph_filter]
        snp_ids = snp_ids[mono_morph_filter]
        snp_stds = snp_stds[mono_morph_filter]

        snp_means = sp.mean(snps, 1)

        if debug:
            if K_thinning_frac < 1:
                print 'Thinning SNPs for kinship calculation'
                thinning_filter = sp.random.random(len(snps)) < K_thinning_frac
                k_snps = snps[thinning_filter]
                k_snp_stds = snp_stds[thinning_filter]

    
            print 'Filtering SNPs with MAF <', maf_thres
            maf_filter = k_snp_stds.flatten() > std_thres
            k_snps = k_snps[maf_filter]
            k_snp_stds = k_snp_stds[maf_filter]
            k_snp_means = sp.mean(k_snps)

            print 'Verifying that the Kinship makes sense'
            norm_snps = (k_snps - k_snp_means[sp.newaxis].T) / k_snp_stds[sp.newaxis].T
            K = sp.dot(norm_snps.T, norm_snps)
            num_snps += len(norm_snps)
            if sp.isclose(sp.sum(sp.diag(K)) / (num_snps * num_eur_indivs), 1.0) and (K - sp.tril(K)).max() < (max_relatedness * 1.5):
                print 'It looks OK!'
            else:
                raise Exception('Kinship looks wrong?')
        

        nts = sp.array([[nt1, nt2] for nt1, nt2 in izip(ref_nts, alt_nts)])

        print 'Writing to disk'
        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps', data=snps)
        cg.create_dataset('snp_means', data=snp_means[sp.newaxis].T)
        cg.create_dataset('snp_stds', data=snp_stds[sp.newaxis].T)
        cg.create_dataset('snp_ids', data=snp_ids)
        cg.create_dataset('positions', data=positions)
        cg.create_dataset('nts', data=nts)
        oh5f.flush()
        print 'Done writing to disk'
        
#         centimorgans = h5f[chrom_str]['centimorgans'][...]
#         cg.create_dataset('centimorgans',data=centimorgans)
#         
#         centimorgan_rates = h5f[chrom_str]['centimorgan_rates'][...]
#         cg.create_dataset('centimorgan_rates',data=centimorgan_rates)
        
    oh5f.close()
    h5f.close()
    print 'Done'
    
    
    
def get_kinship_pca_dict(input_genotype_file, kinship_pca_file, maf_thres, snp_filter_frac):
    if os.path.isfile(kinship_pca_file):
        print ':Loading Kinship and PCA information from %s' % kinship_pca_file
        k_h5f = h5py.File(kinship_pca_file)
        kinship_pca_dict = hu.hdf5_to_dict(k_h5f)
    else:
        kinship_pca_dict = calc_kinship(input_file=input_genotype_file , out_file=kinship_pca_file,
                                                maf_thres=maf_thres, figure_dir=None, snp_filter_frac=snp_filter_frac)
    return kinship_pca_dict

    



def get_genotype_data(in_h5f, chrom_i, maf_thres=0, indiv_filter=None,
                        snp_filter=None, randomize_sign=True, snps_signs=None,
                        return_raw_snps=False, return_snps_info=False,
                        return_normalized_snps=True, debug_filter_frac=1):
        
    chrom_str = 'chr%d' % chrom_i                    
    print 'Loading SNPs'
    snps = in_h5f[chrom_str]['snps'][...]
    
    if return_snps_info:
        positions = in_h5f[chrom_str]['positions'][...]
        snp_ids = in_h5f[chrom_str]['snp_ids'][...]
        nts = in_h5f[chrom_str]['nts'][...]
        
    if indiv_filter is not None:
        snps = snps[:, indiv_filter]
    
    if snp_filter is not None:
        snps = snps[snp_filter]        
        if return_snps_info:
            positions = positions[snp_filter]
            snp_ids = snp_ids[snp_filter]
            nts = nts[snp_filter]

    if debug_filter_frac < 1:
        debug_filter = sp.random.random(len(snps)) < debug_filter_frac
        snps = snps[debug_filter]        
        if return_snps_info:
            positions = positions[debug_filter]
            snp_ids = snp_ids[debug_filter]
            nts = nts[debug_filter]
    
    snp_means = sp.mean(snps, 1)
    snp_means.shape = (len(snp_means), 1)
    snp_stds = sp.std(snps, 1)
    snp_stds.shape = (len(snp_stds), 1)
    
    if maf_thres > 0:
        print 'Filtering SNPs with MAF <', maf_thres
        std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        snp_means = snp_means[maf_filter]
        if return_snps_info:
            positions = positions[maf_filter]
            snp_ids = snp_ids[maf_filter]
            nts = nts[maf_filter]
    
    print '%d SNPs remaining after all filtering steps.' % len(snps)    
    
    
    ret_dict = {'snp_stds':snp_stds, 'snp_means':snp_means}
    
    if return_normalized_snps:
        print 'Normalizing SNPs'
        norm_snps = sp.array((snps - snp_means) / snp_stds)
        ret_dict['norm_snps'] = norm_snps

    if randomize_sign:
        if snps_signs is None:
            snps_signs = 2 * sp.array(sp.random.random(len(norm_snps)) < 0.5, dtype='int8') - 1
            snps_signs.shape = (len(snps_signs), 1)
        norm_snps = norm_snps * snps_signs
        ret_dict['snps_signs'] = snps_signs
    
    if return_raw_snps:
        ret_dict['snps'] = snps
    
    if return_snps_info:
        ret_dict['positions'] = positions
        ret_dict['snp_ids'] = snp_ids
        ret_dict['nts'] = nts

    return ret_dict
    
    
def calc_kinship(input_file='Data/1Kgenomes/1K_genomes_v3.hdf5' , out_file='Data/1Kgenomes/kinship.hdf5',
                  maf_thres=0.01, figure_dir='', figure_fn='', snp_filter_frac=1, indiv_filter_frac=1):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print 'Loading Genotype from '
    in_h5f = h5py.File(input_file)
#     eur_filter = in_h5f['indivs']['continent'][...] == 'EUR'
#     num_indivs = sp.sum(eur_filter)
    indiv_ids = in_h5f['indiv_ids'][...] 
    indiv_filter = None
    if indiv_filter_frac < 1:
        indiv_filter = sp.array(sp.random.random(len(indiv_ids)) < indiv_filter_frac, dtype='bool8')
        indiv_ids = indiv_ids[indiv_filter]
    assert len(sp.unique(indiv_ids)) == len(indiv_ids)
    num_indivs = len(indiv_ids) 
    

    ok_chromosome_dict = {}

    not_done = set(range(1, 23))
    while len(not_done) > 0:
        chromosome_dict = {}
        
        K_all_snps = sp.zeros((num_indivs, num_indivs), dtype='float32')
        num_all_snps = 0
        
        sum_indiv_genotypes_all_chrom = sp.zeros(num_indivs, dtype='float32')
#         snp_cov_all_snps = sp.zeros((num_indivs, num_indivs), dtype='float64')
        
        print 'Calculating kinship'
        
        for chrom in range(1, 23):
            print 'Working on Chromosome %d' % chrom
            chrom_str = 'chr%d' % chrom
            
            snp_filter = None
            if snp_filter_frac < 1:
                snp_filter = sp.random.random(len(in_h5f[chrom_str]['snps'])) < snp_filter_frac

            g_dict = get_genotype_data(in_h5f, chrom, maf_thres, indiv_filter=indiv_filter,
                        snp_filter=snp_filter, randomize_sign=True, snps_signs=None)
            
            norm_snps = g_dict['norm_snps']
            
            sum_indiv_genotypes = sp.sum(g_dict['norm_snps'], 0)
            sum_indiv_genotypes_all_chrom += sum_indiv_genotypes
            
            print 'Calculating chromosome kinship'
            K_unscaled = sp.array(sp.dot(norm_snps.T, norm_snps), dtype='float32')
            assert sp.isclose(sp.sum(sp.diag(K_unscaled)) / (len(norm_snps) * num_indivs), 1.0), '..bug' 
            K_all_snps += K_unscaled
            num_all_snps += len(norm_snps)
    
            print 'SNP-cov normalisation'
            sum_indiv_genotypes = sp.sum(norm_snps, 0)
            sum_indiv_genotypes_all_chrom += sum_indiv_genotypes
            mean_indiv_genotypes = sum_indiv_genotypes / len(norm_snps)
            norm_snps = norm_snps - mean_indiv_genotypes
            
            print 'Calculating SNP covariance unscaled'
            
            snp_cov_unscaled = sp.array(sp.dot(norm_snps.T, norm_snps), dtype='float32')
#             snp_cov_all_snps += snp_cov_unscaled
            
            print 'Storing and updating things'
            chromosome_dict[chrom_str] = {'K_unscaled':K_unscaled, 'num_snps':len(norm_snps),
                                          'sum_indiv_genotypes':sum_indiv_genotypes,
                                          'snp_cov_unscaled':snp_cov_unscaled,
                                          'snps_signs':g_dict['snps_signs']}
            
            if snp_filter_frac < 1:
                chromosome_dict[chrom_str]['snp_filter'] = snp_filter
    
#         snp_cov_all_snps = snp_cov_all_snps / float(num_all_snps)
#         K_all_snps = K_all_snps / float(num_all_snps)
#         print 'K_all_snps.shape: %s' % str(K_all_snps.shape)
#         print 'snp_cov_all_snps.shape: %s' % str(snp_cov_all_snps.shape)
#         print 'sp.diag(snp_cov_all_snps): %s' % str(sp.diag(snp_cov_all_snps))
#         print 'sp.mean(sp.diag(snp_cov_all_snps)_: %s' % str(sp.mean(sp.diag(snp_cov_all_snps)))
        
#         print 'Full kinship and snp-covariance calculation done using %d SNPs\n' % num_all_snps
        
        mean_indiv_genotypes_all_chrom = sum_indiv_genotypes_all_chrom / num_all_snps
        print 'Individual gentoype mean found:'
        print mean_indiv_genotypes_all_chrom
        
        print 'Calculating chromosome-wise SNP-covariance and kinship matrices'
        for chrom in range(1, 23):
            if chrom in not_done:
                print 'Working on Chromosome %d' % chrom
                chrom_str = 'chr%d' % chrom
                
                snp_cov_leave_one_out = sp.zeros((num_indivs, num_indivs), dtype='float32')
                K_leave_one_out = sp.zeros((num_indivs, num_indivs), dtype='float32')
                num_snps_used = 0 
                
                sum_indiv_genotypes = sp.zeros(num_indivs, dtype='float32')
                
                for chrom2 in range(1, 23):
                    chrom2_str = 'chr%d' % chrom2
                    if chrom2 != chrom: 
                        sum_indiv_genotypes += chromosome_dict[chrom2_str]['sum_indiv_genotypes']
                        K_leave_one_out += chromosome_dict[chrom2_str]['K_unscaled']
                        num_snps_used += chromosome_dict[chrom2_str]['num_snps']
                        assert sp.isclose(sp.sum(sp.diag(K_leave_one_out)) / (num_snps_used * num_indivs), 1.0), '..bug' 
        
                mean_indiv_genotypes = sum_indiv_genotypes / num_snps_used
        
                for chrom2 in range(1, 23):
                    chrom2_str = 'chr%d' % chrom2
                    if chrom2 != chrom: 
                        print 'Loading SNPs'
                        snps_signs = chromosome_dict[chrom2_str]['snps_signs']
                        snp_filter = chromosome_dict[chrom2_str]['snp_filter']
                        g_dict = get_genotype_data(in_h5f, chrom2, maf_thres, indiv_filter=indiv_filter,
                                                   snp_filter=snp_filter, randomize_sign=True,
                                                   snps_signs=snps_signs)
                        norm_snps = g_dict['norm_snps']
                        print 'SNP-cov normalisation'
                        norm_snps = norm_snps - mean_indiv_genotypes
                        
                        print 'Calculating SNP covariance unscaled'
                        snp_cov_unscaled = sp.dot(norm_snps.T, norm_snps)
                        snp_cov_leave_one_out += snp_cov_unscaled
                  
                snp_cov_leave_one_out = snp_cov_leave_one_out / num_snps_used
                
                K_leave_one_out = K_leave_one_out / num_snps_used
                assert (K_leave_one_out - sp.diag(K_leave_one_out)).max() < 0.1, '..bug' 
                
                try:
                    cholesky_decomp_inv_snp_cov = linalg.cholesky(linalg.pinv(sp.array(snp_cov_leave_one_out, dtype='float64')))  
                    evals, evecs = linalg.eig(sp.array(K_leave_one_out, dtype='float64')) 
                except:
                    try: 
                        cholesky_decomp_inv_snp_cov = linalg.cholesky(linalg.pinv(sp.array(snp_cov_leave_one_out, dtype='float32')))
                        evals, evecs = linalg.eig(sp.array(K_leave_one_out, dtype='float32')) 
                    except:
                        print 'Failed when obtaining the Cholesky decomposotion or eigen decomposition'
                        print 'Moving on, trying again later.'
                        continue
                
                sort_indices = sp.argsort(evals,)
                ordered_evals = evals[sort_indices]
                print ordered_evals[-10:] / sp.sum(ordered_evals)
                ordered_evecs = evecs[:, sort_indices]
                d = {}
                d['evecs_leave_one_out'] = ordered_evecs
                d['evals_leave_one_out'] = ordered_evals
                d['cholesky_decomp_inv_snp_cov'] = cholesky_decomp_inv_snp_cov
                d['K_leave_one_out'] = K_leave_one_out
                d['K_unscaled'] = chromosome_dict[chrom_str]['K_unscaled']
                d['num_snps'] = chromosome_dict[chrom_str]['num_snps']
                d['snp_cov_leave_one_out'] = snp_cov_leave_one_out
                ok_chromosome_dict[chrom_str] = d
                not_done.remove(chrom)

    # While loop ends here.
    K_all_snps = K_all_snps / float(num_all_snps)
    in_h5f.close()

    assert sp.sum((ok_chromosome_dict['chr1']['K_leave_one_out'] - ok_chromosome_dict['chr2']['K_leave_one_out']) ** 2) != 0 , 'Kinships are probably too similar.'
        
    print 'Calculating PCAs'
    evals, evecs = linalg.eigh(sp.array(K_all_snps, dtype='float64'))  # PCA via eigen decomp
    evals[evals < 0] = 0
    sort_indices = sp.argsort(evals,)[::-1]
    ordered_evals = evals[sort_indices]
    print ordered_evals[:10] / sp.sum(ordered_evals)
    pcs = evecs[:, sort_indices]

    tot = sum(evals)
    var_exp = [(i / tot) * 100 for i in sorted(evals, reverse=True)]
    print 'Total variance explained:', sp.sum(var_exp)

    if figure_dir is not None:
        plt.clf()    
        plt.plot(pcs[:, 0], pcs[:, 1], 'k.')
        plt.title("Overall PCA")
        plt.xlabel('PC1')
        plt.xlabel('PC2')
        plt.tight_layout()
        plt.savefig(figure_dir + '/' + figure_fn, format='pdf')
        plt.clf()
    
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        K_leave_one_out = ok_chromosome_dict[chrom_str]['K_leave_one_out']

    out_h5f = h5py.File(out_file)
    hu.dict_to_hdf5(ok_chromosome_dict, out_h5f)
    out_h5f.close()
    
    return ok_chromosome_dict



def calc_structure_covar():
    pass


def ld_prune_1k_genotypes(in_hdf5_file, out_hdf5_file, local_ld_file_prefix, ld_radius, max_ld=0.2, maf_thres=0.01):
    # Open input and output file
    ih5f = h5py.File(in_hdf5_file)
    oh5f = h5py.File(out_hdf5_file)
    
    for chrom_i in range(1, 23):
        print 'Working on Chromosome %d' % chrom_i
        chrom_str = 'chr%d' % chrom_i
          
        g_dict = get_genotype_data(ih5f, chrom_i, maf_thres=maf_thres, randomize_sign=False, snps_signs=None,
                                   return_raw_snps=True, return_snps_info=True, return_normalized_snps=False)
        snps = g_dict['snps']
        snp_means = g_dict['snp_means']
        snp_stds = g_dict['snp_stds']
        snp_ids = g_dict['snp_ids'] 
        positions = g_dict['positions']
        nts = g_dict['nts']

        local_ld_hdf5_file = '%s_chrom%d_ldradius%d.hdf5' % (local_ld_file_prefix, chrom_i, ld_radius)
        print 'Loading LD information from file: %s' % local_ld_hdf5_file
        ldh5f = h5py.File(local_ld_hdf5_file, 'r')
        ld_dict = hu.hdf5_to_dict(ldh5f)
        ldh5f.close()

        ld_snp_filter = ld.ld_pruning(ld_dict, max_ld=max_ld, verbose=True)
        print ld_snp_filter
        
        assert ld_dict['num_snps'] == len(snps)
        assert ld_dict['num_snps'] == len(ld_snp_filter)
        
        # Pruning SNPs in LD
        snps = snps[ld_snp_filter]
        snp_means = snp_means[ld_snp_filter]
        snp_stds = snp_stds[ld_snp_filter]
        snp_ids = snp_ids[ld_snp_filter]
        positions = positions[ld_snp_filter]
        nts = nts[ld_snp_filter]

        print 'Out of %d SNPs %d were retained' % (ld_dict['num_snps'], len(snps))

        cg = oh5f.create_group(chrom_str)
        cg.create_dataset('snps', data=sp.array(snps, dtype='int8'))
        cg.create_dataset('snp_means', data=snp_means)
        cg.create_dataset('snp_stds', data=snp_stds)
        cg.create_dataset('snp_ids', data=snp_ids)
        cg.create_dataset('positions', data=positions)
        cg.create_dataset('nts', data=nts)

    indiv_ids = ih5f['indiv_ids'][...]
    oh5f.create_dataset('indiv_ids', data=indiv_ids)    
    
    ih5f.close()
    oh5f.close()

        
# For debugging purposes
if __name__ == '__main__':        
    kinship_pca_dict = calc_kinship('/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated_ld_pruned.hdf5',
                         out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1kgenomes_kinship_pca.hdf5',
                         figure_dir='/home/bjarni/tmp', figure_fn='test.pdf',
                         maf_thres=0.01, snp_filter_frac=1, indiv_filter_frac=1)


