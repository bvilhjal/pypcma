"""
"""
import h5py
import scipy as sp
import kgenome
import h5py_util

def pc_score_regression(sum_stats_file, pc_file):

    
    # 1. Get summary stats z-scores
    # 2. Get PC weights
    # 3. Perform regression
    
    pass


def calc_pc_snp_weights(input_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
                        pc_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1kgenomes_kinship_pca_f0.95.hdf5',
                        out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/pc_snp_weights_top20.hdf5',
                        snp_filter_frac=1, maf_thres=0.01, num_pcs=20):
    pcs_h5f = h5py.File(pc_file)
    
    print 'Loading Genotype from '
    in_h5f = h5py.File(input_file)
    out_h5f = h5py.File(out_file, 'w')
#     eur_filter = in_h5f['indivs']['continent'][...] == 'EUR'
#     num_indivs = sp.sum(eur_filter)
    indiv_ids = in_h5f['indiv_ids'][...] 
    indiv_filter = None
    assert len(sp.unique(indiv_ids)) == len(indiv_ids)
    num_indivs = len(indiv_ids) 
    

    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        snp_filter = None
        if snp_filter_frac < 1:
            snp_filter = sp.random.random(len(in_h5f[chrom_str]['snps'])) < snp_filter_frac

        g_dict = kgenome.get_genotype_data(in_h5f, chrom, maf_thres, indiv_filter=indiv_filter,
                    snp_filter=snp_filter, randomize_sign=True, snps_signs=None)
        
        norm_snps = g_dict['norm_snps']
        
        evecs = pcs_h5f[chrom_str]['evecs_leave_one_out'][...]
        evals = pcs_h5f[chrom_str]['evals_leave_one_out'][...]
        sort_indices = sp.argsort(evals,)[::-1]
        ordered_evals = evals[sort_indices]
        pcs_var_expl = sp.array(ordered_evals / sp.sum(ordered_evals), dtype='double')
        pcs = evecs[:, sort_indices]
        pcs = pcs[:, :num_pcs]
        ordered_evals = ordered_evals[:num_pcs]
        norm_pcs = pcs - sp.mean(pcs, axis=1)
        pcs_std = sp.std(norm_pcs, axis=1)
        norm_pcs = norm_pcs / pcs_std
        
        cg = out_h5f.create_group(chrom_str)
        cg.create_dataset('snp_pc_weights', data=sp.dot(norm_snps, norm_pcs) / num_indivs) 
        cg.create_dataset('pcs_var_expl', data=pcs_var_expl)
        out_h5f.flush()
    in_h5f.close()
    out_h5f.close()
            
