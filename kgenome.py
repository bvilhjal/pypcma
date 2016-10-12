"""
Methods for analysing 1000 genomes data.
"""

from itertools import izip
import cPickle
import gzip
import os

from scipy import linalg
import h5py

import scipy as sp


__updated__ = '2016-10-12'

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}
Kg_nt_decoder = {1:'A', 2:'T', 3:'C', 4:'G', }
ok_nts = sp.array(['A', 'T', 'G', 'C'])
 
def dict_to_hdf5(input_dict, hdf5_group):
    """
    Recursively constructs HDF5 file groups from dictionaries.  
    
    Assumes that the data can be undestood by h5py create_dataset function.
    """
    for key in input_dict:
        if isinstance(input_dict[key], dict):
            new_hdf5_group = hdf5_group.create_group(key)
            dict_to_hdf5(input_dict[key], new_hdf5_group)
        else:
            hdf5_group.create_dataset(key, data=input_dict[key])
            
def hdf5_to_dict(hdf5_group):
    output_dict = {}
    for key in hdf5_group.keys():
        if isinstance(hdf5_group[key], h5py._hl.group.Group):
            output_dict[key] = hdf5_to_dict(hdf5_group[key])
        else:
            output_dict[key] = hdf5_group[key][...]
    return output_dict
    
    
def gen_unrelated_eur_1k_data(input_file='/home/bjarni/TheHonestGene/faststorage/1Kgenomes/phase3/1k_genomes_hg.hdf5' ,
                              out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1K_genomes_phase3_EUR_unrelated.hdf5',
                              maf_thres=0.01, max_relatedness=0.05, K_thinning_frac=0.1, debug=False):
    h5f = h5py.File(input_file)
    num_indivs = len(h5f['indivs']['continent'])
    eur_filter = h5f['indivs']['continent'][...] == 'EUR'
    num_eur_indivs = sp.sum(eur_filter)
    print 'Number of European individuals: %d', num_eur_indivs
    K = sp.zeros((num_eur_indivs, num_eur_indivs), dtype='single')
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

        print 'Loading NTs'
        ref_nts = h5f[chrom_str]['variants']['REF'][...]
        alt_nts = h5f[chrom_str]['variants']['ALT'][...]
        
        print 'Filtering multi-allelic SNPs'
        multi_allelic_filter = sp.negative(h5f[chrom_str]['variants']['MULTI_ALLELIC'][...])
        snps = snps[multi_allelic_filter]
        ref_nts = ref_nts[multi_allelic_filter]
        alt_nts = alt_nts[multi_allelic_filter]

        print 'Filter individuals'
        snps = snps[:, indiv_filter]
        
        print 'Filter SNPs with missing NT information'
        nt_filter = sp.in1d(ref_nts, ok_nts)
        nt_filter = nt_filter * sp.in1d(alt_nts, ok_nts)
        if sp.sum(nt_filter) < len(nt_filter):
            snps = snps[nt_filter]
            ref_nts = ref_nts[nt_filter]
            alt_nts = alt_nts[nt_filter]
        
        print 'filter monomorphic SNPs'
        snp_stds = sp.std(snps, 1)
        mono_morph_filter = snp_stds > 0
        snps = snps[mono_morph_filter]
        
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
        
        snp_ids = h5f[chrom_str]['snp_ids'][...]
        snp_ids = snp_ids[mono_morph_filter]
        snp_ids = snp_ids[nt_filter]
        
        positions = h5f[chrom_str]['positions'][...]
        positions = positions[mono_morph_filter]
        positions = positions[nt_filter]

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
        print 'Done'
        
#         centimorgans = h5f[chrom_str]['centimorgans'][...]
#         cg.create_dataset('centimorgans',data=centimorgans)
#         
#         centimorgan_rates = h5f[chrom_str]['centimorgan_rates'][...]
#         cg.create_dataset('centimorgan_rates',data=centimorgan_rates)
        
    oh5f.close()
    h5f.close()
    
    
    
def get_kinship_pca_dict(input_genotype_file, kinship_pca_file, maf_thres, debug_filter):
    if os.path.isfile(kinship_pca_file):
        print ':Loading Kinship and PCA information from %s' % kinship_pca_file
        k_h5f = h5py.File(kinship_pca_file)
        kinship_pca_dict = hdf5_to_dict(k_h5f)
    else:
        kinship_pca_dict = calc_kinship(input_file=input_genotype_file , out_file=kinship_pca_file,
                                                maf_thres=maf_thres, figure_dir=None, debug_filter=debug_filter)
    return kinship_pca_dict

    
def calc_kinship(input_file='Data/1Kgenomes/1K_genomes_v3.hdf5' , out_file='Data/1Kgenomes/kinship.hdf5',
                  maf_thres=0.1, figure_dir='', figure_fn='', debug_filter=1, indiv_filter_frac=0.5):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    in_h5f = h5py.File(input_file)
#     eur_filter = in_h5f['indivs']['continent'][...] == 'EUR'
#     num_indivs = sp.sum(eur_filter)
    indiv_ids = in_h5f['indiv_ids'][...] 
    indiv_filter = sp.random.random(len(indiv_ids)) < indiv_filter_frac
    indiv_filter[0] = False
    indiv_ids = indiv_ids[indiv_filter]
    num_indivs = len(indiv_ids)  # An ugly bug hack!!
    chromosome_dict = {}
    
    std_thres = sp.sqrt(2.0 * (1 - maf_thres) * (maf_thres))

    K_all_snps = sp.zeros((num_indivs, num_indivs), dtype='double')
    num_all_snps = 0
    print 'Calculating kinship'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        
        
        print 'Loading SNPs'
        snps = in_h5f[chrom_str]['snps'][...]
        snps = snps[:, indiv_filter]
        # filter non-europeans.
#         print 'Filtering non-European individuals'
#         snps = snps[:, eur_filter]
#         snp_stds = in_h5f[chrom_str]['snp_stds'][...]
#         snp_means = in_h5f[chrom_str]['snp_means'][...]
        nts = in_h5f[chrom_str]['nts'][...]
#         print 'snps.shape: %s, snp_stds.shape: %s, snp_means.shape: %s' % (str(snps.shape), str(snp_stds.shape), str(snp_means.shape))
        
        if debug_filter < 1:
            debug_snp_filter = sp.random.random(len(snps)) < debug_filter
        snps = snps[debug_snp_filter]
        snp_means = sp.mean(snps, 1)
        snp_means.shape = (len(snp_means), 1)
#         snp_freqs = snp_means / 2
        snp_stds = sp.std(snps, 1)
        snp_stds.shape = (len(snp_stds), 1)
        
#         snp_stds = snp_stds[debug_snp_filter]
#         snp_means = snp_means[debug_snp_filter]
        nts = nts[debug_snp_filter]

        
        print 'Filtering SNPs with MAF <', maf_thres
        maf_filter = snp_stds.flatten() > std_thres
        snps = snps[maf_filter]
        snp_stds = snp_stds[maf_filter]
        snp_means = snp_means[maf_filter]
        nts = nts[maf_filter]
        
#         nt_filter = sp.in1d(nts, ok_nts)
#         if not sp.all(nt_filter):
#             print 'Removing SNPs with missing NT information'
#             snps = snps[nt_filter]
#             snp_stds = snp_stds[nt_filter]
#             snp_means = snp_means[nt_filter]
        
        print '%d SNPs remaining' % len(snps)
        
        print 'Normalizing SNPs'
        
        norm_snps = (snps - snp_means) / snp_stds
        
        print 'Calculating chromosome kinship'
        K_unscaled = sp.dot(norm_snps.T, norm_snps)
        
        assert sp.isclose(sp.sum(sp.diag(K_unscaled)) / (len(norm_snps) * num_indivs), 1.0), '..bug' 
        
        print 'Storing and updating kinships'
        chromosome_dict[chrom_str] = {'K_unscaled':K_unscaled, 'num_snps':len(norm_snps)}
        
        K_all_snps += K_unscaled
        num_all_snps += len(norm_snps)
    K_all_snps = K_all_snps / float(num_all_snps)
    print 'K_all_snps.shape: %s' % str(K_all_snps.shape)
    print 'Full kinship calculation done using %d SNPs\n' % num_all_snps

    in_h5f.close()

    print 'Now generating the leave-one-chromosome-out kinships'
    for chrom in range(1, 23):
        print 'Working on Chromosome %d' % chrom
        chrom_str = 'chr%d' % chrom
        K_leave_one_out = sp.zeros((num_indivs, num_indivs), dtype='double')
        num_snps_used = 0 
        for chrom2 in range(1, 23):
            if not chrom2 == chrom: 
                chrom2_str = 'chr%d' % chrom
                K_leave_one_out += chromosome_dict[chrom2_str]['K_unscaled']
                num_snps_used += chromosome_dict[chrom2_str]['num_snps']
                assert sp.isclose(sp.sum(sp.diag(K_leave_one_out)) / (num_snps_used * num_indivs), 1.0), '..bug' 
        K_leave_one_out = K_leave_one_out / num_snps_used
        chromosome_dict[chrom_str]['K_leave_one_out'] = K_leave_one_out
        evals, evecs = linalg.eig(K_leave_one_out)  # PCA via eigen decomp
        assert sp.all(evals > 0), '...bug'

        evals, evecs = linalg.eig(K_leave_one_out)  # PCA via eigen decomp
        assert sp.all(evals > 0), '...bug'
        sort_indices = sp.argsort(evals,)
        ordered_evals = evals[sort_indices]
        print ordered_evals[-10:] / sp.sum(ordered_evals)
        ordered_evecs = evecs[:, sort_indices]
        chromosome_dict[chrom_str]['evecs_leave_one_out'] = ordered_evecs
        chromosome_dict[chrom_str]['evals_leave_one_out'] = ordered_evals
        
    assert sp.sum((chromosome_dict['chr1']['K_leave_one_out'] - chromosome_dict['chr2']['K_leave_one_out']) ** 2) != 0 , 'Kinships are probably too similar.'
        
    print 'Calculating PCAs'
    evals, evecs = linalg.eig(K_all_snps)  # PCA via eigen decomp
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
        K_leave_one_out = chromosome_dict[chrom_str]['K_leave_one_out']

    
    out_h5f = h5py.File(out_file)
    
    dict_to_hdf5(chromosome_dict, out_h5f)
    out_h5f.close()
    
    return chromosome_dict

    
   
def gen_1k_test_genotypes(kg_file='Data/1Kgenomes/1K_genomes_v3_EUR_unrelated2.hdf5',
                          nt_map_file='tmp/nt_map.pickled', out_prefix='tmp/1k_ind'):
    """
    Generates 1K genotypes in the internal genotype format for validation and other purposes.
    """
    print 'Loading NT map from file: %s' % nt_map_file
    f = open(nt_map_file, 'r')
    snp_map_dict = cPickle.load(f)
    f.close()
    
    h5f = h5py.File(kg_file)
    kg_indivs = h5f['indiv_ids'][...]
    chromosomes = range(1, 23) 
    
    
    # Figure out 1K SNP filter
    chrom_filter_dict = {}
    
    for chrom in chromosomes:
        kg_chrom_str = 'chr%d' % chrom
        cg = h5f[kg_chrom_str]
        sids = cg['snp_ids'][...]
        
        # Get the nucleotides coding map (from 1K genomes project).
        chrom_dict = snp_map_dict[kg_chrom_str]
        ok_sids = chrom_dict['sids']
        snps_filter = sp.in1d(sids, ok_sids)
        chrom_filter_dict[chrom] = snps_filter
        
        reorder_snps = False
        filtered_sids = sids[snps_filter]
        assert len(filtered_sids) == len(ok_sids), '.... bug'
        if not sp.all(filtered_sids == ok_sids):
            sids_indices_dict = {}
            for i, sid in enumerate(sids):
                sids_indices_dict[sid] = i
            snp_order = []
            for sid in ok_sids:
                snp_order.append(sids_indices_dict[sid])
            filtered_sids = filtered_sids[snp_order]
            assert sp.all(filtered_sids == ok_sids), '... bug'
            reorder_snps = True
    
    
    for ind_i in range(10):  # len(kg_indivs)):
        print 'Generating genotype for individual: %d ' % ind_i
        
        # prepare output file
        out_h5fn = out_prefix + '_%d.hdf5' % ind_i
        out_h5f = h5py.File(out_h5fn)
        
        for chrom in chromosomes:
#             print '\nWorking on chromosome %d'%chrom
            kg_chrom_str = 'chr%d' % chrom
            chrom_str = 'Chr%d' % chrom
            cg = h5f[kg_chrom_str]
            
            snps_filter = chrom_filter_dict[chrom]

            sids = (cg['snp_ids'][...])[snps_filter]
            snps = (cg['snps'][...])[:, ind_i]
            snps = snps[snps_filter]
            positions = (cg['positions'][...])[snps_filter]
            nts = (cg['nts'][...])[snps_filter]
            
            if reorder_snps:
                snps = snps[snp_order]
                positions = positions[snp_order]
                nts = nts[snp_order]
                sids = sids[snp_order]
            
            assert len(snps) == len(sids) == len(positions) == len(nts), '..bug'
            # Store information
            ocg = out_h5f.create_group  (chrom_str)
            ocg.create_dataset('snps', data=snps)
            ocg.create_dataset('sids', data=sids)
            ocg.create_dataset('positions', data=positions)
            ocg.create_dataset('nts', data=nts)

        out_h5f.close()
    h5f.close()
        # return genome_dict        



# Coding key
# def prepare_nt_coding_key(K_genomes_snps_map, indiv_genot_file, nt_map_file):
#     """
#     Determines the nucleotide coding for the genotype using the 1K genome  the 10KUK
#     """
#     gf = h5py.File(indiv_genot_file,'r')
#     kgf = h5py.File(K_genomes_snps_map,'r')
#     chromosomes = range(1,23) 
#     snp_map_dict = {}
#     num_snps = 0
#     for chrom in chromosomes:
#         print 'Working on chromosome %d'%chrom
#         kg_chrom_str = 'chr%d'%chrom
#         chrom_str = 'Chr%d'%chrom
#         
#         #Get SNPs from genotype
#         cg = gf[chrom_str]
#         sids = cg['ids'][...]
#         snps = cg['snps'][...]
#         sid_dict = dict(zip(sids, snps))
#         
#         #Get SNP IDs from 1K genomes
#         kcg = kgf[kg_chrom_str]
#         kg_sids = kcg['snp_ids'][...]
#         
#         #Determine overlap between SNPs..
#         kg_filter = sp.in1d(kg_sids,sids)
#         kg_sids = kg_sids[kg_filter]
#         kg_nts = (kcg['nts'][...])[kg_filter]
#         kg_positions = (kcg['positions'][...])[kg_filter]
#         
#         #Check that nt are ok in genotype data, otherwise filter.
#         sid_nt_map = {}
#         positions = []
#         ok_sids = []
#         nts = []
#         snp_i = 0
#         for sid, kg_nt, kg_pos in izip(kg_sids, kg_nts, kg_positions):
#             snp = sid_dict[sid]
#             if tuple(kg_nt) not in ambig_nts:
#                 # All possible (and allowed) nucleotides strings 
#                 ntm = {}
#                 ntm['--']=-9
#                 ntm['-'+kg_nt[0]]=-9
#                 ntm['-'+kg_nt[1]]=-9
#                 ntm[kg_nt[0]+'-']=-9
#                 ntm[kg_nt[1]+'-']=-9
#                 ntm[kg_nt[0]+kg_nt[0]]=0
#                 ntm[kg_nt[1]+kg_nt[0]]=1
#                 ntm[kg_nt[0]+kg_nt[1]]=1
#                 ntm[kg_nt[1]+kg_nt[1]]=2
#                 sid_nt_map[sid]={'ntm':ntm, 'snp_i':snp_i}
#                 positions.append(kg_pos)
#                 nts.append(kg_nt)
#                 ok_sids.append(sid)
#                 snp_i += 1
#         
#         num_snps += len(sid_nt_map)
#         snp_map_dict[kg_chrom_str]={'sid_nt_map':sid_nt_map, 'positions':positions, 'nts':nts, 'sids':ok_sids}
#     
#     print 'Found %d SNPs'%num_snps
#     print 'Writing to file'
#     f = open(nt_map_file, 'wb')
#     cPickle.dump(snp_map_dict, f, protocol=2)
#     f.close()
#     return snp_map_dict
        
# For debugging purposes
if __name__ == '__main__':        
    data_dir = '/home/bjarni/PCMA/faststorage/1Kgenomes/'
    kinship_pca_dict = calc_kinship('/home/bjarni/HeritPartition/faststorage/1Kgenomes_bjarni/1K_genomes_v3_EUR_unrelated.hdf5',
                         out_file='/home/bjarni/PCMA/faststorage/1_DATA/1k_genomes/1kgenomes_kinship_pca.hdf5',
                         figure_dir='/home/bjarni/tmp', figure_fn='test.pdf',
                         maf_thres=0.01, debug_filter=0.05, indiv_filter_frac=0.5)


