"""
"""

from itertools import izip
import cPickle
import gzip
import os
import sys
import getopt
import kgenome
import slurm
import h5py_util as hu

from scipy import linalg
import h5py

import scipy as sp
import time

__updated__ = '2016-10-18'




def get_ld_scores(snps, ld_radius=500, gm=None, gm_ld_radius=None):
    """
    Calculates LD tables, and the LD scores in one go, with or without a genetic map.
    
    Assumes SNPs are standardized.
    """
    
    m, n = snps.shape
    print m, n
    ld_scores = sp.ones(m)
    ret_dict = {}
    for snp_i, snp in enumerate(snps):
        # Calculate D
        start_i = max(0, snp_i - ld_radius)
        stop_i = min(m, snp_i + ld_radius + 1)
        X = snps[start_i: stop_i]
        D_i = sp.dot(snp, X.T) / n
        r2s = D_i ** 2
        lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
        # lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
        ld_scores[snp_i] = lds_i
    ret_dict['ld_scores'] = ld_scores
    
    return ret_dict


def get_ld_table(norm_snps, ld_radius=1000, min_r2=0.2, verbose=True):
    """
    Calculate LD between all SNPs using a sliding LD square
    
    This function only retains r^2 values above the given threshold
    """
    # Normalize SNPs (perhaps not necessary, but cheap)    
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    m, n = norm_snps.shape

    ld_pairs = []
    ld_pair_r2s = []

    ld_snp_indices = set([])
    
    print m, n
    ld_scores = sp.ones(m)
    num_stored = 0
    num_pairs = 0
    for snp_i, snp in enumerate(norm_snps):
        # Calculate D
        if snp_i == m - 1:
            continue
        start_i = max(0, snp_i - ld_radius)
        stop_i = min(m, snp_i + ld_radius + 1)
        X = norm_snps[start_i: stop_i]
        D_i = sp.dot(X, snp.T) / n
        r2s = D_i ** 2
        assert r2s.max() <= 1.00001 and r2s.min() >= -0.00001
        lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
        ld_scores[snp_i] = lds_i

        shift = min(ld_radius, snp_i) + 1
        r2s = r2s[shift:]
        shift_start_i = shift + start_i
        for snp_j in range(shift_start_i, stop_i):
            num_pairs += 1
            ld_vec_i = snp_j - shift_start_i
            if r2s[ld_vec_i] > min_r2:
                ld_pairs.append(sp.array([snp_i, snp_j], dtype='int8')) 
                ld_pair_r2s.append(r2s[ld_vec_i])
                ld_snp_indices.add(snp_i)
                ld_snp_indices.add(snp_j)
                num_stored += 1
            
        if verbose and snp_i % 10000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
        if verbose:
            print 'Correlation between %d pairs was tested' % num_pairs
            print '%d SNPs out of %d had LD-partners' % (len(ld_snp_indices), len(norm_snps))            
            print 'Stored %d (%0.4f%%) correlations that made the cut (r^2>%0.3f).' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2)
        else:
            print '-'
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    
    print 'Average LD score was %f' % sp.mean(ld_scores)
    ld_snp_indices = list(ld_snp_indices)
    ld_snp_indices.sort()
    return {'ld_scores': sp.array(ld_scores, dtype='float32'), 'ld_pairs':sp.array(ld_pairs, dtype='int8'),
            'ld_pair_r2s':sp.array(ld_pair_r2s, dtype='float32'),
            'ld_snp_indices': sp.array(ld_snp_indices, dtype='int8'),
            'num_snps':len(norm_snps)}


def ld_pruning(ld_dict, max_ld=0.5, verbose=False):
    """
    Prunes SNPs in LD, in random order. 
    """
    import random
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    ld_pairs = ld_dict['ld_pairs']
    ld_snp_indices = ld_dict['ld_snp_indices']
    ld_pair_r2s = ld_dict['ld_pair_r2s']
    num_snps = ld_dict['num_snps']
    random_indices = sp.random.permutation(sp.arange(len(ld_pairs)))
    indices_to_keep = set(ld_snp_indices)
    for pair_index in random_indices:
        (i, j) = ld_pairs[pair_index]
        
        if len(indices_to_keep) == 0:
            break
        
        if (i in indices_to_keep and j in indices_to_keep):
            if ld_pair_r2s[pair_index] > max_ld:
                if random.random() > 0.5:
                    indices_to_keep.remove(i)
                else:
                    indices_to_keep.remove(j)
                    
    filter_vector = sp.zeros(num_snps, dtype='bool')
    filter_vector[list(indices_to_keep)] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return filter_vector


def calculate_ld_tables(input_genotype_file, chrom_i, local_ld_hdf5_file, ld_radius,
                        min_r2=0.2, maf_thres=0.01, indiv_filter=None,
                        snp_filter=None, return_void=True, verbose=True):
    """
    Calculate the LD tables for the given radius, and store in the given file.
    """
    if not os.path.isfile(local_ld_hdf5_file):
        h5f = h5py.File(input_genotype_file)
        
        print 'Calculating LD information for chromosome %d w. radius %d' % (chrom_i, ld_radius)

        g_dict = kgenome.get_genotype_data(h5f, chrom_i, maf_thres, indiv_filter=indiv_filter,
                        snp_filter=snp_filter, randomize_sign=False, snps_signs=None)
          
        ld_dict = get_ld_table(g_dict['norm_snps'], ld_radius=ld_radius, min_r2=min_r2, verbose=verbose)
        ld_dict['avg_ld_score'] = sp.mean(ld_dict['ld_scores'])
                    
        print 'Done calculating the LD table and LD score, writing to file:', local_ld_hdf5_file
        oh5f = h5py.File(local_ld_hdf5_file, 'w')
        hu.dict_to_hdf5(ld_dict, oh5f)
        oh5f.close()
        print 'LD information is now stored.'
    else:
        print 'Loading LD information from file: %s' % local_ld_hdf5_file
        oh5f = h5py.File(local_ld_hdf5_file, 'r')
        ld_dict = hu.hdf5_to_dict(oh5f)
        oh5f.close()
    
    if not return_void:
        return ld_dict
        

def get_pruning_filter_dict(input_file='Data/1Kgenomes/1K_genomes_v3.hdf5', max_r2=0.2, maf_thres=0.01):
    """
    LD prune SNPs...
    """
    pass

def _get_sub_run_id_(run_id, chrom_i, ld_radius):
    return 'rid%s_chr%d_ldr%d' % (run_id, chrom_i, ld_radius)

def submit_ld_job(run_id, genotype_file, chrom_i, ld_radius, ld_file_prefix,
                  min_maf=0.01, walltime='12:00:00', queue_id='normal', max_memory='16g', num_cores=4,
                  job_dir='.', script_dir='.', email='bjarni.vilhjalmsson@gmail.com'):
    """
    Submit a command to the cluster
    """
    sub_run_id = _get_sub_run_id_(run_id, chrom_i, ld_radius)
    out_file = '%s.out' % sub_run_id
    err_file = '%s.err' % sub_run_id

    dir_path = os.path.dirname(os.path.realpath(__file__))

    python_command = 'python %s/ld.py' % dir_path
    python_command = '%s --local_ld_file_prefix=%s' % (python_command, ld_file_prefix)
    python_command = '%s --chrom=%d' % (python_command, chrom_i)
    python_command = '%s --ld_radius=%d' % (python_command, ld_radius)
    python_command = '%s --genotype_file=%s' % (python_command, genotype_file)
    python_command = '%s --min_maf=%0.3f' % (python_command, min_maf)

    slurm.submit_job(run_id, python_command, out_file=out_file, err_file=err_file,
               walltime=walltime, queue_id=queue_id, num_cores=num_cores, max_memory=max_memory,
               script_dir=script_dir, job_dir=job_dir, email=email, delete_script_after_submit=False)
    

def submit_all_ld_jobs(genotype_file, ld_file_prefix, run_id, min_maf=0.01, ld_radius=100,
                       job_dir='.', script_dir='.', walltime='6:00:00', num_cores=4,
                       max_memory='16g', email='bjarni.vilhjalmsson@gmail.com'):
    
    for chrom_i in range(1, 23):
        submit_ld_job(run_id, genotype_file, chrom_i, ld_radius, ld_file_prefix,
                      min_maf=min_maf, walltime=walltime, queue_id='normal',
                      max_memory=max_memory, num_cores=num_cores,
                      job_dir=job_dir, script_dir=script_dir, email=email)


def _parse_parameters_():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

    long_options_list = ['ld_radius=', 'local_ld_file_prefix=', 'run_id=',
                         'genotype_file=', 'chrom=', 'min_maf=', 'sub_run_id=', 'h']

    p_dict = {'ld_radius':1000, 'local_ld_file_prefix':'./ld_file', 'run_id':'default_id',
              'genotype_file': None, 'chrom':None, 'min_maf': 0.01, 'sub_run_id':None, }

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
    
        except:
            print "Some problems with parameters.  Please read the usage documentation carefully."
            print "Use the -h option for usage information."
#             traceback.print_exc()
#             print __doc__
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == "-h" or opt == "--h":
                print __doc__
                sys.exit(0)
            elif opt == "--ld_radius": p_dict['ld_radius'] = int(arg)
            elif opt == "--local_ld_file_prefix": p_dict['local_ld_file_prefix'] = arg
            elif opt == "--sub_run_id": p_dict['sub_run_id'] = arg
            elif opt == "--chrom": p_dict['chrom'] = int(arg)
            elif opt == "--genotype_file": p_dict['genotype_file'] = arg
            elif opt == "--min_maf": p_dict['min_maf'] = float(arg)
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


def main():
    p_dict = _parse_parameters_()
    if p_dict['sub_run_id'] is None:
        # Calculate LD...
        local_ld_dict_file = '%s_chrom%d_ldradius%d.hdf5' % (p_dict['local_ld_file_prefix'], p_dict['chrom'], p_dict['ld_radius'])
        calculate_ld_tables(p_dict['genotype_file'], p_dict['chrom'], local_ld_dict_file,
                            p_dict['ld_radius'], min_r2=0.2, maf_thres=p_dict['min_maf'])
    else:
        print 'Nothing happened.  What did you expect?'
    
if __name__ == '__main__':
    main()
    
