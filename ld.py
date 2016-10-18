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

from scipy import linalg
import h5py

import scipy as sp
import time

__updated__ = '2016-10-18'




def get_ld_tables(snps, ld_radius=500, gm=None, gm_ld_radius=None):
    """
    Calculates LD tables, and the LD scores in one go, with or without a genetic map.
    
    Assumes SNPs are standardized.
    """
    
    m, n = snps.shape
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
            lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
            ld_scores[snp_i] = lds_i
        
        avg_window_size = sp.mean(window_sizes)
        print 'Average # of SNPs in LD window was %0.2f' % avg_window_size
        ret_dict['ld_boundaries'] = ld_boundaries
    ret_dict['ld_scores'] = ld_scores
    
    return ret_dict

def ld_pruning(ld_table, max_ld=0.5, verbose=False):
    """
    Prunes SNPs in LD, in random order. 
    """
    import random
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    indices_to_keep = []
    num_snps = len(ld_table)
    indices = sp.random.permutation(num_snps)
    remaining_indices = set(indices)
    for i in indices:
        if len(remaining_indices) == 0:
            break
        elif not (i in remaining_indices):
            continue
        else:
            indices_to_keep.append(i)
            for j in ld_table[i]:
                if ld_table[i][j] > max_ld and j in remaining_indices:
                    remaining_indices.remove(j)
    filter_vector = sp.zeros(num_snps, dtype='bool')
    filter_vector[indices_to_keep] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return filter_vector


def calculate_ld_tables(input_genotype_file, chrom_i, ld_file_prefix, ld_radius,
                        maf_thres, gm_ld_radius=None, indiv_filter=None,
                        snp_filter=None, return_void=True):
    """
    Calculate the LD tables for the given radius, and store in the given file.
    """
    local_ld_dict_file = '%s_lrd%d_chr%d.pickled.gz' % (ld_file_prefix, ld_radius, chrom_i)
    if not os.path.isfile(local_ld_dict_file):
        h5f = h5py.File(input_genotype_file)
        
        chrom_ld_scores_dict = {}
        chrom_ld_dict = {}
        chrom_ref_ld_mats = {}
        if gm_ld_radius is not None:
            chrom_ld_boundaries = {}
        ld_score_sum = 0
        print 'Calculating LD information for chromosome %d w. radius %d' % (chrom_i, ld_radius)

        g_dict = kgenome.get_genotype_data(h5f, chrom_i, maf_thres, indiv_filter=indiv_filter,
                        snp_filter=snp_filter, randomize_sign=False, snps_signs=None)
          
        
        ld_dict = get_ld_tables(g_dict['norm_snps'], ld_radius=ld_radius)
        ld_dict['avg_ld_score'] = sp.mean(ld_dict['ld_scores'])
        ld_dict['num_snps'] = len(g_dict['norm_snps'])
                    
        print 'Done calculating the LD table and LD score, writing to file:', local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'wb')
        cPickle.dump(ld_dict, f, protocol=2)
        f.close()
        print 'LD information is now pickled.'
    else:
        print 'Loading LD information from file: %s' % local_ld_dict_file
        f = gzip.open(local_ld_dict_file, 'r')
        ld_dict = cPickle.load(f)
        f.close()
    
    if not return_void:
        return ld_dict
        

def get_pruning_filter_dict(input_file='Data/1Kgenomes/1K_genomes_v3.hdf5', max_r2=0.2, maf_thres=0.01):
    """
    LD prune SNPs...
    """
    pass

def _get_sub_run_id_(run_id, chrom_i, ld_radius):
    return 'rid%s_chr%i_ldr%i' % (run_id, chrom_i, ld_radius)

def submit_ld_job(run_id, min_maf, genotype_file, chrom_i, ld_radius, ld_file_prefix,
                  walltime='12:00:00', queue_id='normal', max_memory=8000, num_cores=4,
                  job_dir='.', script_dir='.', email='bjarni.vilhjalmsson@gmail.com'):
    """
    Submit a command to the cluster
    """
    sub_run_id = _get_sub_run_id_(run_id, chrom_i, ld_radius)
    out_file = '%s.out' % sub_run_id
    err_file = '%s.err' % sub_run_id

    python_command = 'python ./ld.py'
    python_command = '%s --ld_file_prefix=%s' % (python_command, ld_file_prefix)
    python_command = '%s --chrom=%d' % (python_command, chrom_i)
    python_command = '%s --ld_radius=%d' % (python_command, ld_radius)
    python_command = '%s --genotype_file=%s' % (python_command, genotype_file)
    python_command = '%s --min_maf=%0.3f' % (python_command, min_maf)
    python_command = '%s --ld_file_prefix=%s' % (python_command, ld_file_prefix)

    slurm.submit_job(run_id, python_command, out_file=out_file, err_file=err_file,
               walltime=walltime, queue_id=queue_id, num_cores=num_cores, max_memory=max_memory,
               script_dir=script_dir, job_dir=job_dir, email=email, delete_script_after_submit=False)
    

def submit_all_ld_jobs(genotype_file, ld_file_prefix, run_id, min_maf=0.01, ld_radius=100,
                       job_dir='.', script_dir='.', walltime='4:00:00', num_cores=4,
                       max_memory=8000, email='bjarni.vilhjalmsson@gmail.com'):
    
    for chrom_i in range(1, 23):
        submit_ld_job(run_id, min_maf, genotype_file, chrom_i, ld_radius, ld_file_prefix,
                      walltime=walltime, queue_id='normal', max_memory=max_memory, num_cores=num_cores,
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
              'genotype_file': None, 'chrom':None, 'min_maf': 0.01, 'sub_run_id':None}

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
            elif opt == "--run_id": p_dict['run_id'] = arg
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
    if p_dict['sub_run_id'] is not None:
        # Calculate LD...
        local_ld_dict_file = '%s_ldradius%d.pickled.gz' % (p_dict['local_ld_file_prefix'], p_dict['ld_radius'])
        calculate_ld_tables(p_dict['genotype_file'], p_dict['chrom'], local_ld_dict_file,
                            p_dict['ld_radius'], p_dict['min_maf'])
    else:
        print 'Nothing happened.  What did you expect?'
    
if __name__ == '__main__':
    main()
    
