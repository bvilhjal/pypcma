"""
LDpruning GWAS summary statistics datasets for estimating trait correlations.

Usage: 
python LD_prune_sum_stats.py --coordfile=COORD_FILE  --out=OUT_FILE_PREFIX  --LDradius=LD_RADIUS --LDthres=LD_THRES --KGfile=1KG_FILE  
                         [--nopruning  --wmissing]
 
 - COORD_FILE is a HDF5 with all summary statistics coordinated, i.e. where rsIDs and nucleotides have been coordinated across 
   the different summary statistics. The 1K genomes positions are reported.  This file is in HDF5 format.
 
 - OUT_FILE_PREFIX is the file prefix (w full path) for the Z-score files, one LD pruned and another one without pruning.
 
 - VAL_PLINK_BIM_FILE is a PLINK BIM file which can be used to filter the set of SNPs down to the set of validation SNPs.  To 
   maximize accuracy, it's best to calculate the LDpred weights for the SNPs that are used to calculate the risk scores.
 
 - LD_RADIUS is the maximum distance between SNPs (measured in number of SNPs) for which to calculate LD and prune. 
 
 - LD_THRES is the R2 pruning threshold.
 
 - 1KG_PATH is the full directory path to the 1KG genome impute.legend files
   
 - wmissing: Allow missing data in the coordination step.  Otherwise they are ignored.   

   
 2015 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
          and 
          Hugues Aschard: haschard@hsph.harvard.edu
"""

import scipy as sp
import h5py
import itertools as it
import sys
import time
from sys import argv
import getopt



def calc_ld_table(snps, max_ld_dist=2000, min_r2=0.2, verbose=True, normalize=False):
    """
    Calculate LD between all SNPs using a sliding LD square
    
    This function only retains r^2 values above the given threshold
    """
    # Normalize SNPs (perhaps not necessary, but cheap)
    if normalize:
        snps = snps.T
        snps = (snps - sp.mean(snps, 0)) / sp.std(snps, 0)
        snps = snps.T

    
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    num_snps, num_indivs = snps.shape    
    ld_table = {}
    for i in range(num_snps):
        ld_table[i] = {}

    a = min(max_ld_dist, num_snps)
    num_pairs = (a * (num_snps - 1)) - a * (a + 1) * 0.5
    if verbose:
        print 'Correlation between %d pairs will be tested' % num_pairs
    num_stored = 0
    for i in range(0, num_snps - 1):
        start_i = i + 1
        end_i = min(start_i + max_ld_dist, num_snps)
        ld_vec = sp.dot(snps[i], sp.transpose(snps[start_i:end_i])) / float(num_indivs)
        ld_vec = sp.array(ld_vec).flatten()
        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]
                num_stored += 1
        if verbose:
            if i % 1000 == 0:
                sys.stdout.write('.')
#                 sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(i + 1) / (num_snps - 1)))))
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
        if num_pairs>0:
            print 'Stored %d (%0.4f%%) correlations that made the cut (r^2>%0.3f).' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2)
        else:
            print '-'
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    del snps
    return ld_table



def informed_ld_pruning(scores, ld_table, max_ld=0.5, verbose=False, reverse=False):
    """
    Prunes SNPs in LD, but with smaller scores (p-values or betas)
    
    If using betas, set reversed to True.
    """
    if verbose:
        print 'Performing smart LD pruning'
    t0 = time.time()
    if type(scores) == type([]):
        l = zip(scores, range(len(scores)))
    else:
        l = zip(scores.tolist(), range(len(scores)))
    l.sort(reverse=reverse)
    l = map(list, zip(*l))
    rank_order = l[1]
    indices_to_keep = []
    remaining_indices = set(rank_order)
    for i in rank_order:
        if len(remaining_indices) == 0:
            break
        elif not (i in remaining_indices):
            continue
        else:
            indices_to_keep.append(i)
            for j in ld_table[i]:
                if ld_table[i][j] > max_ld and j in remaining_indices:
                    remaining_indices.remove(j)
    indices_to_keep.sort()
    pruning_vector = sp.zeros(len(scores), dtype='bool8')
    pruning_vector[indices_to_keep] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return pruning_vector             


def ld_pruning(ld_table, max_ld=0.5, verbose=False):
    """
    Prunes SNPs in LD, in random order. 
    """
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


def LD_prune_ss(coord_hdf5_file, out_file, KG_file, r2_thres=0.2, ld_radius=1000, ):
    """
    Uses 1K genome information by default.
    """
    print 'LDpruning SNPs..'
    kg_h5f = h5py.File(KG_file)
    indiv_filter = kg_h5f['indivs']['continent'][...]=='EUR'

    h5f = h5py.File(coord_hdf5_file)
    sums_ids = h5f['sums_ids'][...]
    zs_strings = [s+'_zs' for s in sums_ids]
    ws_strings = [s+'_weights' for s in sums_ids]

    of = open(out_file,'w')
    label_str = 'SID\t'+('\t'.join(zs_strings))+'\t'+('\t'.join(ws_strings))+'\n'
    of.write(label_str)
    for chrom in range(1,23):
        chrom_str = 'chrom_%d' % chrom
        chr_g = h5f[chrom_str]
        sids = chr_g['sids'][...]
        ss_sid_index_map = {}
        for i, sid in enumerate(sids):
            ss_sid_index_map[sid] = i
        zscores = []
        for i, sums_id in enumerate(sums_ids):
            zscores.append((chr_g[sums_id]['zs'][...]).tolist())
        
        zscores = map(list,zip(*zscores))  #transpose
        assert len(sids)==len(zscores), 'WTF?'
    
        kgcg = kg_h5f['chr%d'%chrom]
        kg_sids = kgcg['snp_ids'][...]
        kg_snp_filter = sp.in1d(kg_sids, sids, assume_unique=True)
        
        raw_snps = kgcg['raw_snps'][...]
        raw_snps = raw_snps[kg_snp_filter]
        raw_snps = raw_snps[:,indiv_filter]
        
        ld_table = calc_ld_table(raw_snps, max_ld_dist=ld_radius, min_r2=r2_thres, verbose=True, normalize=True)
        
        filter_vector = ld_pruning(ld_table, max_ld=r2_thres, verbose=False)
        
        filtered_kg_sids = (kg_sids[kg_snp_filter])[filter_vector]
        if 'weights' in chr_g[sums_ids[0]].keys():
            weights = []
            for sums_id in sums_ids:
                weights.append((chr_g[sums_id]['weights'][...]).tolist())        
            weights = map(list,zip(*weights))  #transpose

            assert len(sids)==len(zscores)==len(weights), 'WTF?'
            for sid in filtered_kg_sids:
                zs = zscores[ss_sid_index_map[sid]]
                ws = weights[ss_sid_index_map[sid]]
                out_str = ('%s\t'%(sid))+('\t'.join(map(str,zs)))+'\t'+('\t'.join(map(str,ws)))+'\n'
                of.write(out_str)

        else:
            for sid in filtered_kg_sids:
                zs = zscores[ss_sid_index_map[sid]]
                out_str = ('%s\t'%(sid))+('\t'.join(map(str,zs)))+'\n'
                of.write(out_str)
        
        

def write_out_ss_file(coord_hdf5_file, out_file):
    h5f = h5py.File(coord_hdf5_file)
    sums_ids = h5f['sums_ids'][...]
    zs_strings = [s+'_zs' for s in sums_ids]
    ws_strings = [s+'_weights' for s in sums_ids]
    of = open(out_file,'w')
    label_str = 'Chromosome\tPosition\tSID\tEUR_MAF\t'+('\t'.join(zs_strings))+'\t'+('\t'.join(ws_strings))+'\n'
    of.write(label_str)
    for chrom in range(1,23):
        chrom_str = 'chrom_%d' % chrom
        chr_g = h5f[chrom_str]
        sids = chr_g['sids'][...]
        positions = chr_g['positions'][...]
        chromosomes = sp.repeat(chrom,len(positions))
        eur_mafs = chr_g['eur_mafs'][...]
        zscores = []
        for sums_id in sums_ids:
            zscores.append((chr_g[sums_id]['zs'][...]).tolist())        
        zscores = map(list,zip(*zscores))  #transpose

        if 'weights' in chr_g[sums_ids[0]].keys():
            weights = []
            for sums_id in sums_ids:
                weights.append((chr_g[sums_id]['weights'][...]).tolist())        
            weights = map(list,zip(*weights))  #transpose

            assert len(sids)==len(zscores)==len(weights), 'WTF?'
            for chromosome, position, sid, eur_maf, zs, ws in it.izip(chromosomes,positions,sids,eur_mafs,zscores, weights):
                out_str = ('%s\t%s\t%s\t%0.4f\t'%(chromosome, position, sid, eur_maf))+('\t'.join(map(str,zs)))+'\t'+('\t'.join(map(str,ws)))+'\n'
                of.write(out_str)

        else:
            assert len(sids)==len(zscores), 'WTF?'
            for chromosome, position, sid, eur_maf, zs in it.izip(chromosomes,positions,sids,eur_mafs,zscores):
                out_str = ('%s\t%s\t%s\t%0.4f\t'%(chromosome, position, sid, eur_maf))+('\t'.join(map(str,zs)))+'\n'
                of.write(out_str)
        of.flush() 


def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

                          
    long_options_list = ['coordfile=', 'out=', 'KGfile=', 'LDradius=', 'LDthres=', 'nopruning', 'wmissing']

    p_dict = {'coordfile':None, 'out':None, 'KGfile':'/Users/bjarnivilhjalmsson/data/1Kgenomes/1K_genomes_v3.hdf5', 
              'LDradius':100, 'LDthres':0.2, 'nopruning':False, 'wmissing':False}

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
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt == "--coordfile": p_dict['coordfile'] = arg
            elif opt == "--out": p_dict['out'] = arg
            elif opt == "--KGfile": p_dict['KGfile'] = arg
            elif opt == "--LDradius": p_dict['LDradius'] = int(arg)
            elif opt == "--LDthres": p_dict['LDthres'] = float(arg)
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


if __name__=='__main__':
    p_dict = parse_parameters()
    assert p_dict['coordfile'] is not None, 'Coordinated file is missing.'
    coord_hdf5_file = p_dict['coordfile']

    assert p_dict['out'] is not None, 'Output prefix is missing.'
    zs_file = p_dict['out']+'_zs.txt'
    zs_file_ld_pruned = p_dict['out']+'_zs_ldpruned.txt'

    assert p_dict['KGfile'] is not None, 'The 1K Genomes file is missing.'
    
    
    write_out_ss_file(coord_hdf5_file, zs_file)
    LD_prune_ss(coord_hdf5_file, zs_file_ld_pruned, KG_file=p_dict['KGfile'], r2_thres=p_dict['LDthres'], ld_radius=p_dict['LDradius'])
    
    
