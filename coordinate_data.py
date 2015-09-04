"""
Coordinate summary statistics datasets for performing multiple-trait analysis.

Usage: 
python coordinate_data.py --ssfiles=SUM_STATS_FILE1,SUM_STATS_FILE2,... --combfile=COMB_SS_FILE   
                          --sslabels=LABEL1,LABEL2,... --1KGpath=1KG_PATH  [--coordfile=COORD_FILE --ssf_format=SSF_FORMAT --wmissing]

 - SUM_STATS_FILE1,SUM_STATS_FILE2,...: This is a comma-separated (without space) list of files which should contain a (full path) filename 
   with the GWAS summary statistics.  The STANDARD format is FORMAT1 (see below)

 - COMB_SS_FILE is a HDF5 file that contains all of the GWAS summary statistics, but they have not been coordinated at this stage.
 
 - COORD_FILE is a HDF5 with all summary statistics coordinated, i.e. where rsIDs and nucleotides have been coordinated across 
   the different summary statistics. The 1K genomes positions are reported.  This file is in HDF5 format.
  
 - LABEL1,LABEL2,...: This is a comma-separated (without space) list of labels for the GWAS summary statistics. 
 
 - 1KG_PATH is the full directory path to the 1KG genome impute.legend files
   
 - wmissing: Allow missing data in the coordination step.  Otherwise they are ignored.   
   
 - SSF_FORMAT refers to the summary statistics format. The standard format is the "BASIC" format, which contains of the 
   basic required information, is as follows:
    
    hg19chrc    snpid    a1    a2    bp    or    p       
    chr1    rs4951859    C    G    729679    0.97853    0.2083  
    chr1    rs142557973    T    C    731718    1.01949    0.3298  
    ..
    ..
   
 2015 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
          and 
          Hugues Aschard: haschard@hsph.harvard.edu
"""
import scipy as sp
from scipy import stats
import h5py
import itertools as it
import gzip
import sys
import time
from sys import argv
import random
import getopt


ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
ok_nts = set([('A', 'G'), ('G', 'A'), ('A', 'C'), ('C', 'A'),('G', 'T'), ('T', 'G'),('C', 'T'), ('T', 'C')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

valid_nts = set(['A','T','C','G'])

lc_2_cap_map = {'a':'A', 'c':'C', 'g':'G', 't':'T'}


def get_sid_pos_map(sids, KGenomes_prefix):
    sids = set(sids)
    sid_map = {}
    for chrom_i in range(1,23):
        fn = '%sALL_1000G_phase1integrated_v3_chr%d_impute.legend.gz'%(KGenomes_prefix,chrom_i)
        with gzip.open(fn) as f:
            f.next()
            for line in f:
                l = line.split()
                sid = l[0]
                if sid in sids:
                    sid_map[l[0]]={'pos':int(l[1]), 'chrom':chrom_i}
    return sid_map



# def parse_sum_stats1(filename,
#                     comb_hdf5_file,
#                     ss_id,
#                     N,
#                     bimfile=None):
#     """
#     Input format:
# 
#     hg19chrc    snpid    a1    a2    bp    or    p       
#     chr1    rs4951859    C    G    729679    0.97853    0.2083  
#     chr1    rs142557973    T    C    731718    1.01949    0.3298  
# 
# 
#     ...
#     
#     """
#     h5f = h5py.File(comb_hdf5_file)
#     if bimfile!=None:
#         print 'Parsing SNP list'
#         valid_sids = set()
#         print 'Parsing bim file: %s'%bimfile
#         with open(bimfile) as f:
#             for line in f:
#                 l = line.split()
#                 valid_sids.add(l[1])
#         print len(valid_sids)
#     chrom_dict = {}
# 
# 
#     print 'Parsing the file: %s' % filename
#     with open(filename) as f:
#         print f.next()
#         for line in f:
#             l = (line.strip()).split()
#             chrom_str = l[0]
#             chrom = chrom_str[3:]
#             if chrom.isdigit():
#                 chrom = int(chrom)
#                 pos = int(l[4])
#                 sid = l[1]
#                 if sid in valid_sids:
#                     if not chrom in chrom_dict.keys():
#                         chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[],
#                                              'betas':[], 'nts': [], 'sids': [], 
#                                              'positions': [], 'Ns':[]}
#                     chrom_dict[chrom]['sids'].append(sid)
#                     chrom_dict[chrom]['positions'].append(pos)
#                     pval = float(l[6])
#                     chrom_dict[chrom]['ps'].append(pval)
#                     if random.random()>0.5:
#                         nt = [l[2], l[3]]
#                         raw_beta = sp.log(float(l[5]))
#                     else:
#                         nt = [l[3], l[2]]
#                         raw_beta = sp.log(float(l[5]))
#                     chrom_dict[chrom]['nts'].append(nt)                
#                     chrom_dict[chrom]['log_odds'].append(raw_beta)
#                     beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
#                     chrom_dict[chrom]['betas'].append(beta/sp.sqrt(N))
#                     chrom_dict[chrom]['Ns'].append(N)
#          
#             
# 
#     print 'SS file loaded, now sorting and storing in HDF5 file.'
#     assert not ss_id in h5f.keys(), 'Summary stats with this name are already in the HDF5 file?'
#     ssg = h5f.create_group(ss_id)
#     num_snps = 0
#     for chrom in chrom_dict.keys():
#         print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
#         sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
#                  chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['ps'], chrom_dict[chrom]['Ns'])
#         sl.sort()
#         ps = []
#         betas = []
#         nts = []
#         sids = []
#         positions = []
#         log_odds = []
#         ns = []
#         prev_pos = -1
#         for pos, sid, nt, beta, lo, p, n in sl:
#             if pos == prev_pos:
#                 print 'duplicated position %d' % pos
#                 continue
#             else:
#                 prev_pos = pos
#             ps.append(p)
#             betas.append(beta)
#             nts.append(nt)
#             sids.append(sid)
#             positions.append(pos)
#             log_odds.append(lo)
#             ns.append(n)
#         g = ssg.create_group('chrom_%d' % chrom)
#         g.create_dataset('ps', data=sp.array(ps))
#         g.create_dataset('betas', data=betas)
#         g.create_dataset('log_odds', data=log_odds)
#         num_snps +=len(log_odds)
#         g.create_dataset('nts', data=nts)
#         g.create_dataset('sids', data=sids)
#         g.create_dataset('positions', data=positions)
#         g.create_dataset('Ns', data=ns)
#         h5f.flush()
#     print 'In all, %d SNPs parsed from summary statistics file.'%num_snps
#             


def parse_sum_stats_basic(filename,
                    comb_hdf5_file,
                    ss_id,
                    KGpath,
                    bimfile =None):
    """
    Input format:

    MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       se      p       N
    rs3094315       a       g       0.864   0.0068  0.012   0.57    54178
    rs2905035       a       g       0.106   -0.0075 0.014   0.59    51152
    rs2980319       a       t       0.106   -0.0074 0.014   0.6     51152
    rs4040617       a       g       0.894   0.0082  0.014   0.56    51052

    ...
    
    """
    h5f = h5py.File(comb_hdf5_file)
    if bimfile!=None:
        print 'Parsing SNP list'
        valid_sids = set()
        print 'Parsing bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print len(valid_sids)
    chrom_dict = {}

    print 'Retrieving 1K genomes positions..'
    sids = []
    with open(filename) as f:
        for line in f:
            l = (line.strip()).split()
            sid = l[0]
            sids.append(sid)
    sid_map = get_sid_pos_map(sids,KGpath)


    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        while 1:
            line =  f.next()
            header = line.split()[0]
            if header=='MarkerName' or header=='SNP':
                break
        
        line_i = 0
        for line in f:
            line_i +=1
            l = (line.strip()).split()
            sid = l[0]
            d = sid_map.get(sid,None)
            if d !=None:
                pos = d['pos']
                chrom = d['chrom']
                if not chrom in chrom_dict.keys():
                    chrom_dict[chrom] = {'ps':[], 'zs':[], 'betas':[], 'nts': [], 'sids': [], 
                                         'positions': [], 'Ns':[]}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                pval = float(l[6])
                chrom_dict[chrom]['ps'].append(pval)
                if random.random()>0.5:
                    #nt = [lc_2_cap_map[l[1]], lc_2_cap_map[l[2]]]
                    nt = [l[1], l[2]]
                    raw_beta = float(l[4])
                else:
                    #nt = [lc_2_cap_map[l[2]], lc_2_cap_map[l[1]]]
                    nt = [l[2], l[1]]
                    raw_beta = -float(l[4])

                chrom_dict[chrom]['nts'].append(nt)                
                z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                N = float(l[7])
                if N<=1:
                    print 'N is too small?? %0.2f'%N
                    N = 1
                chrom_dict[chrom]['betas'].append(z / sp.sqrt(N))     
                chrom_dict[chrom]['zs'].append(z)     
                chrom_dict[chrom]['Ns'].append(N)
            if line_i%100000==0:
                print line_i
    
        assert sp.all(sp.isreal(chrom_dict[chrom]['betas'])), 'WTF?'

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not ss_id in h5f.keys(), 'Summary stats with this name are already in the HDF5 file?'
    ssg = h5f.create_group(ss_id)
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['ps'], chrom_dict[chrom]['zs'],
                 chrom_dict[chrom]['Ns'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        zs = []
        ns = []
        prev_pos = -1
        for pos, sid, nt, beta, p, z, n in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            zs.append(z)
            ns.append(n)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        num_snps +=len(betas)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        g.create_dataset('zs', data=zs)
        g.create_dataset('Ns', data=ns)
        h5f.flush()
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps
            

def coordinate_sum_stats(comb_hdf5_file, coord_hdf5_file, filter_ambiguous_nts=True, only_common_snps=True):
    """
    
    """
    h5f = h5py.File(comb_hdf5_file)
    oh5f = h5py.File(coord_hdf5_file)
    sums_ids = h5f.keys()
    sums_ids = [x.encode('UTF8') for x in sums_ids]
    print 'Combining datasets: '+' '.join(sums_ids)
    oh5f.create_dataset('sums_ids', data=sums_ids)
    for chrom in range(1,23):
        chrom_str = 'chrom_%d' % chrom
        if only_common_snps:
            common_sids = set(h5f[sums_ids[0]][chrom_str]['sids'][...])
            for sums_id in sums_ids[1:]:
                chr_g = h5f[sums_id][chrom_str]
                sids = chr_g['sids'][...]
                print len(sids)
                common_sids = common_sids.intersection(sids)
            num_sids = len(common_sids)
            common_sids = sp.array(list(common_sids))
            print 'Found %d SNPs in common on chromosome %d'%(num_sids,chrom)
    
            #Use order and information from first summary stats dataset.
            chr_g = h5f[sums_ids[0]][chrom_str]
            sids = chr_g['sids'][...]
            sids_map = sp.in1d(sids, common_sids)
            filtered_sids = sids[sids_map] 
        
            if filter_ambiguous_nts:
                #Filter SNPs with ambiguous NTs
                ok_sids = set()
                for sums_id in sums_ids:
                    chr_g = h5f[sums_id][chrom_str]
                    sids1 = chr_g['sids'][...]
                    sids_map = sp.in1d(sids1, common_sids)
                    nts1 = chr_g['nts'][...][sids_map]
                    for sid, nt in it.izip(filtered_sids,nts1):
                        nt_tuple = tuple(nt)
                        if  nt_tuple in ok_nts:
                            ok_sids.add(sid)
                
                print "%d SNPs were found with ambiguous or weird nucleotides in some dataset."%(len(sids) - len(ok_sids))
                #If ambiguous SNPs found, then update the SNP map
                if len(ok_sids)< len(common_sids):
                    common_sids = sp.array(list(ok_sids))
                    chr_g = h5f[sums_ids[0]][chrom_str]
                    sids = chr_g['sids'][...]
                    sids_map = sp.in1d(sids, common_sids)
            else:
                pass
        order_sids = chr_g['sids'][...][sids_map]
        positions = chr_g['positions'][...][sids_map]
        nts = chr_g['nts'][...][sids_map]
        out_chr_g = oh5f.create_group(chrom_str)
        out_chr_g.create_dataset('sids', data = order_sids)
        out_chr_g.create_dataset('positions', data = positions)
        out_chr_g.create_dataset('nts', data = nts)
        
        
        #Retrieve effect estimates for other summary statistics.
        for sums_id in sums_ids:
            chr_g = h5f[sums_id][chrom_str]
            sids = chr_g['sids'][...]
            ss_spec_sids_map = sp.in1d(sids, common_sids)
            sids = sids[ss_spec_sids_map]
            betas = chr_g['betas'][...][ss_spec_sids_map]
            zs = chr_g['zs'][...][ss_spec_sids_map]
            pvals = chr_g['ps'][...][ss_spec_sids_map]
            nts2 = chr_g['nts'][...][ss_spec_sids_map]
            
            #Verify order..
            if not sp.all(order_sids==sids):
                print 'Need to re-order SNPs in summary statistics data %s'%sums_id
                snp_info_map = dict(it.izip(sids, it.izip(betas, zs, pvals, nts2)))
                num_snps = len(betas)
                betas = sp.empty(num_snps,dtype='single')
                zs = sp.empty(num_snps,dtype='single')
                pvals = sp.empty(num_snps,dtype='single')
                nts2 = []
                for i, sid in enumerate(order_sids):
                    snp_info = snp_info_map[sid]
                    betas[i]=snp_info[0]
                    zs[i]=snp_info[1]
                    pvals[i]=snp_info[2] 
                    nts2.append(snp_info[3])   
                nts2 = sp.array(nts2)
                sids = order_sids           

            assert sp.all(sp.isreal(betas)), 'WTF?'
            
            
            #Check nucleotide match, try flipping, etc, and perhaps post-filter data...
            for sid_i, sid, nt1, nt2 in it.izip(range(len(betas)), sids, nts, nts2):
                if sp.all(nt1==nt2):
                    continue
                else:
                    os_nt2 = sp.array([opp_strand_dict[nt2[0]], opp_strand_dict[nt2[1]]])
                    if sp.all(nt1 == os_nt2):
                        continue
                    else:
                        flip_nts = (nt1[1] == nt2[0] and nt1[0] == nt2[1]) or (nt1[1] == nt2[0] and nt1[0] == nt2[1])
                        #Try flipping the SS
                        if flip_nts:
                            betas[sid_i] = -betas[sid_i]                        
                            zs[sid_i] = -zs[sid_i]                        
                        else:
                            print "Nucleotides don't match after all? sid_i=%d, sid=%s, nt1=%s, nt2=%s" % (sid_i, sid, str(nt1), str(nt2))
                            betas[sid_i] = 0
                            zs[sid_i] = 0
                            pvals[sid_i] = 1
            
            betas = sp.array(betas, dtype='single')
            out_chr_ss_g = out_chr_g.create_group(sums_id)
            out_chr_ss_g.create_dataset('ps', data=pvals)
            out_chr_ss_g.create_dataset('betas', data=betas)
            out_chr_ss_g.create_dataset('zs', data=zs)
        

def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

                          
    long_options_list = ['ssfiles=', 'combfile=', 'coordfile=', 'sslabels=', '1KGpath=', 'ssf_format=','help']

    p_dict = {'ssfiles':None, 'combfile':None, 'coordfile':None, 'sslabels':None, '1KGpath':'/Users/bjarnivilhjalmsson/data/1Kgenomes/', 
              'ssf_format':'BASIC',}

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
            elif opt =="--ssfiles": p_dict['ssfiles'] = arg.split(',')
            elif opt =="--sslabels": p_dict['sslabels'] = arg.split(',')
            elif opt == "--combfile": p_dict['combfile'] = arg
            elif opt == "--coordfile": p_dict['coordfile'] = arg
            elif opt == "--1KGpath": p_dict['1KGpath'] = arg
            elif opt == "--ssf_format": p_dict['ssf_format'] = arg
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
    assert p_dict['combfile'] is not None, 'Combined SS file is missing.'
    comb_hdf5_file = p_dict['combfile']

    if p_dict['ssfiles'] is not None:                
        ssfiles = p_dict['ssfiles']
        if p_dict['sslabels'] is None:
            ss_labels = ['LABEL%d'%i for i in range(1,1+len(p_dict['ssfiles']))]
        else:
            ss_labels = p_dict['sslabels']

        assert len(ss_labels) == len(ssfiles), "There number of labels doens't match the number of SS files"
        
        assert p_dict['1KGpath'] is not None, 'Path to 1K Genomes is missing.'
        
        for ss_file, ss_label in zip(ssfiles,ss_labels):
            if p_dict['ssf_format']=='BASIC':
                parse_sum_stats_basic(ss_file,comb_hdf5_file,ss_label,KGpath=p_dict['1KGpath'])
            else:
                raise Exception('Unknown GWAS summary statistics file format')

    if p_dict['coordfile'] is not None:
        print 'Coordinating summary statistic datasets'
        coord_hdf5_file = p_dict['coordfile']
    
        coordinate_sum_stats(comb_hdf5_file, coord_hdf5_file)
    
