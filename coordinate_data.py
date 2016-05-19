"""
Coordinate summary statistics datasets for performing multiple-trait analysis.

In order to conduct the multiple-trait analysis the following steps are necessary.

1. Parse the summary statistics and output in a internal hdf5 formatted file.
2. Coordinate the summary statistics and output a coordinated file.
    - At this step several filtering options are possible.
3. LD-prune the summary statistics.
4. Perform the PCMA analysis.

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


"""
Bulick-Sullivan et al., NG 2015:
  
Summary statistic data sets. We selected traits for inclusion in the analysis
for the main text via the following procedure:
1. Begin with all publicly available non-sex-stratified European-only
summary statistics.
2. Remove studies that do not provide signed summary statistics.
3. Remove studies not imputed to at least phase 2 of HapMap.
4. Remove studies that adjust for heritable covariates.
5. Remove all traits with a heritability z score below 4. Genetic correlation
estimates for traits with a heritability z score below 4 are generally
too noisy to report.
6. Prune clusters of correlated phenotypes (for example, obesity
classes 1–3) by picking the trait from each cluster with the highest
heritability z score.
We then applied the following filters (implemented in the script munge_
sumstats.py included with ldsc):
1. For studies that provide a measure of imputation quality, filter to INFO
above 0.9.
2. For studies that provide sample MAF, filter to sample MAF above 1%.
3. To restrict to well-imputed SNPs in studies that do not provide a measure
of imputation quality, filter to SNPs in the HapMap 3 panel61 with a 1000 
Genomes Project EUR (European) MAF above 5%, which tend to be well
imputed in most studies. This step should be skipped if INFO scores are
available for all studies.
4. If the sample size varies from SNP to SNP, remove SNPs with an effective
sample size less than 0.67 times the 90th percentile of sample size.
5. For meta-analyses with specialty chips (for example, the Metabochip),
remove SNPs with a sample size above the maximum GWAS sample size.
6. Remove indels and structural variants.
7. Remove strand-ambiguous SNPs.
8. Remove SNPs whose alleles do not match the alleles in the 1000 Genomes
Project.
"""

#Step 1:  Parse all summary statistics into one HDF5 file.  Do not coordinate.
#Step 2:  For each pair of summary statistics calculate LDscore regression.
#Step 3:  


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
import os 


ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
ok_nts = set([('A', 'G'), ('G', 'A'), ('A', 'C'), ('C', 'A'),('G', 'T'), ('T', 'G'),('C', 'T'), ('T', 'C')])
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

valid_nts = set(['A','T','C','G'])

lc_2_cap_map = {'a':'A', 'c':'C', 'g':'G', 't':'T'}

def parse_1KG_snp_info(KGenomes_prefix='/Users/bjarnivilhjalmsson/data/1Kgenomes/', 
                       outfile='/Users/bjarnivilhjalmsson/data/1Kgenomes/snps2.hdf5',
                       filter_ambiguous=True):
    of = h5py.File(outfile)
    for chrom_i in range(1,23):
        print 'Chromosome %d'%chrom_i
        sids = []
        positions = []
        chromosomes = []
        nts = []
        eur_mafs=[]
        fn = '%sALL_1000G_phase1integrated_v3_chr%d_impute.legend.gz'%(KGenomes_prefix,chrom_i)
        with gzip.open(fn) as f:
            f.next()
            line_i=0
            for line in f:
                l = line.split()
                nt1=l[2]
                nt2=l[3]
                if nt1 not in valid_nts:
                    continue
                if nt2 not in valid_nts:
                    continue
                if filter_ambiguous and (nt1,nt2) in ambig_nts:
                    continue
                line_i +=1
                sids.append(l[0])
                positions.append(int(l[1]))
                chromosomes.append(chrom_i)
                nts.append([nt1,nt2])
                eur_mafs.append(float(l[11]))
                if line_i%100000==0:
                    print line_i
        nts = sp.array(nts)
        print 'Constructed nts'
        g = of.create_group('chrom_%d' % chrom_i)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        g.create_dataset('chromosomes', data=chromosomes)
        g.create_dataset('eur_mafs', data=eur_mafs)
        g.create_dataset('nts', data=nts)
        of.flush()
    of.close()


def get_sid_pos_map(sids, KGenomes_prefix):
    h5fn = '%ssnps.hdf5'%(KGenomes_prefix)
    h5f = h5py.File(h5fn,'r')
    sid_map = {}
    for chrom_i in range(1,23):
        cg = h5f['chrom_%d' % chrom_i]
        sids_1k = cg['sids'][...]
        sids_filter_1k = sp.in1d(sids_1k, sp.array(sids)) 
        common_sids = sids_1k[sids_filter_1k]
        common_positions = cg['positions'][sids_filter_1k]
        eur_mafs = cg['eur_mafs'][sids_filter_1k]
        for sid,pos,eur_maf in it.izip(common_sids,common_positions,eur_mafs):
            sid_map[sid]={'pos':pos, 'chrom':chrom_i, 'eur_maf':eur_maf}
    return sid_map


headers = {'SSGAC1':['MarkerName', 'Effect_Allele', 'Other_Allele', 'EAF', 'Beta', 'SE', 'Pvalue'],
           'SSGAC2':['MarkerName', 'Effect_Allele', 'Other_Allele', 'EAF', 'OR', 'SE', 'Pvalue'],
           'CHIC':['SNP', 'CHR', 'BP', 'A1', 'A2', 'FREQ_A1', 'EFFECT_A1', 'SE', 'P'],
           'GCAN':['chromosome', 'position', 'SNP', 'reference_allele', 'other_allele', 'eaf', 'OR', 
                             'OR_se', 'OR_95L', 'OR_95U', 'z', 'p_sanger', '_-log10_p-value', 'q_statistic', 
                             'q_p-value', 'i2', 'n_studies', 'n_samples', 'effects'],
           'TESLOVICH':['MarkerName', 'Allele1', 'Allele2', 'Weight', 'GC.Zscore', 'GC.Pvalue', 'Overall', 'Direction'],
           'GIANT1':['MarkerName', 'Allele1', 'Allele2', 'FreqAllele1HapMapCEU', 'b', 'se', 'p', 'N'],
           'GIANT1b':['MarkerName', 'Allele1', 'Allele2', 'Freq.Allele1.HapMapCEU', 'b', 'SE', 'p', 'N'],
           'GIANT1c':['MarkerName', 'Chr', 'Pos', 'Allele1', 'Allele2', 'FreqAllele1HapMapCEU', 'b', 'se', 'p', 'N'],
           'GIANT2':['SNP', 'A1', 'A2', 'Freq1.Hapmap', 'b', 'se', 'p', 'N'],
           'MAGIC':['snp', 'effect_allele', 'other_allele', 'maf', 'effect', 'stderr', 'pvalue'],
           'CARDIoGRAM':['SNP', 'chr_pos_(b36)', 'reference_allele', 'other_allele', 'ref_allele_frequency', 'pvalue', 'het_pvalue', 'log_odds', 'log_odds_se', 'N_case', 'N_control', 'model'],
           'DIAGRAM':['SNP', 'CHROMOSOME', 'POSITION', 'RISK_ALLELE', 'OTHER_ALLELE', 'P_VALUE', 'OR', 'OR_95L', 'OR_95U', 'N_CASES', 'N_CONTROLS'],
           'TAG':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A', 'FRQ_U', 'INFO', 'OR', 'SE', 'P'],
           'CD':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_5956', 'FRQ_U_14927', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa'],
           'UC':['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_6968', 'FRQ_U_20464', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa'],
           'IBD': ['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_12882', 'FRQ_U_21770', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa'],
           'GEFOS':['chromosome', 'position', 'rs_number', 'reference_allele', 'other_allele', 'eaf', 'beta', 'se', 'beta_95L', 'beta_95U', 'z', 'p-value', '_-log10_p-value', 'q_statistic', 'q_p-value', 'i2', 'n_studies', 'n_samples', 'effects'],
           'RA':['SNPID','Chr','Position(hg19)','A1','A2','OR(A1)','OR_95%CIlow','OR_95%CIup','P-val'],
           'ASTHMA':['Chr', 'rs', 'position', 'Allele_1', 'Allele_2', 'freq_all_1_min', 'freq_all_1_max', 'OR_fix', 'ORl_fix', 'ORu_fix', 'P_fix'],
           'ICBP': ['ID', 'Analysis', 'ID', 'SNP', 'ID', 'P-value', 'Rank', 'Plot', 'data', 'Chr', 'ID', 'Chr', 'Position', 'Submitted', 'SNP', 'ID', 'ss2rs', 'rs2genome', 'Allele1', 'Allele2', 'Minor', 'allele', 'pHWE', 'Call', 'Rate', 'Effect', 'SE', 'R-Squared', 'Coded', 'Allele', 'Sample', 'size', 'Bin', 'ID'],
           'GLC': ['SNP_hg18', 'SNP_hg19', 'rsid', 'A1', 'A2', 'beta', 'se', 'N', 'P-value', 'Freq.A1.1000G.EUR'],
           }

def parse_sum_stats(filename,
                        comb_hdf5_file,
                        ss_id,
                        KGpath,
                        bimfile=None,):
    """
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
        
        line = f.next()
        header = line.split()        
        if header==['hg19chrc', 'snpid', 'a1', 'a2', 'bp', 'info', 'or', 'se', 'p', 'ngt'] or header==headers['TAG'] or header==headers['CD'] or header==headers['UC'] or header==headers['IBD'] or header==headers['ASTHMA']:
            for line in f:
                l = line.split()
                sids.append(l[1])
        elif header==['Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue'] or header==headers['GCAN'] or header==headers['GEFOS'] or header==headers['ICBP'] or header==headers['GLC']:
#             line_count =0
            for line in f:
#                 line_count +=1
#                 if line_count<100000:
                l = line.split()
                sids.append(l[2])
#                 elif random.random()<0.1:
#                     l = line.split()
#                     sids.append(l[2])
                        
        else:
            for line in f:
                l = line.split()
                sids.append(l[0])
    sid_map = get_sid_pos_map(sids,KGpath)
    assert len(sid_map)>0, 'WTF?'

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        line = f.next()
        header = line.split()
        line_i = 0
        if header== ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'or', 'se', 'pval', 'info', 'ngt', 'CEUaf']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[5]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['hg19chrc', 'snpid', 'a1', 'a2', 'bp', 'info', 'or', 'se', 'p', 'ngt']:    
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[2], l[3]]
                    else:
                        nt = [l[3], l[2]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header== ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'zscore', 'pval', 'CEUmaf']:     
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[5])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['SNP', 'CHR', 'BP', 'A1', 'A2', 'OR', 'SE', 'P', 'INFO', 'EUR_FRQ']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[5]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header==['Chromosome', 'Position', 'MarkerName', 'Effect_allele', 'Non_Effect_allele', 'Beta', 'SE', 'Pvalue']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[7])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[5])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i
        
        elif header ==headers['SSGAC1']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        elif header ==headers['SSGAC2']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[4]))
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        elif header ==headers['CHIC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[],'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[8])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        
        elif header==headers['GCAN']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[11])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = -sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i             
        
        elif header==headers['TESLOVICH']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [lc_2_cap_map[l[1]], lc_2_cap_map[l[2]]]
                    else:
                        nt = [lc_2_cap_map[l[2]], lc_2_cap_map[l[1]]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[3])))
                if line_i%100000==0:
                    print line_i                            
        
        elif header==headers['GIANT1'] or header==headers['GIANT1b'] or header==headers['GIANT2']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [l[1], l[2]]
                    else:
                        nt = [l[2], l[1]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[7])))
                if line_i%100000==0:
                    print line_i                                     
        elif header==headers['GIANT1c']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[6])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    chrom_dict[chrom]['weights'].append(int(float(l[9])))
                if line_i%100000==0:
                    print line_i   
        elif header==headers['MAGIC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                eur_maf = float(l[3])
                if d is not None and eur_maf>0:
                    pos = d['pos']
                    chrom = d['chrom']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[4])
                    if random.random()>0.5:
                        nt = [lc_2_cap_map[l[1]], lc_2_cap_map[l[2]]]
                    else:
                        nt = [lc_2_cap_map[l[2]], lc_2_cap_map[l[1]]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['CARDIoGRAM']: 

            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[7])
                    if random.random()>0.5:
                        nt = [l[2], l[3]]
                    else:
                        nt = [l[3], l[2]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[9]) +float(l[10])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['DIAGRAM']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[5])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[6]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[9]) +float(l[10])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['TAG']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[10])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[8])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['CD'] or header==headers['UC'] or header==headers['IBD']:
            #['CHR', 'SNP', 'BP', 'A1', 'A2', 'FRQ_A_12882', 'FRQ_U_21770', 'INFO', 'OR', 'SE', 'P', 'Direction', 'HetISqt', 'HetPVa']
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                info_score = float(l[7])
                d = sid_map.get(sid,None)
                if d is not None and info_score>0.8:
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[10])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = sp.log(float(l[8]))
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    direction = l[11]
                    num_studies = len(direction)-direction.count('?')                    
                    weight = num_studies
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['GEFOS']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                eur_maf = float(l[5])
                if d is not None and eur_maf>0 and eur_maf<1:
                    pos = d['pos']
                    chrom = d['chrom']
#                     eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[11])
                    chrom_dict[chrom]['ps'].append(pval)
                    raw_beta = float(l[6])
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
#                     weight = (z/raw_beta)**2
#                     weight = float(l[16])  # Number of studies used.
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
#           'RA':['SNPID','Chr','Position(hg19)','A1','A2','OR(A1)','OR_95%CIlow','OR_95%CIup','P-val'],
#           'ASTHMA':['Chr', 'rs', 'position', 'Allele_1', 'Allele_2', 'freq_all_1_min', 'freq_all_1_max', 'OR_fix', 'ORl_fix', 'ORu_fix', 'P_fix'],
        elif header==headers['RA']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[0]
                d = sid_map.get(sid,None)
                if d is not None:
                    raw_beta = sp.log(float(l[5]))
                    if raw_beta==0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['ASTHMA']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[1]
                d = sid_map.get(sid,None)
                if d is not None:
                    raw_beta = sp.log(float(l[7]))
                    pval = float(l[10])
                    if raw_beta==0 or pval == 0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    chrom_dict[chrom]['ps'].append(pval)
                    if random.random()>0.5:
                        nt = [l[3], l[4]]
                    else:
                        nt = [l[4], l[3]]
                        raw_beta = -raw_beta
    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = z**2/((raw_beta**2)*2*eur_maf*(1-eur_maf))
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   
        elif header==headers['ICBP']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                coded_allele = l[16]
                if d is not None and coded_allele in valid_nts:
#                     raw_beta = sp.log(float(l[7]))
                    pval = float(l[3])
                    if pval == 0:
                        continue
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    chrom_dict[chrom]['ps'].append(pval)
#                     if random.random()>0.5:
                    nt = [l[11], l[12]]
                    sign = 1
#                     else:
#                         sign = -1
                    if coded_allele==nt[1] or opp_strand_dict[coded_allele]==nt[1]:
                        nt = [l[12], l[11]]
                        sign = -1#*sign
#                     else:
#                         assert coded_allele==nt[0] or opp_strand_dict[coded_allele]==nt[0]
                    
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sign * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[17])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   

           #'GLC': ['SNP_hg18', 'SNP_hg19', 'rsid', 'A1', 'A2', 'beta', 'se', 'N', 'P-value', 'Freq.A1.1000G.EUR'],
        elif header==headers['GLC']:
            for line in f:
                line_i +=1
                l = line.split()
                sid = l[2]
                d = sid_map.get(sid,None)
                if d is not None :
                    raw_beta = float(l[5])
                    pval = float(l[8])
                    pos = d['pos']
                    chrom = d['chrom']
                    eur_maf = d['eur_maf']
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'zs':[], 'nts': [], 'sids': [], 
                                             'positions': [], 'eur_maf':[], 'weights':[]}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    chrom_dict[chrom]['eur_maf'].append(eur_maf)
                    chrom_dict[chrom]['ps'].append(pval)
                    if random.random()>0.5:
                        nt = [l[3], l[4]]                    
                    else:
                        nt = [l[4], l[3]]                    
                        raw_beta = -raw_beta
                    chrom_dict[chrom]['nts'].append(nt)                
                    z = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['zs'].append(z)     
                    weight = float(l[7])
                    chrom_dict[chrom]['weights'].append(weight)
                if line_i%100000==0:
                    print line_i   

        else:
            raise Exception('Wrong or unknown file format')
    
        assert sp.all(sp.isreal(chrom_dict[1]['zs'])), 'WTF?'

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not ss_id in h5f.keys(), 'Summary stats with this name are already in the HDF5 file?'
    ssg = h5f.create_group(ss_id)
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['ps'], chrom_dict[chrom]['zs'], chrom_dict[chrom]['eur_maf'], 
                 chrom_dict[chrom]['weights'])
        sl.sort()
        ps = []
        nts = []
        sids = []
        positions = []
        zs = []
        eur_mafs = []
        weights = []
        prev_pos = -1
        for pos, sid, nt, p, z, eur_maf, weight in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            zs.append(z)
            eur_mafs.append(eur_maf)
            weights.append(weight)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('eur_mafs', data=eur_mafs)
        g.create_dataset('positions', data=positions)
        g.create_dataset('zs', data=zs)
        g.create_dataset('weights', data=weights)
        num_snps +=len(sids)
        h5f.flush()
    
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps




def coordinate_sum_stats(comb_hdf5_file, coord_hdf5_file, filter_ambiguous_nts=True,
                         ss_labs=None, weight_min=0, weight_max_diff=1, 
                         outlier_thres=0, sd_thres=0, iq_range=None, min_maf=0):
    """
    Coordinate multiple summary statistics
    """
    h5f = h5py.File(comb_hdf5_file)
    oh5f = h5py.File(coord_hdf5_file)
    if ss_labs is None:
        sums_ids = h5f.keys()
        sums_ids = [x.encode('UTF8') for x in sums_ids]
    else:
        all_sums_ids = h5f.keys()
        all_sums_ids = [x.encode('UTF8') for x in all_sums_ids]
        sums_ids = []
        for ss_lab in all_sums_ids:
            if ss_lab in ss_labs:
                sums_ids.append(ss_lab)
    print 'Coordinating summary statistics: '+' '.join(sums_ids)
    oh5f.create_dataset('sums_ids', data=sums_ids)
    for chrom in range(1,23):
        chrom_str = 'chrom_%d' % chrom
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
        #Store the OK sids in the common_sids, which we'll update accordingly..
        
        chr_g = h5f[sums_ids[0]][chrom_str]
        sids = chr_g['sids'][...]
        sids_map = sp.in1d(sids, common_sids)
        filtered_sids = sids[sids_map] 
        assert len(filtered_sids)==len(common_sids),'WTF?'
        common_sids = filtered_sids  #To ensure that they are ordered by position, etc.
    
        if filter_ambiguous_nts:
            #Filter SNPs with ambiguous NTs
            ok_sids = set()
            for sums_id in sums_ids:
                chr_g = h5f[sums_id][chrom_str]
                sids1 = chr_g['sids'][...]
                sids_map = sp.in1d(sids1, common_sids)
                nts1 = chr_g['nts'][...][sids_map]
                sids1 = sids1[sids_map]
                assert len(sids1)==len(common_sids),'WTF?'
                for sid, nt in it.izip(sids1,nts1):
                    nt_tuple = tuple(nt)
                    if  nt_tuple in ok_nts:
                        ok_sids.add(sid)
            
            print "%d SNPs were found with ambiguous or weird nucleotides in some dataset."%(len(sids) - len(ok_sids))
            #If ambiguous SNPs found, then update the SNP map
            if len(ok_sids)< len(common_sids):
                chr_g = h5f[sums_ids[0]][chrom_str]
                sids = chr_g['sids'][...]
                sids_map = sp.in1d(sids, sp.array(list(ok_sids)))
                common_sids = sids[sids_map]  #To ensure that they are ordered by the order in the first sum stats
        
        
        if weight_min>0 or weight_max_diff<1 or outlier_thres>0 or sd_thres>0 or min_maf>0 or iq_range is not None:
            #Filtering SNPs with weight differences.           
            #Calculating the relative weights
            n_snps = len(common_sids)
            n_sums = len(sums_ids)
            rel_weights_mat = sp.zeros((n_snps,n_sums))
            weights_filter =sp.ones((n_snps),dtype='bool')
            for s_i, sums_id in enumerate(sums_ids):
                chr_g = h5f[sums_id][chrom_str]
                sids1 = chr_g['sids'][...]
                sids_map1 = sp.in1d(sids1, common_sids)
                weights=chr_g['weights'][...][sids_map1]
                eur_mafs=chr_g['eur_mafs'][...][sids_map1]
#                 print weights
                sids1 = sids1[sids_map1]
                
                #Verify order..
                if not sp.all(sids1==common_sids):
                    print 'Need to re-order SNPs in summary statistics data %s'%sums_id
                    snp_map = dict(it.izip(common_sids, range(len(sids1))))
                    snp_order  = [snp_map[sid] for sid in sids1]
                    sids1 = sids1[snp_order]
                    weights = weights[snp_order]
                    eur_mafs=eur_mafs[snp_order]
                    
                #Filter mafs
                if min_maf>0:
                    weights_filter = weights_filter*(eur_mafs>min_maf)
                
                max_weight = sp.nanmax(weights)
                min_weight = sp.nanmin(weights)
                print max_weight,min_weight
                weights[sp.isnan(weights)]=min_weight

                if sp.isinf(max_weight):
                    inf_filter = sp.isinf(weights)
                    not_inf_weights = weights[~inf_filter]
                    max_weight = sp.nanmax(not_inf_weights)
                    weights[inf_filter]=max_weight

                median_weight = sp.median(weights)
                print 'Median weight: %0.2f; Minimum weight: %0.2f; Maximum weight: %0.2f'%(median_weight,min_weight,max_weight)
                                    
                #Outlier filter
                if outlier_thres>0:
                    weights_filter = weights_filter*(weights<median_weight+outlier_thres*max_weight)
                    weights_filter = weights_filter*(weights>median_weight-outlier_thres*max_weight)
                elif sd_thres>0:
                    weights_sd = sp.std(weights)
                    print 'Weights SD: ',weights_sd
                    weights_filter = weights_filter*(weights<median_weight+weights_sd*sd_thres)
                    weights_filter = weights_filter*(weights>median_weight-weights_sd*sd_thres)
                elif iq_range is not None:
                    w_min, w_max = sp.percentile(sp.array(weights,dtype='single'), [iq_range[0] ,iq_range[1]])
                    print iq_range
                    print 'Filtering with thresholds:',w_min, w_max
                    weights_filter = weights_filter*(weights>=w_min)
                    weights_filter = weights_filter*(weights<=w_max)
                    print 'Filter ratio:',sp.sum(weights_filter)/float(len(weights_filter))
                    
                    
                
                rel_weights_mat[:,s_i] = weights/float(max_weight)       


            #Calculating the minimum relative weight per SNP
            #Calculating the maximum difference in relative weights.
            min_rel_weights = rel_weights_mat.min(1)
            if weight_min>0:
                median_weight = sp.median(rel_weights_mat)
                max_weight = sp.nanmax(rel_weights_mat)
                min_weight = sp.nanmin(rel_weights_mat)
                min_filter = min_rel_weights>weight_min
                print 'Filter ratio:',sp.sum(weights_filter)/float(len(weights_filter))
                weights_filter = min_filter * weights_filter
                print 'Filter ratio:',sp.sum(weights_filter)/float(len(weights_filter))

                print 'Median weight: %0.2f; Minimum weight: %0.2f; Maximum weight: %0.2f'%(median_weight,min_weight,max_weight)

            max_diffs = sp.absolute(rel_weights_mat.max(1)-min_rel_weights)
            if weight_max_diff<1:
                weights_filter = (max_diffs<weight_max_diff)*weights_filter
                
            num_filtered_snps = len(weights_filter)-sp.sum(weights_filter)
            print 'Filter %d SNPs due to insufficient sample size/weights or to large sample size/weights differences.'%num_filtered_snps
            if num_filtered_snps>0:                    
                #Update SNP map
                common_sids = common_sids[weights_filter]

        chr_g = h5f[sums_ids[0]][chrom_str]
        sids = chr_g['sids'][...]
        sids_map = sp.in1d(sids, common_sids)

        order_sids = chr_g['sids'][...][sids_map]
        positions = chr_g['positions'][...][sids_map]
        eur_mafs = chr_g['eur_mafs'][...][sids_map]
        nts = chr_g['nts'][...][sids_map]
        out_chr_g = oh5f.create_group(chrom_str)
        out_chr_g.create_dataset('sids', data = order_sids)
        out_chr_g.create_dataset('positions', data = positions)
        out_chr_g.create_dataset('eur_mafs', data = eur_mafs)
        out_chr_g.create_dataset('nts', data = nts)
            
            
        
        #Retrieve effect estimates for other summary statistics.
        for sums_id in sums_ids:
            chr_g = h5f[sums_id][chrom_str]
            sids = chr_g['sids'][...]
            ss_spec_sids_map = sp.in1d(sids, common_sids)
            sids = sids[ss_spec_sids_map]
            zs = chr_g['zs'][...][ss_spec_sids_map]
            pvals = chr_g['ps'][...][ss_spec_sids_map]
            nts2 = chr_g['nts'][...][ss_spec_sids_map]
            if 'weights' in chr_g.keys():
                weights = chr_g['weights'][...][ss_spec_sids_map]
            
            #Verify order..
            if not sp.all(order_sids==sids):
                print 'Need to re-order SNPs in summary statistics data %s'%sums_id
                snp_info_map = dict(it.izip(sids, it.izip(zs, pvals, nts2)))
                num_snps = len(zs)
                zs = sp.empty(num_snps,dtype='single')
                pvals = sp.empty(num_snps,dtype='single')
                nts2 = []
                for i, sid in enumerate(order_sids):
                    snp_info = snp_info_map[sid]
                    zs[i]=snp_info[0]
                    pvals[i]=snp_info[1] 
                    nts2.append(snp_info[2])   
                nts2 = sp.array(nts2)
                sids = order_sids           

            assert sp.all(sp.isreal(zs)), 'WTF?'
            
            
            #Check nucleotide match, try flipping, etc, and perhaps post-filter data...
            for sid_i, sid, nt1, nt2 in it.izip(range(len(zs)), sids, nts, nts2):
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
                            zs[sid_i] = -zs[sid_i]                        
                        else:
                            print "Nucleotides don't match after all? sid_i=%d, sid=%s, nt1=%s, nt2=%s" % (sid_i, sid, str(nt1), str(nt2))
                            zs[sid_i] = 0
                            pvals[sid_i] = 1
            
            out_chr_ss_g = out_chr_g.create_group(sums_id)
            out_chr_ss_g.create_dataset('ps', data=pvals)
            out_chr_ss_g.create_dataset('zs', data=zs)
            out_chr_ss_g.create_dataset('weights', data = weights)
        

def coordinate_sum_stats_w_missing(comb_hdf5_file, coord_hdf5_file, KGpath, filter_ambiguous_nts=True, only_common_snps=True):
    """
    
    """
    snps_h5f = h5py.File(KGpath+'snps.hdf5','r')
    h5f = h5py.File(comb_hdf5_file,'r')
    oh5f = h5py.File(coord_hdf5_file)
    sums_ids = h5f.keys()
    sums_ids = [x.encode('UTF8') for x in sums_ids]
    print 'Combining datasets: '+' '.join(sums_ids)
    oh5f.create_dataset('sums_ids', data=sums_ids)
    for chrom in range(1,23):
        chrom_str = 'chrom_%d' % chrom
        snps_chrom_g = snps_h5f[chrom_str]
        all_sids = snps_chrom_g['sids'][...]
        total_sids = set(h5f[sums_ids[0]][chrom_str]['sids'][...])
        for sums_id in sums_ids[1:]:
            chr_g = h5f[sums_id][chrom_str]
            sids = chr_g['sids'][...]
            print len(sids)
            total_sids = total_sids.union(sids)
        num_sids = len(total_sids)
        total_sids = sp.array(list(total_sids))
        print 'Found %d SNPs in common on chromosome %d'%(num_sids,chrom)
    
        #Use order and information from first summary stats dataset.
        chr_g = h5f[sums_ids[0]][chrom_str]
        sids_map = sp.in1d(all_sids, total_sids)
        final_sids = all_sids[sids_map] 
        positions = snps_chrom_g['positions'][...][sids_map]
        nts = snps_chrom_g['nts'][...][sids_map]
        eur_mafs = snps_chrom_g['eur_mafs'][...][sids_map]
        
        num_snps = len(final_sids)
        sid_map = dict(it.izip(final_sids, range(num_snps)))
        
        out_chr_g = oh5f.create_group(chrom_str)
        out_chr_g.create_dataset('sids', data = final_sids)
        out_chr_g.create_dataset('eur_mafs', data = eur_mafs)
        out_chr_g.create_dataset('positions', data = positions)
        out_chr_g.create_dataset('nts', data = nts)
        
        
        #Retrieve effect estimates for other summary statistics.
        for sums_id in sums_ids:
            chr_g = h5f[sums_id][chrom_str]
            sids = chr_g['sids'][...]
            ss_spec_sids_map = sp.in1d(sids, final_sids)
            sids = sids[ss_spec_sids_map]
            zs = chr_g['zs'][...][ss_spec_sids_map]
            pvals = chr_g['ps'][...][ss_spec_sids_map]
            nts2 = chr_g['nts'][...][ss_spec_sids_map]
            if 'weights' in chr_g.keys():
                weights = chr_g['weights'][...][ss_spec_sids_map]
            
            final_zs = sp.empty(num_snps,dtype='single')
            final_zs.fill(sp.nan)
            final_pvals = sp.ones(num_snps,dtype='single')
            final_weights = sp.zeros(num_snps,dtype='single')
            #Check nucleotide match, try flipping, etc, and perhaps post-filter data...
            for i2, sid in enumerate(sids):
                i1 = sid_map[sid]
                nt1 = nts[i1]
                nt2 = nts2[i2]
                if sp.all(nt1==nt2):
                    final_zs[i1]=zs[i2]
                    final_pvals[i1]=pvals[i2]
                    final_weights[i1]=weights[i2]
                    continue
                else:
                    os_nt2 = sp.array([opp_strand_dict[nt2[0]], opp_strand_dict[nt2[1]]])
                    if sp.all(nt1 == os_nt2):
                        final_zs[i1]=zs[i2]
                        final_pvals[i1]=pvals[i2]
                        final_weights[i1]=weights[i2]
                        continue
                    else:
                        flip_nts = (nt1[1] == nt2[0] and nt1[0] == nt2[1]) or (nt1[1] == nt2[0] and nt1[0] == nt2[1])
                        #Try flipping the SS
                        if flip_nts:
                            final_zs[i1]=-zs[i2]
                            final_pvals[i1]=pvals[i2]
                            final_weights[i1]=weights[i2]
                        else:
                            print "Nucleotides don't match after all?  sid=%s, nt1=%s, nt2=%s" % (sid, str(nt1), str(nt2))
            
            out_chr_ss_g = out_chr_g.create_group(sums_id)
            out_chr_ss_g.create_dataset('ps', data=final_pvals)
            out_chr_ss_g.create_dataset('zs', data=final_zs)
            out_chr_ss_g.create_dataset('weights', data = final_weights)


def concatenate_ss_h5files(h5files, outfile, ss_labs = None):
    oh5f = h5py.File(outfile)
    for h5file in h5files:
        ih5f = h5py.File(h5file)
        if ss_labs is None:
            ok_ss_labs = ih5f.keys()
        else:
            ok_ss_labs = []
            for ss_lab in ih5f.keys():
                if ss_lab in ss_labs:
                    ok_ss_labs.append(ss_lab)
        for ss_lab in ok_ss_labs:
            ossg = oh5f.create_group(ss_lab)
            issg = ih5f[ss_lab]
            for chrom in range(1,23):
                chrom_str = 'chrom_%d' % chrom
                ochrom_g = ossg.create_group(chrom_str)
                ichrom_g = issg[chrom_str]
                ochrom_g.create_dataset('ps', data=ichrom_g['ps'][...])
                ochrom_g.create_dataset('nts', data=ichrom_g['nts'][...])
                ochrom_g.create_dataset('sids', data=ichrom_g['sids'][...])
                ochrom_g.create_dataset('eur_mafs', data=ichrom_g['eur_mafs'][...])
                ochrom_g.create_dataset('positions', data=ichrom_g['positions'][...])
                ochrom_g.create_dataset('zs', data=ichrom_g['zs'][...])
                ochrom_g.create_dataset('weights', data=ichrom_g['weights'][...])


def do_analysis():
    """
    Runs coordnations of multiple datasets.
    """
    #PGC..
    
    


def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

                          
    long_options_list = ['ssfiles=', 'combfile=', 'coordfile=', 'sslabels=', '1KGpath=', 'ssf_format=', 'weight_min=', 'weight_max_diff=', 
                         'outlier_thres=', 'sd_thres=', 'iq_range=', 'min_maf=','help', 'wmissing', 'ow']

    p_dict = {'ssfiles':None, 'combfile':None, 'coordfile':None, 'sslabels':None, '1KGpath':'/Users/bjarnivilhjalmsson/data/1Kgenomes/', 
              'ssf_format':'BASIC', 'wmissing':False, 'weight_min': 0.0, 'weight_max_diff': 1, 'outlier_thres':0, 'sd_thres':0, 
              'iq_range':None, 'ow':False, 'min_maf':0.01}

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
            elif opt == "--wmissing": p_dict['wmissing'] = True
            elif opt == "--weight_min": p_dict['weight_min'] = float(arg)
            elif opt == "--weight_max_diff": p_dict['weight_max_diff'] = float(arg)
            elif opt == "--outlier_thres": p_dict['outlier_thres'] = float(arg)
            elif opt == "--sd_thres": p_dict['sd_thres'] = float(arg)
            elif opt == "--iq_range": p_dict['iq_range'] = map(float,arg.split(','))
            elif opt == "--min_maf": p_dict['min_maf'] = float(arg)
            elif opt == "--ow": p_dict['ow'] = True
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict



def parse_all_sum_stats():
    """
    Code for parsing all of the summary stastics on the cluster.
    
    DIAGRAM_T2D    GABRIEL_ASTHMA    GEFOS_BMD-FOREARM    GEFOS_BMD-NECK    GEFOS_BMD-SPINE    GIANT_BMI    GIANT_HEIGHT    GIANT_HIP    GIANT_WC    GIANT_WHR    
    GLG_HDL    GLG_LDL    GLG_TC    GLG_TG    ICBP_DBP    ICBP_MAP    ICBP_PP    ICBP_SBP    
    IIBDGC_CD    IIBDGC_IBD    IIBDGC_UC    
    MAGIC_FAST-GLUCOSE    MAGIC_FAST-INSULIN    MAGIC_HOMA-B    MAGIC_HOMA-IR    
    RA_RA
    """
    diagram_parse_str = '%run coordinate_data --ssfiles=/faststorage/project/PCMA/3_SUMSTATS/DIAGRAMv3.GWAS.T2D/DIAGRAMv3.2012DEC17.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/DIAGRAMv3.GWAS.T2D/DIAGRAM_T2D.hdf5 --sslabels=DIAGRAM_T2D --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
    gabriel_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GABRIEL_asthma/gabriel_asthma.txt --combfile=/project/PCMA/faststorage/3_SUMSTATS/GABRIEL_asthma/GABRIEL_ASTHMA.hdf5 --sslabels=GABRIEL_ASTHMA --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
    gefoss_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/fa2stu.MAF0_.005.pos_.out__0,/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/fn2stu.MAF0_.005.pos_.out_,/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/ls2stu.MAF0_.005.pos_.out_ --combfile=/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_BMD.hdf5 --sslabels=GEFOS_BMD-FOREARM,GEFOS_BMD-NECK,GEFOS_BMD-SPINE --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
    giant_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_HIP_COMBINED_EUR.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_WC_COMBINED_EUR.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_WHR_COMBINED_EUR.txt  --combfile=/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT.hdf5 --sslabels=GIANT_BMI,GIANT_HEIGHT,GIANT_HIP,GIANT_WC,GIANT_WHR --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
    glc_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_HDL.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_LDL.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_TC.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_TG.txt  --combfile=/project/PCMA/faststorage/3_SUMSTATS/GIANT/GLC.hdf5 --sslabels=GLG_HDL,GLG_LDL,GLG_TC,GLG_TG --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
    icbp_parse_str = '%run coordinate_data --ssfiles=/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/dbp_phs000585.pha003589.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/map_phs000585.pha003591.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/pp_phs000585.pha003590.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/sbp_phs000585.pha003588.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/ICBP.hdf5 --sslabels=ICBP_DBP,ICBP_MAP,ICBP_PP,ICBP_SBP --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
    iibdgc_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.CD.gwas.assoc,/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc,/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.UC.gwas.assoc --combfile=/faststorage/project/PCMA/3_SUMSTATS/IBD/IIBDGC.hdf5 --sslabels=IIBDGC_CD,IIBDGC_IBD,IIBDGC_UC --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
    magic_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_FastingGlucose.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_FastingInsulin.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_HOMA-B.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_HOMA-IR.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC.hdf5 --sslabels=MAGIC_FAST-GLUCOSE,MAGIC_FAST-INSULIN,MAGIC_HOMA-B,MAGIC_HOMA-IR --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
    ra_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/RA/RA_GWASmeta_European_v2 --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/RA/RA.hdf5 --sslabels=RA_RA --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
    
    teslovich_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TESLOVICH_comb.hdf5 --sslabels=TAG_cpd,TAG_evrsmk,TAG_former,TAG_logonset --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'

    #Then merge into big file:
    


def process_all_sum_stats():
    #gefoss_coord_str = '%run coordinate_data --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_comb.hdf5  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_coord.hdf5 --iq_range=90,100 --min_maf=0.02 --ow'
    #gefoss_LDprune_str = '%run LD_prune_sum_stats.py  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_coord.hdf5 --LDradius=250 --out=/faststorage/project/PCMA/1_DATA/GEFOS --LDthres=0.2'    
    
    #icbp_coord_str = '%run coordinate_data --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_comb.hdf5  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_coord.hdf5 --iq_range=50,100 --min_maf=0.02 --ow'
    #icbp_LDprune_str = '%run LD_prune_sum_stats.py  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_coord.hdf5 --LDradius=100 --out=/faststorage/project/PCMA/1_DATA/ICBP --LDthres=0.2'
    pass

if __name__=='__main__':
    p_dict = parse_parameters()
    assert p_dict['combfile'] is not None, 'Combined SS file is missing.'
    comb_hdf5_file = p_dict['combfile']

    if p_dict['ssfiles'] is not None:                
        if os.path.isfile(comb_hdf5_file):
            if p_dict['ow']:
                print 'Overwriting the combfile: %s'%comb_hdf5_file
                os.remove(comb_hdf5_file)
            else:
                print 'The combfile %s already exists.  Please use the overwrite parameter.'%comb_hdf5_file
                sys.exit(0)
                
        ssfiles = p_dict['ssfiles']
        if p_dict['sslabels'] is None:
            ss_labels = ['LABEL%d'%i for i in range(1,1+len(p_dict['ssfiles']))]
        else:
            ss_labels = p_dict['sslabels']

        assert len(ss_labels) == len(ssfiles), "There number of labels doens't match the number of SS files"
        
        assert p_dict['1KGpath'] is not None, 'Path to 1K Genomes is missing.'
        
        for ss_file, ss_label in zip(ssfiles,ss_labels):
            parse_sum_stats(ss_file,comb_hdf5_file,ss_label,KGpath=p_dict['1KGpath'])
#             if p_dict['ssf_format']=='BASIC':
#                 parse_sum_stats_basic(ss_file,comb_hdf5_file,ss_label,KGpath=p_dict['1KGpath'])
#             elif p_dict['ssf_format']=='PGC':
#                 parse_sum_stats(ss_file,comb_hdf5_file,ss_label,KGpath=p_dict['1KGpath'])
#             else:
#                 raise Exception('Unknown GWAS summary statistics file format')

    if p_dict['coordfile'] is not None:
        print 'Coordinating summary statistic datasets'
        coord_hdf5_file = p_dict['coordfile']
        if os.path.isfile(coord_hdf5_file):
            if p_dict['ow']:
                print 'Overwriting the coord_hdf5_file: %s'%coord_hdf5_file
                os.remove(coord_hdf5_file)
            else:
                print 'The coord_hdf5_file %s already exists.  Please use the overwrite parameter.'%coord_hdf5_file
                sys.exit(0)

        if p_dict['wmissing']:
            coordinate_sum_stats_w_missing(comb_hdf5_file, coord_hdf5_file, p_dict['1KGpath'])
        else:
            coordinate_sum_stats(comb_hdf5_file, coord_hdf5_file, ss_labs=p_dict['sslabels'], 
                                 weight_min=p_dict['weight_min'], weight_max_diff=p_dict['weight_max_diff'], 
                                 outlier_thres=p_dict['outlier_thres'], sd_thres=p_dict['sd_thres'], iq_range=p_dict['iq_range'],
                                 min_maf=p_dict['min_maf'])
    


