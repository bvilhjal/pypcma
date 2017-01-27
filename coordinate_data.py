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

__updated__ = '2017-01-24'



# Step 1:  Parse all summary statistics into one HDF5 file.  Do not coordinate.
# Step 2:  For each pair of summary statistics calculate LDscore regression.
# Step 3:  

import sum_stats


ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')]) 
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}


lc_2_cap_map = {'a':'A', 'c':'C', 'g':'G', 't':'T'}



def parse_all_sum_stats():
    """
    Code for parsing all of the summary stastics on the cluster.
    
    DIAGRAM_T2D    GABRIEL_ASTHMA    GEFOS_BMD-FOREARM    GEFOS_BMD-NECK    GEFOS_BMD-SPINE    GIANT_BMI    GIANT_HEIGHT    GIANT_HIP    GIANT_WC    GIANT_WHR    
    GLG_HDL    GLG_LDL    GLG_TC    GLG_TG    ICBP_DBP    ICBP_MAP    ICBP_PP    ICBP_SBP    
    IIBDGC_CD    IIBDGC_IBD    IIBDGC_UC    
    MAGIC_FAST-GLUCOSE    MAGIC_FAST-INSULIN    MAGIC_HOMA-B    MAGIC_HOMA-IR    
    RA_RA
    """
#     diagram_parse_str = '%run coordinate_data --ssfiles=/faststorage/project/PCMA/3_SUMSTATS/DIAGRAMv3.GWAS.T2D/DIAGRAMv3.2012DEC17.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/DIAGRAMv3.GWAS.T2D/DIAGRAM_T2D.hdf5 --sslabels=DIAGRAM_T2D --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
#     gabriel_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GABRIEL_asthma/gabriel_asthma.txt --combfile=/project/PCMA/faststorage/3_SUMSTATS/GABRIEL_asthma/GABRIEL_ASTHMA.hdf5 --sslabels=GABRIEL_ASTHMA --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
#     gefoss_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/fa2stu.MAF0_.005.pos_.out__0,/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/fn2stu.MAF0_.005.pos_.out_,/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/ls2stu.MAF0_.005.pos_.out_ --combfile=/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_BMD.hdf5 --sslabels=GEFOS_BMD-FOREARM,GEFOS_BMD-NECK,GEFOS_BMD-SPINE --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
#     giant_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_HIP_COMBINED_EUR.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_WC_COMBINED_EUR.txt,/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT_2015_WHR_COMBINED_EUR.txt  --combfile=/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT.hdf5 --sslabels=GIANT_BMI,GIANT_HEIGHT,GIANT_HIP,GIANT_WC,GIANT_WHR --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
#     glc_parse_str = '%run coordinate_data.py --ssfiles=/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_HDL.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_LDL.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_TC.txt,/project/PCMA/faststorage/3_SUMSTATS/GLC/jointGwasMc_TG.txt  --combfile=/project/PCMA/faststorage/3_SUMSTATS/GIANT/GLC.hdf5 --sslabels=GLG_HDL,GLG_LDL,GLG_TC,GLG_TG --1KGpath=/project/PCMA/faststorage/3_SUMSTATS/1Kgenomes/ --ow'
#     icbp_parse_str = '%run coordinate_data --ssfiles=/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/dbp_phs000585.pha003589.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/map_phs000585.pha003591.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/pp_phs000585.pha003590.txt,/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/sbp_phs000585.pha003588.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/ICBP.hdf5 --sslabels=ICBP_DBP,ICBP_MAP,ICBP_PP,ICBP_SBP --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
#     iibdgc_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.CD.gwas.assoc,/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.IBD.gwas.assoc,/home/bjarni/PCMA/faststorage/3_SUMSTATS/IBD/iibdgc-trans-ancestry-summary-stats/EUR.UC.gwas.assoc --combfile=/faststorage/project/PCMA/3_SUMSTATS/IBD/IIBDGC.hdf5 --sslabels=IIBDGC_CD,IIBDGC_IBD,IIBDGC_UC --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
#     magic_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_FastingGlucose.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_FastingInsulin.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_HOMA-B.txt,/home/bjarni/PCMA/faststorage/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC_ln_HOMA-IR.txt --combfile=/faststorage/project/PCMA/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC.hdf5 --sslabels=MAGIC_FAST-GLUCOSE,MAGIC_FAST-INSULIN,MAGIC_HOMA-B,MAGIC_HOMA-IR --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
#     ra_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/RA/RA_GWASmeta_European_v2.txt --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/RA/RA.hdf5 --sslabels=RA_RA --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'
#     
#     concatenate_ss_h5files(['/faststorage/project/PCMA/3_SUMSTATS/DIAGRAMv3.GWAS.T2D/DIAGRAM_T2D.hdf5',
#                             '/project/PCMA/faststorage/3_SUMSTATS/GABRIEL_asthma/GABRIEL_ASTHMA.hdf5',
#                             '/project/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_BMD.hdf5',
#                             '/project/PCMA/faststorage/3_SUMSTATS/GIANT/GIANT.hdf5',
#                             '/project/PCMA/faststorage/3_SUMSTATS/GIANT/GLC.hdf5',
#                             '/faststorage/project/PCMA/3_SUMSTATS/ICPB_bloodPress/ICBP.hdf5',
#                             '/faststorage/project/PCMA/3_SUMSTATS/IBD/IIBDGC.hdf5',
#                             '/faststorage/project/PCMA/3_SUMSTATS/MAGIC_glycaemic.traits/MAGIC.hdf5',
#                             '/home/bjarni/PCMA/faststorage/3_SUMSTATS/RA/RA.hdf5'], '/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb.hdf5')
#     
#     coordinate_sum_stats_w_missing('/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb.hdf5', '/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb_coord.hdf5', '/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/', 
#                                    only_common_snps=False)
    
    sum_stats.hdf5_coord_file_2_txt('/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb_coord.hdf5', '/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb_coord_zs.txt',
                          '/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb_coord_weights.txt', '/home/bjarni/PCMA/faststorage/3_SUMSTATS/comb_coord_ps.txt',
                          ['DIAGRAM_T2D', 'GABRIEL_ASTHMA', 'GEFOS_BMD-FOREARM', 'GEFOS_BMD-NECK', 'GEFOS_BMD-SPINE', 'GIANT_BMI', 'GIANT_HEIGHT',
                           'GIANT_HIP', 'GIANT_WC', 'GIANT_WHR', 'GLG_HDL', 'GLG_LDL', 'GLG_TC', 'GLG_TG', 'ICBP_DBP', 'ICBP_MAP', 'ICBP_PP',
                           'ICBP_SBP', 'IIBDGC_CD', 'IIBDGC_IBD', 'IIBDGC_UC', 'MAGIC_FAST-GLUCOSE', 'MAGIC_FAST-INSULIN', 'MAGIC_HOMA-B',
                           'MAGIC_HOMA-IR', 'RA_RA'])
    
#     teslovich_parse_str = '%run coordinate_data --ssfiles=/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl,/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TG_ONE_Europeans.tbl --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/TESLOVITCH/TESLOVICH_comb.hdf5 --sslabels=TAG_cpd,TAG_evrsmk,TAG_former,TAG_logonset --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --ow'

    # Then merge into big file:
    


def process_all_sum_stats():
    # gefoss_coord_str = '%run coordinate_data --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_comb.hdf5  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_coord.hdf5 --iq_range=90,100 --min_maf=0.02 --ow'
    # gefoss_LDprune_str = '%run LD_prune_sum_stats.py  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/GEFOS_osteoporosis/GEFOS_coord.hdf5 --LDradius=250 --out=/faststorage/project/PCMA/1_DATA/GEFOS --LDthres=0.2'    
    
    # icbp_coord_str = '%run coordinate_data --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_comb.hdf5  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_coord.hdf5 --iq_range=50,100 --min_maf=0.02 --ow'
    # icbp_LDprune_str = '%run LD_prune_sum_stats.py  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_coord.hdf5 --LDradius=100 --out=/faststorage/project/PCMA/1_DATA/ICBP --LDthres=0.2'
    pass


