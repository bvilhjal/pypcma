import analyze_results as ar
import LD_prune_sum_stats as ldps
import coordinate_data as cd


#GEFOS

#ICBP
%run coordinate_data --combfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_comb.hdf5  --1KGpath=/faststorage/project/PCMA/3_SUMSTATS/1Kgenomes/ --coordfile=/home/bjarni/PCMA/faststorage/3_SUMSTATS/ICPB_bloodPress/ICBP_coord.hdf5 --iq_range=80,100 --min_maf=0.02 --ow
ar.run_all_ts('/home/bjarni/PCMA/faststorage/1_DATA/ICBP_zs_ldpruned_no_weights.txt','/home/bjarni/PCMA/faststorage/1_DATA/ICBP_zs_no_weights.txt', 'ICBP', '/home/bjarni/PCMA/faststorage/2_RESULTS/ICBP')
ar.count_ld_indep_regions('/home/bjarni/PCMA/faststorage/1_DATA/ICBP_zs.txt', '/home/bjarni/PCMA/faststorage/2_RESULTS/PCMA_ICBP_t1.0.txt', ld_reg_map = '/project/PCMA/faststorage/1_DATA/fourier_ls.hdf5')
