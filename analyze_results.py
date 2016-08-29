"""
Analyze results..
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# import matplotlib

import pandas
import scipy as sp
import gzip
import pylab
import itertools as it

import h5py

def get_sid_pos_map(sids):
    sids = set(sids)
    sid_map = {}
    chrom_pos_dict = {}
    for chrom_i in range(1,23):
        chrom_pos_dict[chrom_i] = [] 
    for chrom_i in range(1,23):
        fn = '/project/PCMA/faststorage/1_DATA/1k_genomes/ALL_1000G_phase1integrated_v3_chr%d_impute.legend.gz'%chrom_i
        with gzip.open(fn) as f:
            f.next()
            for line in f:
                l = line.split()
                sid = l[0]
                if sid in sids:
                    pos = int(l[1])
                    sid_map[l[0]]={'pos':pos, 'chrom':chrom_i}
                    chrom_pos_dict[chrom_i].append(pos)
    return {'sid_map':sid_map, 'chrom_pos_dict':chrom_pos_dict}




def plot_manhattan(result_file,fig_filename='/project/PCMA/faststorage/2_RESULTS/figures/manhattan_all.png'):
    """
    Generates a Manhattan plot for the PCMA results...
    """
    res = pandas.read_csv(result_file, delim_whitespace=True)
    sids = list(res.SID)
    print 'Getting SNP positions from 1K Genomes data'
    d = get_sid_pos_map(sids)
    sid_map = d['sid_map']
    chrom_pos_dict = d['chrom_pos_dict']
    print 'Calculating X-axis offsets'
    chrom_offset_dict = {}
    x_tick_pos = []
    x_tick_lab = []
    x_offset = 0
    for chrom_i in range(1,23):
        chrom_offset_dict[chrom_i]=x_offset
        old_x_offset = x_offset
        x_offset += max(chrom_pos_dict[chrom_i])
        x_tick_pos.append((old_x_offset+x_offset)/2.0)
        x_tick_lab.append(str(chrom_i))
    
    print 'Calculating X-axis positions'
    ps = sp.array(res.pvCHI2)
#         ps = sp.array(res.pvCPC)
    x_positions=sp.empty(len(ps))
    chromosomes=sp.empty(len(ps))
    for i, sid in enumerate(sids):
        chrom_i=sid_map[sid]['chrom']
        pos=sid_map[sid]['pos']
        x_offset = chrom_offset_dict[chrom_i]
        x_positions[i]=x_offset+pos
        chromosomes[i]=chrom_i
    
    neg_log_ps = -sp.log10(ps)
    ps_filter = neg_log_ps>3
    filtered_log_ps = neg_log_ps[ps_filter]   
    filtered_pos = x_positions[ps_filter] 
    filtered_chroms = chromosomes[ps_filter]
    
    color_map = {1:{'x_pos':[],'ps':[]}, 2:{'x_pos':[],'ps':[]},
                 3:{'x_pos':[],'ps':[]}, 4:{'x_pos':[],'ps':[]}}
    for lps,pos,chrom in it.izip(filtered_log_ps,filtered_pos,filtered_chroms):
        if chrom%2==0:
            if lps<7.301029:
                color_map[1]['x_pos'].append(pos)
                color_map[1]['ps'].append(lps)
            else:
                color_map[3]['x_pos'].append(pos)
                color_map[3]['ps'].append(lps)
        else:
            if lps<7.301029:
                color_map[2]['x_pos'].append(pos)
                color_map[2]['ps'].append(lps)
            else:
                color_map[4]['x_pos'].append(pos)
                color_map[4]['ps'].append(lps)
        
    print 'Filtering and plotting'
    with pylab.style.context('fivethirtyeight'):
        pylab.figure(figsize=(14,5))
        pylab.plot(color_map[1]['x_pos'],color_map[1]['ps'],'.',color='#1199EE',alpha=0.2)
        pylab.plot(color_map[2]['x_pos'],color_map[2]['ps'],'.',color='#11BB00',alpha=0.2)
        pylab.plot(color_map[3]['x_pos'],color_map[3]['ps'],'.',color='#AA99EE',alpha=0.7)
        pylab.plot(color_map[4]['x_pos'],color_map[4]['ps'],'.',color='#AABB00',alpha=0.7)
        pylab.ylabel('-log(P-value)')
        pylab.xlabel('Chromosomes')
        pylab.xticks(x_tick_pos,x_tick_lab)
        pylab.tight_layout()
        pylab.savefig(fig_filename)
    pylab.clf()


def get_log_quantiles(scores, num_dots=1000, max_val=8):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = sp.array(10 ** ((-sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)) * max_val) \
                * len(scores), dtype='int')
    return -sp.log10(scores[indices])



def _log_qqplot_(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            max_val=5, title=None, text=None, plot_label=None, ax=None, **kwargs):
    storeFig = False
    if ax is None:
        f = pylab.figure(figsize=(5.4, 5))
        storeFig = True
    pylab.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            pylab.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            pylab.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    pylab.ylabel("Observed $-log_{10}(P$-value$)$")
    pylab.xlabel("Expected $-log_{10}(P$-value$)$")
    if title:
        pylab.title(title)
    max_x = max_val
    max_y = max(map(max, quantiles_list))
    pylab.axis([-0.025 * max_x, 1.025 * max_x, -0.025 * max_y, 1.025 * max_y])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        pylab.legend(loc=2, numpoints=2, markerscale=1, prop=fontProp)
    y_min, y_max = pylab.ylim()
    if text:
        f.text(0.05 * max_val, y_max * 0.9, text)
    if plot_label:
        f.text(-0.138 * max_val, y_max * 1.01, plot_label, fontsize=14)
    pylab.tight_layout()
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')


def get_quantiles(scores, num_dots=1000):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = [int(len(scores) * i / (num_dots + 2)) for i in range(1, num_dots + 1)]
    return scores[indices]



def _qqplot_(quantiles_list, png_file=None, pdf_file=None, quantile_labels=None, line_colors=None,
            title=None, text=None, ax=None, plot_label=None, **kwargs):
    storeFig = False
    if ax is None:
        f = pylab.figure(figsize=(5.4, 5))
        storeFig = True
    pylab.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=2.0)
    num_dots = len(quantiles_list[0])
    exp_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)
    for i, quantiles in enumerate(quantiles_list):
        if line_colors:
            c = line_colors[i]
        else:
            c = 'b'
        if quantile_labels:
            pylab.plot(exp_quantiles, quantiles, label=quantile_labels[i], c=c, alpha=0.5, linewidth=2.2)
        else:
            pylab.plot(exp_quantiles, quantiles, c=c, alpha=0.5, linewidth=2.2)
    pylab.ylabel("Observed $P$-value")
    pylab.xlabel("Expected $P$-value")
    if title:
        pylab.title(title)
    pylab.axis([-0.025, 1.025, -0.025, 1.025])
    if quantile_labels:
        fontProp = matplotlib.font_manager.FontProperties(size=10)
        pylab.legend(loc=2, numpoints=2, markerscale=1, prop=fontProp)
    if text:
        f.text(0.05, 0.9, text)
    if plot_label:
        f.text(-0.151, 1.04, plot_label, fontsize=14)
    pylab.tight_layout()
    if storeFig == False:
        return
    if png_file != None:
        f.savefig(png_file)
    if pdf_file != None:
        f.savefig(pdf_file, format='pdf')

def plot_QQ_plots(result_file,png_file_prefix='/Users/bjarnivilhjalmsson/data/tmp/test', num_dots=1000, max_neg_log_val=7,
                  title=''):
    """
Generates a QQ plot of the PCMA results..
    """
    res = pandas.read_table(result_file)
    mvt_ps = sp.array(res.pval)
    combPC_ps = sp.array(res.combPC)
    pvals_list = [mvt_ps, combPC_ps]
    line_colors=['c','m']
    result_labels=['MVT', 'comb. PC']
    qs = []
    log_qs = []
    for pvals in pvals_list:
        qs.append(get_quantiles(pvals, num_dots))
        log_qs.append(get_log_quantiles(pvals, num_dots, max_neg_log_val))
    _qqplot_(qs, png_file_prefix + '_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title)
    _log_qqplot_(log_qs, png_file_prefix + '_log_qq.png', quantile_labels=result_labels,
                line_colors=line_colors, num_dots=num_dots, title=title, max_val=max_neg_log_val)




def plot_overlap_ps(result_file, ss_file='/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt', 
                   fig_filename='/Users/bjarnivilhjalmsson/data/tmp/manhattan_combPC_HGT.png', method='combPC', 
                   ylabel='Comb. PC (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='Height $-log_{10}(P$-value$)$', p_thres = 0.00001):
    #Parse results ans SS file
    res_table = pandas.read_table(result_file)
    ss_table = pandas.read_table(ss_file)
    #Parse 
    res_sids = sp.array(res_table['SNPid'])
    if method=='MVT':
        comb_ps = sp.array(res_table['pval'])
    elif method=='combPC':
        comb_ps = sp.array(res_table['combPC'])
    if 'MarkerName' in ss_table.keys():
        ss_sids = sp.array(ss_table['MarkerName'])
    elif 'SNP' in ss_table.keys():
        ss_sids = sp.array(ss_table['SNP'])
    else:
        raise Exception("Don't know where to look for rs IDs")
    marg_ps = sp.array(ss_table['p'])
    
    # Filtering boring p-values
    res_p_filter = comb_ps<p_thres
    res_sids = res_sids[res_p_filter]
    comb_ps = comb_ps[res_p_filter]
#     ss_p_filter = marg_ps<p_thres
#     ss_sids = ss_sids[ss_p_filter]
#     marg_ps = marg_ps[ss_p_filter]
    
    common_sids = sp.intersect1d(res_sids, ss_sids)
    print 'Found %d SNPs in common'%(len(common_sids))
    ss_filter = sp.in1d(ss_sids, common_sids)
    res_filter = sp.in1d(res_sids, common_sids)
    
    ss_sids = ss_sids[ss_filter]
    res_sids = res_sids[res_filter]
    marg_ps = marg_ps[ss_filter]
    comb_ps = comb_ps[res_filter]
    
    print 'Now sorting'
    ss_index = sp.argsort(ss_sids)
    res_index = sp.argsort(res_sids)
    
    marg_ps=-sp.log10(marg_ps[ss_index])
    comb_ps=-sp.log10(comb_ps[res_index])
    
    with pylab.style.context('fivethirtyeight'):
        pylab.plot(marg_ps,comb_ps,'b.',alpha=0.2)
        (x_min,x_max) = pylab.xlim()
        (y_min,y_max) = pylab.ylim()
        
        pylab.plot([x_min,x_max],[x_min,x_max],'k--',alpha=0.2)
        pylab.ylabel(ylabel)
        pylab.xlabel(xlabel)
        pylab.tight_layout()
        pylab.savefig(fig_filename)
    pylab.clf()
        
        
def plot_overlap_ps(result_file, ss_files=['/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt',
                                           '/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt',
                                           '/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt',
                                           '/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt'], 
                   fig_filename='/Users/bjarnivilhjalmsson/data/tmp/manhattan_combPC_HGT.png', method='combPC', 
                   ylabel='Comb. PC (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='Height $-log_{10}(P$-value$)$', p_thres = 0.00001):
    #Parse results ans SS file
    res_table = pandas.read_table(result_file)
    marg_ps_list = []
    #Parse 
    res_sids = sp.array(res_table['SNPid'])
    if method=='MVT':
        comb_ps = sp.array(res_table['pval'])
    elif method=='combPC':
        comb_ps = sp.array(res_table['combPC'])
    for ss_file in ss_files:
        ss_table = pandas.read_table(ss_file)
        if 'MarkerName' in ss_table.keys():
            ss_sids = sp.array(ss_table['MarkerName'])
        elif 'SNP' in ss_table.keys():
            ss_sids = sp.array(ss_table['SNP'])
        else:
            raise Exception("Don't know where to look for rs IDs")
        marg_ps = sp.array(ss_table['p'])
        marg_ps_list.append(marg_ps)
    min_pvals = sp.minimum(sp.array(marg_ps_list))
    # Filtering boring p-values
    res_p_filter = comb_ps<p_thres
    res_sids = res_sids[res_p_filter]
    comb_ps = comb_ps[res_p_filter]
#     ss_p_filter = marg_ps<p_thres
#     ss_sids = ss_sids[ss_p_filter]
#     marg_ps = marg_ps[ss_p_filter]
    
    common_sids = sp.intersect1d(res_sids, ss_sids)
    print 'Found %d SNPs in common'%(len(common_sids))
    ss_filter = sp.in1d(ss_sids, common_sids)
    res_filter = sp.in1d(res_sids, common_sids)
    
    ss_sids = ss_sids[ss_filter]
    res_sids = res_sids[res_filter]
    marg_ps = marg_ps[ss_filter]
    comb_ps = comb_ps[res_filter]
    
    print 'Now sorting'
    ss_index = sp.argsort(ss_sids)
    res_index = sp.argsort(res_sids)
    
    marg_ps=-sp.log10(marg_ps[ss_index])
    comb_ps=-sp.log10(comb_ps[res_index])
    
    with pylab.style.context('fivethirtyeight'):
        pylab.plot(marg_ps,comb_ps,'b.',alpha=0.2)
        (x_min,x_max) = pylab.xlim()
        (y_min,y_max) = pylab.ylim()
        
        pylab.plot([x_min,x_max],[x_min,x_max],'k--',alpha=0.2)
        pylab.ylabel(ylabel)
        pylab.xlabel(xlabel)
        pylab.tight_layout()
        pylab.savefig(fig_filename)
    pylab.clf()
        

def parse_ldetect_map(file_prefix= '/Users/bjarnivilhjalmsson/REPOS/others/ldetect-data/EUR/fourier_ls-', 
                      out_hdf5_file='/Users/bjarnivilhjalmsson/REPOS/others/ldetect-data/EUR/fourier_ls.hdf5'):
    of = h5py.File(out_hdf5_file)
    for chrom in range(1,23):
        chrom_str = 'chr%d'%chrom 
        map_file = '%s%s.bed'%(file_prefix,chrom_str)
        bin_limits = []
        with open(map_file,'r') as f:
            print f.next()
            for line in f:
                l = line.split()
                bin_limits.append(int(l[1]))
            bin_limits.append(int(l[2]))
        of.create_dataset(chrom_str, data=sp.array(bin_limits))
    of.close()
        
        
def parse_PCMA_results(ss_ps_file, res_file):
# def parse_PCMA_results(ss_ps_file, ss_zs_file, ss_wt_file, res_file):
    #Parse ss file, get various information
    print 'Starting to load data...'
    ss_ps_df = pandas.read_table(ss_ps_file,delim_whitespace=True)
    print 'Parsed p-value file'
#     ss_zs_df = pandas.read_table(ss_zs_file)
#     ss_wts_df = pandas.read_table(ss_wt_file)
#     ss_zs_ids = list(ss_ps_df.columns)[6:]
#     ss_wts_ids = list(ss_ps_df.columns)[6:]
    ss_ps_ids = list(ss_ps_df.columns)[6:]
    num_ss = len(ss_ps_ids)
    res_df = pandas.read_csv(res_file)
    print 'Parsed results file'
    pc_ids = ['pvPC%d'%i for i in range(1,num_ss+1)]
    
    dupl_columns = ss_ps_df.columns[1:6]
    
    #Partition things by chromosome
    chrom_res_dict = {}
    for chrom in range(1,23):
        print 'Working on chromosome %d'%chrom
        res_chrom_df = res_df.loc[res_df['CHR']==chrom]
        ss_ps_chrom_df = ss_ps_df.loc[ss_ps_df['CHR']==chrom]
#         ss_zs_chrom_df = ss_zs_df.loc[ss_zs_df['CHR']==chrom]
#         ss_wts_chrom_df = ss_wts_df.loc[ss_wts_df['CHR']==chrom]
        
        print 'Sub-sampled chromsomes, now merging'
        use_cols = ss_ps_chrom_df.columns - dupl_columns
        merged_df = res_chrom_df.merge(ss_ps_chrom_df[use_cols],on='SID')
        print 'Merge done'
#         use_cols = ss_zs_chrom_df.colums - dupl_columns
#         merged_df = merged_df.merge(ss_zs_chrom_df[use_cols],on='SID')
#         print 'Merge 2 done'
#         use_cols = ss_wts_chrom_df.colums - dupl_columns
#         merged_df = merged_df.merge(ss_wts_chrom_df[use_cols],on='SID')
#         print 'Merge 3 done'

        print merged_df.colums
        chrom_str = 'chr%d'%chrom
        marg_ps = merged_df[ss_zs_ids]
        min_marg_ps = marg_ps.min(1)
        chrom_res_dict[chrom_str] = {'ps':merged_df[ss_ps_ids],'positions':merged_df['POS'], 
                                     'sids':merged_df['SID'], 'maf':merged_df['MAF'],
                                     'marg_ps':marg_ps, 'min_marg_ps':min_marg_ps, 
                                     'comb_ps':merged_df['pvCHI2'], 'pc_ps':merged_df[pc_ids]}
#         chrom_res_dict[chrom_str] = {'zs':merged_df[ss_zs_ids],'weights':merged_df[ss_wts_ids], 'positions':merged_df['POS'], 
#                                      'sids':merged_df['SID'], 'maf':merged_df['MAF'], 'res_df':res_chrom_df, 
#                                      'marg_ps':marg_ps, 'min_marg_ps':min_marg_ps, 'comb_ps':merged_df['pvCHI2'], 'pc_ps':merged_df[pc_ids]}

    
    return chrom_res_dict


def count_ld_indep_regions(ss_file, res_file, ld_reg_map = '/project/PCMA/faststorage/1_DATA/fourier_ls.hdf5'):
    #parse results..
    print 'Parsing PCMA results'
    chrom_res_dict = parse_PCMA_results(ss_file,res_file)
    
    #Filter for good SNPs?
    
    #parse ldetect map
    print 'Loading ldetect map'
    ldr = h5py.File(ld_reg_map,'r')
    
    num_new_hits = 0
    num_comb_hits = 0
    num_marg_hits = 0
    num_shared_hits = 0
    num_missed_hits = 0
        
    chrom_bin_dict = {} 
    
    res_summary_dict = {}
    for chrom in range(1,23):
        print 'Working on chromosome %d'%chrom
        chrom_str = 'chr%d'%chrom
        res_dict = chrom_res_dict[chrom_str]
        chrom_bins = ldr[chrom_str]
        bin_indices = sp.digitize(res_dict['positions'], chrom_bins)
        chrom_bin_dict[chrom_str]={'bin_indices':bin_indices, 'chrom_bins':chrom_bins, 'num_bins':len(chrom_bins)-1}
        
        #Count things..
        print 'Counting hits'
        #assert len(chrom_bins)-1==bin_indices.max()+1, 'WTF?'
        for bin_i in range(bin_indices.max()+1):
            bin_filter = bin_indices==bin_i
            if sp.any(bin_filter):
                min_marg_pv = (res_dict['min_marg_ps'][bin_filter]).min()
                marg_hit = min_marg_pv<5E-8
                comb_ps = res_dict['comb_ps'][bin_filter]
                min_i = comb_ps.idxmin()  
                min_comb_pv = comb_ps[min_i]
                min_sid = res_dict['sids'][min_i]
                
                comb_hit =min_comb_pv<5E-8
    
                if marg_hit:
                    num_marg_hits+=1
                    if comb_hit:
                        num_shared_hits +=1
                        num_comb_hits +=1
                    else:
                        num_missed_hits+=1
                elif comb_hit:
                    num_new_hits+=1
                    num_comb_hits +=1
                    
                    start_pos = chrom_bins[bin_i]
                    if bin_i<len(chrom_bins)-1:
                        end_pos = chrom_bins[bin_i+1]
                    else:
                        end_pos = -1
                    res_summary_dict[bin_i]={'min_marg_pv':min_marg_pv, 'min_comb_pv':min_comb_pv, 
                                             'min_PC_pv': res_dict['pc_ps'].loc[min_i], 'min_sid':min_sid,
                                             'chromsome':chrom, 'positions_bin':(start_pos,end_pos)}
                    #More information on new hits somewhere
        
    print '\nResults summary: \n# new hits: %d \n# missed hits: %d \n# of shared hits: %d \n# multivar. hits: %d \n# marg. hits: %d \n'%(num_new_hits, num_missed_hits, num_shared_hits, num_comb_hits, num_marg_hits)
    print res_summary_dict
                
                
def run_all_ts(pruned_file, ss_file, name, out_prefix, ts=[0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4]):
    """  
    """
    import os

    with open(ss_file,'r') as f:
        header = f.next()
        l = header.split()
        ss_ids = l[1:]
    weights_fn = out_prefix+'_ss_weights.txt'
    with open(weights_fn,'w') as f:
        f.write('Study    Weight\n')
        for ss_id in ss_ids:
            f.write('%s    %d\n'%(ss_id,1))
    
    
    for t in ts:
        print 'Working on t=%0.2f'%t
        run_id = name+('_t%0.1f'%t)
        out_file = out_prefix+('_t%0.1f'%t)+'.out'
        command_str = '/home/bjarni/PCMA/0_PROGS/PCMA/Debug/PCMA -p %s -i %s -n %s -t %0.1f -w %s --v --f > %s'%(pruned_file, ss_file, run_id, t, weights_fn, out_file)
        print command_str
        os.system(command_str)
        command_str = 'mv PCMA_%s.txt /home/bjarni/PCMA/faststorage/2_RESULTS/'%run_id
        print command_str
        os.system(command_str)

def run_t(pruned_file, ss_file, name, out_prefix, t=1):
    """  
    """
    import os

    with open(ss_file,'r') as f:
        header = f.next()
        l = header.split()
        ss_ids = l[1:]
    weights_fn = out_prefix+'_ss_weights.txt'
    with open(weights_fn,'w') as f:
        f.write('Study    Weight\n')
        for ss_id in ss_ids:
            f.write('%s    %d\n'%(ss_id,1))
    
    
    print 'Working on t=%0.2f'%t
    run_id = name+('_t%0.1f'%t)
    out_file = out_prefix+('_t%0.1f'%t)+'.out'
    command_str = '/home/bjarni/PCMA/0_PROGS/PCMA/Debug/PCMA -p %s -i %s -n %s -t %0.1f -w %s --v --f > %s'%(pruned_file, ss_file, run_id, t, weights_fn, out_file)
    print command_str
    os.system(command_str)
    command_str = 'mv PCMA_%s.txt /home/bjarni/PCMA/faststorage/2_RESULTS/'%run_id
    print command_str
    os.system(command_str)


def parse_corr_matrices(ss_file, res_prefix, ts=[0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4]):
    """
    """
    import os
    with open(ss_file,'r') as f:
        header = f.next()
        l = header.split()
        ss_ids = l[1:]
    
    print ss_ids
    num_ss = len(ss_ids)
    
    res_dict = {}
    
    for t in ts:
        zz_corr_mat = sp.empty((num_ss,num_ss))
        out_file = res_prefix+('_t%0.1f'%t)+'.out'
        try:
            with open(out_file,'r') as f:
                while not (f.next()).startswith(' *** Zscore correlation matrix:'):
                    pass
                for i in range(num_ss):
                    line = f.next()
                    l = line.split()
                    zz_corr_mat[i]=sp.array(map(float,l))
            res_dict['t%0.1f'%t] = zz_corr_mat
        except Exception:
            pass
    print res_dict
    return res_dict

def plot_corr_mat_convergence(corr_mat_dict, title_str=''):
    
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    # set ticks and tick labels
    ax.set_xlim((0, 2.4))
    ax.set_xticks([0.0, 0.8,1.6,2.4])
    ax.set_xticklabels(['0', '0.8', '1.6','2.4'])
    ax.set_ylim((-1.3, 1.3))
    ax.set_yticks([-1, 0, 1])
    
    ax.spines['left'].set_bounds(-1, 1)
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    
    t_strs = corr_mat_dict.keys()
    t_strs.sort()
    corr_mat = corr_mat_dict[t_strs[0]]
    num_ss = len(corr_mat)
    res_list = [[] for i in range((num_ss)*(num_ss-1)/2)]
    ts = []
    for t_str in t_strs:
        t = float(t_str[1:])
        corr_mat = corr_mat_dict[t_str]
        l_i = 0
        for i in range(num_ss):
            for j in range(i):
                v = corr_mat[i,j]
                res_list[l_i].append(v)
                l_i +=1
        ts.append(t)
        
    for l in res_list:
        plt.plot(ts,l,alpha=0.3,color='dodgerblue', linewidth=2)
    plt.ylabel('Correlation estimate')
    plt.xlabel('Z-score threshold')
    
    plt.title(title_str)
    
    # Finally, save the figure as a PNG.
    plt.savefig('test.png', bbox_inches='tight')    






if __name__=='__main__':
#     plot_manhattan('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt',fig_filename='/Users/bjarnivilhjalmsson/data/tmp/manhattan_combPC.png',method='combPC')
#     plot_manhattan('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt',fig_filename='/Users/bjarnivilhjalmsson/data/tmp/manhattan_MVT.png',method='MVT')
#     plot_overlap_ps('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt', ss_file='/Users/bjarnivilhjalmsson/data/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt', 
#                    fig_filename='/Users/bjarnivilhjalmsson/data/tmp/ps_MVT_HGT.png', method='MVT', 
#                    ylabel='MVT (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='Height $-log_{10}(P$-value$)$')
#     plot_overlap_ps('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt', ss_file='/Users/bjarnivilhjalmsson/data/GIANT/GIANT_2015_HIP_COMBINED_EUR.txt', 
#                    fig_filename='/Users/bjarnivilhjalmsson/data/tmp/ps_MVT_HIP.png', method='MVT', 
#                    ylabel='MVT (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='HIP $-log_{10}(P$-value$)$')
#     plot_overlap_ps('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt', ss_file='/Users/bjarnivilhjalmsson/data/GIANT/SNP_gwas_mc_merge_nogc.tbl.uniq', 
#                    fig_filename='/Users/bjarnivilhjalmsson/data/tmp/ps_MVT_BMI.png', method='MVT', 
#                    ylabel='MVT (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='BMI $-log_{10}(P$-value$)$')
#     plot_overlap_ps('/Users/bjarnivilhjalmsson/REPOS/pcma/Debug/PCMA_test.txt', ss_file='/Users/bjarnivilhjalmsson/data/GIANT/GIANT_2015_WC_COMBINED_EUR.txt', 
#                    fig_filename='/Users/bjarnivilhjalmsson/data/tmp/ps_MVT_WC.png', method='MVT', 
#                    ylabel='MVT (HIP,WC,HGT,BMI) $-log_{10}(P$-value$)$', xlabel='WC $-log_{10}(P$-value$)$')

    run_all_ts('/home/bjarni/PCMA/faststorage/1_DATA/IMMUNE_REL_4_zs_ldprunedno_weights.txt','/home/bjarni/PCMA/faststorage/1_DATA/IMMUNE_REL_4_zsno_weights.txt', 'IMMUNE_REL4', '/home/bjarni/PCMA/faststorage/2_RESULTS/IMMUNE_REL4')
    count_ld_indep_regions('/home/bjarni/PCMA/faststorage/1_DATA/IMMUNE_REL_4_zs.txt', '/home/bjarni/PCMA/faststorage/2_RESULTS/PCMA_IMMUNE_REL4_t0.8.txt', ld_reg_map = '/project/PCMA/faststorage/1_DATA/fourier_ls.hdf5')
