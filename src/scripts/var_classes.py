"""funcs to classifiy variant status."""
import pandas, csv

def write_pandas_header(in_file, out_file, new_fields, rm_fields, delimiter):
    with open(in_file) as f, open(out_file, 'w') as fout:
        fields = f.readline().strip().split(delimiter) + new_fields
        fields_use = [field for field in fields if not field in rm_fields]
        print(delimiter.join(fields_use), file=fout)

def flag_pass_1per_and_exac_filter(df):
    oneKgs = ['af_1kg_%s' % (x,) for x in ('afr', 'amr', 'eas', 'eur', 'sas')]
    crit = df.apply(lambda row: max([row[x] for x in oneKgs])<.001
                    and ( str(row['ExACfilterMatchRefAlt']) == 'nan'
                          or str(row['ExACfilterMatchRefAlt']) == 'PASS'), axis=1)
    return crit

def apply_point_one_percent(in_file, out_file):
    """This is just a filter."""
    delimiter = '\t'
    chunk_size = 5000
    write_pandas_header(in_file, out_file, [], [], delimiter)
    for df in pandas.read_csv(in_file, delimiter=delimiter, chunksize=chunk_size):
        crit = flag_pass_1per_and_exac_filter(df)
        df[crit].to_csv(out_file, index=False, sep=delimiter, header=False, mode='a')

def check_tumor(in_file, out_file):
    delimitelr = '\t'
    chunk_size = 5000
    write_pandas_header(in_file, out_file, [], [], delimiter)

# dbnsfp_radialsvm_pred has the full list. This wasn't broken apart by vt
# dbnsfp_radialsvm summarizes into one call
# same w/ dbnsfp_lr 
def is_lp(row):
    """Picks up a lof of splice_region_variant with MED impact."""
    c = -1
    if isinstance(row['cadd_phred'], str):
        if row['cadd_phred'].strip() and row['cadd_phred'].strip() != 'None':
            c = float(row['cadd_phred'])
    else:
        c = row['cadd_phred']
    return ( (row['dbnsfp_lr'] == 'D' or row['dbnsfp_radialsvm'] == 'D')
             and c>=20.0) or ((row['eff_indel_splice'] == 1 or row['is_splicing'] == 1)
                              and row['impact_severity'] in ('HIGH','MED'))

#lall is_lof are HIGH
def is_p(row):
    """Picks up splice_donor_variant and splice_acceptor_variant not marked as lof."""
    h = row['hgmd_phen'].strip()
    has_hgmd = h != 'None' and h != ''
    return has_hgmd or row['impact_severity'] == 'HIGH'

def apply_p_or_lp(in_file, out_file, p_or_lp_func, p_or_lp_label):
    delimiter = '\t'
    chunk_size = 5000
    write_pandas_header(in_file, out_file, [p_or_lp_label], [], delimiter)
    for df in pandas.read_csv(in_file, delimiter=delimiter, chunksize=chunk_size):
        df.loc[:, p_or_lp_label] = df.apply(p_or_lp_func, axis=1)
        df.to_csv(out_file, index=False, sep=delimiter, header=False, mode='a')

# def merge_p_and_lp(p_file, lp_file, out_file):
#     with open(p_file) as p, open(lp_file) as lp, open(out_file, 'w') as fout:
#         header = p_file.readline().strip()
#         print(header + ['lp']

def fix_exac_cov(row):
    if isinstance(row['totexaccov_median'], float):
        return float(row['totexaccov_median'])
    if isinstance(row['totexaccov_median'], str):
        if row['totexaccov_median'] != 'None':
            return max( [float(x) for x in row['totexaccov_median'].split(',')] )
        return 0

def pad_decimal(d):
    if len(d.split('.')[-1]) == 2:
        return d
    return d + '0'

def fix_cadd(row, col):
    """cadd_phred"""
    cf = row[col]
    if isinstance(cf, float):
        return pad_decimal('%.2f' % (cf,))
    if '.' in cf:
        return pad_decimal('%.2f' % (float(cf),))
    return cf

def fix_tot_exac_cov(row):
    cf = row['totexaccov_10']
    if isinstance(cf, float):
        return pad_decimal('%.2f' % (cf,))
    if '.' in cf:
        if ',' in cf:
            m = max( [float(x) for x in cf.split(',')] )
            return pad_decimal('%.2f' % (m, ))
        else:
            return pad_decimal('%.2f' % (float(cf), ))
    return cf

def flag_vqsr(vqsr_filter, var_type):
    """For cap and wxs, is vqsr passed?"""
    if var_type == 'SNV':
        return str('None' == vqsr_filter)
    else:
        return str(vqsr_filter in ('VQSRTrancheINDEL99.90to100.00',
                                   'None'))

var_class_ignore_impact = ('intron_variant', 'upstream_gene_variant',
                           'downstream_gene_variant', 'synonymous_variant',
                           'intergenic_variant', 'intragenic_variant',
                           'intergenic_region')

def clean_cgi(in_file, out_file):
    with open(in_file) as f, open(out_file, 'w') as fout:
        reader = csv.DictReader(f, delimiter='\t')
        fields = reader.fieldnames
        print('\t'.join(fields), file=fout)
        for row in reader:
            if not row['impact'] in var_class_ignore_impact:
                ls = [row[x] for x in fields]
                print('\t'.join(ls), file=fout)

def getSrrToTarget(srrFile, popFile):
    cols = ['Run_s', 'Sample_Name_s', 'histological_type_s']
    sampleDf = pandas.read_csv(srrFile, delimiter='\t')[cols]
    popDf = pandas.read_csv(popFile, delimiter='\t')[ ['Run_s_pre', 'pop'] ]
    popDf.loc[:, 'Run_s'] = popDf.apply(lambda row: row['Run_s_pre'].split('/')[0], axis=1)
    sampleDf.loc[:, "perryHist"] = sampleDf.apply(lambda row: 'WT' if 'Kidney' in row['histological_type_s'] else row['histological_type_s'], axis=1)
    df = pandas.merge(sampleDf, popDf, on='Run_s', how='left').fillna('??')
    newCols = ['Run_s', 'Sample_Name_s', 'perryHist', 'pop']
    srrToId = { srr:'_'.join( (id,hist,pop) ) for srr,id,hist,pop in
                list(df[newCols].values) }
#    print('debug', srrToId['SRR1767937'])
    return srrToId

# Normal_SRR331718_PAMDAL
def transpose_pindel(in_file, sample_file, ancestry_file,
                     out_file, var_type, platform_label):
    """for wxs and capture"""
    srrToTgt = getSrrToTarget(sample_file, ancestry_file)
    rm_cols = ['vqslod', 'culprit', 'variant_id', 'totexaccov_median', 'cadd_phred',
               'totexaccov_10', 'dels', 'cadd_raw', 'vcf_id', 'qual.1', 'qual',
               'ru', 'set', 'excesshet', 'end']
    new_cols = ['sample', 'histotype', 'pop', 'run_num', 'TorN', 'refReads',
                'altReads', 'zygosity', 'depth',
                'totexaccov_median_use', 'CADD_phred', 'CADD_raw',
                'tot_exac_cov_frac_10', 'st_bed', 'end_bed', 'platform',
                'pass_vqsr', 'in_tumor', 'tumor_ref_reads', 'tumor_alt_reads']
    infoLs = ('gt_ref_depths', 'gt_alt_depths', 'gt_types', 'gt_depths')
    with open(in_file) as f, open(out_file, 'w') as fout:
        # header
        s = f.readline().strip().split('\t')
        init_cols = s
        f2i = {field:idx for idx,field in enumerate(s)}
        fields = f2i

        samples = set( ['.'.join( _.split('.')[1:] ) for _ in fields
                        if 'Normal' in _ or 'Tumor' in _] )
        normal_samples = [n for n in samples if 'Normal' in n]
        tumor_samples = [n for n in samples if 'Tumor' in n]
        df_normal = pandas.DataFrame({'normal_sample_name':normal_samples,
                                      'tgt':[x.split('_')[-1]
                                             for x in normal_samples]})
        df_tumor = pandas.DataFrame({'tumor_sample_name':tumor_samples,
                                     'tgt':[x.split('_')[-1]
                                            for x in tumor_samples]})
        df_names = pandas.merge(df_normal, df_tumor, on='tgt', how='outer')
        tmp_names = ['normal_sample_name', 'tumor_sample_name', 'tgt']
        normal_to_tumor = {n:t for n,t,_ in df_names[tmp_names].values}
                                     
        head = new_cols + [_ for _ in s if not 'Normal' in _
                           and not 'Tumor' in _ and not _ in rm_cols]
        print('\t'.join(head), file=fout)
        for line in f:
            sp = line.strip().split('\t')
            row = {}
            for field in f2i:
                row[field] = sp[ f2i[field] ]
            row['st_bed'] = row['start']
            row['end_bed'] = row['end']
            # correct position for hg19 checks
            row['start'] = str( int(row['start']) + 1 )
            for sample_pre in samples:
                row['gt_depths.' + sample_pre] = sum([int(row[x + sample_pre])
                                                      for x in ['gt_ref_depths.',
                                                                'gt_alt_depths.']])
                # genotype is not ref 0 or unknown 2
                if (
                    int( float(row['gt_types.' + sample_pre]) ) not in (1, 2)
                    and int(row['gt_depths.' + sample_pre]) >= 10
                    and not row['impact'] in var_class_ignore_impact
                    ):
                    #Normal_SRR331718_PAMDAL
                    t_or_n, srr, simple_tgt = sample_pre.split('_')
                    row['in_tumor'] = 'NA'
                    row['tumor_ref_reads'] = 'NA'
                    row['tumor_alt_reads'] = 'NA'
                    if sample_pre in normal_samples:
                        tumor_sample = normal_to_tumor[sample_pre]
                        row['tumor_ref_reads'] = row['gt_ref_depths.' + tumor_sample]
                        row['tumor_alt_reads'] = row['gt_alt_depths.' + tumor_sample]
                        row['in_tumor'] = int(float(row['gt_types.' + tumor_sample])) not in (1,2)
                    sample = srrToTgt[srr]
                    sampleName = sample.split('_')[0]
                    row['sample'] = sampleName
                    row['TorN'] = t_or_n[0]
                    row['run_num'] = srr
                    _, row['histotype'], row['pop'] = sample.split('-')[-1].split('_')
                    (row['refReads'],
                     row['altReads'], row['zygosity'],
                     row['depth']) = [ row[_ + '.' + sample_pre] for _ in infoLs ]
                    row['totexaccov_median_use'] = str(fix_exac_cov(row))
                    row['CADD_phred'] = str(fix_cadd(row, 'cadd_phred'))
                    row['CADD_raw'] = str(fix_cadd(row, 'cadd_raw'))
                    row['tot_exac_cov_frac_10'] = fix_tot_exac_cov(row)
                    row['platform'] = platform_label
                    row['pass_vqsr'] = flag_vqsr(str(row['filter']), var_type)
                    ls = [str(row[x]) for x in head]
                    print('\t'.join(ls), file=fout)

def transpose(in_file, out_file, var_type, platform_label):
    """for wxs and capture"""
    rm_cols = ['vqslod', 'culprit', 'variant_id', 'totexaccov_median', 'cadd_phred',
               'totexaccov_10', 'dels', 'cadd_raw', 'vcf_id', 'qual.1',
               'ru', 'set', 'excesshet', 'end']
    new_cols = ['sample', 'histotype', 'pop', 'run_num', 'qual', 'refReads',
                'altReads', 'zygosity', 'depth',
                'totexaccov_median_use', 'CADD_phred', 'CADD_raw',
                'tot_exac_cov_frac_10', 'st_bed', 'end_bed', 'platform',
                'pass_vqsr']
    infoLs = ('gt_quals', 'gt_ref_depths', 'gt_alt_depths', 'gt_types', 'gt_depths')
    with open(in_file) as f, open(out_file, 'w') as fout:
        s = f.readline().strip().split('\t')
        f2i = {field:idx for idx,field in enumerate(s)}
        fields = f2i
        samples = set( ['.'.join( _.split('.')[1:] ) for _ in fields if 'TARGET' in _] )
        print(samples)
        print(len(samples))
        head = new_cols + [_ for _ in fields if not 'TARGET' in _ and not _ in rm_cols]
        print('\t'.join(head), file=fout)
        for line in f:
            sp = line.strip().split('\t')
            row = {}
            for field in f2i:
                row[field] = sp[ f2i[field] ]
            row['st_bed'] = row['start']
            row['end_bed'] = row['end']
            # correct position for hg19 checks
            row['start'] = str( int(row['start']) + 1 )
            for sample in samples:
                if (
                    int( float(row['gt_quals.' + sample]) ) != -1
                    and int(row['gt_depths.' + sample]) >= 10
                    and not row['impact'] in var_class_ignore_impact
                    ):
                    sampleName = '-'.join(sample.split('_')[:5])
                    row['sample'] = sampleName
                    row['histotype'], row['pop'], row['run_num'] = sample.split('_')[5:]
                    (row['qual'], row['refReads'],
                     row['altReads'], row['zygosity'],
                     row['depth']) = [ row[_ + '.' + sample] for _ in infoLs ]
                    row['totexaccov_median_use'] = str(fix_exac_cov(row))
                    row['CADD_phred'] = str(fix_cadd(row, 'cadd_phred'))
                    row['CADD_raw'] = str(fix_cadd(row, 'cadd_raw'))
                    row['tot_exac_cov_frac_10'] = fix_tot_exac_cov(row)
                    row['platform'] = platform_label
                    row['pass_vqsr'] = flag_vqsr(str(row['filter']), var_type)
                    ls = [row[x] for x in head]
                    print('\t'.join(ls), file=fout)

def transpose_capture(in_file, out_file, var_type):
    """for capture. some indel fields are too big for csv"""
    rm_cols = ['vqslod', 'culprit', 'variant_id', 'totexaccov_median', 'cadd_phred',
               'totexaccov_10', 'dels', 'cadd_raw', 'vcf_id', 'qual.1',
               'ru', 'set', 'excesshet', 'end']
    new_cols = ['sample', 'histotype', 'pop', 'TorN', 'qual', 'refReads',
                'altReads', 'zygosity', 'depth',
                'totexaccov_median_use', 'CADD_phred', 'CADD_raw',
                'tot_exac_cov_frac_10', 'st_bed', 'end_bed', 'platform',
                'pass_vqsr']
    infoLs = ('gt_quals', 'gt_ref_depths', 'gt_alt_depths', 'gt_types', 'gt_depths')
    with open(in_file) as f, open(out_file, 'w') as fout:
        s = f.readline().strip().split('\t')
        f2i = {field:idx for idx,field in enumerate(s)}
        fields = f2i
        samples = set( ['.'.join( _.split('.')[1:] ) for _ in fields if 'TARGET' in _] )
        head = new_cols + [_ for _ in fields if not 'TARGET' in _ and not _ in rm_cols]
        print('\t'.join(head), file=fout)
        for line in f:
            sp = line.strip().split('\t')
            row = {}
            for field in f2i:
                row[field] = sp[ f2i[field] ]

            row['st_bed'] = row['start']
            row['end_bed'] = row['end']
            # correct position for hg19 checks
            row['start'] = str( int(row['start']) + 1 )
            for sample in samples:
                if (
                    int( float(row['gt_quals.' + sample]) ) != -1
                    and int(row['gt_depths.' + sample]) >= 10
                    and not row['impact'] in var_class_ignore_impact
                    ):
                    sampleName = '-'.join(sample.split('_')[:5])
                    row['sample'] = sampleName
                    row['histotype'], row['TorN'] = sample.split('_')[5:]
                    row['pop'] = '???'
                    (row['qual'], row['refReads'],
                     row['altReads'], row['zygosity'],
                     row['depth']) = [ row[_ + '.' + sample] for _ in infoLs ]
                    row['totexaccov_median_use'] = str(fix_exac_cov(row))
                    row['CADD_phred'] = str(fix_cadd(row, 'cadd_phred'))
                    row['CADD_raw'] = str(fix_cadd(row, 'cadd_raw'))
                    row['tot_exac_cov_frac_10'] = fix_tot_exac_cov(row)
                    row['platform'] = 'CAP'
                    row['pass_vqsr'] = flag_vqsr(str(row['filter']), var_type)
                    ls = [row[x] for x in head]
                    print('\t'.join(ls), file=fout)
                    
def is_p_fix(row):
    h = row['hgmd_phen'].strip()
    has_hgmd = h != 'None' and h != '' and 'DM' in str(row['hgmd_class_ls']).split(';')
    return has_hgmd or row['impact_severity'] == 'HIGH'

    # h = row['hgmd_phen'].strip()
    # has_hgmd = h != 'None' and h != '' and 'DM' in str(row['hgmd_class_ls']).split(';')
    # return has_hgmd or row['is_lof'] == 1 or ( (row['eff_indel_splice'] == '1' or
    #                                             row['is_splicing'] == '1') and row['impact_severity'] == 'HIGH' )

def is_p_fix_again(row):
    h = row['hgmd_phen'].strip()
    has_hgmd = h != 'None' and h != '' and 'DM' in str(row['hgmd_class_ls_use']).split(';')
    return has_hgmd or row['impact_severity'] == 'HIGH'

def is_lp_fix(row):
    """Cannot be p and lp"""
    if row['is_pathogenic']:
        return False
    return row['lp']

def is_lp_fix_again(row):
    """Cannot be p and lp"""
    if row['is_pathogenic_final']:
        return False
    return row['lp']

def update_p_and_lp(in_file, hgmd, out_file):
    delimiter = '\t'
    chunk_size = 2500
    h = pandas.read_csv(hgmd, sep='\t').rename(index=str, columns={'start':'start_tmp'})
    #I think it is breaking because it has start_tmp values as strings not ints so trying to fix that
    h['chrom'] = h['chrom'].astype(str)
    h['ref'] = h['ref'].astype(str)
    h['alt'] = h['alt'].astype(str)
    h['start_tmp'] = h['start_tmp'].astype(int)
    hgmd_cols = list(h.columns.values)[4:]
    write_pandas_header(in_file, out_file,
                        hgmd_cols + ['is_pathogenic', 'is_likely_pathogenic'],
                        [], delimiter)
    
    key = ('chrom', 'start_tmp', 'ref', 'alt')
    for df in pandas.read_csv(in_file, delimiter=delimiter, chunksize=chunk_size):
        df.loc[:, 'start_tmp'] = df.apply(lambda row: row['start']+1, axis=1)
        #Again making sure the int columns are of int type
        df['chrom'] = df['chrom'].astype(str)
        df['start_tmp'] = df['start_tmp'].astype(int)
        m = pandas.merge(df, h, on=key, how='left')
        m.loc[:, 'is_pathogenic'] = m.apply(is_p_fix, axis=1)
        m.loc[:, 'is_likely_pathogenic'] = m.apply(is_lp_fix, axis=1)
        m.to_csv(out_file, index=False, sep=delimiter, header=False, mode='a')

def update_p_and_lp_again(in_file, hgmd, out_file):
    delimiter = '\t'
    chunk_size = 2500
    h = pandas.read_csv(hgmd, sep='\t').rename(index=str, columns={'hgmd_class_ls':'hgmd_class_ls_use'})
    hgmd_cols = list(h.columns.values)[4:]
    write_pandas_header(in_file, out_file,
                        hgmd_cols + ['is_pathogenic_final', 'is_likely_pathogenic_final'],
                        [], delimiter)
    
    key = ('chrom', 'start', 'ref', 'alt')
    for df in pandas.read_csv(in_file, delimiter=delimiter, chunksize=chunk_size):
        m = pandas.merge(df, h, on=key, how='left').fillna('None')
        m.loc[:, 'is_pathogenic_final'] = m.apply(is_p_fix_again, axis=1)
        m.loc[:, 'is_likely_pathogenic_final'] = m.apply(is_lp_fix_again, axis=1)
        m.to_csv(out_file, index=False, sep=delimiter, header=False, mode='a')

def filter_vqsr(in_file, out_file):
    """This is just a filter."""
    delimiter = '\t'
    chunk_size = 5000
    write_pandas_header(in_file, out_file, [], [], delimiter)
    dtypes = {'CADD_phred':str, 'CADD_raw':str, 'tot_exac_cov_frac_10':str}
    for df in pandas.read_csv(in_file, delimiter=delimiter, chunksize=chunk_size, dtype=dtypes):
        keep_cols = [x for x in list(df.columns.values)]
        df[keep_cols].to_csv(out_file, index=False, sep=delimiter, header=False, mode='a')
    
# or ( (row['eff_indel_splice'] == 1 or
#                           row['is_splicing'] == 1) and row['impact_severity'] == 'HIGH' )

