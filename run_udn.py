#! /usr/bin/env python3

"""
input should be the project path
the folder should contains
0. the settings file (yaml), the config file
1. the VCF files
2. the phenotype/symptom txt file
3. query_terms, in this file, include the terms in short, to filter on the AMELIE result, each term per line


output:
anno_exp = f'{pw}/{sample_id}_{v[1]}.annotated.txt'
anno_filter = f'{pw}/{sample_id}_{v[1]}.filtered.txt
proband_genelist = f'{pw}/{proband}_{v[1]}.filtered.genelist
anno_extract = f'{pw}/{sample_id}_{v[1]}.extracted.txt
anno_tally = f'{pw}/{lb}.for_tally.pdict'

hpo = f'{proband}_terms_hpo.txt
hpo_pure = f'{proband}_terms_pure_hpo.txt
amelie_lite = f'{proband}.amelie.lite.txt
amelie_final = f'{proband}.amelie.matched_terms.txt

result_raw = f'{proband}.final_result.all.sorted.xlsx
result_selected = f'{proband}.final_result.selected.xlsx
report = f'{proband}.report.xlsx'

"""

import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('config', help="""the config file, yaml format""")
args = ps.parse_args()

import os
import sys
import re
import yaml
import logging
import logging.config

# get the config
cfg = yaml.safe_load(open(args.config).read())
vcf_filter = 'PASS'
pw_code = os.path.dirname(os.path.realpath(__file__))


##########
prj = cfg['prj']
sv_type_convert = {'<DEL>': 'deletion', '<DUP>': 'duplication'}
col_keep_new_name = [
    "anno_id", "chr_", "pos_s", "pos_e", "sv_type", "filter", "QUAL", "data", "gn", "sv_len", "exon_span",
    "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
    "annot_ranking"]


col_keep_raw_name = ['AnnotSV ID',
                     'SV chrom',
                     'SV start',
                     'SV end',
                     'SV type',
                     'FILTER',
                     'QUAL',
                     'FORMAT',
                     'Gene name',
                     'tx length',
                     'location',
                     'DGV_GAIN_Frequency',
                     'DGV_LOSS_Frequency',
                     'GD_AF',
                     'DDD_mode',
                     'DDD_disease',
                     'Mim Number',
                     'Phenotypes',
                     'Inheritance',
                     'AnnotSV ranking']


# config the logging
def get_logger(prj):

    fl_log_conf = f'{pw_code}/logging_setting.yaml'
    log_prefix = prj
    log_prefix = log_prefix + '_' if log_prefix else ''
    with open(fl_log_conf) as f:
        cfg = yaml.safe_load(f.read())
        cfg['handlers']['console']['level'] = 'INFO'
        cfg['handlers']['file']['mode'] = 'w'
        cfg['handlers']['file']['filename'] = log_prefix + cfg['handlers']['file']['filename']
        cfg['handlers']['error']['filename'] = log_prefix + cfg['handlers']['error']['filename']
        logging.config.dictConfig(cfg)
    return logging.getLogger('main')


logger = get_logger(prj)


def get_annotsv(annotsv):
    if annotsv and os.path.exists(annotsv):
        return annotsv

    import platform
    machine_platform = platform.system()
    if machine_platform == 'Darwin':
        home = os.path.expanduser('~')
        dock_path = f'{home}/dock/annotsv.sif'
    else:
        dock_path = '/scratch/cqs/chenh19/dock/annotsv.sif'
    # the annot_sv_rank and gnomad_AF filtering would only take effect on the proband, would not filter on the parent

    if os.path.exists(dock_path):
        return f'singularity exec {dock_path}  /tools/AnnotSV/bin/AnnotSV'
    logger.warning('No annotSV specified')
    return None


def verify_config(cfg):
    path = cfg['path']
    prj = cfg['prj']
    pheno_file = cfg['files']['pheno_file'].strip()
    pheno_file_select = cfg['files']['pheno_file_select'].strip()

    proband_id = cfg['IDs']['proband_id']
    proband_sex = cfg['IDs']['proband_sex']
    father_id = cfg['IDs']['fater_id']
    mother_id = cfg['IDs']['mother_id']
    male_sibling_id = cfg['IDs']['male_sibling_id']
    female_sibling_id = cfg['IDs']['female_sibling_id']

    annotsv_config = cfg['default']['annotsv']  # default=None
    thres_annot_sv_ranking = cfg['default']['thres_annot_sv_ranking']   # default is 2
    thres_gnomad_maf = cfg['default']['thres_gnomad_maf']   # default is 2
    vcf_file_path = cfg['default']['vcf_file_path']
    # default is the same as path
    vcf_file_path = vcf_file_path or path

    valid_family_mem = {}
    for i, sex, type_, type_short in zip(
        [proband_id, father_id, mother_id, male_sibling_id, female_sibling_id],
        [proband_sex, 1, 2, 1, 2],
        ["Proband", 'Father', 'Mother', 'Male_sibling', 'Female_sibbing'],
            "PFMSS"):
        if not i:
            continue
        i = i.strip().upper()
        if re.match(r'(UDN)?\d+$', i):
            # element = type_, type_short, sex, file_vcf, file_anno, anno_filtered, anno_extracted, anno_pdict
            sample_id = "UDN" + i.replace("UDN", '')
            valid_family_mem[sample_id] = [type_, type_short, sex, 0, 0, 0, 0, 0]
        else:
            logger.warning(f'invalid {type_} ID:  {i}')

    if len(valid_family_mem) < 1:
        logger.warning('Only proband are specified, no parental/sibling data available')

    if int(proband_sex) not in [1, 2]:
        logger.error(f'wrong proband sex, should be 1 or 2: ({proband_sex})')
        sys.exit(1)

    # valid file
    if not os.path.exists(pheno_file):
        logger.error(f'phenotype File not exist: {pheno_file}')
        sys.exit(1)
    pheno_pure_fn = pheno_file.rsplit('/', 1)[-1]
    os.symlink(pheno_file, f'{prj}/{pheno_pure_fn}')

    thres_annot_sv_ranking = thres_annot_sv_ranking or 2
    thres_gnomad_maf = thres_gnomad_maf or 0.01

    assert isinstance(thres_gnomad_maf, float), 'thres_gnomad_maf in config file should be float'
    assert isinstance(thres_annot_sv_ranking, int), 'thres_annot_sv_ranking in config file should be int'

    return [path, prj, pheno_file, pheno_file_select, vcf_file_path, valid_family_mem, annotsv_config, thres_annot_sv_ranking, thres_gnomad_maf]


from amelie_api import main
# from backports import csv  # for writing unicode circle into the csv file


def parse_pw(pw, vcf_file_path, valid_family_mem):
    proband = 'proband'
    for sample_id, v in valid_family_mem.items():
        if v[0] == 'Proband':
            proband = sample_id
        # get vcf file
        vcf = os.popen(
            f'find {vcf_file_path} -type f -iname "*{sample_id}*.vcf.gz" -o -iname "*{sample_id}*.vcf"').read().split('\n')
        vcf = [_.strip() for _ in vcf if _.strip()]
        if len(vcf) == 0:
            logger.error(f'{v[0]}  {sample_id} vcf file not found! ')
            sys.exit(1)
        elif len(vcf) == 1:
            valid_family_mem[sample_id][3] = vcf[0]
        else:
            logger.error(
                f'multiple({len(vcf)}) VCF file for {v[0]}  {sample_id} under {vcf_file_path} found, please check ')
            sys.exit(1)

        # get annotated file
        anno_exp = f'{pw}/{sample_id}_{v[1]}.annotated.txt'
        if not os.path.exists(anno_exp):
            valid_family_mem[sample_id][4] = anno_exp

    file_exp = {}
    file_exp['hpo'] = f'{pw}/{proband}_terms_hpo.txt'
    file_exp['amelie'] = f'{pw}/{proband}_amelie_query_terms.txt'

    return valid_family_mem, file_exp, proband


def parse_hpo_ref(hpo_ref, pw_code):
    """
    input file should have 2 columns,
    col1 = HPO ID, col2 = phenotype
    sep by tab
    """
    try:
        with open(hpo_ref) as fp:
            data = [_.strip().split('\t') for _ in fp]
            data = {k.strip(): v.strip() for v, k in data}
            with open(f'{pw_code}/hpo_ref.pdict', 'w') as out:
                print(data, file=out)
            return data
    except:
        return 0


def get_hpo_ref_file():
    # build the ref file
    hpo_ref = args.hpo or f'{pw_code}/hpo_ref.pdict'
    logger.debug(hpo_ref)

    if not os.path.exists(hpo_ref):
        if os.path.exists(f'{pw_code}/hpo_ref.txt'):
            logger.debug('rebuild hpo_ref.pdict')
            hpo_db = parse_hpo_ref(f'{pw_code}/hpo_ref.txt', pw_code)
            if not hpo_db:
                logger.error('fail to build hpo_ref.pdict')
                sys.exit(1)
        else:
            logger.error(f'hpo_ref file not found')
            sys.exit(1)
    elif hpo_ref.find('hpo_ref.pdict') < 0:
        logger.debug('rebuild hpo_ref.pdict')
        hpo_db = parse_hpo_ref(hpo_ref, pw_code)
        if not hpo_db:
            logger.error(f'fail to build hpo_ref.pdict from  {hpo_ref}')
            sys.exit(1)
    else:
        hpo_db = eval(open(hpo_ref).read())

    return hpo_db


def get_hpo_id(pw, proband, pheno, hpo_db):
    fl_result = f'{pw}/{proband}_terms_hpo.txt'
    fl_result1 = f'{pw}/{proband}_terms_pure_hpo.txt'
    if os.path.exists(fl_result):
        logger.info(f'the term to HPO ID file already done! {fl_result}')
        return 0

    logger.info('Map HPO terms to IDs')
    res = []
    res_pure_hpo = []
    bad = []  # terms not exactly found in hpo db
    candidate = []

    for i in open(pheno):
        i = i.strip()
        try:
            res.append(f'{i}\t{hpo_db[i]}')
            res_pure_hpo.append(hpo_db[i])
        except:
            iquery = set(re.split(r'\s+', i))
            len1 = len(iquery) * 0.8
            for ipheno in hpo_db:
                ipheno = set(re.split(r'\s+', ipheno))
                if len(ipheno & iquery) > len1:
                    candidate.append('{i}\t{hpo_db[ipheno]}\t{ipheno}')
            if len(candidate) == 0:
                bad.append(i)

    if len(bad) > 0:
        logger.debug(f'{len(bad)} terms are not converted to HPO ID')
        with open(f'{pw}/error_hpo_terms_not_converted.txt', 'w') as out:
            print('\n'.join(bad), file=out)
    else:
        logger.debug('All terms are converted!')

    if len(candidate) > 0:
        logger.debug(f'{len(candidate)} terms are ambiguous, please check {pw}/{proband}_terms_ambiguous.txt')
        with open(f'{pw}/{proband}_terms_ambiguous.txt', 'w') as out:
            print('UDN_symptom\tHPO_term\tHPO_ID', file=out)
            print('\n'.join(candidate), file=out)

    with open(fl_result, 'w') as out:
        print('\n'.join(res), file=out)
    with open(fl_result1, 'w') as out:
        print('\n'.join(res_pure_hpo), file=out)
    return 0


def get_annot_col(ianno):
    tmp = os.popen(f"head -1 {ianno}").read().strip()
    cols = tmp.split('\t')
    logger.info(f'annotation file columns = {len(cols)}')

    col_keep = {}

    for new, raw in zip(col_keep_new_name, col_keep_raw_name):
        try:
            col_keep[new] = cols.index(raw)
        except:
            logger.error('column not found in annot file: {raw} (for {new} ) ')
            return 'err'
    col_keep['data'] += 1
    return col_keep


def anno_filter(lb, col_keep, type_, thres_gnomad_maf, thres_annot_sv_ranking):
    f_anno_exp = f'{pw}/{lb}.annotated.txt'
    f_anno_filter = f'{pw}/{lb}.filtered.txt'
    # type_ = the type of this file, P/M/F/S
# ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "QUAL", "data", "gn", "sv_len", "exon_span",
#              "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
#              "annot_ranking"]
    # treat the annot_SV_ranking column diffently between proband and parents
    # for the proband, we demand the rank >=2
    # for the parents, the rank have no limit

    col_gnomad_AF = col_keep['af_gnomad'] + 1
    col_annot_sv_ranking = col_keep['annot_ranking'] + 1
    col_filter = col_keep['filter'] + 1
    col_gn = col_keep['gn'] + 1

    if type_ == 'P':  # proband
        extra_filter = f' && ${col_gnomad_AF}<{thres_gnomad_maf} && ${col_annot_sv_ranking}>={thres_annot_sv_ranking}'
    else:
        extra_filter = ''

    logger.info(f'{lb}:  filter criteria = FILTER==PASS {extra_filter}')

    cmd = f"""head -1 {f_anno_exp} > {f_anno_filter};awk -F $'\\t'  '${col_filter}=="PASS" {extra_filter}' {f_anno_exp} >> {f_anno_filter} """
    # logger.info(cmd)
    os.system(cmd)

    # if proband, extract the genelist
    if type_ == 'P':
        cmd = f"""cut -d $'\\t' -f {col_gn} {f_anno_filter}|sed '1d' |sort|uniq > {lb}.filtered.genelist"""
        # logger.info(cmd)
        os.system(cmd)

    lines = os.popen(f'wc -l < {lb}.filtered.genelist').read().strip()
    try:
        lines = int(lines)
    except:
        logger.error(f'fail to get the line number of {lb}.filtered.genelist')

    if lines == 0:
        logger.error(f'{lb}:  no gene was found')


def extract_anno(lb, col_keep, sex):
    f_anno_filter = f'{pw}/{lb}.filtered.txt'
    f_anno_extract = f'{pw}/{lb}.extracted.txt'

    file_start = 1
    out = open(f_anno_extract, 'w')
    print(
        '\t'.join(
            ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "qual", "exon_span_tag", "gn", "sv_len", "exon_span",
             "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
             "annot_ranking", 'copy_number']),
        file=out)

    with open(f_anno_filter) as fp:
        for i in fp:
            a = i.strip().split('\t')

            # exclude the first header line
            if file_start:
                if i.find('AnnotSV') > -1:
                    file_start = 0
                    continue

            anno_id, chr_, pos_s, pos_e, sv_type, _filter, qual, data, gn, \
                sv_len, exon_span, af_dgv_gain, af_dgv_loss, af_gnomad, \
                ddd_mode, ddd_disease, omim, phenotype, inheritance, annot_ranking = [a[col_keep[_]] for _ in col_keep_new_name]

            # data = expon_span_tag/FORMAT
            # sv_len = tx length (overlap of the SV and transcript)
            # exon_span = location

            chr_ = 'chr' + chr_.upper()
            try:
                sv_len = f'{int(sv_len)/1000:.1f}kbp'
            except:
                logger.warning('wrong sv_len format: sv_len={sv_len}  : {gn}  {anno_id}')

            # copy number
            try:
                data = int(data.split(':')[1])  # CN
                if data == 0:
                    copy_number = '-'
                else:
                    copy_number = 'o' * data
            except:
                copy_number = 'NA'
                logger.warning('wrong copy number format: copy number={data}  : {gn}  {anno_id}')

            if chr_ == 'chrX' and sex == 1:
                copy_number += 'y'
            elif chr_ == 'chrY' and sex == 2 and copy_number != '-' and copy_number != 'NA':
                copy_number = 'invalid: ' + copy_number
            elif chr_ == 'chrY' and sex == 1:
                copy_number = 'x' + copy_number

            # sv type
            try:
                sv_type = sv_type_convert[sv_type]
            except:
                sv_type = 'NA'
                logger.warning('wrong sv_type format: sv_type={sv_type}  : {gn}  {anno_id}')

            # gnomad
            try:
                af_gnomad = float(af_gnomad)
                if af_gnomad == -1:
                    af_gnomad = "'-"
                elif af_gnomad >= 0.01:
                    af_gnomad = f'{af_gnomad:.2f}'
                else:
                    af_gnomad = f'{af_gnomad:.2e}'
            except:
                af_gnomad = 'NA'
                logger.warning('wrong af_gnomad format: af_gnomad={af_gnomad}  : {gn}  {anno_id}')

            # dgv gain
            try:
                af_dgv_gain = float(af_dgv_gain)
                if af_dgv_gain == -1:
                    af_dgv_gain = '-'
                elif af_dgv_gain >= 0.01:
                    af_dgv_gain = f'{af_dgv_gain:.2f}'
                else:
                    af_dgv_gain = f'{af_dgv_gain:.2e}'
            except:
                af_dgv_gain = 'NA'
                logger.warning('wrong af_dgv_gain format: af_dgv_gain={af_dgv_gain}  : {gn}  {anno_id}')

            # dgv loss
            try:
                af_dgv_loss = float(af_dgv_loss)
                if af_dgv_loss == -1:
                    af_dgv_loss = '-'
                elif af_dgv_loss >= 0.01:
                    af_dgv_loss = f'{af_dgv_loss:.2f}'
                else:
                    af_dgv_loss = f'{af_dgv_loss:.2e}'
            except:
                af_dgv_loss = 'NA'
                logger.warning('wrong af_dgv_loss format: af_dgv_loss={af_dgv_loss}  : {gn}  {anno_id}')

            # exons_span_tag
            exon_span_tag = 'please_check,this_is_init_state'
            m = re.match(r'^(intron|exon|)\d+$', exon_span)
            if m:
                ele_type = m.group(1)
                if ele_type == 'exon':
                    exon_span_tag = 'cover one exon'
                else:
                    exon_span_tag = ''
            else:
                m1 = re.match(r'^(intron|exon|txStart)(\d*)-(intron|exon|txEnd)(\d*)$', exon_span)

                if m1:
                    head, head_num, tail, tail_num = m1.groups()
                    head = head.strip()
                    tail = tail.strip()

                    if tail == 'txEnd':
                        if head == 'txStart':
                            exon_span_tag = 'covered multiple exons'
                        else:
                            exon_span_tag = 'please check the exon_span'
                    elif head == 'txStart':
                        try:
                            tail_num = int(tail_num)
                        except:
                            logger.warning(f'error exon_span format: exon_span={exon_span}, {gn}  {anno_id}')
                            continue
                        if tail_num == 1:
                            exon_span_tag = 'cover one exon'
                        elif tail_num > 1:
                            exon_span_tag = 'covered multiple exons'
                        else:
                            logger.warning(f'error exon_span format: exon_span={exon_span}, {gn}  {anno_id}')
                            continue
                    else:
                        try:
                            head_num = int(head_num)
                            tail_num = int(tail_num)
                        except:
                            logger.warning(f'error exon_span format: exon_span={exon_span}, {gn}  {anno_id}')
                            continue

                        if tail_num - head_num == 0 and 'exon' in [head, tail]:
                            exon_span_tag = 'cover one exon'

                        elif tail_num - head_num == 0 and head == 'intron' and tail == 'intron':
                            exon_span_tag = ''

                        elif tail_num - head_num == 1 and 'intron' in [head, tail]:
                            exon_span_tag = 'cover one exon'
                        else:
                            exon_span_tag = 'covered multiple exons'
                else:
                    exon_span_tag = 'error, no pattern matched'

            print('\t'.join([anno_id, chr_, pos_s, pos_e, sv_type, qual, exon_span_tag, gn,
                             sv_len, exon_span, af_dgv_gain, af_dgv_loss, af_gnomad,
                             ddd_mode, ddd_disease, omim, phenotype, inheritance, annot_ranking, copy_number]), file=out)
    out.close()


def build_parent_dict(suffix):

    # add the parent information
    # first , we read the father and mother, build a dict
    # {chr1: {
    #   1:[[starts], [ends]],
    #   2:[[starts], [ends]],
    #   3:[[starts], [ends]],
    #   }
    # chr2: {
    #   1:[[starts], [ends]],
    #   2:[[starts], [ends]],
    #   3:[[starts], [ends]],
    #   }
    # }

    # the 2nd layer key, is the  int(start/100000)  i.e. 100k bin,
    # this is to avoid the proband site iter over all the sites
    # however, the proband - father/mother site is not always matched,
    # on the proband search step, the each site would search the upper and lower bin, e.g. the proband site is 123,456,789,  the bin is 123, but, I would search both 122, 123, 124 bin in the father/mother dict
    # the list , starts, ends, they should be sorted.

    extra = ''
    if suffix not in 'MFS':
        logger.error(f'parent should be M / D / F/ S, your input={suffix}')
        return 0
    elif suffix == 'F':
        extra = f'-o -iname "*_D.extracted_fields.txt" '

    parent = os.popen(f'find {pw} -type f -iname "*_{suffix}.extracted_fields.txt" {extra}').read().split('\n')
    parent = [_.strip() for _ in parent if _.strip()]

    if len(parent) == 0:
        logger.error(f'Extracted fields file for parent(suffix={suffix}) is not found')
        return 0
        # sys.exit(1)
    elif len(parent) > 1:
        logger.error(f'Extracted fields file for parent(suffix={suffix}): expect=1, {len(parent)} found.')
        sys.exit(1)

    d_parent = {}
    with open(parent[0]) as fp:
        for i in fp:
            a = i.strip().split('\t')
            # chr_, s, e, sv_type, gn, cn = [a[_] for _ in [1, 2, 3, 4, 7, -1]]
            chr_, s, e, sv_type, gn, cn = [a[_] for _ in [1, 2, 3, 4, 7, -1]]
            try:
                s = int(s)
                e = int(e)
            except:
                continue
            bin_ = int(s / 100000)
            try:
                d_parent[chr_][bin_].append([s, e, sv_type, gn, cn])
            except:
                try:
                    d_parent[chr_][bin_] = [[s, e, sv_type, gn, cn]]
                except:
                    d_parent[chr_] = {}
                    d_parent[chr_][bin_] = [[s, e, sv_type, gn, cn]]

    # sort the list
    for k1, v1 in d_parent.items():
        for k2, v2 in v1.items():
            v3 = sorted(v2, key=lambda x: x[0])
            d_parent[k1][k2] = v3

    return d_parent


def get_sites_in_bin(n, chr_, bin_, site_id, task_parent, dict_parent, parent_type):
    def get_sites_in_bin_single_run(bin_new, phase):
        # phase =  -1, 0, or 1, which means the position relative to the query_bin
        try:
            task_parent.extend(dict_parent[chr_][bin_new])
        except:
            if chr_ not in dict_parent:
                logger.debug(f'n={n} - phase={phase}:  chr not found in {parent_type} dict : query_id = {site_id}')
            elif not dict_parent[chr_].get(bin_ - 1):
                logger.debug(f'n={n} - phase={phase}:  100k bin not found in {parent_type} dict : query_id = {site_id}')

    get_sites_in_bin_single_run(bin_ - 1, 'bin-1')
    get_sites_in_bin_single_run(bin_, 'bin')
    get_sites_in_bin_single_run(bin_ + 1, 'bin+1')


def run_proband_match(proband_final, d_amelie, trio):

    # trio is a list,  2 element
    # indicate the sample type beside proband
    # could be m, f   or m, s (mother and sibling)  or s, s(2 sibling)

    trio1, trio2 = trio
    lb_map = {'m': 'mother', 'd': 'father', 'f': 'father', 's': 'sibling'}
    dict_map = {'m': dict_mother, 'd': dict_father, 'f': dict_father, 's': dict_sibling}

    trio_lb1 = lb_map[trio1]
    trio_lb2 = lb_map[trio2]

    final = open(proband_final, 'w')
    # logger.debug(proband_final)
    with open(file_proband[0]) as fp:
        task1 = {}
        task2 = {}
        n = 0
        for i in fp:
            i = i.strip()
            n += 1
            if n == 1:
                print(
                    "AMELIE\t" + i + '\t' +
                    f'{trio_lb1}_copy_number,{trio_lb2}_copy_number,{trio_lb1}_gn,{trio_lb2}_gn,{trio_lb1}_sv_type,{trio_lb1}_s,{trio_lb1}_e,{trio_lb1}_overlap,{trio_lb2}_sv_type,{trio_lb2}_s,{trio_lb2}_e,{trio_lb2}_overlap'.
                    replace(',', '\t'),
                    file=final)
                continue
            a = i.split('\t')
            site_id, chr_, s, e, gn = [a[_] for _ in [0, 1, 2, 3, 7]]

            # logger.info(f'{gn}\t{site_id}')

            s = int(s)
            e = int(e)
            bin_ = int(s / 100000)

            # info is the final appended info result
            info = {}
            try:
                info[amelie] = d_amelie[gn]
            except:
                logger.debug(f'n={n}: gene not found in AMELIE list, due to the score is too low gene = {gn}')
                continue

            if bin_ not in task1:
                task1[bin_] = []
                task2[bin_] = []

                # father
                get_sites_in_bin(n, chr_, bin_, site_id, task1[bin_], dict_map[trio1], trio_lb1)

                # mother
                get_sites_in_bin(n, chr_, bin_, site_id, task2[bin_], dict_map[trio2], trio_lb2)

            for task, lb in zip([task1, task2], [trio_lb1, trio_lb2]):
                matched = ''
                for isite in task[bin_]:
                    # d_parent[chr_][bin_] = [[s, e, sv_type, gn, cn]]
                    s_q, e_q = isite[0], isite[1]

                    # no overlap
                    # if s_q > e or e_q < s:
                    #     pass
                    if s_q <= s and e_q >= s and e_q <= e:
                        overlap = e_q - s
                        if not matched or matched[0] < overlap:
                            matched = [overlap] + isite
                    elif s_q >= s and e_q <= e:
                        overlap = e_q - s_q
                        if not matched or matched[0] < overlap:
                            matched = [overlap] + isite
                    elif s_q <= s and e_q >= e:
                        overlap = e - s
                        if not matched or matched[0] < overlap:
                            matched = [overlap] + isite

                    elif s_q >= s and s_q <= e and e_q >= e:
                        overlap = e - e_q
                        if not matched or matched[0] < overlap:
                            matched = [overlap] + isite
                if matched:
                    # [[overlap, s, e, sv_type, gn, cn]]
                    matched = [str(_) for _ in matched]

                    # new order = # copy_number, gn, sv_type, s, e, overlap
                    new = [matched[_] for _ in [-1, -2, -3, 1, 2, 0]]
                    logger.debug(f'n={n}: {lb}: {new}')
                    info[lb] = new
                else:
                    logger.debug(f'n={n}: {lb}: no match found')
                    if chr_.lower().find('chrx') > -1:
                        info[lb] = ['oy', 'not_found'] + [''] * 4
                    elif chr_.lower().find('chry') > -1:
                        info[lb] = ['wild_type', 'not_found'] + [''] * 4
                    else:
                        info[lb] = ['oo', 'not_found'] + [''] * 4
            res = [info[trio_lb1][0], info[trio_lb2][0], info[trio_lb1][1],
                   info[trio_lb2][1]] + info[trio_lb1][2:] + info[trio_lb2][2:]
            print(info[amelie] + '\t' + i + "\t" + '\t'.join(res), file=final)
    final.close()


def get_amelie_dict(amelie):
    d_amelie = [_.strip().split('\t') for _ in open(amelie)]
    d_amelie = {k: v for k, v in d_amelie}
    logger.info(f'AMELIE count = {len(d_amelie)}')
    return d_amelie


if __name__ == "__main__":
    pw, prj, pheno_file, pheno_file_select, \
        vcf_file_path, valid_family_mem, annotsv_config, \
        thres_annot_sv_ranking, thres_gnomad_maf \
        = verify_config(cfg)

    # valid_family_mem, a dict, key=ID, value=list [type_, type_short, sex, file_vcf=3, file_anno=4]

    # could be real path or None
    annot_sv_app = get_annotsv(annotsv_config)
    os.chdir(pw)

    # get file
    valid_family_mem, file_exp, proband = parse_pw(pw, vcf_file_path, valid_family_mem)

    #  get_hpo_id
    if not os.path.exists(file_exp['hpo']):
        logger.info('step 1:  get HPO id for the phenotype')
        hpo_db = get_hpo_ref_file()
        get_hpo_id(pw, proband, pheno_file, hpo_db)

    # run annotSV
    for sample_id, v in valid_family_mem.items():
        vcf = v[3]
        anno = v[4]
        lb = f'{sample_id}_{v[1]}'
        anno_exp = f'{pw}/{lb}.annotated.txt'

        if anno:
            logger.info(f'annotation for {v[0]} : {sample_id} already done')
            continue

        if annot_sv_app:
            cmd = f'{annot_sv_app} -genomeBuild GRCh38 -typeOfAnnotation split -outputDir AnnotSV -SVinputFile {vcf}'
            logger.debug(cmd)
            os.system(cmd)
        else:
            logger.error('annotSV is not specified, also, annotation file is not found')
            sys.exit(1)

        # rename/hard-link the annotation result
        vcf_short = vcf.rsplit('/', 1)[-1].replace('.vcf', '').replace('.gz', '')
        os.system(f'ln {pw}/AnnotSV/{vcf_short}.annotated.tsv {anno_exp}')

        if not os.path.exists(anno_exp):
            logger.error(f'Annotation failed: {anno_exp}')
            sys.exit(1)
        v[4] = anno_exp

    # get column number
    col_keep = get_annot_col(valid_family_mem[proband][4])

    # filter the annotation file
    for sample_id, v in valid_family_mem.items():
        lb = f'{sample_id}_{v[1]}'
        type_ = v[1]  # type_short
        sex = v[2]  # 1=male, 2=female
        f_anno_exp = f'{pw}/{lb}.annotated.txt'
        f_anno_filter = f'{pw}/{lb}.filtered.txt'
        f_anno_extract = f'{pw}/{lb}.extracted.txt'
        f_anno_tally = f'{pw}/{lb}.for_tally.pdict'

        # filter annotation file
        if os.path.exists(f_anno_filter):
            logger.info(f'anno filter already done {f_anno_filter}')
        else:
            logger.info(f'step 2: filter the annotation file : {lb}')
            anno_filter(lb, col_keep, type_, thres_gnomad_maf, thres_annot_sv_ranking)

        # extract fields
        if os.path.exists(f_anno_extract):
            logger.info(f'anno extracted fields already done {f_anno_extract}')
        elif os.path.exists(f_anno_filter):
            logger.info(f'step 3: extract the annotation file : {lb}')
            extract_anno(lb, col_keep, sex)
        else:
            logger.error('filtered annotation result not found:  {lb}')


    # get the amelie result

    amelie = os.popen(f'find {pw} -type f -iname "*.amelie.lite.txt"').read().split('\n')
    amelie = [_.strip() for _ in amelie if _.strip()]

    if len(amelie) == 0:
        logger.info('getting information from AMELIE')
        genelist = os.popen(f'find {pw} -iname "*_P.filtered.genelist"').read().strip().split('\n')
        genelist = [_.strip() for _ in genelist if _.strip()]
        if len(genelist) == 0:
            logger.error('No proband genelist found, name pattern = *_P.filtered.genelist')
            sys.exit(1)
        main(prj, genelist[0], hpo_pure_id)
        amelie = f'{prj}.amelie.lite.txt'
        d_amelie = get_amelie_dict(amelie)

    else:
        # add the amelie information, combine the father and mother copy number
        # the amelie result should ends with .amelie.txt
        amelie = amelie[0]
        d_amelie = get_amelie_dict(amelie)

    # add the parent information
    # first , we read the father and mother, build a dict
    # {chr1: {
    #   1:[[starts], [ends]],
    #   2:[[starts], [ends]],
    #   3:[[starts], [ends]],
    #   }
    # chr2: {
    #   1:[[starts], [ends]],
    #   2:[[starts], [ends]],
    #   3:[[starts], [ends]],
    #   }
    # }

    # the 2nd layer key, is the  int(start/100000)  i.e. 100k bin,
    # this is to avoid the proband site iter over all the sites
    # however, the proband - father/mother site is not always matched,
    # on the proband search step, the each site would search the upper and lower bin, e.g. the proband site is 123,456,789,  the bin is 123, but, I would search both 122, 123, 124 bin in the father/mother/sibling dict
    # the list , starts, ends, they should be sorted.

    trio = []
    # Father
    suffix = 'F'
    file_dict_parent = os.popen(f'find {pw} -type f -iname "father.for_tally.pdict"').read().split('\n')
    file_dict_parent = [_.strip() for _ in file_dict_parent if _.strip()]
    if len(file_dict_parent) == 0:
        logger.info(f'building parent dict (suffix = {suffix})')
        dict_father = build_parent_dict(suffix)
        with open('father.for_tally.pdict', 'w') as out:
            print(dict_father, file=out)
    else:
        logger.info('father dict already done')
        dict_father = eval(open(file_dict_parent[0]).read())

    if dict_father:
        trio.append('f')

    # mother
    suffix = 'M'
    file_dict_parent = os.popen(f'find {pw} -type f -iname "mother.for_tally.pdict"').read().split('\n')
    file_dict_parent = [_.strip() for _ in file_dict_parent if _.strip()]
    if len(file_dict_parent) == 0:
        logger.info(f'building parent dict (suffix = {suffix})')
        dict_mother = build_parent_dict(suffix)
        with open('mother.for_tally.pdict', 'w') as out:
            print(dict_mother, file=out)
    else:
        logger.info('mother dict already done')
        dict_mother = eval(open(file_dict_parent[0]).read())

    if dict_mother:
        trio.append('m')

    # sibling
    suffix = 'S'
    file_dict_parent = os.popen(f'find {pw} -type f -iname "sibling.for_tally.pdict"').read().split('\n')
    file_dict_parent = [_.strip() for _ in file_dict_parent if _.strip()]
    if len(file_dict_parent) == 0:
        logger.info(f'building parent dict (suffix = {suffix})')
        dict_sibling = build_parent_dict(suffix)
        with open('sibling.for_tally.pdict', 'w') as out:
            print(dict_sibling, file=out)
    else:
        logger.info('sibling dict already done')
        dict_sibling = eval(open(file_dict_parent[0]).read())

    if dict_sibling:
        trio.append('s')

    # search the proband
    # d_parent[chr_][bin_] = [[s, e, sv_type, gn, cn]]
    file_proband = os.popen(f'find {pw} -type f -iname "*_P.extracted_fields.txt"').read().split('\n')
    file_proband = [_.strip() for _ in file_proband if _.strip()]

    if len(file_proband) == 0:
        logger.error('Extracted fields file for proband is not found')
        sys.exit(1)

    logger.info(f'start to match father, mother /sibling information with proband')
    logger.info(f'the probind file used for matching is {file_proband[0]}')
    logger.debug(f'trio={trio}')

    proband_final = f'{prj}_final_result.tsv'
    run_proband_match(proband_final, d_amelie, trio)

    # amelie score is column 21
    final_sorted = proband_final.replace('.tsv', '.sorted.tsv')
    os.system(
        f"head -1 {proband_final} > {final_sorted};sed '1d' {proband_final} |sort -t $'\\t' -k 1nr  >> {final_sorted}")
