#! /usr/bin/env python3

"""
input should be the project path
the folder should contains
1. the VCF files
2. the phenotype/symptom txt file, ends with "_terms.txt"

"""
import platform
import os
import sys
import re
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('prj', help="""project name""")
ps.add_argument('pw', help="""optioinal, project path, default is current folder""", nargs='?', default=os.getcwd())
ps.add_argument('-annotsv', help="""path for annotSV software, default is using the annotSV in ACCRE dock""", nargs='?', default='')
ps.add_argument(
    '-hpo', help="""designate the HPO refence file, optional. 2columns, col1=HPO ID, col2=phenotype, sep by tab""",
    nargs='?', default=None)
args = ps.parse_args()

# load logging
import logging
import logging.config
import yaml
from amelie_api import main
# from backports import csv  # for writing unicode circle into the csv file
home = os.path.expanduser('~')
machine_platform = platform.system()

pw_code = os.path.dirname(os.path.realpath(__file__))
prj = args.prj


fl_log_conf = f'{pw_code}/logging_setting.yaml'
log_prefix = 'UDN'
log_prefix = log_prefix + '_' if log_prefix else ''
with open(fl_log_conf) as f:
    cfg = yaml.safe_load(f.read())
    cfg['handlers']['console']['level'] = 'INFO'
    cfg['handlers']['file']['mode'] = 'w'
    cfg['handlers']['file']['filename'] = log_prefix + cfg['handlers']['file']['filename']
    cfg['handlers']['error']['filename'] = log_prefix + cfg['handlers']['error']['filename']
    logging.config.dictConfig(cfg)
logger = logging.getLogger('main')


#######################

# annotSV output file filter settings
thres_annot_sv_ranking = 2  # equal or largher than this
thres_gnomad_maf = 0.01   # less than this value
vcf_filter = 'PASS'  # the filter filed should  use this value

if machine_platform == 'Darwin':
    dock_path = f'{home}/dock/annotsv.sif'
else:
    dock_path = '/scratch/cqs/chenh19/dock/annotsv.sif'
# the annot_sv_rank and gnomad_AF filtering would only take effect on the proband, would not filter on the parent

if os.path.exists(args.annotsv):
    annot_sv_app = args.annotsv
elif os.path.exists(dock_path):
    annot_sv_app = args.annotsv or f'singularity exec {dock_path}  /tools/AnnotSV/bin/AnnotSV'
else:
    logger.warning('No annotSV specified')
    annot_sv_app = ''

###################
pw = args.pw
os.chdir(pw)

sv_type_convert = {'<DEL>': 'deletion', '<DUP>': 'duplication'}
col_keep_new_name = ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "filter", "QUAL", "data", "gn", "sv_len", "exon_span",
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




def parse_pw(pw):
    vcf = os.popen(f'find {pw} -type f -iname "*.vcf.gz" -o -iname "*.vcf"').read().split('\n')
    vcf = [_.strip() for _ in vcf if _.strip()]
    logger.debug(f'vcf file list: {vcf}')

    anno = os.popen(f'find {pw} -type f -iname "*annotated.txt" -o  -iname "*annotated.tsv" ').read().split('\n')
    anno = [_.strip() for _ in anno if _.strip()]

    pheno = os.popen(f'find {pw} -type f -iname "*_terms.txt"').read().split('\n')
    pheno = [_.strip() for _ in pheno if _.strip()]
    logger.debug(f'phenotype file: {pheno}')

    run_anno = 1

    if len(vcf) == 0:
        if len(anno) > 0:
            logger.warning('VCF file not found, annotation file found. would skip the annotSV step')
            run_anno = 0
        else:
            logger.error(f'no VCF file found under {pw}')
            sys.exit(1)
    else:
        if len(anno) < len(vcf):
            pass
        elif len(anno) > len(vcf):
            logger.warning(
                f'VCF file count({len(vcf)}) different from annotation file count({len(anno)}), would re-run the annotation')
            run_anno = 0
        else:
            logger.info('annotation file already exist, would skip the annotSV step')
            run_anno = 0

    if len(pheno) == 0:
        logger.error(f'no phenotype found under {pw}')
        sys.exit(1)
    elif len(pheno) > 1:
        logger.error(f'{len(pheno)} was found, please check file ends with _terms.txt')
        sys.exit(1)
    return vcf, pheno[0], anno, run_anno

def refine_anno_name(anno):
    """
    refine the annotation filename
    input = anno name list
    this would rename it, and return the new name
    """
    p = re.compile(r'.*?(?:\d+[_-])?(?:UDN)?(\d{4,})[_.-]([PMSDF])[._].*annotated\.(?:tsv|txt)$')
    new_anno = []
    for ianno in anno:
        m = re.match(p, ianno)

        if m:
            udn_id, trio_type = m.groups()
        else:
            logger.error(f'annotatioin file pattern not match:  {ianno}')
            return 'error'

        new_name = f'{pw}/UDN-{udn_id}_{trio_type}_annotated.txt'
        logger.debug(f'rename: {ianno} {new_name}')
        os.system(f'ln -s {ianno} {new_name} 2>/dev/null')
        new_anno.append(new_name)

    return new_anno



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


def get_hpo_id(pw, pheno, hpo_db):

    lb = os.path.basename(pheno).rsplit('.', 1)[0]
    fl_result = f'{pw}/{lb}_hpo.txt'
    fl_result1 = f'{pw}/{lb}_pure_hpo.txt'
    if os.path.exists(fl_result):
        logger.info(f'the term to HPO ID file already done! {fl_result}')
        return fl_result1

    logger.info('get HPO id')
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
        logger.debug(f'{len(candidate)} terms are ambiguous, please check {pw}/ambiguous_hpo_terms.txt')
        with open(f'{pw}/ambiguous_hpo_terms.txt', 'w') as out:
            print('UDN_symptom\tHPO_term\tHPO_ID', file=out)
            print('\n'.join(candidate), file=out)

    with open(fl_result, 'w') as out:
        print('\n'.join(res), file=out)
    with open(fl_result1, 'w') as out:
        print('\n'.join(res_pure_hpo), file=out)
    return fl_result1


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


def anno_filter(ianno, col_keep):
    anno_lb = ianno.replace(".annotated.txt", "").replace("_annotated.txt", "")
    filtered = f'{anno_lb}_annotated.filtered.txt'

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

    if ianno.find("_P_annotated.txt") > 0:
        extra_filter = f' && ${col_gnomad_AF}<{thres_gnomad_maf} && ${col_annot_sv_ranking}>={thres_annot_sv_ranking}'
    else:
        extra_filter = ''

    logger.info(f'{anno_lb}:  filter criteria = FILTER==PASS {extra_filter}')


    cmd = f"""head -1 {ianno} > {filtered};awk -F $'\\t'  '${col_filter}=="{vcf_filter}" {extra_filter}' {ianno} >> {filtered} """
    # logger.info(cmd)
    os.system(cmd)


    cmd = f"""cut -d $'\\t' -f {col_gn} {filtered}|sed '1d' |sort|uniq > {anno_lb}.filtered.genelist"""
    # logger.info(cmd)
    os.system(cmd)


    lines = os.popen(f'wc -l < {anno_lb}.filtered.genelist').read().strip()
    try:
        lines = int(lines)
    except:
        logger.error(f'fail to get the line number of {anno_lb}.filtered.genelist')

    if lines == 0:
        logger.error(f'{anno_lb}:  no gene was found')

    else:
        os.system(f'echo finished > {anno_lb}.filter.done')


def extract_anno(ianno_filter, anno_lb_lite, col_keep):
    # logger.info(f'file enter into extraction : {ianno_filter} ')

    file_start = 1
    out = open(f'{anno_lb_lite}.extracted_fields.txt', 'w')
    print(
        '\t'.join(
            ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "qual", "exon_span_tag", "gn", "sv_len", "exon_span",
             "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
             "annot_ranking", 'copy_number']),
        file=out)

    with open(ianno_filter) as fp:
        for i in fp:
            a = i.strip().split('\t')

            if file_start:
                if i.find('AnnotSV') > -1:
                    file_start = 0
                    continue

            anno_id, chr_, pos_s, pos_e, sv_type, _filter, qual, data, gn, \
                sv_len, exon_span, af_dgv_gain, af_dgv_loss, af_gnomad, \
                ddd_mode, ddd_disease, omim, phenotype, inheritance, annot_ranking = [a[col_keep[_]] for _ in col_keep_new_name]

            # ddd_mode  such as monoallelic

            chr_ = 'chr' + chr_
            try:
                sv_len = f'{int(sv_len)/1000:.1f}kbp'
            except:
                logger.warning('wrong sv_len format: sv_len={sv_len}  : {gn}  {anno_id}')

            # copy number
            try:
                data = int(data.split(':')[1])
                if data == 0:
                    copy_number = '-'
                else:
                    copy_number = 'o' * data
            except:
                logger.warning('wrong copy number format: copy number={data}  : {gn}  {anno_id}')

            try:
                sv_type = sv_type_convert[sv_type]
            except:
                logger.warning('wrong sv_type format: sv_type={sv_type}  : {gn}  {anno_id}')

            try:
                af_gnomad = float(af_gnomad)
                if af_gnomad == -1:
                    af_gnomad = '-'
                elif af_gnomad >= 0.01:
                    af_gnomad = f'{af_gnomad:.2f}'
                else:
                    af_gnomad = f'{af_gnomad:.2e}'
            except:
                logger.warning('wrong af_gnomad format: af_gnomad={af_gnomad}  : {gn}  {anno_id}')

            try:
                af_dgv_gain = float(af_dgv_gain)
                if af_dgv_gain == -1:
                    af_dgv_gain = '-'
                elif af_dgv_gain >= 0.01:
                    af_dgv_gain = f'{af_dgv_gain:.2f}'
                else:
                    af_dgv_gain = f'{af_dgv_gain:.2e}'
            except:
                logger.warning('wrong af_dgv_gain format: af_dgv_gain={af_dgv_gain}  : {gn}  {anno_id}')

            try:
                af_dgv_loss = float(af_dgv_loss)
                if af_dgv_loss == -1:
                    af_dgv_loss = '-'
                elif af_dgv_loss >= 0.01:
                    af_dgv_loss = f'{af_dgv_loss:.2f}'
                else:
                    af_dgv_loss = f'{af_dgv_loss:.2e}'
            except:
                logger.warning('wrong af_dgv_loss format: af_dgv_loss={af_dgv_loss}  : {gn}  {anno_id}')

            exon_span_tag = 'please_check_init_state'
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
            res = [info[trio_lb1][0], info[trio_lb2][0], info[trio_lb1][1], info[trio_lb2][1]] + info[trio_lb1][2:] + info[trio_lb2][2:]
            print(info[amelie] + '\t' + i + "\t" + '\t'.join(res), file=final)
    final.close()

def get_amelie_dict(amelie):
        d_amelie = [_.strip().split('\t') for _ in open(amelie)]
        d_amelie = {k: v for k, v in d_amelie}
        logger.info(f'AMELIE count = {len(d_amelie)}')
        return d_amelie


if __name__ == "__main__":

    hpo_db = get_hpo_ref_file()
    # logger.info(type(hpo_db))

    # get file
    vcf, pheno, anno, run_anno = parse_pw(pw)
    #  get_hpo_id
    logger.info('step 1:  get HPO id for the phenotype')
    hpo_pure_id = get_hpo_id(pw, pheno, hpo_db)

    # run annotSV
    if run_anno:
        if annot_sv_app:
            for ivcf in vcf:
                ivcf = os.path.realpath(ivcf)
                cmd = f'{annot_sv_app} -genomeBuild GRCh38 -typeOfAnnotation split -SVinputFile {ivcf}'
                logger.debug(cmd)
                os.system(cmd)
            vcf, pheno, anno, run_anno = parse_pw(pw)

            if len(anno) == 0:
                logger.error('no anno file found after annotSV, please check! ')
                sys.exit(1)
        else:
            logger.error('annotSV is not specified, you need to run the annotation first')
            sys.exit(1)

    # refine the anno filename
    anno = refine_anno_name(anno)
    if anno == 'error':
        logger.error('fail to refine anno name')
        sys.exit(1)
    logger.info(anno)

    # get column number
    col_keep = get_annot_col(anno[0])
    # logger.info(col_keep)

    # filter / extract information the annotation file
    for ianno in anno:
        anno_lb = ianno.replace(".annotated.txt", "").replace("_annotated.txt", "")
        anno_lb_lite = anno_lb.rsplit('/', 1)[-1]

        if os.path.exists(f'{anno_lb}.filter.done'):
            logger.info(f'anno filter already done {anno_lb_lite}')
        else:
            logger.info(f'step 2: filter the annotation file : {anno_lb_lite}')
            anno_filter(ianno, col_keep)

        if os.path.exists(f'{anno_lb}.filter.done'):
            if os.path.exists(f'{anno_lb}.extract.done'):
                logger.info(f'anno extract already done {anno_lb_lite}')
            else:
                logger.info(f'step 3: extract the annotation file : {anno_lb_lite}')
                extract_anno(f'{anno_lb}_annotated.filtered.txt', anno_lb_lite, col_keep)
        else:
            logger.error('filtered annotation result not found:  {anno_lb_lite}')

        try:
            if os.path.getsize(f'{anno_lb_lite}.extracted_fields.txt') > 100:
                os.system(f'echo extract done > {anno_lb}.extract.done')
            else:
                logger.warning(f'extract_field file, filesize too small : {anno_lb_lite}.extracted_fields.txt')
        except:
            logger.warning('extracted file not found')

    # this step should be done manually, go to amelie website, and get the AMELIE score
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
