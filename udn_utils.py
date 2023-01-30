"""
2020-08-18
1. add the UDN gateway API, allow to download the VCF file/ phenotip
2. add the encrypt / decrypt module for the API token
3. pheno.keywords.txt allow comment using #

to-do
add the downloading data / uploading data to AWS server


2021-11-30
1. fix the firefox path issue
2. change SVminSize to 1 for AnnotSV
3. change the extraction step, if sv_len < 50 and the SV is not BND and the SV is not covering exon, exclude it
4. get the gene list in the extract step (prev in filter step)
5. for vcf entry like below, the segment is super long, but it's not actually an SV, exlude it
chr1    817859  DRAGEN:REF:chr1:817859-2649577  N       .       94      PASS    END=2649577;REFLEN=1831719      GT:SM:CN:BC:PE  ./.:1.00865:2:1428:7,138


2021-12-21
1. updated the annotSV to 3.1.1 version , the output file header and content changed a lot
-typeOfAnnotation -> -annotationMode
# -overlap , Minimum overlap (%) between user features (User BED) and the annotated SV to be reported Range values: [0-100], default = 100

new column Exon_count, type=int

The gain/loss af were split to  P_  and B_, P=pathogenic, B=benign. the AF thres can be set by parameter


"""

import os
import sys
import re
import json
import yaml
import time
import logging
import logging.config
import requests

from bs4 import BeautifulSoup as bs
import pickle
# import pandas as pd
from selenium.webdriver.support.ui import WebDriverWait as wait
# basic data
pw_code = os.path.dirname(os.path.realpath(__file__))
sys.path.append(pw_code)
from . import amelie_api

platform = sys.platform


redundant_words = set('and,am,is,in,the,are,to,for,of,a,an,one,two,three,four,at,on,type,activity,increased,decreased,low,level,high'.split(','))

sv_type_convert = {'<DEL>': 'deletion', '<DUP>': 'duplication'}

# this one, excluded the data, data_end_flag and format
col_keep_final = [
    "anno_id", "chr_", "pos_s", "pos_e", "sv_len", "sv_type", "filter", "QUAL", "gn", "exon_span",
    "gain_b_source", 'gain_b_af',
    "loss_b_source", 'loss_b_af',
    "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
    "annot_ranking"]

col_keep_new_name = [
    "anno_id", "chr_", "pos_s", "pos_e", "sv_len", "sv_type", "filter", "QUAL", "format", "data", 'data_end_flag', "gn", "exon_span", 'exon_count',
    "gain_b_source", 'gain_b_af',  "loss_b_source", 'loss_b_af',  "inv_b_source", 'inv_b_af',  "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
    "annot_ranking", 'info']

col_keep_raw_name = ['AnnotSV_ID',
                     'SV_chrom',
                     'SV_start',
                     'SV_end',
                     'SV_length',
                     'SV_type',
                     'FILTER',
                     'QUAL',
                     'FORMAT',
                     'FORMAT',
                     'Annotation_mode',
                     'Gene_name',
                     'Location',
                     'Exon_count',
                     'B_gain_source',
                     'B_gain_AFmax',
                     'B_loss_source',
                     'B_loss_AFmax',
                     'B_inv_source',
                     'B_inv_AFmax',
                     'DDD_mode',
                     'DDD_disease',
                     'OMIM_ID',
                     'OMIM_phenotype',
                     'OMIM_inheritance',
                     'AnnotSV_ranking_score',
                     'INFO',
                     ]

# config the logging
def get_logger(pw, prj):
    prefix = prj
    fn_log = f'{pw}/{prefix}.runudn.log'
    fn_err = f'{pw}/{prefix}.runudn.err'

    fmt = logging.Formatter('%(asctime)s  %(levelname)-6s %(funcName)-20s  line: %(lineno)-5s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    console.setLevel('INFO')

    fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
    fh_err = logging.FileHandler(fn_err, mode='a', encoding='utf8')

    fh_file.setLevel('DEBUG')
    fh_file.setFormatter(fmt)

    fh_err.setLevel('ERROR')
    fh_err.setFormatter(fmt)

    try:
        logger = logging.getLogger(__file__)
    except:
        logger = logging.getLogger('terminal')
    logger.setLevel('DEBUG')
    logger.addHandler(console)
    logger.addHandler(fh_file)
    logger.addHandler(fh_err)

    return logger

def line_count(fn):
    n = 0
    for _ in open(fn):
        n += 1
    return n


intermediate_folder = 'intermediate'
thres_quality_score = 40



def test_pheno_pattern(l_pheno, l_pheno_raw, logger):
    patterns = []
    for ipheno, ipheno_raw in zip(l_pheno, l_pheno_raw):
        tmp = []
        for _ in ipheno:
            _ = _.lower()
            if _[0] == '-':  # negative match, exclude
                word = f'exclude: {_[1:]}'
            else:
                if _[0] == '@' or _.find('|') > 0:
                    extra_flag = r'\b'
                    if _[0] == '@':
                        word = _[1:]
                    else:
                        word = _

                    if word.find('|') > -1 and word.find('(') < 0:
                        word = re.sub(r'\b((?:\w+\|)+\w+)', r'(\g<1>)', word)
                else:
                    word = _
                    extra_flag = ''
                word = fr'.*{extra_flag}({word}[a-z]*){extra_flag}'
            tmp.append(word)
        patterns.append([' '.join(ipheno), tmp])

    tmp = json.dumps(patterns, indent=4)
    logger.info(tmp)




class UDN_case():
    def __init__(self, config_file):
        cfg = yaml.safe_load(open(config_file).read())
        self._cfg = cfg
        self.prj = cfg['prj']
        self.pw = cfg['path'] or os.getcwd()
        self.proband_id = cfg['IDs']['proband_id']
        logger = get_logger(self.pw, self.prj)
        self.logger = logger

        self.header_prefix = 'coord,note,match_count,mut_type,AMELIE'.split(',')

        # verify the config
        cfg_tested = self.verify_config()
        # logger.info(cfg_tested)


        self.path = cfg_tested['path']
        self.sv_caller = cfg_tested['sv_caller']
        self.pheno_file = cfg_tested['pheno_file']
        self.pheno_file_select = cfg_tested['pheno_file_select']
        self.family = cfg_tested['family']
        self.vcf_file_path = cfg_tested['vcf_file_path']
        self.thres_annot_sv_ranking = cfg_tested['thres_annot_sv_ranking']
        self.thres_gnomad_maf = cfg_tested['thres_gnomad_maf']
        self.amelie_dict = {}
        self.udn_match = cfg_tested['udn_match']

        # copy the config file
        os.system(f'cp {config_file} {self.pw}/{self.prj}_settings.yaml 2>/dev/null')
        os.system(f'mkdir -p {self.pw}/{intermediate_folder} 2>/dev/null')

        # update self.family, add the vcf_file
        self.parse_pw()

        # convert hpo terms to hpo id
        self.get_hpo_id()
        os.system(f'ln {self.pw}/{intermediate_folder}/{self.prj}.terms_hpo.txt {self.pw}/{self.prj}.terms_hpo.txt 2>/dev/null')

        # get annotSV path
        self.annotsv_app = self.get_annotsv()

        # logger.info(self.annotsv_app)

        if not self.annotsv_app:
            logger.warning('No annotSV specified...')
                # sys.exit(1)

        # verfify if the automatic process is already done
        self.done_phase1 = 0
        self.done_phase2 = 0
        self.done_phase3 = 0
        # logger.info(cfg)
        # logger.info(f'pw={self.pw}, prj={self.prj}, vcfpath={self.vcf_file_path}')
        fn_final = f'{self.pw}/merged.sorted.{self.prj}.xlsx'
        if os.path.exists(fn_final):
            f_size = os.path.getsize(fn_final)
            if f_size < 10000:
                logger.warning(f'the merged anno result already exist: {fn_final}, but the size is only {f_size/1000:.2f}')
            else:
                self.done_phase1 = 1

        # verify if the selected gene list is done
        pw = self.pw
        prj = self.prj
        for fn_tmp in [f'{pw}/{prj}.selected.genes.xlsx', f'{pw}/selected.genes.xlsx', f'{pw}/{prj}.selected.genes.txt', f'{pw}/selected.genes.txt']:
            if fn_tmp and os.path.exists(fn_tmp):
                fn_selected = fn_tmp
                f_size = os.path.getsize(fn_selected)
                if f_size < 200:
                    logger.warning(f'selected gene list already exist: {fn_selected}, but the size is only {f_size/1000:.2f}')
                else:
                    self.done_phase2 = 1
                    break

        # verify if the final report is done
        fn_report = f'{self.pw}/report.{self.prj}.xlsx'
        if os.path.exists(fn_report):
            f_size = os.path.getsize(fn_report)
            if f_size < 200:
                logger.warning(f'report already exist {fn_report}, but the size is only {f_size/1000:.2f}')
            else:
                self.done_phase3 = 1
        try:
            self.get_annot_col()
        except:
            pass



    def omim_query(self, force_rerun=False):
        """
        the input is the merged.sorted.tsv file  and the phenotype keyword file
        1. build the gene comment dict based on the previous manual input
        2. read the local OMIM disease discription, if not found, would run the online scrapy
        3. query the keyword file for each gene in the merged.sorted.tsv file, and rate(hit number) as reference for selecting

        the phenotype file, accept keyword in 3 mode
        word1+word2  this 2 word come together
        @word1 match exact word, would add \b at both side
        -word  do not include this word
        """
        fn = f'{self.pw}/merged.sorted.{self.prj}.tsv'
        fn_pheno = f'{self.pw}/pheno.keywords.txt'
        udn = self.prj
        pw = self.pw
        logger = self.logger
        driver=None


        try:
            cols = self.cols
        except:
            self.get_annot_col()
            cols = self.cols

        if not force_rerun and os.path.exists(f'{pw}/omim_match_gene.{udn}.md'):
            logger.debug(f'OMIM query report already done')
            return 0

        if not os.path.exists(fn_pheno) or os.path.getsize(fn_pheno) < 10:
            fn_pheno = f'{self.pw}/origin/{self.prj}.terms.txt'
            if not os.path.exists(fn_pheno):
                logger.error(f'No phenotype file found for query')
                return 1
            logger.warning(f'OMIM pheno: pheno.keywords.txt not exist, use {fn_pheno} instead')

        # get the gene list with OMIM desease
        d_gene = {}
        n_gene_total = 0
        genes_without_omim_id = set()
        with open(fn) as fp:
            header = fp.readline()
            header = header.strip().split('\t')

            # get the column number for proband cn and last cn of the last family number
            idx = {}
            idx['proband'] = header.index('annot_ranking') + 1

            l = ['qual', 'gn', 'omim', 'sv_type', 'sv_len', 'exon_span_tag', 'AMELIE']
            for _ in l:
                idx[_] = header.index(_)

            for n, i in enumerate(header):
                if i.lower().startswith('gn_'):
                    idx['end_idx'] = n
                    break
            else:
                idx['end_idx'] = len(header)

            # for n, i in enumerate(header):
            #     if idx['proband'] == 'na' and i.lower().startswith('proband'):
            #         idx['proband'] = n
            #     if idx['end_idx'] == 'na' and i.lower().startswith('gn_'):
            #         idx['end_idx'] = n
            try:
                sum(idx.values())
                valid_cn = 1
                cn_header = header[idx['proband']: idx['end_idx']]
            except:
                valid_cn = 0
                logger.error('invalid CN header found: {fn}')

            # logger.info(f'index# : {idx}\nheader={header}')

            header_other_info = ['coord', 'sv_type', 'gain_af_max', 'loss_af_max', 'filter']
            for _ in header_other_info:
                idx[_] = header.index(_)

            n_gene_proband_2copy = 0
            multiple_hit = 0
            for i in fp:
                a = i.split('\t')
                gn = a[idx['gn']].strip()
                omim_id = a[idx['omim']].strip()
                sv_type = a[idx['sv_type']]
                sv_len = a[idx['sv_len']]
                qual_proband = a[idx['qual']]
                n_gene_total += 1

                other_info = ['',]
                for _ in header_other_info:
                    sep = '\t' if len(_) >= 7 else '\t\t'
                    other_info.append(f'{_}:{sep}{a[idx[_]]}')
                try:
                    amelie_score = float(a[idx['AMELIE']])
                except:
                    logger.warning(f'wrong amelie score: {a}')
                    amelie_score = -1

                cover_exon_flag = 1 if a[idx['exon_span_tag']].strip() else 0
                copy_number = ''
                if valid_cn:
                    copy_number = a[idx['proband'] :idx['end_idx']]

                    # if the copy number for proband is 2, then unable to interpret
                    if copy_number[0][:3] == 'oo@':
                        n_gene_proband_2copy += 1
                        continue

                    # copy_number[0] = f'{copy_number[0]}@{qual_proband}'
                    copy_number = [_.strip() for _ in copy_number if _.strip]

                    copy_number = [' = '.join(_) for _ in zip(cn_header, copy_number)] + other_info
                    copy_number = '\n\t\t'.join(copy_number)
                    copy_number = f'SV={sv_type}, svlen={sv_len} @@\n\t\t' + copy_number

                if omim_id:
                    omim_id = str(int(float(omim_id.split(';')[0])))
                else:
                    omim_id = 'NA'
                    genes_without_omim_id.add(gn)

                if gn not in d_gene:
                    d_gene[gn] = [omim_id, cover_exon_flag, amelie_score, '\n\t' + copy_number]
                else:
                    d_gene[gn][3] += '\n\t' + copy_number
                    multiple_hit += 1

                d_gene[gn][1] = max(cover_exon_flag, d_gene[gn][1])


        logger.info(f'd_gene gene count = {len(d_gene)}, genes with 2 copies = {n_gene_proband_2copy}, multiple_hit={multiple_hit},  expected_total = {len(d_gene) + n_gene_proband_2copy + multiple_hit}, n_gene_total = {n_gene_total}')

        # build the predefined gene list

        # logger.warning(f'POU6F2:{d_gene["POU6F2"]}')

        # fn_gene_description_omim_pkl = f'{pw_code}/omim.gene_description_omim.pkl'

        fn_omim_pkl = f'{pw_code}/omim_description.scrapy.pkl'
        # k1 = hgnc_gn, v = [omim_gene_id_, []]
        # ele2 = [pheno, pheno_id, full_pheno_title, description]
        # usually, we only need the full_pheno_title and description
        # or "no_linked_pheno"

        dump_omim_pkl = 0
        d_gene_comment_scrapy = {}

        # logger.warning(d_gene['NALCN'])
        # sys.exit(1)

        if os.path.exists(fn_omim_pkl):
            logger.info('load omim info from pickle')
            with open(fn_omim_pkl, 'rb') as f:
                d_gene_comment_scrapy = pickle.load(f)
        else:
            logger.info(f'file not found: {fn_omim_pkl}')

        omim_gn_total = set()
        non_omim_gn_total = set()
        for gn, v in d_gene_comment_scrapy.items():
            if isinstance(v[1], list) and len(v[1]) > 0:
                omim_gn_total.add(gn)
            elif v[1] == 'no_linked_pheno':
                non_omim_gn_total.add(gn)

        logger.info(f'local total gene number with OMIM description={len(omim_gn_total)};  gene without linked omim = {len(non_omim_gn_total)}')

        logger.info(fn_omim_pkl)
        genename_list = set(d_gene)
        gene_with_omim = genename_list & omim_gn_total
        genes_without_linked_pheno = genename_list & non_omim_gn_total
        genes_need_scrapy = genename_list - set(d_gene_comment_scrapy) - genes_without_omim_id
        
        logger.info(f'gene count in merged.sorted.tsv: total = {n_gene_total}, genes in d_gene = {len(d_gene)}')
        logger.info(f'\n\tgene count in with OMIM description: {len(gene_with_omim)}\n\tgenes without linked pheno = {len(genes_without_linked_pheno)},  \n\tgenes need scrapy OMIM = {len(genes_need_scrapy)}')
        

        # build the phenotype keyword list
        # if 2 words link together, use +, if match the exact word, add a prefix @
        # e.g. lung+cancer  @short+stature
        l_pheno = []
        l_pheno_raw = []

        force_exact_match = {'mental', }  # the words in this set will force exact match

        with open(fn_pheno) as fp:
            for i in fp:
                new = i.lower().split('#', 1)[0].strip()
                if not new:  # pure comment line
                    continue

                l_pheno_raw.append(i.strip())

                if new[0] == '!':  # regex pattern
                    l_pheno.append([new])
                    continue

                a = re.split(r'\s+', new)
                a = [_.strip() for _ in a if _.strip()]
                a = [f'@{_}' if _ in force_exact_match else _ for _ in a]

                if len(a) == 0:
                    continue
                l_pheno.append([_.replace('+', ' ') for _ in a])

        logger.warning(l_pheno)


        # refine the l_pheno,

        # res , key=gene, value = [0, 0]  first int = count of exact match phenotype, second int = count of partial match phenotype
        res = {}

        # logger.info(f'l_pheno={l_pheno}')

        def refine_comment(s, symbol_to_word):
            for symbol, word in symbol_to_word.items():
                s = re.sub(symbol, f'**{word}**', s)
            s = re.sub(r'\*\*\*\*', '**', s)
            s = re.sub(r'\*\*\s+\*\*', ' ', s)
            # # **developmental **delay**,
            # s = re.sub(r'\*\*(\w[^*]*?\s)\*\*(\w[^*]*?)\w\*\*', r' **\g<1>\g<2>**', s)
            # # s = re.sub(r'\*\*\s+\*\*', ' ', s)
            # s = re.sub(r'\*\*([^* ][^*]*?[^* ])\*\*([^* ][^*]*?[^* ])\*\*', r'**\g<1>\g<2>**', s)
            # s = re.sub(r'\*\*([^* ][^*]*?[^* ])\*\*([^* ][^*]*?[^* ])\*\*', r'**\g<1>\g<2>**', s)
            return s

        # debug = 1
        n_no_omim = 0
        for gn, v in d_gene.items():
            omim_id, cover_exon_flag, amelie_score, copy_number = v
            if omim_id == 'NA':
                continue
            if gn in res:
                continue  # avoid the duplicate running
            res[gn] = [[], [], '', cover_exon_flag, amelie_score, copy_number]  # match pheno, partial match pheno, comment, cover_exon_flag, amelie rank
            comment = []

            if gn in omim_gn_total:
                omim_pheno_list = d_gene_comment_scrapy[gn][1]
            else:
                # run the scrapy
                # comment is a dict, key=pheno_id, v=pheno description
                logger.warning(f'run scrapy: {gn}')
                driver, ires = run_omim_scrapy(gn, omim_id, logger, driver=driver)
                if ires:
                    omim_pheno_list = ires[1]
                    d_gene_comment_scrapy[gn] = ires
                    dump_omim_pkl = 1
                else:
                    n_no_omim += 1
                    # logger.info(f'No pheno description found on OMIM: gn={gn}')
                    continue

            total_omim = len(omim_pheno_list)
            for sn_omim, iomim_pheno in enumerate(omim_pheno_list):
                try:
                    omim_title, omim_desc = iomim_pheno[2:]
                except:
                    logger.warning(f'invalid omim pickle item, gene = {gn}, v = {iomim_pheno}')
                    continue
                comment.append('\n' + '@' * 50 + f'\n\n\n### ({sn_omim + 1} / {total_omim}) **' + omim_title + '**')

                omim_desc = omim_desc.replace('Descriptioin\n', '#### Descriptioin \n').replace('Clinical\n', '#### Clinical \n').replace('Molecular\n', '#### Molecular \n')

                comment.append(f'\n{omim_desc}')

            comment_compact = '\n'.join(comment).lower().replace('\n', '')

            # query the phenotype
            highlighted_words = set() | redundant_words
            symbol_to_word = {}
            # print(highlighted_words)
            n_symbol = 0
            for ipheno, ipheno_raw in zip(l_pheno, l_pheno_raw):
                n_word_match_meaning = 0  # for not exporting the items match of, is and to export into partial match file
                n_matched_word = 0
                matched_word = []
                n_total_word = len(ipheno)
                for _ in ipheno:
                    _ = _.lower()

                    if _[0] == '-':  # negative match, exclude
                        word = _[1:]
                        if comment_compact.lower().find(word) < 0:
                            n_matched_word += 1
                            matched_word.append(_)
                    else:
                        if _[0] == '@' or _.find('|') > 0:
                            extra_flag = r'\b'
                            if _[0] == '@':
                                word = _[1:]
                            else:
                                word = _
                            if word.find('|') > -1 and word.find('(') < 0:
                                word = re.sub(r'\b((?:\w+\|)+\w+)', r'(\g<1>)', word)
                            pattern_suffix = r'[a-z]*'

                        elif _[0] == '!':  # the input is the regex pattern
                            word = _[1:]
                            extra_flag = ''
                            pattern_suffix = ''

                        else:
                            word = _
                            extra_flag = ''
                            pattern_suffix = r'[a-z]*'

                        #  {word}s  means that, the word could contain an extra s
                        m = re.match(fr'.*{extra_flag}({word}{pattern_suffix}){extra_flag}', comment_compact)

                        # if debug:
                        #     logger.info(m)
                        #     logger.info(fr'.*\b({word})\b')
                        # debug = 0
                        if m:
                            word = m.group(1)
                            n_matched_word += 1
                            matched_word.append(word)
                            # avoid the highligh of meaningless single letter or pure number
                            if len(word) == 1 or re.match(r'^\d+$', word):
                                continue

                            if word not in redundant_words:
                                n_word_match_meaning += 1
                            else:
                                continue

                            if word not in highlighted_words:
                                n_symbol += 1
                                symbol = f'@@{n_symbol}@@'
                                # logger.info(f'{symbol}, {word} {gn}')
                                highlighted_words.add(word)
                                symbol_to_word[symbol] = word
                                comment = [re.sub(word, symbol, icomment, flags=re.I) if icomment[:3] != '###' and icomment[-2:] != '**' else icomment for icomment in comment]

                if n_matched_word == n_total_word:
                    res[gn][0].append(f'full_match: {ipheno_raw}')
                elif n_word_match_meaning > 0:
                    res[gn][1].append(f'partial: {ipheno_raw}:{matched_word}')


            if n_word_match_meaning > 0:
                logger.debug(f'{gn}\t{highlighted_words - redundant_words}')

            comment = [refine_comment(_, symbol_to_word) for _ in comment]
            res[gn][2] = '\n'.join(comment)

        # logger.info(f'genes not found in OMIM: {n_no_omim} / {len(d_gene)}')

        # for gn, v1 in res.items():
        #     # [[], [], '', cover_exon_flag, amelie_score]  # match pheno, partial match pheno, comment, cover_exon_flag, amelie rank
        #     print(gn, [v1[_] for _ in [0, 1]])

        if dump_omim_pkl and len(d_gene_comment_scrapy) > 5000:
            logger.info('dumping OMIM records')
            with open(fn_omim_pkl, 'wb') as o:
                pickle.dump(d_gene_comment_scrapy, o)

        # logger.info(f'OMIM query count={len(res)}')


        # amelie result
        try:
            with open(f'{pw}/{udn}.amelie.matched_query.pkl', 'rb') as f:
                tmp = pickle.load(f)
            matched_amelie = {}
            for gn, v in tmp.items():
                amelie_str = v[2]
                ires = ['### amelie match']
                for amelie_pheno, amelie_v1 in amelie_str.items():
                    amelie_v1_dedup = []
                    dedup = set()
                    for i in amelie_v1:
                        if i[0] not in dedup:
                            dedup.add(i[0])
                            amelie_v1_dedup.append(i)


                    ires.append(f'#### {amelie_pheno} n articles = {len(amelie_v1_dedup)}')
                    for amelie_v2 in amelie_v1_dedup:
                        ires.append('- ' + '\t'.join(amelie_v2[:4]))

                amelie_str = '\n'.join(ires)
                matched_amelie[gn] = amelie_str

            logger.info(f'amelie matched genes = {len(matched_amelie)}')
        except:
            matched_amelie = {}



        # tally the result
        out_full_match = open(f'{pw}/omim_match_gene.{udn}.md', 'w')
        out_partial_match = open(f'{pw}/omim_partial_match_gene.{udn}.md', 'w')
        out_all_genes = open(f'{pw}/omim_all_genes.{udn}.md', 'w')
        # test = next(iter(res.items()))
        # logger.info(f'res_entry_len={len(test[1])}')
        # res 0 = matched phenotype, 1=partial matched phenotype, 2 = comment. 3 = exon_flag, 4=amelie score
        res1 = sorted(res.items(), key=lambda _: (_[1][3], len(_[1][0]), _[1][4]), reverse=True)
        res2 = sorted(res.items(), key=lambda _: (_[1][3], len(_[1][1]), _[1][4]), reverse=True)

        with open(f'{pw}/intermediate/omim_match_result.pkl', 'wb') as o:
            pickle.dump(res, o)

        n1 = 0
        n2 = 0
        n3 = 0
        gene_match = set()
        
        for gn, v in res.items():
            match, partial_match, comment, cover_exon_flag, amelie_score, copy_number = v
            if comment.strip() == '':
                logger.debug(f'gene with no OMIM description: {gn}')
                continue
            n3 += 1
            print(f'## {n3}:\t{gn} : {amelie_score} \tcover_exon={cover_exon_flag}\n{match}\n{partial_match}\n{copy_number}\n\n### main\n\n', file=out_all_genes)  # all genes
            print(comment, file=out_all_genes)
            print('#' * 50 + '\n\n\n', file=out_all_genes)

        s_amelie = set(matched_amelie)
        s_omim = set([_[0] for _ in res1 if len(_[1][0]) > 0])
        gn_not_included = s_amelie - s_omim
        total_full_match = len(s_amelie | s_omim)
        for v in res1:
            gn = v[0]
            amelie_str = matched_amelie.get(gn) or ''

            match, partial_match, comment, cover_exon_flag, amelie_score, copy_number = v[1]
            partial_match = '\n'.join(partial_match)
            match = '\n'.join(match)
            if len(match) > 0:
                gene_match.add(gn)
                n1 += 1
                print(f'## {n1}/{total_full_match}:\t{gn} : {amelie_score}\tcover_exon={cover_exon_flag}\n{match}\n{partial_match}\n{copy_number}\n{amelie_str}\n\n### main\n\n', file=out_full_match)
                print(comment, file=out_full_match)
                print('#' * 50 + '\n\n\n', file=out_full_match)


        not_in_d_gene = []
        for gn in gn_not_included:
            try:
                omim_id, cover_exon_flag, amelie_score, copy_number = d_gene[gn]
            except:
                not_in_d_gene.append(gn)
                continue
            
            try:
                comment = res[gn][2]
            except:
                comment = ''

            n1 += 1
            amelie_str = matched_amelie.get(gn) or ''
            print(f'## {n1}/{total_full_match}:\t{gn} amelie match only\nAmelie score = {amelie_score}\n{copy_number}\n\n{amelie_str}\n\n{comment}', file=out_full_match)
            print('#' * 50 + '\n\n\n', file=out_full_match)

        if len(not_in_d_gene) > 0:
            logger.warning(f'the following genes are not found in d_gene: \n\t' + '\n\t'.join(not_in_d_gene))


        total_partial_match = len([v[0] for v in res1 if len(v[1][0]) == 0 and len(v[1][1]) > 0])
        for v in res2:
            gn = v[0]
            if gn in gene_match:
                continue
            match, partial_match, comment, cover_exon_flag, amelie_score, copy_number = v[1]
            partial_match = '\n'.join(partial_match)
            if len(partial_match) > 0:
                n2 += 1
                print(f'## {n2}/{total_partial_match}:\t{gn} : {amelie_score}, exon = {cover_exon_flag}\n{copy_number}\tcover_exon={cover_exon_flag}\n{partial_match}{copy_number}\n\n### main\n\n', file=out_partial_match)
                print(comment, file=out_partial_match)
                print('#' * 50 + '\n\n\n', file=out_partial_match)
        out_full_match.close()
        out_partial_match.close()
        out_all_genes.close()
        logger.info(f'gene with OMIM match= {n1}')
        logger.info(f'gene with OMIM partial match= {n2}')


    def verify_config(self):
        cfg = self._cfg
        logger = self.logger
        path = cfg['path']
        pw = path or os.getcwd()
        # prj = cfg['prj']
        pheno_file = cfg['pheno_file'].strip()

        fn_udn_match = cfg.get('udn_match_gene_file')
        udn_match = {}

        if fn_udn_match:
            tmp = fn_udn_match.rsplit('/', 1)
            if len(tmp) == 1:
                pwtmp = pw
                fn_udn_match = tmp[0]
            else:
                pwtmp, fn_udn_match = tmp
            fn_found = os.popen(f'find {pwtmp} -iname "{fn_udn_match}"').read().strip().split('\n')
            fn_found = [_ for _ in fn_found if _.strip()]

            if fn_found:
                try:
                    with open(fn_found[0]) as fp:
                        fp.readline() # skip the header
                        for i in fp:
                            a = i.strip().split('\t')
                            gn = a[0]
                            ct_hits = int(a[2])
                            udn_match[gn] = ct_hits
                except FileNotFoundError:
                    logger.error(f'file for UDN match specified in config file, but not exist: {fn_udn_match}')
                    sys.exit(1)
                except ValueError:
                    logger.warning(f'wrong wrong format for the match count: {i.strip()}')

            # logger.info(f'udn_match_dict size = {len(udn_match)} fn_found={fn_found}')

        pheno_file_select = cfg['default']['pheno_file_select']

        # valid pheno file
        if pheno_file.find('/') < 0:
            tmp = os.popen(f'find {pw} -iname {pheno_file}').read().strip().split('\n')
            if len(tmp) == 0:
                logger.error(f'phenotype File not exist: {pheno_file}')
                sys.exit(1)
            else:
                pheno_file = tmp[0]
        elif not os.path.exists(pheno_file):
            logger.error(f'phenotype File not exist: {pheno_file}')
            sys.exit(1)


        if pheno_file_select:
            tmp = pheno_file_select.rsplit('/', 1)
            if len(tmp) == 1:
                pwtmp = pw
                pheno_file_select = tmp[0]
            else:
                pwtmp, pheno_file_select = tmp
            fn_found = os.popen(f'find {pwtmp} -iname "{pheno_file_select}"').read().strip().split('\n')
            fn_found = [_ for _ in fn_found if _.strip()]

        if fn_found:
            pheno_file_select = fn_found[0]
        else:
            pheno_file_select = pheno_file

        proband_id = cfg['IDs']['proband_id']
        proband_sex = cfg['IDs']['proband_sex']

        if proband_sex not in [1, 2]:
            logger.error('Wrong proband sex, valid=1 / 2, where 1=male, 2=female')
            sys.exit(1)

        thres_annot_sv_ranking = cfg['default']['thres_annot_sv_ranking']   # default is 2
        thres_gnomad_maf = cfg['default']['thres_gnomad_maf']   # default is 2
        vcf_file_path = cfg['default']['vcf_file_path']
        # default is the same as path
        vcf_file_path = vcf_file_path or path or os.getcwd()
        pw_origin = f'{vcf_file_path}/origin'
        if os.path.exists(pw_origin):
            vcf_file_path = pw_origin


        sv_caller = cfg.get('sv_caller') or 'dragen'

        family = {}
        family[proband_id] = {
                'lb': f'{proband_id}_P',
                'rel': 'proband',
                'type': 'Proband',
                'sex': proband_sex,
                'type_short': 'P',
                'vcf': None,
                'aff_state': 1,
        }
        for k, v in cfg['IDs'].items():
            if not v:
                continue
            if k.lower().startswith('proband'):
                continue

            try:
                sample_id, aff_state = v[:2]
            except:
                logger.error(f'Wrong value specified in IDs scope of config file, key={k}, value={v}.  value shoudl have 2 elements')
                sys.exit(1)

            if len(v) == 3:
                rel_to_proband = re.sub(r'\s+', '_', v[2].strip())
            else:
                rel_to_proband = None

            if aff_state not in [0, 1, -1]:
                logger.error('Wrong affect state sepcified, key={k}, aff_state = {aff_state}.  valid = 0 / 1 ,  -1 for unkown')
                sys.exit(1)

            sample_id = sample_id.upper()
            m = re.match(r'(UDN)?\d+$', sample_id)
            if not m:
                logger.error(f'Wrong UDN ID for {k} : value={sample_id}')
                sys.exit(1)


            if k == 'father':
                family[sample_id] = {
                'lb': f'{sample_id}_F',
                'rel': k,
                'type': 'Father',
                'sex': 1,
                'type_short': 'F',
                'vcf': None,
                'aff_state': aff_state,
                }
            elif k == 'mother':
                family[sample_id] = {
                'lb': f'{sample_id}_M',
                'rel': k,
                'type': 'Mother',
                'sex': 2,
                'type_short': 'M',
                'vcf': None,
                'aff_state': aff_state,
                }
            elif k.lower().startswith('male'):
                family[sample_id] = {
                'lb': f'{sample_id}_S',
                'rel': 'male',
                'type': rel_to_proband or k.title(),
                'sex': 1,
                'type_short': 'S',
                'vcf': None,
                'aff_state': aff_state,
                }
            elif k.lower().startswith('female'):
                family[sample_id] = {
                'lb': f'{sample_id}_S',
                'rel': 'female',
                'type': rel_to_proband or k.title(),
                'sex': 2,
                'type_short': 'S',
                'vcf': None,
                'aff_state': aff_state,
                }
            else:
                logger.error(f'Invalid key found under IDs scope in config file: key={k}, value={v}')
                sys.exit(1)

        if len(family) < 2:
            logger.warning('Only proband are specified, no parental/sibling data available')

        pheno_pure_fn = pheno_file.rsplit('/', 1)[-1]
        try:
            os.symlink(pheno_file, f'{pw}/{intermediate_folder}/{pheno_pure_fn}')
        except:
            pass
        thres_annot_sv_ranking = thres_annot_sv_ranking or 2
        thres_gnomad_maf = thres_gnomad_maf or 0.01

        assert isinstance(thres_gnomad_maf, float), 'thres_gnomad_maf in config file should be float'
        assert isinstance(thres_annot_sv_ranking, int), 'thres_annot_sv_ranking in config file should be int'

        res = {'path': path,
               'pheno_file': pheno_file,
               'pheno_file_select': pheno_file_select,
               'family': family,
               'sv_caller': sv_caller,
               'vcf_file_path': vcf_file_path,
               'thres_annot_sv_ranking': thres_annot_sv_ranking,
               'thres_gnomad_maf': thres_gnomad_maf,
               'udn_match': udn_match if len(udn_match) > 0 else None
               }
        return res

    def get_annotsv(self) -> "None or the annotsv app path":
        """
        get the path for annotSV
        """
        annotsv = self._cfg['default']['annotsv']
        # logger = self.logger
        if annotsv and os.path.exists(annotsv):
            return annotsv

        import platform
        machine_platform = platform.system()
        machine_node = platform.node()
        if machine_platform == 'Darwin':
            dock_path = f'/Users/files/dock/annotsv.sif'
            mount = '-B /Users/files'
        elif machine_node.find('viccbiostat120') > -1:
            # personal desktop
            dock_path = '/mnt/d/dock/annotSV.sif'
            mount = '-B /mnt/d -B /data'
        else:
            # dock_path = '/scratch/cqs/chenh19/dock/annotsv.sif'
            dock_path = '/data/cqs/chenh19/dock/annotsv.sif'
            mount = ''
        # the annot_sv_rank and gnomad_AF filtering would only take effect on the proband, would not filter on the parent

        if os.path.exists(dock_path):
            return f'singularity exec {mount} {dock_path} /tools/AnnotSV/bin/AnnotSV'
        # logger.warning('No annotSV specified')
        return None

    def parse_pw(self):
        """
        set the vcf file path, annot file path
        """
        logger = self.logger
        pw = self.pw
        valid_family_mem = self.family
        vcf_file_path = self.vcf_file_path
        sv_caller = self.sv_caller
        proband_id = self.proband_id


        if sv_caller == 'pacbio':
            vcf = os.popen(f'find {pw} -type f -iname "*.pbsv.vcf"  -o -iname "*.pbsv.vcf.gz" 2>/dev/null').read().strip().split('\n')
            vcf = [_ for _ in vcf if _.strip()]
            v = valid_family_mem[proband_id]

            if len(vcf) != 1:
                logger.warning(f'Pacbio SV file: *.pbsv.vcf.gz not found')
                sys.exit(1)
            v['vcf'] = vcf[0]
            v['lb'] = v['lb'] + '.pacbio'
            return 0

        # if the vcf not exist, try to download it
        vcfs = os.popen(f'find {vcf_file_path} -type f -iname "*.cnv.vcf"  -o -iname "*.cnv.vcf.gz" 2>/dev/null').read().strip().split('\n')
        if len(vcfs) == 0:
            try:
                os.system('find . -iname "*.download_cnv.sh" -exec bash {} \\; 2>/dev.null')
            except:
                pass



        exit_flag = 0
        for sample_id, v in valid_family_mem.items():
            # get vcf file
            vcf = os.popen(
                f'find {vcf_file_path} -type f -iname "*{sample_id}*.vcf.gz" -o -iname "*{sample_id}*.vcf"').read().split('\n')
            vcf = [_.strip() for _ in vcf if _.strip() and not _.rsplit('/', 1)[-1].startswith('new.')]
            if len(vcf) == 0:
                logger.warning(f'CNV file not downloaded yet, now downloading')
                os.system(f'bash {pw}/download.cnv.*.sh 2>/dev/null')
                time.sleep(2)
                vcf = os.popen(
                    f'find {vcf_file_path} -type f -iname "*{sample_id}*.vcf.gz" -o -iname "*{sample_id}*.vcf"').read().split('\n')
                vcf = [_.strip() for _ in vcf if _.strip() and not _.rsplit('/', 1)[-1].startswith('new.')]

            if len(vcf) == 0:
                logger.error(f'{v["type"]}  {sample_id} vcf file not found! ')
                exit_flag = 1
            elif len(vcf) == 1:
                v['vcf'] = vcf[0]
            else:

                vcf1 = [_ for _ in vcf if _.find('update') > -1]
                if len(vcf1) == 1:
                    v['vcf'] = vcf[0]
                    logger.warning(f'multiple vcf ({len(vcf)}) found, would use the file named as "updated" ')
                else:
                    logger.error(
                    f'multiple({len(vcf)}) VCF file for {v}, UDN={sample_id} under {vcf_file_path} found, please check \n\t{vcf} ')
                    exit_flag = 1
        if exit_flag:
            sys.exit(1)

    def build_hpo_dict(self, hpo_ref=None):
        """
        incase the hpo_ref.pkl is not present, would build it
        input file should have 2 columns,
        col1 = HPO ID, col2 = phenotype
        sep by tab
        """
        hpo_ref = hpo_ref or f'{pw_code}/hpo_ref.txt'
        try:
            with open(hpo_ref) as fp:
                data = [_.strip().split('\t') for _ in fp]
                data = {k.strip(): v.strip() for v, k in data}

            with open(f'{pw_code}/hpo_ref.pkl', 'wb') as out:
                pickle.dump(data, out)
            return data
        except:
            return 0

    def get_hpo_db(self):
        """
        get he hpo ref dict
        """
        hpo_pickle = f'{pw_code}/hpo_ref.pkl'

        if not os.path.exists(hpo_pickle):
            self.build_hpo_dict()
        else:
            hpo_db = pickle.load(open(hpo_pickle, 'rb'))
        return hpo_db

    def get_hpo_id(self) -> 'write to the {pw}/{intermediate_folder}/{prj}_terms_hpo.txt':
        """
        convert the HPO terms to HPO ID
        """
        pw = self.pw
        prj = self.prj
        logger = self.logger
        pheno = self.pheno_file

        fl_result = f'{pw}/{intermediate_folder}/{prj}.terms_hpo.txt'
        fl_result1 = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'

        tmp1 = os.popen(f'find {pw} -iname "{prj}.terms_hpo.txt"').read().strip().split('\n')
        tmp1 = [_.strip() for _ in tmp1 if _.strip()]

        tmp2 = os.popen(f'find {pw} -iname "{prj}.terms_pure_hpo.txt"').read().strip().split('\n')
        tmp2 = [_.strip() for _ in tmp2 if _.strip()]




        if len(tmp1) > 0:
            # logger.debug('HPO file already exist')
            if not os.path.exists(fl_result):
                os.system(f'ln -sf {tmp1[0]} {fl_result}')
            if len(tmp2) == 0:
                os.system(f'cut -f2 "{tmp1[0]}" > {fl_result1}')
            elif not os.path.exists(fl_result1):
                os.system(f'ln -sf {tmp2[0]} {fl_result1}')
            return 0


        if os.path.exists(fl_result):
            # logger.info(f'HOP id file already exists {fl_result}')
            return 0

        logger.info('Map HPO terms to IDs')
        res = []
        res_pure_hpo = []
        bad = []  # terms not exactly found in hpo db
        candidate = []

        # get the hpo_db
        hpo_db = self.get_hpo_db()

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
            fl_ambiguous = f'{pw}/{prj}_terms_ambiguous.txt'
            logger.debug(f'{len(candidate)} terms are ambiguous, please check {fl_ambiguous}')
            with open(fl_ambiguous, 'w') as out:
                print('UDN_symptom\tHPO_term\tHPO_ID', file=out)
                print('\n'.join(candidate), file=out)

        with open(fl_result, 'w') as out:
            print('\n'.join(res), file=out)
        with open(fl_result1, 'w') as out:
            print('\n'.join(res_pure_hpo), file=out)
        return 0

    def get_annot_col(self):
        """
        get the column index for the kept columns
        """
        logger = self.logger
        pw = self.pw
        logger = self.logger
        proband_id = self.proband_id
        lb = self.family[proband_id]['lb']

        for v in self.family.values():
            lb = v['lb']

            try:
                f_anno_exp = f'{pw}/{intermediate_folder}/{lb}.annotated.tsv'
                with open(f_anno_exp) as f:
                    anno_header = f.readline().strip()
                cols = anno_header.split('\t')
                break
            except:
                pass
        else:
            logger.error('Annotation file header not found, use default')
            cols = 'AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;Samples_ID;ID;REF;ALT;QUAL;FILTER;INFO;FORMAT;955567-UDN222476;Annotation_mode;CytoBand;Gene_name;Gene_count;Tx;Tx_start;Tx_end;Overlapped_tx_length;Overlapped_CDS_length;Overlapped_CDS_percent;Frameshift;Exon_count;Location;Location2;Dist_nearest_SS;Nearest_SS_type;Intersect_start;Intersect_end;RE_gene;P_gain_phen;P_gain_hpo;P_gain_source;P_gain_coord;P_loss_phen;P_loss_hpo;P_loss_source;P_loss_coord;P_ins_phen;P_ins_hpo;P_ins_source;P_ins_coord;P_snvindel_nb;P_snvindel_phen;B_gain_source;B_gain_coord;B_gain_AFmax;B_loss_source;B_loss_coord;B_loss_AFmax;B_ins_source;B_ins_coord;B_ins_AFmax;B_inv_source;B_inv_coord;B_inv_AFmax;TAD_coordinate;ENCODE_experiment;GC_content_left;GC_content_right;Repeat_coord_left;Repeat_type_left;Repeat_coord_right;Repeat_type_right;Gap_left;Gap_right;SegDup_left;SegDup_right;ENCODE_blacklist_left;ENCODE_blacklist_characteristics_left;ENCODE_blacklist_right;ENCODE_blacklist_characteristics_right;ACMG;HI;TS;DDD_HI_percent;DDD_status;DDD_mode;DDD_consequence;DDD_disease;DDD_pmid;ExAC_delZ;ExAC_dupZ;ExAC_cnvZ;ExAC_synZ;ExAC_misZ;GenCC_disease;GenCC_moi;GenCC_classification;GenCC_pmid;OMIM_ID;OMIM_phenotype;OMIM_inheritance;OMIM_morbid;OMIM_morbid_candidate;LOEUF_bin;GnomAD_pLI;ExAC_pLI;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class'.split(';')

        col_keep = {}
        cols = [_.lower() for _ in cols]

        for new, raw in zip(col_keep_new_name, col_keep_raw_name):
            try:
                col_keep[new] = cols.index(raw.lower())
            except:
                logger.error(f'column not found in annot file: {raw} (for {new} ) header = {cols}')
                sys.exit(1)
            # logger.info(f'{new} {raw} {cols.index(raw)}')
        # the actual "data" is the next column for "FORMAT"
        col_keep['data'] += 1
        self.cols = col_keep


    def annotate(self, sample_id) -> '{pw}/{intermediate_folder}/{lb}.annotated.txt':
        """
        annotate the VCF files using annotSV
        """
        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        vcf = self.family[sample_id]['vcf']

        f_anno_exp = f'{pw}/{intermediate_folder}/{lb}.annotated'
        if os.path.exists(f_anno_exp + '.tsv'):
            run_annot = 0
        else:
            run_annot = 1

        if not run_annot:
            logger.info(f'annotation already done: {lb}')
            self.get_annot_col()
            return 0

        # refine the vcf file
        ext = vcf.rsplit('.', 1)[-1]
        if ext == 'gz':
            zcat = 'gunzip -c '
        elif ext == 'vcf':
            zcat = 'cat '
        tmp =  vcf.rsplit('/', 1)
        if len(tmp) == 1:
            vcf_new = 'new.' + vcf
        else:
            vcf_new = tmp[0] + '/new.' + tmp[1]

        if not os.path.exists(vcf_new):
            cmd = f"""{zcat} {vcf}|awk '!match($3, /DRAGEN:REF:/)' |gzip > {vcf_new} """
            os.system(cmd)

        size_new = os.path.getsize(vcf_new)
        size_old = os.path.getsize(vcf)

        logger.info(f'AnnotSV: {sample_id}: vcf new ={vcf_new}, size old = {size_old}, size new={size_new}')

        # new options
        # -overlap , Minimum overlap (%) between user features (User BED) and the annotated SV to be reported Range values: [0-100], default = 100

        cmd = f'{self.annotsv_app} -genomeBuild GRCh38 -annotationMode split -outputFile {f_anno_exp} -SVinputFile {vcf_new} -overlap 40 -SVminSize 1 -benignAF 0.0001 >{pw}/{intermediate_folder}/{lb}.annotsv.log 2>&1'

        logger.debug(cmd)
        os.system(cmd)
        self.get_annot_col()

    def anno_filter(self, sample_id) ->'{pw}/{intermediate_folder}/{lb}.filtered.txt':
        """filter the raw annotSV result"""
        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        sv_caller = self.sv_caller
        type_short = self.family[sample_id]['type_short']
        col_keep = self.cols


        f_anno_exp = f'{pw}/{intermediate_folder}/{lb}.annotated.tsv'
        f_anno_filter = f'{pw}/{intermediate_folder}/{lb}.filtered.txt'

        # type_short = the type of this file, P/M/F/S

        # treat the annot_SV_ranking column diffently between proband and parents
        # for the proband, we demand the rank >=2
        # for the parents, the rank have no limit

        if os.path.exists(f_anno_filter) and os.path.getsize(f_anno_filter) > 1000:
            logger.info(f'filter file already exists: {lb}')
            return 0
        else:
            # col_filter = col_keep['filter'] + 1

            # the col_filter is disabled, because for dragen,  the SV with FILTER as "cnvLenght" can still be valid (because the thres is 10k bp, which is too long)

            col_quality = col_keep['QUAL'] + 1

            if type_short == 'P':  # proband
                # extra_filter = f' && ${col_gnomad_AF}<{thres_gnomad_maf}'
                extra_filter = ''
                if sv_caller == 'dragen':
                    extra_filter += f'${col_quality}>{thres_quality_score}'
                    # cmd = f"""head -1 {f_anno_exp} > {f_anno_filter};awk -F $'\\t'  '{extra_filter}' {f_anno_exp} >> {f_anno_filter} """
                    cmd = f"""awk -F $'\\t'  '{extra_filter}' {f_anno_exp} >> {f_anno_filter} """

                    logger.info(f'{lb}:  filter criteria = {extra_filter}')
                    logger.debug(cmd)
                else:
                    cmd = f'ln -s {f_anno_exp} {f_anno_filter}'
            else:
                extra_filter = ''
                cmd = f'ln -s {f_anno_exp} {f_anno_filter}'

            # cmd = f"""head -1 {f_anno_exp} > {f_anno_filter};awk -F $'\\t'  '${col_filter}=="PASS" {extra_filter}' {f_anno_exp} >> {f_anno_filter} """
            # logger.info(cmd)
            os.system(cmd)


    def get_copy_number_dragen(self, data, chr_, sex, idx_cn):
        """
        get the copy number from dragen vcf file
        """
        logger = self.logger
        # GT:SM:CN:BC:PE
        # SM:CN:BC:PE

        try:
            data = int(data.split(':')[idx_cn])  # CN
            if data == 0:
                copy_number = "-"
            elif data < 6:
                copy_number = 'o' * data
            else:
                copy_number = f'cn={data}'
        except:
            logger.warning(f'wrong copy number format: copy number={data}')
            return 'NA'


        if chr_.lower().find('chrx') > -1 and sex == 1:
            copy_number += 'y'
        elif chr_.lower().find('chry') > -1 and sex == 2 and copy_number != '-' and copy_number != 'NA':
            copy_number = 'invalid: ' + copy_number
        elif chr_.lower().find('chry') > -1 and sex == 1:
            copy_number = 'x' + copy_number

        return copy_number


    def get_copy_number_pacbio(self, data, sv_type, format_info):
        """
        format_info could be CN or GT:AD:DP, it would define how the data would be parased
        SV could be
    'DEL': 22449,
    'INS': 25805,
    'DUP': 2830,
    'BND': 78,
    'INV': 95,
    'cnv': 5

    format_info and data is like
     GT:AD:DP        0/1:5,9:14
        """
        solid_circle = chr(9679)
        logger = self.logger
        if format_info == 'CN':
            try:
                data = int(data)
            except:
                copy_number = 'NA'
            else:
                if data == 0:
                    copy_number = "-"
                else:
                    copy_number = 'o' * data
            return copy_number

        idx = {k.strip(): n for n, k in enumerate(format_info.split(':'))}
        # k = GT, AD, DP

        data = data.split(':')
        try:
            gt = data[idx['GT']]
        except:
            # invalid GT
            return 'NA'

        try:
            allele_depth =  '@' + data[idx['AD']]
        except:
            allele_depth = ''

        alleles = sorted(re.split(r'[/|]', gt))
        if alleles == ['0', '0']:
            return 'oo' + allele_depth

        if alleles == ['0', '1']:
            return solid_circle + 'o' + allele_depth

        if alleles == ['1', '1']:
            return solid_circle * 2 + allele_depth
        else:
            return gt + allele_depth

    def anno_extract(self, sample_id) ->'{pw}/{intermediate_folder}/{lb}.extracted.txt':
        prj = self.prj
        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        type_short = self.family[sample_id]['type_short']
        sex = self.family[sample_id]['sex']  # 1 or 2 for male and female
        col_keep = self.cols
        rel = self.family[sample_id]['type'].lower()
        sv_caller = self.sv_caller

        f_anno_filter = f'{pw}/{intermediate_folder}/{lb}.filtered.txt'
        f_anno_extract = f'{pw}/{intermediate_folder}/{lb}.extracted.txt'

        # validate existence
        if not os.path.exists(f_anno_filter):
            logger.error(f'filtered annot file not exist: {lb}')
            return 'err'

        if os.path.exists(f_anno_extract):
            lines = os.popen(f'cat {f_anno_extract}|wc -l').read().strip()
            lines = int(lines)
            logger.info(f'extracted file lines = {lines}')
            if lines > 0:
                logger.info(f'extracted anno file already exists: {lb}')
                return 0

        out = open(f_anno_extract, 'w')
        header = ["filter", "anno_id", "chr_", "pos_s", "pos_e", "sv_type", "qual", "exon_span_tag", "gn", "sv_len", "exon_span",
                "gain_af_max", "loss_af_max", 'gain_source', 'loss_source', "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
                "annot_ranking"]

        svtype1 = {}
        gene_list = set()


        with open(f_anno_filter) as fp:
            first_line = fp.readline()
            if first_line.find('AnnotSV') < 0:
                # return to file start
                fp.seek(0)
                header += ['copy_number']
            else:
                if sv_caller == 'dragen':
                    header += ['copy_number']
                else:
                    tmp = first_line.split('\t')[col_keep['data']: col_keep['data_end_flag']]
                    header_extra = []
                    for i in tmp:
                        m = re.match(r'.*(UDN\d+)', i)
                        if m:
                            udn_id = m.group(1)
                            try:
                                aff_state = self.family[udn_id]['aff_state']
                                rel_to_proband = self.family[udn_id]['type']
                                aff_state = {0: 'unaff', 1:'affected', -1: 'unknown'}[aff_state]
                            except:
                                aff_state = '@please_check_aff@'
                                rel_to_proband = '@check_rel_to_proband@'
                            header_extra.append(f'{rel_to_proband} {udn_id} {aff_state}')
                        else:
                            header_extra.append(i)

                    header += header_extra
                    logger.info(f'family_ID in header={header_extra}')
            print('\t'.join(header), file=out)

            for i in fp:
                a = i.strip().split('\t')
                anno_id, chr_, pos_s, pos_e, sv_len, sv_type, \
                _filter, qual, gn, exon_span, \
                gain_b_source, gain_b_af, \
                loss_b_source, loss_b_af, \
                ddd_mode, ddd_disease, omim, \
                phenotype, inheritance,annot_ranking = [a[col_keep[_]] for _ in col_keep_final]

                # data = expon_span_tag/FORMAT
                # sv_len = tx length (overlap of the SV and transcript)
                # exon_span = location

                # copy number
                # the the cols for the GT of each family member
                data_all = a[col_keep['data']: col_keep['data_end_flag']]
                format_info_raw = a[col_keep['format']]
                format_info = format_info_raw.split(':')
                # GT:SM:CN:BC:PE
                # SM:CN:BC:PE

                if sv_caller == 'dragen':
                    if 'CN' not in format_info and gn != "Gene name" :
                        logger.warning(f'CN not in FORMAT column: {format_info}, line={a}')
                        continue
                    idx_cn = format_info.index('CN')


                copy_number_all = []
                for data in data_all:
                    if sv_caller == 'dragen':
                        copy_number = self.get_copy_number_dragen(data, chr_, sex, idx_cn)
                    else:
                        copy_number = self.get_copy_number_pacbio(data, sv_type, format_info)

                    copy_number_all.append(copy_number)

                # sv type
                try:
                    sv_type = sv_type_convert[sv_type]
                except:
                    pass
                    # sv_type = 'NA'
                    # logger.warning(f'wrong sv_type format: sv_type={sv_type}  : {gn}  {anno_id}')

                # dgv gain
                try:
                    gain_af_max = float(gain_b_af) if gain_b_af else -1
                    if gain_af_max == -1:
                        gain_af_max = '-'
                    elif gain_af_max >= 0.01:
                        gain_af_max = f'{gain_af_max:.2f}'
                    else:
                        gain_af_max = f'{gain_af_max:.4g}'

                    # if gain_af_max != '-':
                    #     gain_af_max = f'{gain_af_max}; {gain_b_source}'
                except:
                    gain_af_max = 'NA'
                    logger.warning('wrong gain_af_max format: gain_af_max={gain_af_max}  : {gn}  {anno_id}')

                # dgv loss
                try:
                    loss_af_max = float(loss_b_af)  if loss_b_af else -1
                    if loss_af_max == -1:
                        loss_af_max = '-'
                    elif loss_af_max >= 0.01:
                        loss_af_max = f'{loss_af_max:.2f}'
                    else:
                        loss_af_max = f'{loss_af_max:.4g}'

                    # if loss_af_max != '-':
                    #     loss_af_max = f'{loss_af_max}; {loss_b_source}'
                except:
                    loss_af_max = 'NA'
                    logger.warning('wrong loss_af_max format: loss_af_max={loss_af_max}  : {gn}  {anno_id}')

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
                                exon_span_tag = 'covered multiple exons'
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

                # SV len
                sv_len1 = int(pos_e) - int(pos_s)
                try:
                    sv_len = abs(int(sv_len))
                except:
                    sv_len = sv_len1

                if sv_len < 50 :
                    try:
                        svtype1[sv_type] += 1
                    except:
                        svtype1[sv_type] = 1
                    if sv_type != 'BND' and exon_span_tag.find('exon') < 0:  # small indel, and not BND, not in exon
                        continue

                if sv_type == 'BND':
                    info = a[col_keep['info']]
                    info = {k: v for k, v in [_.split('=') for _ in info.split(';')]}
                    mate_id = info.get('MATEID') or ''
                    mate_dist = info.get('MATEDIST') or ''
                    sv_len = f'{mate_id};dist={mate_dist}'

                chr_ = 'chr' + chr_.lower()
                try:
                    if sv_len > 1000:
                        sv_len = f'{int(sv_len)/1000:.2f}kbp'
                    else:
                        sv_len = f'{sv_len} bp'
                except:
                    # logger.warning(f'wrong sv_len format: sv_len={sv_len}  : {gn}  {anno_id}')
                    # continue
                    pass

                if type_short == 'P':
                    gene_list.add(gn)

                # copy number
                if sv_caller == 'dragen':
                    copy_number_all = [f'{copy_number}@QUAL={qual}' for copy_number in copy_number_all]
                try:
                    copy_number_all = '\t'.join(copy_number_all)
                except:
                    logger.error(copy_number_all)
                    sys.exit(1)
                print('\t'.join([_filter, anno_id, chr_, pos_s, pos_e, sv_type, qual, exon_span_tag, gn,
                                sv_len, exon_span, gain_af_max, loss_af_max, gain_b_source, loss_b_source,
                                ddd_mode, ddd_disease, omim, phenotype, inheritance, annot_ranking, copy_number_all]), file=out)
        out.close()
        if len(svtype1) > 0:
            logger.info(f'SV with len < 50: count = {svtype1}')

        # write the genelist
        if type_short == 'P':
            n_genes = len(gene_list)
            if n_genes == 0:
                logger.error(f'no gene list found !')
                sys.exit(1)
            elif n_genes > 3000:
                logger.error(f'too many genes in the gene list ( {n_genes}) ,please check!')
                sys.exit(1)

            fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
            logger.info(f'gene list count = {n_genes}')
            with open(fn_genelist, 'w') as o:
                print('\n'.join(gene_list), file=o)

    def group_sv_into_bin(self, sample_id) -> 'family[sample_id]["sv_dict"]':
        # read the father and mother annotation result, build a dict
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
        # 1st layer = chr
        # the 2nd layer key = the bin# (100kb as a bin), is the  int(start/100000)  i.e. 100k bin,
        # this is to avoid the proband site iter over all the sites
        # however, the proband - father/mother site is not always matched,
        # on the proband search step, the each site would search the upper and lower bin, e.g. the proband site is 12,345,678,  the bin is 123, but, I would search both 122, 123, 124 bin in the father/mother dict
        # the list , starts, ends, they should be sorted.
        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        # sex = self.family[sample_id]['sex']  # 1 or 2 for male and female

        fn = f'{pw}/{intermediate_folder}/{lb}.extracted.txt'
        fnout = f'{pw}/{intermediate_folder}/{lb}.sv_to_bin.pkl'

        if not os.path.exists(fn):
            logger.error(f'extracted fields file for {lb} not found!\nexiting...')
            sys.exit(1)

        d_parent = {}
        n = 0
        with open(fn) as fp:
            header = fp.readline().strip().split('\t')
            idx = [header.index(_) for _ in ['chr_', 'pos_s', 'pos_e', 'sv_type', 'gn', 'copy_number']]
            fp.seek(0)

            for i in fp:
                n+=1
                a = i.strip().split('\t')

                try:
                    chr_, s, e, sv_type, gn, cn = [a[_] for _ in idx]
                except:
                    logger.error(fn)
                    logger.error(a)
                    logger.error(n)
                    sys.exit(1)
                chr_ = chr_.lower()
                try:
                    s = int(s)
                    e = int(e)
                except:
                    if s != 'pos_s':
                        logger.debug(f'wrong positon: {lb}: chr={chr_}, s={s}, e={e}')
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

        with open(fnout, 'wb') as out:
            pickle.dump(d_parent, out)

        return d_parent


    def check_amelie_missing(self):
        pw = self.pw
        prj = self.prj
        logger = self.logger

        fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
        fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'

        total_genes = line_count(fn_genelist)
        total_hpo = line_count(fn_hop_pure_id)

        if total_genes > 1000:
            logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), max input is 1000')

        if os.path.exists(f'{pw}/err_amelie.{prj}.genes_not_added.txt'):
            bad_genes = line_count(f'{pw}/err_amelie.{prj}.genes_not_added.txt')
            logger.warning(f'AMELIE missing genes: {bad_genes}/{total_genes} not included in AMELIE query')


    def query_amelie(self, force=False) -> 'build all the amelie result files':
        """
        run the amelie API
        """
        pw = self.pw
        prj = self.prj
        logger = self.logger

        fn = f'{pw}/{prj}.amelie.lite.txt'

        if os.path.exists(fn) and not force:
            try:
                self.amelie_dict = self.get_amelie_dict()
            except:
                logger.info(f'rerun amelie query')
                os.system(f'rm {fn} 2>/dev/null')
                force=True
            else:
                # self.check_amelie_missing()
                return 0
        logger.info('getting information from AMELIE')

        fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
        fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'
        fn_pheno = f'{pw}/pheno.keywords.txt'

        total_genes = line_count(fn_genelist)
        # total_hpo = line_count(fn_hop_pure_id)

        if total_genes > 1000:
            logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), would split the genelist')

        # main(prj, genelist, pheno_hpo_list, pw=None, pw_main=None, force=False)
        # amelie_api.main(prj, fn_genelist, fn_hop_pure_id, pw=pw, force=force)
        amelie_api.main_old(prj, fn_genelist, fn_hop_pure_id, pheno_for_match=fn_pheno,pw=pw, force=force)
        logger.info('amelie done')
        # self.check_amelie_missing()

        # # link the final amelie files
        # os.system(f'ln {pw}/{intermediate_folder}/{prj}.amelie.matched_query.txt {pw} 2>/dev/null')
        # os.system(f'ln {pw}/{intermediate_folder}/{prj}.amelie.lite.txt {pw} 2>/dev/null')

    def get_amelie_dict(self):
        """
        after querying the amelie API
        convert {prj}.amelie.lite.txt(only the final score of gene) to dict
        """
        pw = self.pw
        prj = self.prj
        logger = self.logger
        fn = f'{pw}/{prj}.amelie.lite.txt'

        try:
            d_amelie = [_.strip().split('\t') for _ in open(fn)]
            d_amelie = {k: v for k, v in d_amelie}
            logger.info(f'AMELIE count = {len(d_amelie)}')
        except:
            if not os.path.exists(fn):
                logger.error(f'amelie result not found: {fn}, would use dummy amelie score')
                fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
                return {gn.strip(): 0 for gn in open(fn_genelist) if gn.strip()}

        return d_amelie

    def get_family_sv_in_bin(self, sample_id, chr_, bin_):
        """
        sample_id  -> self.family[sample_id]['sv_dict']
        then, get all the sv_id in  bin_ -1, bin_, and bin_ + 1 list

        each element is [s, e, sv_type, gn, cn]
        """
        try:
            d_sv = self.family[sample_id]['sv_dict'][chr_]
        except:
            # the whole chr_ is not found in the corresponding dict
            logger = self.logger
            sample_type = self.family[sample_id]['type']
            logger.debug(f'{chr_} not found in sv dict for {sample_type}')
            return []

        res = []
        for i in [bin_ -1, bin_, bin_ +1]:
            if i < 0:
                continue
            if i not in d_sv:
                continue
            res.extend(d_sv[i])

        return res


    def run_proband_sv_match_pacbio(self) -> "{pw}/merged.sorted.{prj}.tsv":
        pw = self.pw
        proband_id = self.proband_id
        family = self.family
        prj = self.prj
        logger = self.logger
        lb = self.family[proband_id]['lb']
        fn_proband = f'{pw}/{intermediate_folder}/{lb}.extracted.txt'
        merged_table = f'{pw}/{intermediate_folder}/{prj}.merged.tsv'
        final = open(merged_table, 'w')
        proband_id = self.proband_id.lower()

        d_amelie = self.amelie_dict
        header_prefix = self.header_prefix

        with open(fn_proband) as f:
            header_raw = f.readline()
            tmp = header_raw.split('\t')
            idx_proband = [n for n, _ in enumerate(tmp) if _.lower().find(proband_id) > -1]
            idx = [tmp.index(_) for _ in ['chr_', 'pos_s', 'pos_e', 'sv_type', 'gn']]


            if len(idx_proband) == 1:
                idx_proband = idx_proband[0]
            else:
                logger.error(f'wrong header for {fn_proband}:\n{tmp}')
                sys.exit(1)

            header_prefix.append(header_raw)
            header = '\t'.join(header_prefix)
            print(header, file=final)
            for i in f:
                i = i.strip()
                a = i.split('\t')
                chr_, s, e, gn = [a[_] for _ in idx]
                chr_ = chr_.lower()
                proband_cn = a[idx_proband]

                if proband_cn.split('@')[0] == 'oo':
                    # the proband is wt
                    continue

                s = int(s)
                e = int(e)
                coord = f'{chr_}:{s}-{e}'

                # info stores the extra data need to be added to this line
                try:
                    amelie = d_amelie[gn]
                except:
                    logger.debug(f'gene not found in AMELIE list, gene = {gn}')
                    amelie = -99

                line = [coord, '', 'na', '', str(amelie), i]
                print('\t'.join(line), file=final)
        final.close()

        # sort the result by amelie score
        sorted_file = f'{pw}/merged.sorted.{prj}.tsv'
        sorted_excel = f'{pw}/merged.sorted.{prj}..xlsx'

        import pandas as pd
        df = pd.read_csv(merged_table, sep='\t')
        # df['match_strong'] =  df['chr_'] + ':' + df['pos_s'].astype(str) + "-" + df['pos_e'].astype(str)
        # df.drop(['match_strong', 'match_weak', 'match_count_udn'], axis=1, inplace=True)
        df.sort_values('AMELIE', ascending=False, inplace=True)
        df.to_csv(sorted_file, index=False, na_rep='', sep='\t')
        df.to_excel(sorted_excel, index=False, na_rep='')


    def run_proband_sv_match(self) -> "{pw}/merged.sorted.{prj}.tsv":
        """
        based on the proband extracted.txt,
        add the copy number/genename, sv_type, overlap_len from other family member
        the order is father, mother, brother, sister
        """
        pw = self.pw
        proband_id = self.proband_id
        family = self.family
        prj = self.prj
        logger = self.logger
        fn_proband = f'{pw}/{intermediate_folder}/{proband_id}_P.extracted.txt'
        merged_table = f'{pw}/{intermediate_folder}/{prj}.merged.tsv'
        final = open(merged_table, 'w')

        # logger.info(self.family[proband_id].keys())
        # logger.info(self.family[proband_id]['sv_dict'])

        d_amelie = self.amelie_dict

        # determine the order of non-proband samples
        trio_order = []
        for sample_id, v in family.items():
            t = v['type']
            if t == 'Proband':
                continue
            if t == 'Father':
                trio_order.append([sample_id, 1])
            elif t == 'Mother':
                trio_order.append([sample_id, 2])
            else:
                m = re.findall(r'\w+(\d+)$', t)
                if len(m) > 0:
                    sn = int(m[0])
                else:
                    sn = 0
                trio_order.append([sample_id, 10+sn])

        trio_order = sorted(trio_order, key=lambda x: x[1])
        trio_order = [_[0] for _ in trio_order]

        proband_id_lite = proband_id.replace('UDN', '')
        header_suffix_copy_number = [f'Proband {proband_id_lite}']
        header_suffix_gn = []
        header_suffix_sv_type = []
        header_suffix_overlap = []

        for i in trio_order:
            aff_state = family[i]['aff_state']
            aff_state = {0: 'unaff', 1:'affected', -1: 'unknown'}[aff_state]
            type_ = family[i]['type']
            sample_id_lite = i.replace('UDN', '')

            header_suffix_copy_number.append(f'{type_} ({aff_state}) {sample_id_lite}')

            header_suffix_gn.append(f'gn_{type_}')
            header_suffix_sv_type.append(f'sv_type_{type_}')
            header_suffix_overlap.append(f'overlap_{type_}')


        header_suffix = header_suffix_copy_number + header_suffix_gn + header_suffix_sv_type + header_suffix_overlap

        header_suffix = '\t'.join(header_suffix)

        with open(fn_proband) as fp:
            # header = next(fp).strip().replace('copy_number', '')
            header_raw = fp.readline().strip()
            header = header_raw.replace('copy_number', '')
            header_prefix = '\t'.join(self.header_prefix)
            header = f'{header_prefix}\t{header}{header_suffix}'

            header_list = header_raw.split('\t')
            idx = [header_list.index(_) for _ in ['chr_', 'pos_s', 'pos_e', 'gn']]
            print(header, file=final)

            family_sv_cache = {}
            bad_chr_all = {}
            # key1 = chr_
            # key2 = bin_
            # key3 = family id
            # v3 = list of all sv info in bin_-1, bin_, bin +1
            # sv_info = [s, e, sv_type, gn, cn]

            for i in fp:
                i = i.strip()
                a = i.split('\t')
                chr_, s, e, gn = [a[_] for _ in idx]
                chr_ = chr_.lower()

                s = int(s)
                e = int(e)
                bin_ = int(s / 100000)
                coord = f'{chr_}:{s}-{e}'

                # info stores the extra data need to be added to this line
                info = {}
                try:
                    amelie = d_amelie[gn]
                except:
                    logger.debug(f'gene not found in AMELIE list, gene = {gn}')
                    amelie = -99

                # get the same bin in the same family annotation results
                for sample_id in trio_order:
                    # logger.info(f'{sample_id}, {self.family[sample_id]["rel"]}')
                    # logger.info(self.family[sample_id].keys())
                    sex = family[sample_id]['sex']
                    sv_dict = self.family[sample_id]['sv_dict']
                    sample_type = self.family[sample_id]['type']
                    try:
                        bad_chr = bad_chr_all[sample_id]
                    except:
                        bad_chr = set()
                        bad_chr_all[sample_id] = bad_chr

                    if chr_ in bad_chr:
                        sv_in_bin = []
                    elif chr_ not in sv_dict:
                        if not (sex == 2 and chr_ == 'chry'):
                            bad_chr.add(chr_)
                            logger.warning(f'{chr_} not found in sv dict for {sample_type}')
                            sv_in_bin = []
                    else:
                        # use try except
                        try:
                            sv_in_bin = family_sv_cache[chr_][bin_][sample_id]
                        except:
                            sv_in_bin = self.get_family_sv_in_bin(sample_id, chr_, bin_)

                            try:
                                family_sv_cache[chr_][bin_][sample_id] = sv_in_bin
                            except:
                                try:
                                    family_sv_cache[chr_][bin_] = {sample_id: sv_in_bin}
                                except:
                                    family_sv_cache[chr_] = {bin_: {sample_id: sv_in_bin}}

                    # compare to find overlap
                    matched = ''
                    for isite in sv_in_bin:
                        # [s, e, sv_type, gn, cn]
                        s_q, e_q = isite[0], isite[1]

                        # no overlap
                        # if s_q > e or e_q < s:
                        #     pass
                        if s_q <= s <= e_q and e_q <= e:
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
                        info[sample_id] = new
                    else:  # the sv is not found in family, suppose to be normal
                        if chr_.lower().find('chrx') > -1 and sex == 1:
                            cn = 'oy'
                        elif chr_.lower().find('chry') > -1 and sex == 1:
                            cn = 'xo'
                        elif chr_.lower().find('chry') > -1 and sex == 2:
                            cn = 'na'
                        else:
                            cn = 'oo'
                        info[sample_id] = [cn, 'not_found'] + [''] * 4

                # if self.udn_match:  # a gene list provided by UDN gateway, usually not exist or not used
                #     i_udn_match = self.udn_match.get(gn) or 'na'
                # else:
                #     i_udn_match = 'na'

                # if not found in both amelie and udn_match, would skip this gene
                # UDN match means that the gene in the UDN report from Baylor, force include these genes


                # if amelie == -99 and i_udn_match == 'na':
                #     # logger.info(f'{gn}:not found in both amelie and udn_match')
                #     continue
                # elif amelie == -99:  # not found in amelie but foudn in udn_match
                #     amelie = 999
                res = f'{coord}\tna\tna\t\t{amelie}\t{i}'
                res_suffix_cn = []
                res_suffix_gn = []
                res_suffix_sv_type = []
                res_suffix_overlap = []

                for sample_id in trio_order:
                    res_suffix_cn.append(info[sample_id][0])
                    res_suffix_gn.append(info[sample_id][1])
                    res_suffix_sv_type.append(info[sample_id][2])
                    res_suffix_overlap.append(info[sample_id][-1])

                res_suffix = res_suffix_cn + res_suffix_gn + res_suffix_sv_type + res_suffix_overlap

                res = [res] + res_suffix
                res = '\t'.join(res)
                print(res, file=final)
        final.close()

        # sort the result by amelie score
        sorted_file = f'{pw}/merged.sorted.{prj}.tsv'
        sorted_excel = f'{pw}/merged.sorted.{prj}.xlsx'

        import pandas as pd
        df = pd.read_csv(merged_table, sep='\t')

        # try:
        #     df['coord'] =  df['chr_'] + ':' + df['pos_s'].astype(str) + "-" + df['pos_e'].astype(str)
        # except:
        #     print(df.head())
        #     sys.exit(1)
        # df.drop(['match_strong', 'match_weak', 'match_count_udn'], axis=1, inplace=True)
        df.sort_values('AMELIE', ascending=False, inplace=True)
        df.to_csv(sorted_file, index=False, na_rep='', sep='\t')
        df.to_excel(sorted_excel, index=False, na_rep='')

def get_omim_map(logger):
    """
    input = genemap2.txt
    return = {pheno_id1: pheno_desc1}
    """
    if os.path.exists(f'{pw_code}/omim_map.pkl'):
        try:
            d_omim_map = pickle.load(open(f'{pw_code}/omim_map.pkl', 'rb'))
        except:
            return None
        else:
            return d_omim_map
    elif os.path.exists(f'{pw_code}/genemap2.txt'):
        d_omim_map = {}
        with open(f'{pw_code}/genemap2.txt') as fp:
            for i in fp:
                if i[0] == '#':
                    continue

                a = i.split('\t')
                gn_id = a[5]
                pheno = a[12].strip()
                if not pheno:
                    d_omim_map[gn_id] = {}
                else:
                    pheno = pheno.split(';')
                    d_omim_map[gn_id] = {}
                    for ipheno in pheno:
                        ipheno = ipheno.replace('Autosomal recessive', 'AR').replace('Autosomal dominant', 'AD').replace('X-linked recessive', 'XLR').replace('X-linked dominant', 'XLD')
                        m = re.match(r".*(?:, )?(\d{6})(?: \(\d+\))", ipheno)
                        if m:
                            d_omim_map[gn_id][int(m.group(1))] = ipheno.strip()
                        else:
                            logger.warning(f'wrong OMIM phenotype format: ipheno={ipheno},gn={a[8]}')
        with open(f'{pw_code}/omim_map.pkl', 'wb') as out:
            pickle.dump(d_omim_map, out)
        return d_omim_map
    return None


global_flag = {'omim_blocked': 0}


def get_driver(driver, logger):
    try:
        driver.current_url
    except:
        from selenium import webdriver
        from selenium.webdriver.firefox.options import Options
        options = Options()
        options.headless = True
        platform = sys.platform
        if platform == 'darwin':
            exe_firefox = '/Users/files/work/package/firefox/geckodriver'
        else:
            exe_firefox = f'/home/chenh19/tools/geckodriver'
        driver = webdriver.Firefox(options=options, executable_path=exe_firefox)
        logger.warning('firefox initiated')
        driver.get('https://www.google.com/')
        driver.save_screenshot('firefox_init.png')
    return driver


def run_selenium(driver, url_omim, gn, gn_omim_id, logger):

    driver = get_driver(driver, logger)
    try:
        driver.current_url
    except:
        logger.error(f'fail to load selenium')
        return 0, 0

    gn = re.sub(r'\W+', '_', gn)
    driver.get(url_omim)

    # skip the donation
    try:
        cancel_donation = wait(driver, 10).until(lambda _: _.find_element_by_xpath('//button[@id="donationPopupCancel"]'))
        cancel_donation.click()
    except:
        pass

    html = driver.page_source
    if html.find('Phenotype-Gene') < 0:
        logger.warning(f'seems an empty html returned: len={len(html)}')
        driver.save_screenshot(f'OMIM_{gn}_{gn_omim_id}.png')

    return driver, html



def gn_to_pheno(driver, gn, gene_id, logger):
    """
    return
    None - error, page not loaded
    'empty',  no linked phenotype for this gene
    list, element = [pheno_id, pheno_name, merged_title]
    """

    gene_id = str(int(float(gene_id)))
    url = f'https://www.omim.org/entry/{gene_id}'
    try:
        driver.get(url)
    except:
        driver = get_driver(None, logger)
        driver.get(url)
    try:
        # test if the page is ready
        wait(driver, 5).until(lambda _: _.find_element_by_xpath('//a[@id="preferredTitle"]'))
    except:
        try:
            if driver.find_element_by_xpath('//p[contains(text(), "has been blocked because it was identified as a crawler")]'):
                logger.error('fail to scrapy because blocked by OMIM website')
                return driver, 'blocked'
        except:
            pass

        logger.warning(f'page not loaded for {gn}, omim gene id = {gene_id}')
        return driver, None

    try:
        ele = driver.find_element_by_xpath('//div[./span/strong[text()="Gene-Phenotype Relationships"]]')
    except:
        return driver, 'empty'

    try:
        tb = ele.find_element_by_xpath('./following-sibling::div/table/tbody')
        rows = tb.find_elements_by_tag_name('tr')
        rows = [[cell.get_attribute('innerText') for cell in _.find_elements_by_tag_name('td')] for _ in rows]
    except:
        logger.error(f'gene-pheno section found, but fail to get table: {gn} - {gene_id}')
        raise
        # return driver,  None

    res = []
    dedup = set()
    for i in rows:
        # ['6q25.3', '?ACAT2 deficiency', '614055', 'IC', '1']
        # location, pheno, pheno_id, inher, pheno_mapping_key
        # the row can be 5 columns or 4 columns
        try:
            pheno, pheno_id, inher, k = i[-4:]
        except:
            logger.warning(f'wrong gene-pheno map row gene={gn}, gene_id = {gene_id}, row =  {i}')
            continue
        # like Intellectual developmental disorder and retinitis pigmentosa, 618195 (3), AR
        if not pheno.strip():
            continue
        title = f'{pheno}, {pheno_id} ({k}), {inher}'
        if pheno_id not in dedup:
            res.append([pheno, pheno_id, title])
            dedup.add(pheno_id)
    if len(res) == 0:
        return driver, 'empty'

    logger.info(f'\t    pheno count: nrow = {len(rows)}, dedup = {len(res)}')
    return driver, res


def get_pheno_detail(driver, pheno_id, pheno_title, logger):
    """
    return
    None, there is error to get the description
    str,  the whole description
    """
    url = f'https://www.omim.org/entry/{pheno_id}'
    try:
        driver.get(url)
    except:
        driver = get_driver(None, logger)
        driver.get(url)
    try:
        # test if the page is ready
        wait(driver, 5).until(lambda _: _.find_element_by_xpath('//a[@id="preferredTitle"]'))
    except:
        logger.warning(f'page not loaded for {pheno_title}, pheno ID = {pheno_id}')
        return driver, None

    try:
        desc = driver.find_element_by_xpath('//div[@id="descriptionFold"]').get_attribute('innerText')
        if desc.strip():
            desc = 'Descriptioin\n' + desc
    except:
        desc =  ''

    try:
        clin = driver.find_element_by_xpath('//div[@id="clinicalFeaturesFold"]').get_attribute('innerText')
        if clin.strip():
            clin = 'Clinical\n' + clin
    except:
        clin =  ''

    try:
        mol = driver.find_element_by_xpath('//div[@id="molecularGeneticsFold"]').get_attribute('innerText')
        if mol.strip():
            mol = 'Molecular\n' + mol
    except:
        mol =  ''

    res = '\n'.join([desc, clin, mol])

    return driver, res

def cleanup_memory(logger, driver):
    try:
        logger.debug('now killing the chrome process')
    except:
        print('\tnow killing the chrome process')
    try:
        driver.close()
        driver.quit()
    except:
        pass
    os.system("""ps -ef|grep -i "chrome"|awk '{print $2}'|sort|uniq|parallel 'kill -9 {} 2>/dev/null' """)
    os.system("""ps -ef|grep -i "firefox"|awk '{print $2}'|sort|uniq|parallel 'kill -9 {} 2>/dev/null' """)



def run_omim_scrapy(gn, gene_id, logger, driver=None):
    start = time.time()
    v = [gene_id, []]
    gene_id = str(int(float(gene_id)))
    driver, tmp = gn_to_pheno(driver, gn, gene_id, logger)
    time_sleep = 5
    time.sleep(time_sleep)

    if tmp == 'blocked':
        end = time.time()
        dure = int(end - start)
        cleanup_memory(logger, driver)
        logger.error(f'IP blocked by OMIM: gn = {gn} total duration = {dure}')
        sys.exit(1)

    if tmp is None:
        v[1] = 'err'
        logger.warning(f'{gn} - {gene_id}:  fail to get pheno')
        return driver, None
    elif tmp == 'empty':
        # no linked pheno
        logger.debug(f'{gn} - {gene_id}:  no linked pheno')
        v[1] = 'no_linked_pheno'
    else:
        n_pheno_tmp = 0
        for pheno_info in tmp:
            n_pheno_tmp += 1
            pheno, pheno_id, title = pheno_info
            if not pheno_id.strip():
                logger.debug(f'\tempty pheno ID: gn = {gn}, - {gene_id} pheno = {title}')
                ires = [pheno, pheno_id, title, '']
                v[1].append(ires)
                continue

            driver, desc = get_pheno_detail(driver, pheno_id, title, logger)
            time.sleep(5)

            if desc is not None:
                ires = [pheno, pheno_id, title, desc]
                try:
                    v[1].append(ires)
                except:
                    logger.error(f'error, gn = {gn}, v = {v} fail to append')
                    raise
            else:
                v[1] = 'err'
                logger.info(f'\t{gn} - {gene_id}fail get pheno detail')
                break
        else:
            logger.debug(f'{gn} - {gene_id}: ok, n pheno = {len(v[1])}')

    return driver, v


def run_omim_scrapy_old(gn, gn_omim_id, logger, res_prev=None, driver=None):
    """
    get the omim disease discription
    return key=pheno_id, value=the description and clinicalfeatures text
    """
    headers = {'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.67 Safari/537.36 Edg/87.0.664.52', 'sec-fetch-site': 'cross-site', 'origin': 'https://www.omim.org', 'referer': 'https://fonts.googleapis.com/', 'accept': '*/*', 'accept-encoding': 'gzip, deflate, br'}

    # gene to phenotype
    # first try the d_omim_map, key=gene id, v = omim phenotype list
    d_omim_map = get_omim_map(logger)
    gene_page_scrapy = 1
    omim_id_list = set(gn_omim_id.split(';'))
    if len(omim_id_list) > 1:
        logger.warning(f'multiple OMIM id find, gn={gn}, omim_id={gn_omim_id}')
        res = {}
        for omim_id_tmp  in omim_id_list:
            res_tmp, driver = run_omim_scrapy_old(gn, omim_id_tmp.strip().split('.', 1)[0], logger, res_prev=res, driver=driver)
            if res_tmp:
                res.update(res_tmp)
        return res, driver
    else:
        gn_omim_id = list(omim_id_list)[0].strip().split('.', 1)[0]

    # if gn == 'PKHD1':
    #     logger.error(f'gn=PKHD1, omim_id_list={omim_id_list}, gn_omim_id={gn_omim_id}')

    url_omim = f'https://www.omim.org/entry/{gn_omim_id}'
    if d_omim_map:
        try:
            pheno_list = d_omim_map[gn_omim_id]
            if len(pheno_list) > 1:
                gene_page_scrapy = 0
        except:
            logger.warning(f'gene not found in d_omim_map dict: {gn}  gn_omim_id={gn_omim_id} url={url_omim}')

    if gene_page_scrapy:
        try:
            r = requests.request('GET', url_omim, headers=headers).text
        except:
            logger.warning(f'getting phenon list: OMIM blocked (request), would try selenium: {gn}')
            global_flag['omim_blocked'] = 1

        if global_flag['omim_blocked']:
            driver, r = run_selenium(driver, url_omim, gn, gn_omim_id, logger)
            if not r:
                logger.warning(f'fail to run selenium: {gn}')
                return 1
        r = bs(r, features='lxml')

        gn_web = r.find('em', text=re.compile('HGNC Approved Gene Symbol'))

        try:
            gn_web = gn_web.text.rsplit(' ', 1)[-1]
        except:
            logger.warning(f'gene name not found on website, OMIM_id={gn_omim_id}, gn={gn}, html={gn_web}')
            return 0

        if gn_web != gn:
            logger.warning(f'gene name not match, excel={gn}, web={gn_web}, omim_id={gn_omim_id}')
            return 0

        try:
            rows = r.find('table').findAll('tr')
        except:
            logger.warning(f'no phenotype table found, gn={gn}')
            return res_prev if res_prev else 0, driver
        rows = [_.findAll('td') for _ in rows]
        rows = [_ for _ in rows if _]

        pheno_list = {}
        start = 1
        for i in rows:
            if start:
                idx_pheno_id = 2
                idx_phenotype = 1
                start = 0
            else:
                idx_pheno_id = 1
                idx_phenotype = 0
            try:
                pheno_id = i[idx_pheno_id].text.strip()
            except:
                tmp = [_.text.strip() for _ in i]
                logger.warning('wrong table format: gn={gn}, pheno row={tmp}')
            else:
                info = [_.text.strip() for _ in i[idx_phenotype:]]
                pheno_desc = f'{info[0]}, {info[1]} ({info[3]}) {", " + info[2] if info[2] else ""}'

                try:
                    pheno_list[int(pheno_id)] = pheno_desc
                except:
                    logger.warning(f'wrong OMIM phenotype ID: {pheno_id}, info={info}')
        logger.info(f'pheno list for {gn}: {pheno_list}')

    # scrapy the pheno page
    # get the descriptioin and clinical features
    res = {}
    for pheno_id, pheno_desc in pheno_list.items():
        url_omim = f'https://www.omim.org/entry/{pheno_id}'

        if global_flag['omim_blocked']:
            logger.info(f'get OMIM using selenium')
            driver, r = run_selenium(driver, url_omim, gn, gn_omim_id, logger)
            if not r:
                logger.warning(f'fail to run selenium: {gn}')
                return 1
        else:
            try:
                r = requests.request('GET', url_omim, headers=headers).text
            except:
                logger.warning(f'OMIM blocked (request), would try selenium: {gn}')
                global_flag['omim_blocked'] = 1
                logger.info(f'get OMIM using selenium')
                driver, r = run_selenium(driver, url_omim, gn, gn_omim_id, logger)
                if not r:
                    logger.warning(f'fail to run selenium: {gn}')
                    return 1

        r = bs(r, features='lxml')
        try:
            desc = r.find('div', attrs={'id': 'descriptionFold'})
            desc = '\n'.join([_.text.strip() for _ in desc.find_all('p')])
        except:
            desc = ''

        try:
            clin = r.find('div', attrs={'id': 'clinicalFeaturesFold'})
            clin = '\n'.join([_.text.strip() for _ in clin.find_all('p')])
        except:
            clin = ''
        tmp = f'\n\n**{pheno_desc}**\n' + '\n'.join([desc, clin]).strip()
        tmp = tmp.strip()
        res[pheno_id] = tmp
    tmp = res_prev.update(res) if res_prev else res

    return tmp, driver
