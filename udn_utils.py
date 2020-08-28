"""
2020-08-18
1. add the UDN gateway API, allow to download the VCF file/ phenotip
2. add the encrypt / decrypt module for the API token

to-do
add the downloading data / uploading data to AWS server

"""

import os
import sys
import re
import yaml
import logging
import logging.config

from . import amelie_api

# basic data
pw_code = os.path.dirname(os.path.realpath(__file__))


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
def get_logger(pw, prj):
    fl_log_conf = f'{pw_code}/logging_setting.yaml'
    log_prefix = prj
    log_prefix = pw + '/' + log_prefix + '_' if log_prefix else ''
    with open(fl_log_conf) as f:
        cfg = yaml.safe_load(f.read())
        cfg['handlers']['console']['level'] = 'INFO'
        cfg['handlers']['file']['mode'] = 'w'
        cfg['handlers']['file']['filename'] = log_prefix + cfg['handlers']['file']['filename']
        cfg['handlers']['error']['filename'] = log_prefix + cfg['handlers']['error']['filename']
        logging.config.dictConfig(cfg)
    return logging.getLogger('main')

def line_count(fn):
    n = 0
    for _ in open(fn):
        n += 1
    return n


intermediate_folder = 'intermediate'

class UDN_case():
    def __init__(self, config_file):
        cfg = yaml.safe_load(open(config_file).read())
        self._cfg = cfg
        self.prj = cfg['prj']
        self.pw = cfg['path']
        self.proband_id = cfg['IDs']['proband_id']
        logger = get_logger(self.pw, self.prj)
        self.logger = logger

        # verify the config
        cfg_tested = self.verify_config()
        self.path = cfg_tested['path']
        self.pheno_file = cfg_tested['pheno_file']
        self.pheno_file_select = cfg_tested['pheno_file_select']
        self.family = cfg_tested['family']
        self.vcf_file_path = cfg_tested['vcf_file_path']
        self.thres_annot_sv_ranking = cfg_tested['thres_annot_sv_ranking']
        self.thres_gnomad_maf = cfg_tested['thres_gnomad_maf']
        self.amelie_dict = {}

        # copy the config file
        os.system(f'cp {config_file} {self.pw}/{self.prj}_settings.yaml 2>/dev/null')
        os.system(f'mkdir -p {self.pw}/{intermediate_folder} 2>/dev/null')

        # update self.family, add the vcf_file and run_annot
        self.parse_pw()

        # convert hpo terms to hpo id
        self.get_hpo_id()

        # get annotSV path
        self.annotsv_app = self.get_annotsv()
        if not self.annotsv_app:
            if self.family[self.proband_id]['run_annot']:
                logger.error('No annotSV specified, exiting...')
                sys.exit(1)

        # get the colum number for the kept fields
        # the dict was stored in self.cols
        self.cols = self.get_annot_col()

        # verfify if the automatic process is already done
        self.done_phase1 = 0
        self.done_phase2 = 0
        self.done_phase3 = 0
        fn_final = f'{self.pw}/{self.prj}.merged.sorted.tsv'
        if os.path.exists(fn_final):
            f_size = os.path.getsize(fn_final)
            if f_size < 10000:
                logger.warning(f'the merged anno result already exist: {fn_final}, but the size is only {f_size/1000:.2f}')
            else:
                self.done_phase1 = 1

        # verify if the selected gene list is done
        fn_selected = f'{self.pw}/{self.prj}.selected.genes.txt'
        if os.path.exists(fn_selected):
            f_size = os.path.getsize(fn_selected)
            if f_size < 200:
                logger.warning(f'selected gene list already exist: {fn_selected}, but the size is only {f_size/1000:.2f}')
            else:
                self.done_phase2 = 1

        # verify if the final report is done
        fn_report = f'{self.pw}/{self.prj}.report.xlsx'
        if os.path.exists(fn_report):
            f_size = os.path.getsize(fn_report)
            if f_size < 200:
                logger.warning(f'report already exist {fn_report}, but the size is only {f_size/1000:.2f}')
            else:
                self.done_phase3 = 1


    def verify_config(self):
        cfg = self._cfg
        logger = self.logger
        path = cfg['path']
        pw = path
        # prj = cfg['prj']
        pheno_file = cfg['pheno_file'].strip()
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
        pheno_file_select = pheno_file_select or pheno_file


        proband_id = cfg['IDs']['proband_id']
        proband_sex = cfg['IDs']['proband_sex']

        if proband_sex not in [1, 2]:
            logger.error('Wrong proband sex, valid=1 / 2, where 1=male, 2=female')
            sys.exit(1)

        thres_annot_sv_ranking = cfg['default']['thres_annot_sv_ranking']   # default is 2
        thres_gnomad_maf = cfg['default']['thres_gnomad_maf']   # default is 2
        vcf_file_path = cfg['default']['vcf_file_path']
        # default is the same as path
        vcf_file_path = vcf_file_path or path

        family = {}
        family[proband_id] = {
                'lb': f'{proband_id}_P',
                'type': 'Proband',
                'sex': proband_sex,
                'type_short': 'P',
                'vcf': None,
                'run_annot': 1,
                'aff_state': 1,
        }
        for k, v in cfg['IDs'].items():
            if not v:
                continue
            if k.startswith('proband'):
                continue

            try:
                sample_id, aff_state = v
            except:
                logger.error(f'Wrong value specified in IDs scope of config file, key={k}, value={v}.  value shoudl have 2 elements')
                sys.exit(1)

            if aff_state not in [0, 1]:
                logger.error('Wrong affect state sepcified, key={k}, aff_state = {aff_state}.  valid = 0 / 1')
                sys.exit(1)

            sample_id = sample_id.upper()
            m = re.match(r'(UDN)?\d+$', sample_id)
            if not m:
                logger.error(f'Wrong UDN ID for {k} : value={sample_id}')
                sys.exit(1)


            if k == 'father':
                family[sample_id] = {
                'lb': f'{sample_id}_F',
                'type': 'Father',
                'sex': 1,
                'type_short': 'F',
                'vcf': None,
                'run_annot': 1,
                'aff_state': aff_state,
                }
            elif k == 'mother':
                family[sample_id] = {
                'lb': f'{sample_id}_M',
                'type': 'Mother',
                'sex': 2,
                'type_short': 'M',
                'vcf': None,
                'run_annot': 1,
                'aff_state': aff_state,
                }
            elif k.startswith('brother'):
                family[sample_id] = {
                'lb': f'{sample_id}_S',
                'type': k.title(),
                'sex': 1,
                'type_short': 'S',
                'vcf': None,
                'run_annot': 1,
                'aff_state': aff_state,
                }
            elif k.startswith('sister'):
                family[sample_id] = {
                'lb': f'{sample_id}_S',
                'type': k.title(),
                'sex': 2,
                'type_short': 'S',
                'vcf': None,
                'run_annot': 1,
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
               'vcf_file_path': vcf_file_path,
               'thres_annot_sv_ranking': thres_annot_sv_ranking,
               'thres_gnomad_maf': thres_gnomad_maf,
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
        if machine_platform == 'Darwin':
            home = os.path.expanduser('~')
            dock_path = f'{home}/dock/annotsv.sif'
        else:
            dock_path = '/scratch/cqs/chenh19/dock/annotsv.sif'
        # the annot_sv_rank and gnomad_AF filtering would only take effect on the proband, would not filter on the parent

        if os.path.exists(dock_path):
            return f'singularity exec {dock_path}  /tools/AnnotSV/bin/AnnotSV'
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

        for sample_id, v in valid_family_mem.items():
            # get vcf file
            vcf = os.popen(
                f'find {vcf_file_path} -type f -iname "*{sample_id}*.vcf.gz" -o -iname "*{sample_id}*.vcf"').read().split('\n')
            vcf = [_.strip() for _ in vcf if _.strip()]
            if len(vcf) == 0:
                logger.error(f'{v["type"]}  {sample_id} vcf file not found! ')
                sys.exit(1)
            elif len(vcf) == 1:
                v['vcf'] = vcf[0]
            else:
                logger.error(
                    f'multiple({len(vcf)}) VCF file for {v["lb"]}  {sample_id} under {vcf_file_path} found, please check ')
                sys.exit(1)

            # get annotated file
            anno_exp = f'{pw}/{intermediate_folder}/{v["lb"]}.annotated.tsv'
            if os.path.exists(anno_exp):
                v['run_annot'] = 0

    def build_hpo_dict(self, hpo_ref=None):
        """
        incase the hpo_ref.pkl is not present, would build it
        input file should have 2 columns,
        col1 = HPO ID, col2 = phenotype
        sep by tab
        """
        import pickle
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
            import pickle
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
        header = open(f'{pw_code}/annotsv.header.txt').read().strip()
        cols = header.split('\t')

        col_keep = {}

        for new, raw in zip(col_keep_new_name, col_keep_raw_name):
            try:
                col_keep[new] = cols.index(raw)
            except:
                logger.error('column not found in annot file: {raw} (for {new} ) ')
                sys.exit(1)
        # the actual "data" is the next column for "FORMAT"
        col_keep['data'] += 1
        return col_keep

    def annotate(self, sample_id) -> '{pw}/{intermediate_folder}/{lb}.annotated.txt':
        """
        annotate the VCF files using annotSV
        """
        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        vcf = self.family[sample_id]['vcf']
        run_annot = self.family[sample_id]['run_annot']


        if not run_annot:
            logger.info(f'annotation already done: {lb}')
            return 0

        f_anno_exp = f'{pw}/{intermediate_folder}/{lb}.annotated'

        cmd = f'{self.annotsv_app} -genomeBuild GRCh38 -typeOfAnnotation split -outputFile {f_anno_exp} -SVinputFile {vcf} >{pw}/{intermediate_folder}/{lb}.annotsv.log 2>&1'
        logger.debug(cmd)
        os.system(cmd)

    def anno_filter(self, sample_id) ->'{pw}/{intermediate_folder}/{lb}.filtered.txt':
        """filter the raw annotSV result"""
        pw = self.pw
        prj = self.prj
        logger = self.logger
        lb = self.family[sample_id]['lb']
        type_short = self.family[sample_id]['type_short']
        thres_gnomad_maf = self.thres_gnomad_maf
        thres_annot_sv_ranking = self.thres_annot_sv_ranking
        col_keep = self.cols

        f_anno_exp = f'{pw}/{intermediate_folder}/{lb}.annotated.tsv'
        f_anno_filter = f'{pw}/{intermediate_folder}/{lb}.filtered.txt'

        # type_ = the type of this file, P/M/F/S

        # ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "QUAL", "data", "gn", "sv_len", "exon_span",
        #              "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
        #              "annot_ranking"]
        # treat the annot_SV_ranking column diffently between proband and parents
        # for the proband, we demand the rank >=2
        # for the parents, the rank have no limit

        if os.path.exists(f_anno_filter) and os.path.getsize(f_anno_filter) > 1000:
            logger.info(f'filter file already exists: {lb}')
            return 0
        else:
            col_gnomad_AF = col_keep['af_gnomad'] + 1
            col_annot_sv_ranking = col_keep['annot_ranking'] + 1
            col_filter = col_keep['filter'] + 1
            col_gn = col_keep['gn'] + 1

            if type_short == 'P':  # proband
                extra_filter = f' && ${col_gnomad_AF}<{thres_gnomad_maf} && ${col_annot_sv_ranking}>={thres_annot_sv_ranking}'
            else:
                extra_filter = ''

            logger.info(f'{lb}:  filter criteria = FILTER==PASS {extra_filter}')

            cmd = f"""head -1 {f_anno_exp} > {f_anno_filter};awk -F $'\\t'  '${col_filter}=="PASS" {extra_filter}' {f_anno_exp} >> {f_anno_filter} """
            # logger.info(cmd)
            os.system(cmd)

        # if proband, extract the genelist
        if type_short == 'P':
            cmd = f"""cut -d $'\\t' -f {col_gn} {f_anno_filter}|sed '1d' |sort|uniq > {pw}/{intermediate_folder}/{prj}.genelist"""
            # logger.info(cmd)
            os.system(cmd)
            lines = os.popen(f'wc -l < {pw}/{intermediate_folder}/{prj}.genelist').read().strip()
            try:
                lines = int(lines)
            except:
                logger.error(f'fail to get the line number of {prj}.genelist')
                sys.exit(1)

            if lines == 0:
                logger.error(f'{lb}:  no gene was found')
                sys.exit(1)

    def anno_extract(self, sample_id) ->'{pw}/{intermediate_folder}/{lb}.extracted.txt':

        pw = self.pw
        logger = self.logger
        lb = self.family[sample_id]['lb']
        sex = self.family[sample_id]['sex']  # 1 or 2 for male and female
        col_keep = self.cols

        f_anno_filter = f'{pw}/{intermediate_folder}/{lb}.filtered.txt'
        f_anno_extract = f'{pw}/{intermediate_folder}/{lb}.extracted.txt'

        # validate existence
        if not os.path.exists(f_anno_filter):
            logger.error(f'filtered annot file not exist: {lb}')
            return 'err'

        if os.path.exists(f_anno_extract):
            logger.info(f'extracted anno file already exists: {lb}')
            return 0
        else:
            out = open(f_anno_extract, 'w')
            print(
                '\t'.join(
                    ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "qual", "exon_span_tag", "gn", "sv_len", "exon_span",
                    "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
                    "annot_ranking", 'copy_number']),
                file=out)

            with open(f_anno_filter) as fp:
                first_line = fp.readline()
                if first_line.find('AnnotSV') < 0:
                    # return to file start
                    fp.seek(0)
                for i in fp:
                    a = i.strip().split('\t')
                    anno_id, chr_, pos_s, pos_e, sv_type, _filter, qual, data, gn, \
                        sv_len, exon_span, af_dgv_gain, af_dgv_loss, af_gnomad, \
                        ddd_mode, ddd_disease, omim, phenotype, inheritance, annot_ranking = [a[col_keep[_]] for _ in col_keep_new_name]

                    # data = expon_span_tag/FORMAT
                    # sv_len = tx length (overlap of the SV and transcript)
                    # exon_span = location

                    chr_ = 'chr' + chr_.lower()
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
                            af_gnomad = "-"
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
        import pickle
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
        with open(fn) as fp:
            for i in fp:
                a = i.strip().split('\t')
                # ["anno_id", "chr_", "pos_s", "pos_e", "sv_type", "qual", "exon_span_tag", "gn", "sv_len", "exon_span",
                #     "af_dgv_gain", "af_dgv_loss", "af_gnomad", "ddd_mode", "ddd_disease", "omim", "phenotype", "inheritance",
                #     "annot_ranking", 'copy_number']
                chr_, s, e, sv_type, gn, cn = [a[_] for _ in [1, 2, 3, 4, 7, -1]]
                chr_ = chr_.lower()
                try:
                    s = int(s)
                    e = int(e)
                except:
                    if s !='pos_s':
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

        if os.path.exists(f'{pw}/err_amelie.{prj}.HPO_not_added.txt'):
            bad_hpo = line_count(f'{pw}/err_amelie.{prj}.HPO_not_added.txt')
            logger.warning(f'AMELIE missing HPOs: {bad_hpo}/{total_hpo} not included in AMELIE query')


    def query_amelie(self, force=False) -> 'build all the amelie result files':
        """
        run the amelie API
        """
        pw = self.pw
        prj = self.prj
        logger = self.logger

        fn = f'{pw}/{prj}.amelie.lite.txt'
        if os.path.exists(fn) and not force:
            logger.info('amelie query already done')
            self.check_amelie_missing()
            return 0
        logger.info('getting information from AMELIE')

        fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
        fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'

        total_genes = line_count(fn_genelist)
        total_hpo = line_count(fn_hop_pure_id)

        if total_genes > 1000:
            logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), would split the genelist, not deployed yet...')

        amelie_api.main(prj, fn_genelist, fn_hop_pure_id, self.pheno_file_select, f'{pw}/{intermediate_folder}', pw_main=pw, force=force)
        self.check_amelie_missing()

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

        if not os.path.exists(fn):
            logger.error(f'amelie result not found: {fn}')
            sys.exit(1)

        d_amelie = [_.strip().split('\t') for _ in open(fn)]
        d_amelie = {k: v for k, v in d_amelie}
        logger.info(f'AMELIE count = {len(d_amelie)}')

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

    def run_proband_sv_match(self) -> "{pw}/{prj}.merged.sorted.tsv":
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
            elif t.startswith('Brother'):
                sn = int(t.replace('Brother', ''))
                trio_order.append([sample_id, 10+sn])
            elif t.startswith('Sister'):
                sn = int(t.replace('Sister', ''))
                trio_order.append([sample_id, 100+sn])

        trio_order = sorted(trio_order, key=lambda x: x[1])
        trio_order = [_[0] for _ in trio_order]

        proband_id_lite = proband_id.replace('UDN', '')
        header_suffix_copy_number = [f'Proband {proband_id_lite}']
        header_suffix_gn = []
        header_suffix_sv_type = []
        header_suffix_overlap = []

        for i in trio_order:
            aff_state = family[i]['aff_state']
            aff_state = {0: 'unaff', 1:'affected'}[aff_state]
            type_ = family[i]['type']
            sample_id_lite = i.replace('UDN', '')

            header_suffix_copy_number.append(f'{type_} ({aff_state}) {sample_id_lite}')

            header_suffix_gn.append(f'gn_{type_}')
            header_suffix_sv_type.append(f'sv_type_{type_}')
            header_suffix_overlap.append(f'overlap_{type_}')


        header_suffix = header_suffix_copy_number + header_suffix_gn + header_suffix_sv_type + header_suffix_overlap

        header_suffix = '\t'.join(header_suffix)

        with open(fn_proband) as fp:
            header = next(fp).strip().replace('copy_number', '')
            header = f'mut_type\tAMELIE\t{header}{header_suffix}'
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
                _, chr_, s, e, gn = [a[_] for _ in [0, 1, 2, 3, 7]]
                chr_ = chr_.lower()

                s = int(s)
                e = int(e)
                bin_ = int(s / 100000)

                # info stores the extra data need to be added to this line
                info = {}
                try:
                    amelie = d_amelie[gn]
                except:
                    logger.debug(f'gene not found in AMELIE list, gene = {gn}')
                    continue

                # get the same bin in the same family annotation results
                for sample_id in trio_order:
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
                    else:
                        cn = 'oo'
                        if chr_.lower().find('chrx') > -1 and sex == 1:
                            cn = 'oy'
                        elif chr_.lower().find('chry') > -1 and sex == 1:
                            cn = 'xo'
                        elif chr_.lower().find('chry') > -1 and sex == 2:
                            cn = 'na'
                        else:
                            cn = 'oo'
                        info[sample_id] = [cn, 'not_found'] + [''] * 4

                res = f'\t{amelie}\t{i}\t'
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

                res_suffix = '\t'.join(res_suffix)

                res += res_suffix
                print(res, file=final)
        final.close()

        # sort the result by amelie score
        sorted_file = f'{pw}/{prj}.merged.sorted.tsv'
        os.system(f"head -1 {merged_table} > {sorted_file};sed '1d' {merged_table} |sort -t $'\\t' -k 1nr  >> {sorted_file}")
