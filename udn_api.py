#! /usr/bin/env python3
"""
this is to get the files from UDN gateway
https://documenter.getpostman.com/view/1367615/RVu1JBWH#0ab3e9c3-792c-b292-6615-761a505df5bb
https://nbviewer.jupyter.org/github/hms-dbmi/udn-gateway-api/blob/master/Sequence%20files%20download%20guide.ipynb

update:
the res dict, the key may be duplicate, such as sister, brother, before changing, if a proband have multiple brothers, then only the last one would be recorded into the dict, the previous one would be overwritten

update 2020-11-03
add chromedriver path when running under ACCRE
update 2021-03-24  add uploading to ACCRE option
"""
import sys
# import pprint
import os
import time
import re
import glob
import json
import pickle
import logging
import requests
# from selenium import webdriver
# import selenium.common.exceptions as sel_exp
# from selenium.webdriver.support.ui import WebDriverWait as wait
# from selenium.webdriver.common.keys import Keys
import udn_igv
from termcolor import colored
from datetime import datetime
fromtimestamp = datetime.fromtimestamp

import hashlib
from getpass import getpass
import base64
from cryptography.fernet import Fernet
import cryptography.fernet as fernet

from bs4 import BeautifulSoup as bs
# basic settings
base_url = 'https://gateway.undiagnosed.hms.harvard.edu/api/2.0'
platform = sys.platform.lower()

ft_convert = {'bai': 'bam', 'cnv.vcf': 'cnv', 'gvcf': 'gvcf', 'fq': 'fastq', 'vcf': 'vcf'}
ft_convert.update({_: _ for _ in ft_convert.values()})

class Credential():
    def __init__(self, password, fn_pickle):
        """
        if fn_pickle == None, this could be used as dump the new credential
        """
        self.key = self.create_encrypt_key(password)
        self.fnout = fn_pickle
        self.dump = self.encrypt_file
        self.load = self.decrypt_file

    def create_encrypt_key(self, password):
        """
        create a key for cryptography package
        """
        keymd5 = hashlib.md5(password.encode()).hexdigest().encode()
        key = base64.urlsafe_b64encode(keymd5)
        return Fernet(key)

    def encrypt_file(self, data, fnout=None):
        """
        the data should be a byte obj, otherwise, would convert to byte obj using pickle
        """
        fnout = fnout or self.fnout
        if not isinstance(data, bytes):
            data = pickle.dumps(data)
        data_en = self.key.encrypt(data)
        # dump the encrypted data
        with open(fnout, 'wb') as out:
            pickle.dump(data_en, out)

    def decrypt_file(self):
        # load and decode the file
        try:
            with open(self.fnout, 'rb') as fp:
                data_en = pickle.load(fp)
        except pickle.UnpicklingError:
            print('ERROR: file is not a pickle file')
            return 0
        try:
            res = pickle.loads(self.key.decrypt(data_en))
        except fernet.InvalidToken:
            print(f'Passcode for the credential is not correct !: {self.fnout}')
            sys.exit(1)
        return res

def get_driver(driver, browser='firefox', headless=True):
    try:
        driver.current_url
    except:
        from selenium import webdriver
        from selenium.webdriver.firefox.options import Options
        options = Options()
        options.headless = headless
        platform = sys.platform
        if platform == 'darwin':
            exe_firefox = '/Users/files/work/package/firefox/geckodriver'
        else:
            exe_firefox = f'/home/chenh19/tools/geckodriver'
        driver = webdriver.Firefox(options=options, executable_path=exe_firefox)
        logger.warning('firefox initiated')

    return driver


def get_json(action=None, payload=None, url=None, header=None):

    headers = header or header1
    if not action and not url:
        logger.error(f'Neither action nor URL was specified: {action}')
        return 0
    url = url or f'{base_url}/{action}'

    payload = payload or {}
    r = requests.get(url, headers=headers, data=payload)
    if r.status_code != 200:
        logger.error(f'action={action} - status_code = {r.status_code}')
        return 0
    try:
        return r.json()
    except:
        logger.error(f'action={action} - fail to get json')
        return 0

def dump_json(obj, fn):
    with open(fn, 'wb') as out:
        pickle.dump(obj, out)


def get_file_extension(fn):
    m = re.match(r'.*?\.(bam|bai|cnv|fastq|fq|gvcf|vcf)\b(\.gz)?', fn.lower())
    if m:
        try:
            return ft_convert[m.group(1)], m.group(2)
        except:
            return m.group(1), m.group(2)
    else:
        logger.warning(f'unclear file type: {fn}')
        return None, None

def format_comment(comment):
    tmp = [_.strip() for _ in comment.split('\n') if _.strip()]
    tmp = [f'#### {_}' if len(_) < 50 and _[-1] == ':' else _ for _ in tmp]
    return '\n\n'.join(tmp)


def ts(timestamp):
    try:
        return fromtimestamp(timestamp/1000).strftime('%Y-%m-%d')
    except:
        return 'NA'

def get_age(timestamp):
    try:
        dob = fromtimestamp(timestamp/1000)
        today = datetime.today()
        return f'{(today - dob).days / 365:.1f} yo'
    except:
        return 'NA'


def download_report(reports):
    for rep_id, rep_name in reports:
        if not os.path.exists(f'intermed/{rep_name}'):
            logger.info(f'\tdownloading seq report {rep_name}')
            try:
                rep_url = get_json(f'participants/{udn}/sequencing/reports/{rep_id}')['downloadLink']
                os.system(f'wget "{rep_url}" -O intermed/{rep_name} > /dev/null 2>&1')
            except:
                logger.warning(f'\tfail to download report: {rep_name}')
        else:
            logger.info(f'\treport  already downloaded: {rep_name}')

def parse_seq_files(info, rel_to_proband):
    reports = []
    files = []
    udn = info['udnId']

    for i in info['requests']:
        for i1 in i['reports']:
            try:
                rep_id = i1['id']
                rep_name = f'{rel_to_proband}-{i1["filename"]}'
            except:
                logger.warning(f'invalid report:  \n{i1}')
                raise
                continue
            reports.append((rep_id, rep_name))

        for ifl in i['files']:
            ires = {}
            try:
                ires['url'] = ifl['storageLocation']
            except:
                ires['url'] = 'NA'

            try:
                ires['file_id'] = ifl['id']
            except:
                ires['file_id'] = 'NA'
            try:
                ires['size'] = ifl['filesize']
            except:
                ires['size'] = 'NA'
            try:
                ires['uploaded'] = ts(ifl['uploaded'])
            except:
                ires['uploaded'] = 'NA'


            ires['download'] = 'NA'

            for k1, k2 in [('sequenceid', 'seq_id'), ('assembly', 'build'), ('md5sum', 'md5'), ('sequencing_type', 'seq_type'), ('filename', 'fn')]:
                try:
                    ires[k2] = ifl['metadata'][k1]
                except:
                    ires[k2] = 'NA'

            if ires['fn'] == 'NA':
                logger.error(f'invalid file info: \n{ifl}')
                sys.exit(1)

            files.append(ires)


    return reports, files


def get_download_link(udn, fl_info, get_aws_ft, gzip_only=None, force_update=False):
    # return, if 0, still valid,
    # if 1, error getting link
    # if 2, the file type is not included in the get amazon url list
    # else, the actual link

    try:
        fn = fl_info['fn']
    except:
        logger.error(f'invalid file info: udn = {udn}, file info = {fl_info}')
        sys.exit(1)

    if not isinstance(gzip_only, set):
        gzip_only = set()

    # which include the amazon presigned URL
    ext, gz = get_file_extension(fn)
    try:
        ext = ft_convert[ext]
    except:
        logger.error(f'invalid file extension: {fn}: ext = {ext}')
        return 1

    if not gz and ext in gzip_only:
        logger.warning(f'file skipped due to gzip file only: {fn}')
        return 1

    if ext not in get_aws_ft:
        logger.debug(f'skip update amazon link due to file type limit: {fn} ext={ext}, valid ft={get_aws_ft}')
        return 2

    if ext == 'cnv' and os.path.exists(f'intermed/download.cnv.{udn}.done'):
        logger.debug(f'cnv already done, skipped')

    # check if the url is still valid
    get_url = force_update
    if not get_url:
        prev_url = fl_info['download']
        if prev_url == 'NA':
            get_url = True
        else:
            err_line_count = validate_download_link(prev_url)
            if err_line_count > 0:
                get_url = True
            else:
                return 0

    if get_url:
        try:
            file_id = fl_info['file_id']
        except:
            logger.error(f'file id not found: info = \n{fl_info}')
        try:
            download_json = get_json(f'participants/{udn}/sequencing/files/{file_id}')
            return download_json['downloadLink']
        except:
            logger.warning(f'fail to get download link json: fn = {fn}, file_id = {file_id}')
            return 1


def red(msg):
    return colored(msg, 'red')

def green(msg):
    return colored(msg, 'green')



def get_all_info(udn, res_all=None, get_aws_ft=None, udn_raw=None, valid_family=None,  udn_proband=None, gzip_only=None, sv_caller='dragen', newname_prefix=None):
    """
    res_all,  when running specific relative, use this dict to accumulate among the family member
    gzip only, a set, if the ext is found in this set, then, the file must be gziped to get involved. e.g. if vcf in gzip_only, then,  the file named .vcf$ would not be included
    """

    udn_raw = udn_raw or udn

    res_all = res_all or {}

    res = {}

    # get the basic info
    basic_info = get_json(f'participants/{udn}')
    sequence_info = get_json(f'participants/{udn}/sequencing')
    is_proband = False

    get_report = True if 'report' in get_aws_ft else False
    res['simpleid'] = udn


    relation = basic_info['relation']
    if relation is None:
        rel_to_proband = 'Proband'
        res['rel_to_proband'] = rel_to_proband
        res['rel_to_proband_short'] = 'P'
        is_proband = True
    else:
        rel_to_proband = relation['name']
        res['rel_to_proband'] = rel_to_proband
        res['rel_to_proband_short'] = relation['shortName']

    res['reports'], res['files'] = parse_seq_files(sequence_info, rel_to_proband)

    # add the download link
    uploaded_date_ready = 0
    res['seq_uploaded'] = 'NA'

    for ifl in res['files']:
        if uploaded_date_ready == 0:
            if ifl['fn'].find('fastq') > -1:
                tmp = ifl['uploaded']
                if tmp != 'NA':
                    uploaded_date_ready = 1
                    res['seq_uploaded'] = tmp

        ifn = ifl['fn']
        download_link = get_download_link(udn, ifl, get_aws_ft, gzip_only=None)
        if download_link == 0:
            logger.info(colored(f'\tamazon link still valid: {ifl["fn"]}', 'green'))
        elif download_link == 1:
            logger.info(colored(f'\tfail to get amazon link {ifl["fn"]}', 'red'))
        elif download_link == 2:
            continue
        else:
            logger.info(green(f'\turl updated: {ifn}'))
            ifl['download'] = download_link

    if get_report:
        download_report(res['reports'])
    # else:
    #     logger.warning(red(f'\treport downloading skipped: rel_to_proband'))

    check_items = [
        ('firstname',['nameFirst'], basic_info),
        ('lastname',['nameLast'], basic_info),
        ('dob',['dateOfBirth'], basic_info, ts),
        ('age',['dateOfBirth'], basic_info, get_age),
        ('race',['races', 'name', 0], basic_info),
        ('gender',['birthAssignedSex', 'name'], basic_info),
        ('alive',['deceased'], basic_info, lambda _: int(not _)),
        ('affected',['affected', 'name'], basic_info),
        ('phenotip',['phenoTipsId'], basic_info),
        ]

    if is_proband:
        family_members = get_json(f'participants/{udn}/family')
        phenotip_info = get_json(f'participants/{udn}/phenotips')['phenotips']
        application_id = basic_info['application']['id']
        application = get_json(f'applications/{application_id}')
        baylor_candidate_raw = basic_info['diagnoses']


        res_family = []
        baylor_candidate = []
        res_pheno = []
        for i in family_members['familyMembers']:
            try:
                name = f'{i["nameFirst"]} {i["nameLast"]}'
            except:
                name = 'NA'

            try:
                rel_short = i['relation']['shortname']
            except:
                rel_short = 'NA'

            try:
                rel_full = i['relation']['name']
            except:
                rel_full = 'NA'


            try:
                fam_udn = i['udnId']
            except:
                fam_udn = 'NA'

            ires = [rel_full, name, fam_udn]


            res_all = get_all_info(fam_udn, res_all=res_all, get_aws_ft=get_aws_ft, udn_raw=udn_raw, valid_family=valid_family,  udn_proband=udn, gzip_only=gzip_only, sv_caller=sv_caller, newname_prefix=newname_prefix)
            res_family_member = res_all[rel_full]
            try:
                ires.append(res_family_member['gender'])
                ires.append(res_family_member['affected'])
                ires.append(res_family_member['seq_uploaded'])
            except:
                logger.error(f'missing keys for family member res:\n{res_family_member}')
                raise

            # | Relative | Name | UDN | Gender | Affected| Sequence_Uploaded
            res_family.append(ires)

        res['family'] = res_family

        for i in baylor_candidate_raw:
            ires = []
            k = 'OMIM_gene'
            try:
                v = i['geneName']
            except:
                v = 'NA'
            ires.append((k, v))

            k = 'OMIM_gene_id'
            try:
                v = i['geneMIM']
            except:
                v = 'NA'
            ires.append((k, v))

            k = 'certainty'
            try:
                v = i['certainty']['shortDescription']
            except:
                v = 'NA'
            ires.append((k, v))

            k = 'phenothype_explain'
            try:
                v = i['majorPhenotypeDegreeExplain']
            except:
                v = 'NA'
            ires.append((k, v))
            baylor_candidate.append(ires)

        res['baylor_candidate'] = baylor_candidate


        for i in phenotip_info['features']:
            if i['observed'] != 'yes':
                continue

            try:
                res_pheno.append((i['id'], i['label']))
            except:
                logger.warning(f'invalid pheno, participants/{udn}/phenotip\n{i}')

        res['symp'] = res_pheno
        try:
            res['family_history'] = phenotip_info['notes']['family_history']
        except:
            res['family_history'] = 'NA'

        try:
            res['medical_history'] = phenotip_info['notes']['medical_history']
        except:
            res['medical_history'] = 'NA'


        check_items.extend(
            [
                ('similar_symp',['familyWithSimilarSymptoms'], application),
                ('summary',['review', 'narrativeSummary', ], application),
                ('evaluation',['review', 'pertinentPriorEvaluations', ], application),
                ('comment',['review', 'comments'], application),
            ]
        )


    for i in check_items:
        k, json_keys, v = i[:3]
        try:
            func = i[3]
        except:
            func = None
        try:
            for ktmp in json_keys:
                v = v[ktmp]
            if func:
                v = func(v)
        except:
            v = 'NA'
        res[k] = v

    res_all[rel_to_proband] = res

    if is_proband:
        try:
            parse_api_res(res_all, udn_raw=udn_raw, valid_family=valid_family, sv_caller=sv_caller, newname_prefix=newname_prefix)
        except:
            logger.error(f'fail to parse the result dict, try again')
            raise

    return res_all


def get_logger(prefix=None, level='INFO'):
    prefix = prefix or 'udn_api_query'
    fn_log = f'{prefix}.log'
    fn_err = f'{prefix}.err'

    fmt = logging.Formatter(
        '%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    console.setLevel(level)

    fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
    fh_err = logging.FileHandler(fn_err, mode='a', encoding='utf8')

    fh_file.setLevel('DEBUG')
    fh_file.setFormatter(fmt)

    fh_err.setLevel('ERROR')
    fh_err.setFormatter(fmt)

    logger = logging.getLogger(prefix)
    logger.setLevel('DEBUG')
    logger.addHandler(console)
    logger.addHandler(fh_file)
    # logger.addHandler(fh_err)
    return logger


def create_cred():
    token = getpass(prompt='Your UDN token: ')
    email = getpass(prompt='Your email to login to UDN gateway: ')
    password = getpass(prompt='Your password to login to UDN gateway: ')
    vunetID = getpass(prompt='Your VUNetID to login to UDN gateway: ')

    if not token:
        logger.error('Please input the UDN token')
        return 0
    if not email:
        logger.error('Please input the UDN email')
        return 0
    if not password:
        logger.error('Please input the UDN login password')
        return 0

    if not vunetID:
        logger.error('Please input the UDN login vunetID')
        return 0

    cred = {'token': token, 'email': email, 'password': password, 'vunetID': vunetID}
    key.dump(cred)
    return cred

def validate_download_link(url):
    """
    Does the url contain a downloadable resource
    return = error lines, so if == 0, means the url is still valid
    """
    error_code = os.popen(f'curl "{url}" 2>&1|head -c 10000|egrep -i "(Error|AccessDenied)"|wc -l').read().strip()

    error_code = int(error_code)
    return error_code


    # try:
    #     h = requests.head(url, allow_redirects=True)
    # except:
    #     return False
    # header = h.headers
    # content_type = header.get('content-type')
    # if 'text' in content_type.lower():
    #     return False
    # if 'html' in content_type.lower():
    #     return False
    # if 'xml' in content_type.lower():
    #     return False

    # return True

def parse_api_res(res, renew_amazon_link=False, update_aws_ft=None, pkl_fn=None, udn_raw=None, valid_family=None, gzip_only=None, sv_caller='dragen', newname_prefix=None):
    """
    res = {proband: {}, relative: [{family1}, {family2}]}
    1. build a report for the proband
    2. extract the phenotype info HPO
    3. build the download script for fastq and vcf files
    4. build the udn config file based on the relative info
    5. IGV session

    file extension = .bam, .bai,
    .fastq.gz
    .cnv.vcf.gz

    un-used = .vcf.gz, .joint.vcf,  .gvcf.gz,
    """
    if not isinstance(gzip_only, set):
        gzip_only = set()
    update_aws_ft = set([ft_convert[_] for _ in update_aws_ft]) if update_aws_ft else set(['fastq', 'cnv'])


    logger.debug(f'file type to update url = {update_aws_ft}')
    # pw_raw = os.getcwd().rsplit('/', 1)[-1]
    # pw = f'/data/cqs/chenh19/udn/{pw_raw}'
    if not os.path.isdir('origin'):
        os.system('mkdir origin 2>/dev/null')

    logger.info(f'keys for the result: {res.keys()}')
    proband_id = res['Proband']['simpleid']

    # filter on the valid_family
    if valid_family is not None:
        for irel in list(res.keys()):
            if irel == 'Proband':
                continue
            iudn = res[irel]['simpleid']
            if iudn not in valid_family:
                logger.warning(f'family member not in select family list, excluded: {irel} - {iudn}')
                del res[irel]

    # update the URL
    if renew_amazon_link:
        for irel, v1 in res.items():
            iudn = v1['simpleid']

            for ifl in v1['files']:
                download_link = get_download_link(iudn, ifl, get_aws_ft=update_aws_ft, gzip_only=gzip_only)

                if download_link == 0:
                    logger.info(colored(f'\tamazon link still valid: {ifl["fn"]}', 'green'))
                elif download_link == 1:
                    logger.info(colored(f'\tfail to get amazon link {ifl["fn"]}', 'red'))
                elif download_link == 2:
                    continue
                else:
                    logger.info(green(f'\turl updated: {ifl["fn"]}'))
                    ifl['download'] = download_link

    # report for proband
    proband = res['Proband']
    udn = proband['simpleid']
    hpo = proband['symp']



    fn_pdf = f'basic_info.{udn}.pdf'
    fn_md = f'basic_info.{udn}.md'
    force_new_pdf = False
    if force_new_pdf or not os.path.exists(fn_pdf):
        with open(fn_md, 'w') as out:
            out.write(f'# {udn}\n## 1. Basic Info\n- UDN\t{udn}\n')
            for k, v in zip('firstname,lastname,gender,race,dob,age,alive,affected,phenotip,seq_uploaded'.split(','), 'FirstName,LastName,Gender,Race,DateBirth,Age,Alive,Affect_state,Phenotip_ID,Seq_uploaded'.split(',')):
                out.write(f'- {v}: - \t{proband.get(k) or "NA"}\n')
            out.write('\n\n')

            # family members

            out.write('## 2. Family Members\n\n')
            out.write('| Relative | Name | UDN | Gender | Affected|Seq_uploaded \n')
            out.write('| :----: | :----: | :----: | :----: | :----: | :----: |\n')
            for ifam in proband["family"]:
                out.write('|' + '|'.join(ifam) + '|\n')

            out.write('\n\n')

            # HPO terms
            out.write(f'## 3. HPO terms\n')

            for k, v in hpo:
                out.write(f'- {k}\t{v}\n')
            out.write('\n\n')


            # symp,uuid_case
            out.write('## 4. Summary and diagnosis\n')
            for k, k_formal in zip('similar_symp,family_history,comment,medical_history,summary,evaluation'.split(','), 'Similar_Sympt,Family_history,Medical_history,Comment,Summary,Evaluation'.split(',')):
                v = proband.get(k) or 'NA'
                out.write(f'\n\n### {k_formal}\n\n{v}\n')
            out.write('\n\n')

        logger.info('converting basic info to pdf')
        os.system(f"pandoc  -t pdf {fn_md} --pdf-engine xelatex -o {fn_pdf}")
        try:
            if os.path.getsize(fn_pdf) > 1000:
                os.unlink(fn_md)
        except:
            pass


    # build the HPO terms file
    with open(f'origin/{udn}.terms.txt', 'w') as out:
        for k, v in hpo:
            out.write(v+'\n')

    if not os.path.exists(f'pheno.keywords.txt'):
        with open(f'origin/{udn}.terms.txt') as f:
            content = f.read()

        with open('pheno.keywords.txt', 'w') as o:
            print('# if need to exact match the word, add @ prefix\n#if this line is a regex pattern, use ! as prefix\n\n', file=o)
            o.write(content)

        os.system(f'cp origin/{udn}.terms.txt pheno.keywords.txt')


    with open(f'origin/{udn}.terms_hpo.txt', 'w') as out:
        for k, v in hpo:
            out.write(f'{v}\t{k}\n')

    with open(f'origin/{udn}.terms_pure_hpo.txt', 'w') as out:
        for k, v in hpo:
            out.write(k+'\n')

    # build the udn config file
    sex_convert = {'Female' : 2, 'Male': 1}
    affect_convert_to_01 = {'affected': 1, 'unaffected': 0}

    cfg_info = {'father': '', 'mother':'', 'male':[], 'female':[], 'other': []}

    for rel_to_proband, irel in res.items():
        if rel_to_proband == 'Proband':
            continue
        try:
            rel_udn = irel['simpleid']
            rel_gender = irel['gender']
        except:
            tmp_set = set(['simpleid', 'affected', 'gender'])
            logger.info(f'{rel_to_proband}: key not found: {[_ for _ in tmp_set if _ not in irel]}, keys={irel.keys()}')
            continue

        if not irel.get('files') and sv_caller != 'pacbio':
            logger.warning(red(f'{rel_to_proband}: no files available, skip... '))
            continue

        try:
            rel_aff = affect_convert_to_01[irel['affected'].lower()]
        except:
            logger.warning(f"{rel_to_proband}: unkown affect state: {irel['affected']}")
            continue

        if rel_to_proband.lower() == 'father':
            cfg_info['father'] = [rel_udn, rel_aff]

        elif rel_to_proband.lower() == 'mother':
            cfg_info['mother'] = [rel_udn, rel_aff]

        elif rel_gender.lower() == 'female':
            cfg_info['female'].append([rel_udn, rel_aff, rel_to_proband])

        elif rel_gender.lower() == 'male':
            cfg_info['male'].append([rel_udn, rel_aff, rel_to_proband])

        else:
            cfg_info['other'].append([rel_udn, rel_aff, rel_to_proband])

    logger.debug(cfg_info)

    male = '\n  '.join([f'male{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['male'])])
    female = '\n  '.join([f'female{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['female'])])
    other = '\n  '.join([f'other{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['other'])])


    cfg = f"""prj: {udn}  # the project name
path:
# if the file is under the path above, could use the filename only
sv_caller: {sv_caller}
pheno_file: {udn}.terms.txt  # include the phenotype terms
udn_match_gene_file:  # if available, specify here, this file contains the Number of occurrences in the 3rd column, gene symbol as the 1st column
# if not available, leave blank


# please keep the indentation below
IDs:
  # with or wihout "UND" as prefix. if not applicable, use "NA" or leave blank
  proband_id: {udn}
  proband_sex: {sex_convert[proband["gender"]]}  # 1=male, 2=female

  father: {cfg_info["father"]}  # [udn, affect_state, relative_to_proband(optional)], affected state 0=unaffected, 1=affected  for unaffected and affected respectively
  mother: {cfg_info["mother"]}
  {male}
  {female}
  {other}

default:
  # if not sure, leave this as blank (using default settings)
  pheno_file_select: pheno.select.txt # optional, default= leave blankinclude the terms to select the matched terms from amelie results could be different from above file
  vcf_file_path:  # default is same as "path" above (leave blank here), if the vcf file is not under the path above, specify here
  annotsv: # path for annotSV software, default is using the annotSV in ACCRE dock
  thres_annot_sv_ranking: 2  # annotsv ranking equal or largher than this
  thres_gnomad_maf: 0.01   # less than this value

"""
    with open(f'{udn}_settings.yaml', 'w') as out:
        out.write(cfg)


    # download file script
    out_info = open(f'download.info.{udn}.txt', 'w')
    out_md5 = open(f'download.{udn}.md5', 'w')

    if 'cnv' in update_aws_ft:
        fn_cnv_sh = f'intermed/download.cnv.{udn}.sh'
        out_cnv = open(fn_cnv_sh, 'w')
        out_cnv.write('missing=0\n\n')

    tmp = udn_raw or udn
    remote_pw = re.sub(r'^\d+_', '', tmp, 1)
    m = re.match(r'(.*_)?(UDN\d+)(_.*)?', remote_pw).groups()
    if not m:
        logger.warning(f'wrong project name fromat: {remote_pw}')
    else:
        if m[0] is not None:
            p1 = re.sub('^_*', '', m[0])
            p1 = re.sub('_*$', '', p1)
        else:
            p1 = None

        if m[2] is not None:
            p2 = re.sub('^_*', '', m[2])
            p2 = re.sub('_*$', '', p2)
        else:
            p2 = None

        p_trailing = '_'.join([_ for _ in (p1, p2) if _])
        p_trailing = '_' + p_trailing if p_trailing else ''
        remote_pw = m[1] + p_trailing

    out_info.write('rel_to_proband\tfn\trename\tsample_list\turl\tudn\tseq_type\tsize\tbuild\tmd5\turl_s3\tupload_type\tdownload_type\tremote_pw\n')

    check_dup = {}
    for rel_to_proband, irel in res.items():
        if not irel.get('files'):
            continue
        udn = irel['simpleid']
        for ifl in irel['files']:
            url = ifl['download']
            fn = ifl['fn']
            url_s3 = ifl['url']
            size = ifl['size']
            build = ifl['build']
            seq_type = ifl.get('seq_type') or 'NA'
            md5 = ifl['md5']
            newname = ''
            newname_base = ''

            ext, gz = get_file_extension(fn)

            if ext not in update_aws_ft:
                # logger.info(f'skipped due to not target file type: {ext}')
                continue

            if newname_prefix and re.match(r'.+\.vcf', fn.lower()) and not re.match(r'(cnv|joint)', fn.lower()):
                # 971146-UDN131930-P_971147-UDN771313-M_reheadered.joint.repeats.merged.vcf.gz
                newname_base = f'{newname_prefix}_{rel_to_proband}'

                fn_prune = fn.rsplit('/', 1)[-1].replace('.vcf', '').replace('.gz', '').replace('.tbi', '')
                m = re.match(r'.*UDN\d+-[a-zA-Z]+[\W_]+(.*)', fn_prune)
                if m:
                    newname = f'{newname_base}_{m.group(1)}'
                else:
                    try:
                        n_prev = check_dup[newname_base]
                    except:
                        check_dup[newname_base] = 0
                    else:
                        check_dup[newname_base] += 1
                        newname = f'{newname_base}_{n_prev + 1}'


            if re.match(r'.+\.cnv\.vcf(\.gz)?$', fn) and 'cnv' in update_aws_ft and fn.find('joint') < 0:
                if not os.path.exists(f'origin{fn}'):
                    if url.lower() == 'na':
                        logger.error(f'invalid CNV file url')
                    else:
                        out_cnv.write(f'size=$(stat -c %s origin/{fn} 2>/dev/null);if [[ "$size" -lt 1000 ]];then wget "{url}" -c -O "origin/{fn}";\nelse echo already exist: "origin/{fn}";fi\nsize=$(stat -c %s origin/{fn} 2>/dev/null);\nif [[ "$size" -lt 1000 ]];then missing=1;fi\n\n')


            out_info.write(f'{rel_to_proband}\t{fn}\t{newname}\t{newname_base}\t{url}\t{udn}\t{seq_type}\t{size}\t{build}\t{md5}\t{url_s3}\n')
            out_md5.write(f'{md5}  {fn}\n')


    if 'cnv' in update_aws_ft:
        out_cnv.write(f'if [[ $missing -eq 0 ]];\n    then touch intermed/download.cnv.{udn}.done;\nelse \n    rm intermed/download.cnv.{udn}.done 2>/dev/null;\nfi')
        out_cnv.close()

    out_md5.close()
    out_info.close()

    if 'bam' in update_aws_ft:
        logger.info('building IGV script')
        udn_igv.main('.', udn, logger)


def resolve_udn_id(s):
    """
    extract the UDN number from the string
    """
    s = s.upper()
    return re.findall(r'(UDN\d+)', s)


def resolve_upload_list():
    """
    just a demo of a task
    """
    s = """
    UDN117150_Proband, UDN482847_Father, UDN976077_Mother (Case 123)
    UDN508581_Proband (Case 218)
    UDN079252_Proband, UDN258649_Father, UDN260402_Mother (Case 50)
    UDN991897_Proband UDN187507_Mother, UDN464491_Father (Case 206)
    UDN527443_Proband, UDN414554_Father, UDN822987_Mother (case solved no need to reupload)
    UDN320766_Proband, UDN230088_Father, UDN710646_Mother (Case 201)
    UDN956557_Proband, UDN367189_Mother, UDN648275_Father, UDN036509_Afected Brother (Case 85)
    UDN335667_Proband, UDN388209_Mother, UDN757864_Father (Case 124)
    UDN279217_Proband, UDN197070_Father, UDN281869_Mother (Case 228)
    UDN327145_Proband, UDN316269_Mother, UDN430357_Father (Case 216)
    UDN179293_Proband, UDN216310_Mother, UDN027382_Father (Case 222)"""
    s = s.split('\n')
    res = []
    for i in s:
        i = i.strip()
        udns = re.findall(r'UDN\d+', i)
        case = re.search(r'\(Case (\d+)\)', i)
        if not case:
            continue
        proband = f'{udns[0]}_case{case.group(1)}'
        res.append('@'.join([proband, *udns[1:]]))
    return ' '.join(res)


def upload_files(nodename, pw_accre_data, pw_accre_scratch, udn_raw, rename=False, udn_raw_old=None, initial=False, upload_only=False):

    if rename:
        os.system(f"""sftp {nodename} <<< $'rename {pw_accre_data}/{udn_raw_old} {pw_accre_data}/{udn_raw};rename {pw_accre_scratch}/{udn_raw_old} {pw_accre_scratch}/{udn_raw};'
        """)
    # upload files
    if initial:
        logger.info(f'initial uploading to ACCRE: {udn_raw}')
        fls = []
        for ifl in upload_file_list:
            fls.append(f'put {ifl} {pw_accre_scratch}/{udn_raw}/')
        fls = '\\n'.join(fls)
        sftp_cmd = ''
        if not upload_only:
            sftp_cmd += f"""mkdir {pw_accre_data}/{udn_raw}/origin/\\nput -r * {pw_accre_data}/{udn_raw}/"""

        sftp_cmd += f"mkdir {pw_accre_scratch}/{udn_raw}/\\n{fls}"

        cmd = f"sftp {nodename} <<< $'{sftp_cmd}' 2>/dev/null >&2;touch {initial_upload_flag}"

        os.system(cmd)
    else:
        fls = []  # file name can't contain space
        if not upload_only:
            fls.append(f'mkdir {pw_accre_data}/{udn_raw}/')
        fls.append(f'mkdir {pw_accre_scratch}/{udn_raw}/')

        for ifl in upload_file_list:
            if not upload_only:
                fls.append(f'put {ifl} {pw_accre_data}/{udn_raw}/')

            fls.append(f'put {ifl} {pw_accre_scratch}/{udn_raw}/')
        fls = '\\n'.join(fls)
        cmd = f"""sftp {nodename} <<< $'{fls}' >/dev/null 2>/dev/null"""
        logger.info(f'update {upload_file_list} to  {nodename}: {udn_raw},  {pw_accre_data}/{udn_raw}/, {pw_accre_scratch}/{udn_raw}/')
        # logger.info(f'command = {cmd}')
        os.system(cmd)


def parse_phillips_map_file(fn):
    if not os.path.exists(fn):
        logger.error(f'file not exist: {fn}')
        sys.exit(1)
    invalid_lines = []
    udn_list = []
    rename_list = []
    with open(fn) as f:
        for i in f:
            i = i.strip()
            if not i:
                continue
            try:
                udn, casename = re.split(r'\s+', i)
            except:
                invalid_lines.append(i)
                continue
            udn = udn.upper()
            if not re.match(r'UDN\d+$', udn):
                invalid_lines.append(i)
                continue

            if casename.lower().find('case') < 0:
                casename = f'case_{casename}'


            udn_list.append(udn)
            rename_list.append(casename)
    if len(invalid_lines) > 0:
        logger.error(f'invalid lines found for {fn}:\n' + '\n'.join(invalid_lines))
        sys.exit(1)

    return udn_list, rename_list


if __name__ == "__main__":

    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('udn', help="""one or multiple UDN ID, sep by space, could also contain the prefix, eg. 11_UDN123456.  if only selected family number are needed, append with @, eg. UDN123@UDN456_father@UDN_789@mother, in this case, other family member would not be resolved""", nargs='*')
    ps.add_argument('-pw', help='working path, default is udn/cases, if use current folder, specify as . or pwd', default=None)
    ps.add_argument(
        '-fn_cred', '-cred',
        help="""filename for the credential file, must include the UDN token, login email, login username(vunetID), login_password""")
    ps.add_argument('-create', '-encry', '-enc', help="""enter into create credential mode""", action='store_true')
    ps.add_argument('-sv_caller', '-svcaller', '-svsoftware', help="""the software used for sv calling, default=dragen, alt=pacbio""", default='dragen', choices=['pacbio', 'dragen'])
    ps.add_argument('-update_cred', '-uc',  help="""update the credential info""", action='store_true')
    ps.add_argument('-url', '-urlonly', '-ckfls', help='get the files URL only (for each family member) and quit', action='store_true'),
    ps.add_argument('-demo', help="""don't resolve the aws download URL""", action='store_true')
    ps.add_argument('-ft', help="""could be multiple, if specified, would only update the amazon URL for these file types, if gzip only, add the .gz suffix. e.g. vcf.gz  would not match .vcf$""", nargs='*', default=None)
    ps.add_argument('-renew', '-new', '-update', help="""flag, update the amazon shared link. default is not renew""", action='store_true')
    ps.add_argument('-lite', help="""flag, donot download the cnv files and the bayler report""", action='store_true')
    ps.add_argument('-noupload', '-noup', help="""don't upload the files to ACCRE""", action='store_true')
    ps.add_argument('-force_upload', '-fu', help="""force upload all files to ACCRE""", action='store_true')
    ps.add_argument('-force_rerun', '-force', help="""force rerun the API query, ignore the local pikcle file""", action='store_true')
    ps.add_argument('-upload_only', '-uo', help="""only upload the files to upload folder, not to the data analysis folder""", action='store_true')
    ps.add_argument('-show', '-chrome', '-noheadless', help="""disable headless mode""", action='store_true')
    ps.add_argument('-showcred', help="""show the credential content and exit""", action='store_true')
    ps.add_argument('-reltable', '-rel', help="""force save the relative table, even in the parse_api_result mode""", action='store_true')
    ps.add_argument('-v', '-level', help="""set the logger level, default is info""", default='INFO', choices=['INFO', 'WARN', 'DEBUG', 'ERROR'])
    ps.add_argument('-nodename', '-node', '-remote', help="remote node name, default is va", choices=['va', 'pc'])
    ps.add_argument('-rename', '-newname', help="rename the case name to this, only valid for vcf file, number should match with prj number", nargs='*')
    ps.add_argument('-phillips', '-philips', '-phil', help="run the tasks for Phillips, need to rename the case. argument = filename for UDN-case number map")


    args = ps.parse_args()
    case_prefix = ''
    pw_accre_scratch_suffix = ''
    fn_phillips = args.phillips
    if fn_phillips:
        args.ft = ['vcf']
        args.pw = '.'
        args.upload_only = True
        args.lite = True
        # args.renew = True
        args.udn, args.rename = parse_phillips_map_file(fn_phillips)
        pw_accre_scratch_suffix = '/upload_for_phillips'
        case_prefix = fn_phillips.rsplit('/', 1)[-1].rsplit('.', 1)[0].replace('_rerun', '') + '_'

    print(args)

    force = args.force_rerun
    udn_list = args.udn or [os.getcwd().rsplit('/')[-1]]
    newname_prefix_l = args.rename
    if newname_prefix_l:
        newname_prefix_l = [re.sub('case[\W_]+', '', _) for _ in newname_prefix_l]

    ft_input = args.ft

    force_upload = args.force_upload
    sv_caller = args.sv_caller



    save_rel_table = args.reltable
    headless = not args.show
    if platform == 'darwin':
        root = args.pw or '.'
        if root.lower() in ['.', 'pwd']:
            root = os.getcwd()

        if root.find('/udn/cases/done/') > -1:
            root = root.rsplit('/done/', 1)[0] + '/done'
        elif root.find('/udn/cases/') > -1:
            root = root.rsplit('/udn/cases/', 1)[0] + '/udn/cases'
        os.chdir(root)
    else:
        logger.info('platform= {platform}')
        root = os.getcwd()
    logger = get_logger(level=args.v)

    pw_script = os.path.realpath(__file__).rsplit('/', 1)[0]
    lite = args.lite
    upload = not args.noupload or force_upload

    if not args.fn_cred and os.path.exists('udn.credential.pkl.encrypt'):
        fn_cred = 'udn.credential.pkl.encrypt'
    else:
        fn_cred = f'{pw_script}/udn.credential.pkl.encrypt'

    demo = args.demo
    pw_script = os.path.realpath(__file__).rsplit('/', 1)[0]
    urlonly = args.url

    upload_only = args.upload_only


    nodename = args.nodename or 'va'
    pw_accre_data = '/data/cqs/chenh19/udn'
    pw_accre_scratch = '/fs0/members/chenh19/tmp/upload' if nodename == 'va' else '/scratch/h_vangard_1/chenh19/udn/upload'


    pw_accre_scratch += pw_accre_scratch_suffix
    # pw_accre_upload = pw_accre_scratch

    upload_file_list = ['pheno.keywords.txt', 'download.*.txt', 'download.*.md5']  # upload these files to scratch and data of ACCRE

    gzip_only = set()
    update_aws_ft = set()

    if ft_input is None:
        update_aws_ft = ['fastq']
    else:
        err = 0
        ft_input = [_ for _ in ft_input if _.strip()]
        for i in ft_input:
            ext = i.replace('.gz', '')
            if i[-2:] == 'gz':
                gzip_only.add(ext)
            try:
                update_aws_ft.add(ft_convert[ext])
            except:
                logger.error(f'invalid filetype: {i}')
                err = 0
        if err:
            sys.exit(1)


    passcode = getpass('Input the passcode for the encrypted credential file: ')
    if not fn_cred:
        fn_cred = input('Specify the credential file name, if not exist, would create it: ')

    if not os.path.exists(fn_cred):
        fn_cred = f'{pw_script}/udn.credential.pkl.encrypt'
    key = Credential(passcode, fn_cred)

    logger.info(fn_cred)

    if args.create:
        cred = create_cred()
    elif not os.path.exists(fn_cred):
        logger.warning(f'The credential file not found: {fn_cred}, would create it')
        cred = create_cred()
    else:
        cred = key.load()

    # unpack the cred

    # print(cred)
    # sys.exit(1)
    update_cred = args.update_cred

    if update_cred:
        cred_key = input(f'Choose which key do you need to modify: {cred.keys()}:   ')

        if cred_key in {'password', 'pw'}:
            cred['password'] = getpass(prompt='Your password to login to UDN gateway: ')
        else:
            cred[cred_key] = input(f'new value for {cred_key}:   ')
        update_cred = True

    try:
        token = cred['token']
    except KeyError:
        logger.error('Token not found in credential file')
        token = input(prompt='You UDN token: ')

        cred['token'] = token
        update_cred = True

    try:
        email = cred['email']
    except KeyError:
        logger.error('email not found in credential file')
        email = input(prompt='Your email to login to UDN gateway: ')
        cred['email'] = email
        update_cred = True

    try:
        password = cred['password']
    except KeyError:
        logger.error('password not found in credential file')
        password = getpass(prompt='Your password to login to UDN gateway: ')
        cred['password'] = password
        update_cred = True

    try:
        vunetID = cred['vunetID']
    except KeyError:
        logger.error('vunetID not found in credential file')
        vunetID = getpass(prompt='Your VUNetID to login to UDN gateway: ')
        cred['vunetID'] = vunetID
        update_cred = True

    if update_cred:
        logger.warning('now dumping the new cred')
        key.dump(cred)
        logger.info('now exiting')
        sys.exit(1)

    api_token_gateway = cred['token']
    api_token_fileservice = cred['fs_token']

    if args.showcred:
        logger.info(cred)
        sys.exit(1)


    udn_list = [_.strip() for _ in udn_list if _.strip()]

    # validate udn
    p = re.compile(r'(?:.*)?(UDN\d+)')
    validated = []
    for n, i in enumerate(udn_list):
        i = [_.strip() for _ in i.split('@')]
        try:
            case_sn = f'case{newname_prefix_l[n]}_'
        except:
            case_sn = ''

        case_prefix_tmp = case_prefix.replace(case_sn, '')

        udn_raw = case_prefix_tmp + case_sn + i[0]

        if len(i) == 1:
            valid_family = None
        else:
            valid_family = i[1:]
            valid_family = [resolve_udn_id(_) for _ in valid_family]
            valid_family = [_[0] for _ in valid_family if _]
            logger.info(f'valid family={valid_family}')

        m = re.match(p, udn_raw)
        if not m:
            logger.warning(f'Invalid UDN_ID: {udn_raw}')
            continue
        udn = m.group(1)
        validated.append([udn_raw, udn, valid_family])

    if len(validated) == 0:
        logger.error('No UDN case passed in, please check')
        sys.exit(1)


    # sys.exit(1)
    logger.info('\n\n\n'+'#' *30)

    # driver, cookie_token = get_cookie(None, headless=headless)
    # logger.info(f'cookie_token=\n{cookie_token}')
    # driver = None
    # cookie_token = None

    header1 = {
        'Content-Type': 'application/json',
        'Authorization': f'Token {token}',
    }


    # get the current UDN list in current folder
    udn_exist = {}  # key = pure_udn_id, v = folder name

    folders = os.popen(f'ls -d */').read().split('\n')
    for i in folders:
        i = re.sub(r'\/$', '', i).rsplit('/', 1)[-1]
        m = re.match(r'.*(UDN\d+)(?:\W|$)', i)
        if m:
            udn_exist[m.group(1).upper()] = i

    if newname_prefix_l:
        if len(newname_prefix_l) != len(validated):
            logger.error(f'number of newname label ({len(newname_prefix_l)}) should match number of cases({len(validated)}), quit...')
            sys.exit(1)
    else:
        newname_prefix_l = [None for _ in validated]

    for tmp, newname_prefix in zip(validated, newname_prefix_l):
        udn_raw, udn, valid_family = tmp
        rename_remote = False
        udn_raw_old = ''
        if udn in udn_exist:
            udn_raw_old = udn_exist[udn].rsplit('/', 1)[-1]
            # print(f'mv {root}/{udn_raw_old} {root}/{udn_raw}')
            if udn_raw != udn_raw_old:
                logger.warning(f'rename the folder, old = {udn_raw_old} , new = {udn_raw}')
                os.system(f'mv {root}/{udn_raw_old} {root}/{udn_raw}')
                rename_remote = True
        else:
            try:
                os.makedirs(f'{root}/{udn_raw}/origin', exist_ok=True)
            except:
                pass
            try:
                os.makedirs(f'{root}/{udn_raw}/intermed', exist_ok=True)
            except:
                pass


        logger.info(f'now running {udn_raw}')
        pw_case = f'{root}/{udn_raw}'
        os.chdir(pw_case)

        logger = get_logger(f'{root}/{udn_raw}/{udn_raw}', level=args.v)
        logger.info(os.getcwd() + f', rename to {newname_prefix}' if newname_prefix else '')

        fn_udn_api_pkl = f'{root}/{udn_raw}/intermed/{udn}.udn_api_query.pkl'

        if os.path.exists(fn_udn_api_pkl) and not force:
            logger.info(green('\tdirectly load API query result from pickle file'))
            res = pickle.load(open(fn_udn_api_pkl, 'rb'))
            if demo:
                renew_amazon_link = False
            else:
                renew_amazon_link = args.renew

            # regular
            parse_api_res(res, update_aws_ft=update_aws_ft, renew_amazon_link=renew_amazon_link, udn_raw=udn_raw, valid_family=valid_family, gzip_only=gzip_only, sv_caller=sv_caller, newname_prefix=newname_prefix)

            with open(fn_udn_api_pkl, 'wb') as out:
                pickle.dump(res, out)
        else:
            res = get_all_info(udn, get_aws_ft=update_aws_ft, udn_raw=udn_raw, valid_family=valid_family, udn_proband=udn, gzip_only=gzip_only, sv_caller=sv_caller, newname_prefix=newname_prefix)
            with open(fn_udn_api_pkl, 'wb') as out:
                pickle.dump(res, out)


        if not lite and 'cnv' in update_aws_ft:
            if os.path.exists(f'intermed/download.cnv.{udn}.done'):
                logger.info(green('\tCNV file already downloaded'))
            else:
                logger.info(f'\tdownloading CNV vcf file')
                os.system(f'bash intermed/download.cnv.{udn}.sh >/dev/null 2>&1')
        try:
            os.unlink('geckodriver.log')
        except:
            pass

        initial_upload_flag = 'origin/initial_uploaded'
        if upload and os.path.exists(fn_udn_api_pkl) and not force_upload and os.path.exists(initial_upload_flag):
            initail_upload_toggle = False
        else:
            initail_upload_toggle = True


        if upload:
            upload_files(nodename, pw_accre_data, pw_accre_scratch, udn_raw, rename=rename_remote, udn_raw_old=udn_raw_old, initial=initail_upload_toggle, upload_only=upload_only)

        print('\n########################\n\n')
