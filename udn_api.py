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
import hashlib
from getpass import getpass
import base64
from cryptography.fernet import Fernet
import cryptography.fernet as fernet
import logging
import requests
from selenium import webdriver
import selenium.common.exceptions as sel_exp
from selenium.webdriver.support.ui import WebDriverWait as wait
import udn_igv

from bs4 import BeautifulSoup as bs
# basic settings
base_url = 'https://gateway.undiagnosed.hms.harvard.edu/api'
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
            logger.error(f'Passcode for the credential is not correct !: {self.fnout}')
            sys.exit(1)
        return res

# def get_driver(driver, headless=True):
#     try:
#         driver.current_url
#     except:
#         pass
#     else:
#         return driver
#     if platform == 'darwin':
#         path_chromedriver = '/Users/files/work/package/chrome_v90/chromedriver'
#     else:
#         path_chromedriver = '/home/chenh19/tools/chromedriver2.35'
#     option = webdriver.ChromeOptions()
#     if headless:
#         option.add_argument('--headless')
#     option.add_argument('--no-sandbox')
#     option.add_argument(f"--window-size=1080,720")
#     if platform == 'linux':
#         option.add_argument('--disable-dev-shm-usage')
#         option.binary_location = '/home/chenh19/tools/chrome/chrome'
#     driver = webdriver.Chrome(executable_path=path_chromedriver, options=option)
#     return driver


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



def login(driver, headless=False, source=None):
    """
    get cookie using selenium
    """
    timeout = 30

    try:
        driver.current_url
    except:
        logger.warning(f'driver not available: {driver}, source={source}')
        driver = get_driver(driver, headless=headless)

    current_url = driver.current_url
    if current_url.find('https://gateway.undiagnosed.hms.harvard.edu') == 0:
        logger.debug(f'session still active')
        return driver
    url = 'https://gateway.undiagnosed.hms.harvard.edu/login/dashboard/'
    driver.get(url)
    if current_url.find('hms-dbmi.auth0.com') > -1:
        # logged out, maybe due to timeout
        # <div class="auth0-lock-social-button-text">hua-chang.chen@vumc.org</div>
        logger.warning(f'you got logged out, try login again')

    # sys.exit(1)
    # if the user name already exist
    run_sso = True
    type_email = False
    try:
        email = cred["email"]
        button = wait(driver, 5).until(lambda x: x.find_element_by_xpath(f'//div[text()="{email}"]'))
    except:
        type_email = True
    else:
        # email saved in previous session'
        for _ in range(10):
            try:
                button.click()
            except:
                time.sleep(1)
                logger.debug(f'click button before sso: - {_}')
            else:
                break
        else:
            logger.error(f'fail to enter into SSO page')
            # driver.save_screenshot('fail_to_enter_sso')
            return None

    if type_email:
        email = wait(driver, timeout).until(lambda x: x.find_element_by_xpath('//input[@type="email" and @name="email"]'))
        password = driver.find_element_by_xpath('//input[@type="password"]')
        button = driver.find_element_by_xpath('//button[@class="auth0-lock-submit"]')

        for _ in range(10):
            if button.is_displayed():
                break
            time.sleep(1)
        else:
            logger.error('fail to display to UDN gateway interface')

        # driver.save_screenshot('udn_login_pre_sso.png')
        logger.debug('driver is set')

        for _ in range(10):
            try:
                email.click()
                email.send_keys(cred['email'])
            except:
                time.sleep(1)
            else:
                break
        else:
            logger.error(f'fail to input email, error')
        # driver.save_screenshot('udn_login_step1.png')

        # driver.save_screenshot('email_input.png')
        time.sleep(2)
        try:
            password.click()
            password.send_keys(cred['password'])
            run_sso = False
        except webdriver.remote.errorhandler.ElementNotVisibleException:
            run_sso = True
        except sel_exp.ElementNotInteractableException:
            run_sso = True
        button.click()

    current_url = driver.current_url
    if current_url.find('https://gateway.undiagnosed.hms.harvard.edu') == 0:
        logger.info(f'session still active')
        return driver

    if run_sso:
        # driver.save_screenshot('udn_sso.png')
        username = wait(driver, timeout).until(lambda x: x.find_element_by_xpath('//input[@id="username"]'))
        password = driver.find_element_by_xpath('//input[@id="password"]')
        button = driver.find_element_by_xpath('//a[@class="ping-button normal allow"]')

        username.click()
        username.send_keys(cred['vunetID'])

        password.click()
        password.send_keys(cred['password'])

        # driver.save_screenshot('udn_login_sso.png')
        button.click()
    for _ in range(10):
        if driver.current_url.find('gateway.undiagnosed.hms.harvard.edu') > -1:
            break
        time.sleep(1)
    else:
        logger.error(f'fail to login to UDN gateway, url=\n{driver.current_url}')
        return 0
    # url = "https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/6145/files/"
    # driver.get(url)
    # driver.save_screenshot('udn_gateway.png')
    return driver

def get_cookie(driver, headless=True):
    driver = login(driver, headless=headless, source='get_cookie')
    cookies = driver.get_cookies()
    cookie_token = {_['name']: _['value'] for _ in cookies}

    # csrftoken=2DIZXSDSd8ffFuwxfOcO3VKCC3W0a4v0
    # sessionid=0lkbrirzd0mz92o8h1ygpc2y3c7daf21
    # mod_auth_openidc_session=e36bf6ee-f612-4b5f-8e55-6f6ea7f3156c
    return driver, cookie_token


def cleanup_memory(driver):
    # try:
    #     print('now killing the chrome process')
    # except:
    #     print('\tnow killing the chrome process')
    try:
        driver.close()
        driver.quit()
    except:
        pass
    os.system("""ps -ef|grep -i "chrome"|awk '{print $2}'|sort|uniq|parallel 'kill -9 {} 2>/dev/null' """)


def get_json(udn=None, action=None, payload=None, files=None, url=None, header=None):
    udn = udn or 'na'

    headers = header or header1
    if not action and not url:
        logger.error(f'Neither action nor URL was specified: {udn}')
        return 0
    url = url or f'{base_url}/{action}/{udn}/'   # the trailing / is as must

    payload = payload or {}
    files = files or {}
    r = requests.request("GET", url, headers=headers, data=payload, files=files)
    if r.status_code != 200:
        logger.error(f'{udn} : action={action} - status_code = {r.status_code}')
        return 0
    try:
        return r.json()
    except:
        logger.error(f'{udn} : action={action} - fail to get json')
        return 0


def get_response(udn=None, action=None, payload=None, files=None, url=None, header=None):
    udn = udn or 'na'
    headers = header or header1
    if not action and not url:
        logger.error(f'Neither action nor URL was specified: {udn}')
        return 0
    url = url or f'{base_url}/{action}/{udn}/'   # the trailing / is as must

    payload = payload or {}
    files = files or {}
    r = requests.request("GET", url, headers=headers, data=payload, files=files)
    if r.status_code != 200:
        logger.error(f'{udn} : action={action} - status_code = {r.status_code}')
        return 0
    return r

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


def get_all_info(udn, cookie_token, rel_to_proband=None, res_all=None, info_passed_in=None, sequence_type_desired=None, demo=False, get_aws_ft=None, udn_raw=None, valid_family=None, lite_mode=None, udn_proband=None, gzip_only=None, driver=None, sv_caller='dragen'):
    """
    rel_to_proband, if run as proband, this is None, otherwize, this func is called during geting the relatives of proband, and specified
    info_passed_in,  {'affected': rel_aff, 'seq_status': have_seq}   include the affect state, sequenced_state, only valid when called during geting the relatives of proband, and specified
    res_all,  when running specific relative, use this dict to accumulate among the family member
    sequence_type_desired : valid = all, wgs, chip, rna, wes ,  is a list, default is [wgs]
    demo: do not resolve the amazon s3 link
    gzip only, a set, if the ext is found in this set, then, the file must be gziped to get involved. e.g. if vcf in gzip_only, then,  the file named .vcf$ would not be included
    """

    if not isinstance(gzip_only, set):
        gzip_only = set()
    udn_raw = udn_raw or udn
    set_seq_type = set('wgs,wes,rna,all,chip'.split(','))
    sequence_type_desired = sequence_type_desired or ['wgs']
    tmp = set(sequence_type_desired) - set_seq_type
    if len(tmp) > 0:
        logger.error(f'wrong sequence type desired found: valid is all, wgs, chip, rna, wes:  unmatch={list(tmp)}')
        return 0

    aff_convert = {'Unaffected': 'Unaffected', 'Affected': 'Affected', '3': 'Unaffected', '2': 'Affected', 'Unknown': 'Unknown', 3: 'Unaffected', 2: 'Affected'}
    res_all = res_all or {}
    rel_to_proband = rel_to_proband or 'Proband'
    info_passed_in = info_passed_in or {}

    res = {}
    res['rel_to_proband'] = rel_to_proband
    res['affect_from_proband_list'] = info_passed_in.get('affected') or 'Unknown'
    res['seq_status'] = info_passed_in.get('seq_status') or 'Unkown'

    if not cookie_token:
        driver, cookie_token = get_cookie(driver)

    if not cookie_token:
        logger.error('fail to get cookie, can\'t get the amazon s3 presigned URL, exit')
        return 0


    header_cookie = get_header_cookie(cookie_token)

    # prepare to write each json result to the pickle
    fn_json_pk = f'{udn}.json.pkl'
    fn_json_pk_bk = f'{udn}.json.pkl'
    d_json = {}
    if os.path.exists(fn_json_pk):
        size1 = os.path.getsize(fn_json_pk)
        try:
            size_bk = os.path.getsize(fn_json_pk_bk)
        except:
            size_bk = 0

        if size1 > size_bk:
            os.system(f'mv {fn_json_pk} {fn_json_pk_bk}')
        elif size1 < size_bk:
            logger.warning(f'json.pickle file not substituted due to size of new file is smaller than previous backup: {udn}')

    # get followup files
    action = 'followup/attachment'
    json = get_json(udn, action)
    # some case don't have followup files
    if json:
        res['followup'] = []
        if isinstance(json, list):
            for ijson in json:
                try:
                    res['followup'].append(ijson['file'])
                except KeyError:
                    logger.warning(f'fail to get followup file: {udn} - json={ijson}')
        else:
            try:
                res['followup'].append(json['file'])
            except KeyError:
                logger.warning(f'fail to get followup file: {udn} - json={json}')
        if not lite_mode:
            for n, ifl in enumerate(res['followup']):
                logger.info(f'downloading followup file: {udn}-{n+1}')
                r = requests.request('GET', ifl, headers=header_cookie)
                with open(f'origin/{udn}.followup-{n}.pdf', 'wb') as out:
                    out.write(r.content)

    # basic patient info
    action = 'patientrecord'
    json = get_json(udn, action)

    d_json['basic_info'] = json
    dump_json(d_json, fn_json_pk)

    # pickle.dump(json, open(f'{udn}.patientrecord.json.pkl', 'wb'))

    if json:
        try:
            res['gender'] = json['patient']['gender']['name']
        except:
            res['gender'] = json['patient']['sex']['name']

        res['firstname'] = json['patient']['patientfirst']
        res['lastname'] = json['patient']['patientlast']

        try:
            res['race'] = json['patient']['race'][0]['name']  # 1=male
        except:
            res['race'] = json['patient']['race']
        res['uuid'] = json['patient']['uuid']  #
        res['uuid_case'] = json['uuid']  #
        res['dob'] = json['patient']['dob']  # birthday

        try:
            affect_state_raw = json['patient']['affected'] # 2 = affected, 3=unaffected
            affect_state_raw = aff_convert[affect_state_raw]
        except:
            affect_state_raw = 'Unknown'

        if res['affect_from_proband_list'] == affect_state_raw:
            res['affect_final'] = affect_state_raw
        elif rel_to_proband == 'Proband':
            res['affect_final'] = affect_state_raw
        else:
            res['affect_final'] = f'{affect_state_raw} / {res["affect_from_proband_list"]}'

        res['phenotip'] = json['patient']['phenotipsid']
        res['alive'] = not json['deceased']

        res['simpleid'] = json['patient']['simpleid']
        try:
            res['similar_symp'] = format_comment(json['patient'].get('similarsymptomsexplain'))
        except:
            res['similar_symp'] = None

        try:
            res['comment'] = format_comment(json['comment'])
        except:
            res['comment'] = None


    # patien info
    # url = "https://gateway.undiagnosed.hms.harvard.edu/api/patient/{uuid}/"
    # the result is just like above patientrecord

    # proband only code
    if rel_to_proband.lower() == 'proband':
        # download the phenotip PDF file
        url_pdf = f'https://gateway.undiagnosed.hms.harvard.edu/phenotips/export/data/{res["phenotip"]}?format=pdf&pdfcover=0&pdftoc=0&pdftemplate=PhenoTips.PatientSheetCode'

        r = requests.request('GET', url_pdf, headers=header_cookie)

        with open(f'{udn}.phenotip.pdf', 'wb') as out:
            out.write(r.content)

        # get relatives
        uuid_case = res['uuid_case']
        url_relative = f'https://gateway.undiagnosed.hms.harvard.edu/patient/relatives/{uuid_case}/'
        if not uuid_case:
            logger.error(f'case UUID not found: patient record={res}')
            return 0

        try:
            driver.get(url_relative)
        except:
            logger.error(f'fail to get relatives page using selenium')
            try:
                response_relative = requests.request('GET', url_relative, headers=header_cookie).text
            except:
                logger.error(f'fail to get relatives info: {udn}')
                return 0
        else:
            try:
                _tmp_ele = wait(driver, 20).until(lambda _:_.find_element_by_xpath(f'//a[text()="Sequence uploaded"]'))
                # logger.warning(f'type={type(response_relative)}, len={len(response_relative)}')
                driver.save_screenshot(f'{udn}.relatives_table.png')
            except:
                logger.warning(f'{udn} : This case seems donnot have any relatives info !, check screenshot')
                driver.save_screenshot(f'{udn}.no_relative_found.relatives_table.png')
            response_relative = driver.page_source

        r = bs(response_relative, features='lxml')
        try:
            tb = r.find('tbody').find_all('tr')
        except:
            tb = []
        if len(tb) == 0:
            logger.error('Relative list is empty')

        # for this part, even the sequence state is noted as "-"  not uploaded, when clicking the actual case, sometimes, the sequence file still exist
        for i in tb:
            try:
                seq_status = i.find(attrs={'class': 'sequploaded'})
                have_seq = seq_status.text
            except:
                logger.warning(f'fail to get the sequence uploaded info from relative: {i}')
                continue

            rel_udn = i.find(attrs={'class': 'patientsimpleid'}).text
            rel_aff = i.find(attrs={'class': 'affected'}).text.strip()
            rel_aff = aff_convert[rel_aff]
            rel_to_proband_tmp = i.find(attrs={'class': 'relative'}).text

            if valid_family is not None:
                if rel_udn not in valid_family:
                    logger.warning(f'family member skipped due to not in valid family list: {rel_to_proband_tmp} - @{rel_udn}@ {valid_family}')
                    logger.info('\n\n****************')

                    continue

            # avoid the overwritten of res_all keys, such as multipel sisters /brothers
            try:
                rel_in_prev_run = [_.split('#')[-1] for _ in res_all if re.match(rel_to_proband_tmp, _)]
            except:
                logger.error(res_all)
                sys.exit(1)

            # incase multiple sister, multiple brother ...
            rel_seq = 0
            for _ in rel_in_prev_run:
                try:
                    n = int(_)
                except:
                    n = 1
                rel_seq = n if n > rel_seq else rel_seq

            if rel_seq > 0:
                # if no previous key was found, rel_seq shoul be zero
                rel_seq += 1
                rel_to_proband_tmp = f'{rel_to_proband_tmp}#{rel_seq}'


            res_all, driver = get_all_info(rel_udn, cookie_token, rel_to_proband=rel_to_proband_tmp, res_all=res_all, info_passed_in={'affected': rel_aff, 'seq_status': have_seq}, demo=demo, get_aws_ft=get_aws_ft, valid_family=valid_family, lite_mode=lite_mode, udn_proband=udn_proband, driver=driver)
            # logger.info(res_all.keys())

        # patient dignosis
        action = 'application'
        json = get_json(udn, action)

        d_json['doctor_review'] = json
        dump_json(d_json, fn_json_pk)

        if json:
            review = json[0]['application_review_set']
            if len(review) > 0:
                res['summary'] = format_comment(review[0]['summary'])
                res['evaluation'] = format_comment(review[0]['evaluations'])
            else:
                res['summary'] = res['evaluation'] = None


        # pehnotip
        try:
            uuid = res['uuid']
        except:
            logger.warning(f'{udn}: fail to get phenotip, UUID not available')
        else:
            url = f'{base_url}/phenotips/{uuid}/'
            json = get_json(udn, url=url)
            if json:
                symp_raw = json['phenotips']['features']
                symp = [[_['id'], _['label']] for _ in symp_raw if _['observed'] == 'yes']
                res['symp'] = symp

    # sequencing
    # https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/5366/
    action = 'sequences'
    json = get_json(udn, action)

    d_json['files'] = json
    dump_json(d_json, fn_json_pk)

    print('\n\n****************')
    logger.info(f'\tnow getting information of {rel_to_proband} : {udn}')
    # if rel_to_proband.find('#') > -1:
    #     logger.warning(f'Family member with same relative to proband: {rel_to_proband}')

    json_all = json
    if len(json_all) > 1:
        logger.warning(f'multi: {len(json_all)}  records found for sequences API, would use {sequence_type_desired}')
    elif not json or len(json_all) == 0:
        logger.error(f'\t{rel_to_proband}: fail to get sequencing json')

    # json_all is a list
    # each element is a dict, which is a sequencing request like
    # {'id': 1984,
#   'patient': {'simpleid': 'UDN139529'},
#   'sequencingfiles': [...]}
    seq_type_convert = {2: 'wes', 3: 'wgs', 4: 'rna', 5: 'reanalysis', 1: 'targeted variant'}
    seq_type_convert.update({str(k): v for k, v in seq_type_convert.items()})
    files = []
    res['bayler_report'] = []

    for json in json_all:
        seq_id = json["id"]
        sequence_type_raw = json['sequencing_type']
        try:
            sequence_type = seq_type_convert[sequence_type_raw]
        except:
            logger.warning(f'{rel_to_proband} seq_id = {seq_id}: unkown sequence type: type code={sequence_type_raw}')
            continue

        if sequence_type in sequence_type_desired  or'all' in  sequence_type_desired or sequence_type == 'reanalysis':
            logger.info(f'desired seq type found: {sequence_type}')
        else:
            logger.warning(f'{rel_to_proband}: skipped due to unmatch sequence type, seq_id = {seq_id}, desired={sequence_type_desired} found={sequence_type}')
            continue

        res['bayler_report'] += [
            f'https://gateway.undiagnosed.hms.harvard.edu/patient/sequencing-report/{_}/'
            for _ in json['sequencereports']]

        files_raw = json['sequencingfiles']
        logger.info(f'\t{rel_to_proband}: {sequence_type} seq_id={seq_id}  total files={len(files_raw)} , ft_to_download={get_aws_ft}')

        for n_fl, ifl in enumerate(files_raw):
            # each file is like
            # demo_format = {'udn_id': ['UDN525121'],
            # 'sequence_id': [6345],
            # 'uuid': '81c137f0-0a1d-403b-befb-83c40a24b81e',
            # 'fileserviceuuid': 'a09300e9-879a-4309-af82-b52e4aea9228',
            # 'sequencing_type': ['genome'],
            # 'filename': '939336-UDN525121-P.bam.bai',
            # 'sequencing_site': 'baylorseq',
            # 'complete': True,
            # 'file_data': {'creationdate': '2021-03-26T00:09:06-04:00',
            # 'modifydate': '2021-03-26T00:09:06-04:00',
            # 'description': '',
            # 'tags': [],
            # 'locations': [{'url': 'S3://udnarchive/8471503a-18a0-4a20-ae95-88e3681628ff/939336-UDN525121-P.bam.bai',
            #     'storagetype': 's3',
            #     'uploadComplete': '2021-03-26T00:09:08-04:00',
            #     'id': 933980,
            #     'filesize': 9481464}],
            # 'uuid': 'a09300e9-879a-4309-af82-b52e4aea9228',
            # 'filename': '939336-UDN525121-P.bam.bai',
            # 'expirationdate': '2021-10-12',
            # 'owner': {'email': 'imachol@bcm.edu'},
            # 'permissions': ['udn'],
            # 'id': 80272,
            # 'metadata': {'sequenceid': '6345',
            # 'assembly': 'hg38',
            # 'description': 'Binary Index File.',
            # 'md5sum': 'd21c930dbcbb160d6ecdde4069aa11c5',
            # 'patientid': '9fd8379b-a074-4457-83e5-b6f84e4fadf8',
            # 'sequencing_type': 'genome',
            # 'filename': '939336-UDN525121-P.bam.bai',
            # 'sequencing_location': 'baylor',
            # 'permissions': 'udn',
            # 'sequencing_sampleid': '939336-UDN525121'}}}

    # udn_id - list
    # sequence_id - list
    # sequencing_type - list
    # fileserviceuuid * useful
    # filename
    # complete

    # file_data
    #         uuid # same as fileserviceuuid in above layer
    #         filename
    #         metadata
    #                 sequenceid
    #                 assembly
    #                 md5sum
    #                 sequencing_type
    #                 filename
    #         locations  # is a list, ele=dict
    #                 url = S3 url
    #                 filesize

            completed = ifl['complete']

            try:
                fl_loc = ifl['file_data']['locations'][0]
                fl_meta = ifl['file_data']['metadata']
            except:
                logger.error(f"wrong file json: {ifl}")
                sys.exit(1)

            fn = ifl['file_data']['filename']

            # use fileserviceuuid key, rather than the uuid key
            file_uuid = ifl['fileserviceuuid']

            # the above file_uuid is used for api query, but the file_uuid used for the js query is ifl['uuid']
            file_uuid_js = ifl['uuid']
            file_uuid1 = ifl['file_data']['uuid']
            if file_uuid != file_uuid1:
                demo_d = {'fileserviceuuid': file_uuid, 'file_data': {'uuid': file_uuid1}}
                logger.warning(f'fs_uuid not match within json: {demo_d}')

            try:
                fl_url = fl_loc['url']
                fl_size = fl_loc['filesize']
            except:
                logger.warning(f'unexpected fl_info {fn}: {fl_loc}')
                fl_url = fl_size = 'NA'

    #         metadata
    #                 sequenceid
    #                 assembly
    #                 md5sum
    #                 sequencing_type
    #                 filename
            try:
                seq_id = fl_meta['sequenceid']
            except:
                pass
            try:
                fl_assembly = fl_meta.get('assembly')
                fl_md5 = fl_meta.get('md5sum') or fl_meta.get('md5')
            except:
                logger.warning(f'unexpected fl meta info  {fn}: {fl_meta}')
                fl_assembly = fl_md5 = 'NA'

            if re.match(r'.+\.bai$', fn):
                pass
            elif re.match(r'.+\.bam$', fn):
                pass
            elif re.match(r'.+\.g?vcf', fn):
                pass
            elif re.match(r'.+\.fastq.gz$', fn):
                pass
            else:
                logger.warning(f'unkown filetype: {sequence_type}: {fn}')
                continue

            # to get the actual amazon presined URL, there is a trick.
            # on each file URL, the click action is binded to
            # <a href="#" onclick="sequencing_file.show_file_download_modal(event, this)" file_uuid="c9ecb8c4-2760-4cdb-a3fe-d8fa61521da2" filename="921192-UDN830680-P.bam">921192-UDN830680-P.bam</a>
            # the action is sequencing_file.show_file_download_modal
            # after checking https://gateway.undiagnosed.hms.harvard.edu/static/js/patient/sequencing_file.js

            # there is a download link function
            # download_url = "/patient/downloadsequenceurl/" + file_uuid + "/";

            # request this URL would responde a json object
            # {"messages": [], "download_url": "https://udnarchive.s3.amazonaws.com:443/db488479-8923-4a60-9779-48f58b2f1dad/921192-UDN830680-P.bam?Signature=iZaOBprDz69DFzCZ%2BF0Pn9qEjDk%3D&Expires=1600299340&AWSAccessKeyId=AKIAZNIVXXYEAPHA7BGD", "success": true}

            # which include the amazon presigned URL
            ext, gz = get_file_extension(fn)
            try:
                ext = ft_convert[ext]
            except:
                logger.error(f'invalid file extension: {fn}')
                sys.exit(1)


            if not gz and ext in gzip_only:
                logger.debug(f'file skipped due to gzip file only: {fn}')
                download = 'NA'

            if get_aws_ft and len(get_aws_ft) > 0 and ext not in get_aws_ft:
                logger.debug(f'skip update amazon link due to file type limit: {fn} ext={ext}, valid ft={get_aws_ft}')
                download = 'NA'
            else:
                time.sleep(2)
                cookie_token, download, res_source = get_amazon_download_link(fn, file_uuid, cookie_token, seq_id=seq_id, demo=demo, driver=driver, file_uuid_js=file_uuid_js)
                if not download:
                    logger.debug(f'error when getting_amazon_download_link cookie_token = {cookie_token} \nfn="{fn}"\nfile_uuid="{file_uuid}"')
                    download = 'NA'
                else:
                    logger.info(f'{sequence_type} seq_id={seq_id}: {fn}, amazon link resolved ({res_source})')

            files.append({'download': download, 'relative': rel_to_proband, 'file_uuid': file_uuid, 'file_uuid_js': file_uuid_js, 'fn': fn, 'url': fl_url, 'complete': completed, 'size': fl_size, 'build': fl_assembly, 'assembly': fl_assembly, 'md5': fl_md5, 'seq_type': sequence_type, 'seq_id': seq_id})
    res['files'] = files
    res_all[rel_to_proband] = res

    if rel_to_proband == 'Proband':
        try:
            parse_api_res(res_all, cookie_token=cookie_token, udn_raw=udn_raw, valid_family=valid_family, driver=driver, sv_caller=sv_caller)
        except:
            logger.error(f'fail to parse the result dict, try again')
            raise
    return res_all, driver

def get_download_link_by_selinium(driver, fn, seq_id, logger, n=0):
    driver = login(driver, source='get_download_link_by_selinium')
    time.sleep(5)
    url = f'https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/{seq_id}/files/'
    # logger.debug(f'getting url {url}')
    if re.sub('/#', '', driver.current_url) != url:
        driver.get(url)
    try:
        ele_file = wait(driver, 5).until(lambda _:_.find_element_by_xpath(f'//a[@filename="{fn}"]'))
        ele_file.click()
        logger.debug(f'clicked: {fn}')
        download_link = wait(driver, 5).until(lambda _: _.find_element_by_xpath('//span[@class="sequencing_file_url"]'))
    except Exception as e:
        logger.debug(f'try - {n} - error={e}')
    else:
        if download_link.text.strip().find('https') > -1:
            logger.debug(f'download link got from selenium: {download_link.get_attribute("outerHTML")}')
            return download_link.text
    if n > 2:
        logger.debug(f'    fail to get download link for {fn}')
        return None
    n += 1
    get_download_link_by_selinium(driver, fn, seq_id, logger, n=n)

def get_amazon_download_link_api(file_uuid):
    # # get all files
    # udn_id = 'UDN665600'
    # gw_headers = {
    #     'Content-Type': 'application/json',
    #     'Authorization': f'Token {api_token_gateway}',
    #     'FSAuthorization': f'FSToken {api_token_fileservice}'
    # }
    # gw_host = 'gateway.undiagnosed.hms.harvard.edu'
    # url = f'https://{gw_host}/api/sequences/{udn_id}/'
    # r = requests.get(url, headers=gw_headers)

    # FileSerivce API needs a separate set of headers
    # file_uuid = 'b3f47f61-a45f-4e6e-a32d-39a79d39a65f'
    fs_headers = {
        'Content-Type': 'application/json; charset=UTF-8',
        'Authorization': f'Token {api_token_fileservice}',
        # 'Authorization': f'Token {api_token_gateway}'
    }
    # setup the host to easily switch between systems
    # fs_host = 'fileservicedev.aws.dbmi.hms.harvard.edu'
    fs_host = 'fileservice.dbmi.hms.harvard.edu'
    # fs_host = 'fileservice-ci.dbmi.hms.harvard.edu'
    url = f'https://{fs_host}/filemaster/api/file/{file_uuid}/download/'
    #
    r = requests.get(url, headers=fs_headers)
    try:
        download_url = r.json()['url']
    except:
        logger.debug(f'fail to get the downloadlink json: {r}')
        return None
    return download_url


def get_amazon_download_link(fn, file_uuid, cookie_token,  seq_id, driver, n=0, demo=False, file_uuid_js=None):
    # file_uuid_js, used for js query
    download_link = get_amazon_download_link_api(file_uuid)
    if download_link:
        return cookie_token, download_link, 'api'

    if demo:
        logger.debug(f'demo mode for {fn}')
        return cookie_token, None, None

    if not cookie_token:
        driver, cookie_token = get_cookie(driver)

    header_cookie = get_header_cookie(cookie_token)
    try:
        csrftoken = cookie_token['csrftoken']
        # csrftoken = 'CWdClEbLTBFuIw7NIDpt2KB8WeG82xaV'
    except:
        logger.error(f'csrftoken not found in cookie')
        logger.error(cookie_token)
        logger.error(header_cookie)
        sys.exit(1)

    if not file_uuid_js:
        logger.warning(f'file_uuid_js is not defined, use file_uuid instead')
        file_uuid_js = file_uuid

    url_download_info = f'https://gateway.undiagnosed.hms.harvard.edu/patient/downloadsequenceurl/{file_uuid_js}/?csrfmiddlewaretoken={csrftoken}'
    response_download_link = requests.request('GET', url_download_info, headers=header_cookie)
    if n > 2:
        return cookie_token, None, None

    # https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/6145/files/
    try:
        json = response_download_link.json()
    except Exception as e:
        if str(e.__class__).find('JSONDecodeError') > -1:
            # logger.info(f'not a json response,  url_download_info= "{url_download_info}", headers={header_cookie}')
            n += 1
            time.sleep(2)
            # logger.info(header_cookie)
            download_link = get_download_link_by_selinium(driver, fn, seq_id, logger=logger)
            res_source = 'selenium'
            if not download_link:
                cookie_token, download_link = get_amazon_download_link(fn, file_uuid, cookie_token, seq_id=seq_id, n=n, demo=demo, driver=driver, file_uuid_js=file_uuid_js)
                res_source = 'JS'
            if download_link:
                return cookie_token, download_link, res_source
        else:
            raise
    else:
        try:
            return cookie_token, json['download_url'], 'JS'
        except:
            # logger.debug(f'{fn}: \niter - {n} json does not contain the download link:\n {json}\nurl=  {url_download_info}\ncookie_token={cookie_token}\n\n')
            logger.debug(f'{fn}: \niter - {n} json does not contain the download link')
            n += 1
            time.sleep(2)
            download_link = get_download_link_by_selinium(driver, fn, seq_id, logger=logger)
            res_source = 'selenium'
            if not download_link:
                cookie_token, download_link = get_amazon_download_link(fn, file_uuid, cookie_token, seq_id=seq_id, n=n, demo=demo, driver=driver, file_uuid_js=file_uuid_js)
                res_source = 'JS'
            if download_link:
                return cookie_token, download_link, res_source

            return cookie_token, None, None


def get_header_cookie(cookie_token):
    # csrftoken=2DIZXSDSd8ffFuwxfOcO3VKCC3W0a4v0
    # sessionid=0lkbrirzd0mz92o8h1ygpc2y3c7daf21
    # mod_auth_openidc_session=e36bf6ee-f612-4b5f-8e55-6f6ea7f3156c

    cookies_list = [f"{k}={v}" for k, v in cookie_token.items()]
    cookie = '; '.join(cookies_list)
    return {'Connection': 'keep-alive', 'Pragma': 'no-cache', 'Cache-Control': 'no-cache', 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.162 Safari/537.36 Edg/80.0.361.109', 'Cookie': cookie}


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
    """
    error_code = os.popen(f'curl "{url}" 2>&1|head -c 10000|egrep -i "Error|Could not resolve host"|wc -l').read().strip()

    error_code = int(error_code)
    return not error_code


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

def parse_api_res(res, cookie_token=None, renew_amazon_link=False, update_aws_ft=None, pkl_fn=None, demo=False, udn_raw=None, valid_family=None, gzip_only=None, driver=None, sv_caller='dragen'):
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

    # get the screenshot of relative table
    fn_reltive_png = f'{proband_id}.relatives_table.png'
    if save_rel_table or not os.path.exists(fn_reltive_png):
        logger.info(f'now save the screenshot for relatives table')
        uuid_case = res['Proband']['uuid_case']
        url_relative = f'https://gateway.undiagnosed.hms.harvard.edu/patient/relatives/{uuid_case}/'
        if not uuid_case:
            logger.error(f'case UUID not found: patient record={res["Proband"]}')
            return 0

        driver = login(driver, headless=headless, source='parse_api_res, get relative table')

        try:
            driver.get(url_relative)
        except:
            logger.error(f'fail to get relatives page using selenium')
        else:
            driver.save_screenshot('browser.init.png')
            try:
                _tmp_ele = wait(driver, 10).until(lambda _:_.find_element_by_xpath(f'//a[text()="Sequence uploaded"]'))

                driver.save_screenshot(fn_reltive_png)
                logger.info(f'Screenshot: relative table saved!')
            except:
                driver.save_screenshot(f'{proband_id}.no_relatives_found.table.png')
                try:
                    os.symlink(f'{proband_id}.no_relatives_found.table.png', fn_reltive_png)
                except:
                    pass
                logger.warning(f'no relatives found for {proband_id}')

    # if not cookie_token:
    #     driver, cookie_token = get_cookie(driver)
    # header_cookie = get_header_cookie(cookie_token)
    # if not cookie_token:
    #     logger.error('fail to get cookie, can\'t get the amazon s3 presigned URL, exit')
    #     return 0

    # export the sequencing html table

    # filter on the valid_family
    if valid_family is not None:
        for irel in list(res.keys()):
            if irel == 'Proband':
                continue
            iudn = res[irel]['simpleid']
            if iudn not in valid_family:
                logger.warning(f'family member not in select family list, excluded: {irel} - {iudn}')
                del res[irel]


    affect_convert_to_01 = {'Affected': 1, 'Unaffected': 0, '3 / Unaffected': 0, '2 / Affected': 1, 'Unaffected / 3': 0, 'Affected / 2': 1, 'Unknown': -1}


    if renew_amazon_link:
        for rel_to_proband, v in res.items():
            logger.info(f'renew amazon download link for {rel_to_proband}')
            if 'files' not in v:
                logger.info(f'{rel_to_proband}: no files found')
                continue
            try:
                seq_id = v['seq_id'][0][0]
                logger.info(f'seq_id found in the outer of res dict ')
            except:
                # logger.warning(f'seq ID not found in dict: {rel_to_proband}')
                pass
                # continue
            seq_id_source = 'outer'
            for ifl in v['files']:
                fn = ifl['fn']
                file_uuid = ifl['file_uuid']
                file_uuid_js = ifl.get('file_uuid_js')
                old_download_link = ifl['download']
                try:
                    seq_id = ifl['seq_id']
                    seq_id_source = 'from_file'
                except:
                    pass

                try:
                    seq_id = int(seq_id)
                except:
                    logger.error(f'wrong seq ID: {seq_id} ,  source = {seq_id_source}')
                    sys.exit(1)
                logger.debug(f'seq_id, source = {seq_id_source}')

                if update_aws_ft and len(update_aws_ft) > 0:
                    ext, gz = get_file_extension(fn)
                    logger.debug(f'      {fn}: ext={ext}, gz={gz}')
                    if not gz and ext in gzip_only:
                        logger.debug(f'        file skipped due to gzip file only: {fn}')
                        continue
                    if ext not in update_aws_ft:
                        logger.debug(f'      skipp update amazon link due to file type limit: {fn} - ext={ext}, selected file type={update_aws_ft}')
                        continue

                if validate_download_link(old_download_link):
                    logger.info(f'      {fn}: href still valid')
                    # print(old_download_link)
                    continue

                time.sleep(2)
                cookie_token, download, res_source = get_amazon_download_link(fn, file_uuid, cookie_token, seq_id=seq_id, demo=demo, driver=driver, file_uuid_js=file_uuid_js)
                if not download:
                    download = 'NA'
                    logger.warning(f'      {fn}:fail to get amazon link')
                    continue
                logger.info(f'      {fn}: update amazon_link done ({res_source})')
                ifl['download'] = download

        pkl_fn = pkl_fn or f'{proband_id}.udn_api_query.pkl'
        if os.path.exists(pkl_fn) and os.path.getsize(pkl_fn) > 5000:
            os.system(f'mv {pkl_fn} {pkl_fn}.bk')
        pickle.dump(res, open(pkl_fn, 'wb'))

    # report for proband


    proband = res['Proband']
    udn = proband['simpleid']
    hpo = proband['symp']
    # logger.info(f'Proband files = {proband["files"]}')


    with open(f'{udn}.basic_info.md', 'w') as out:
        out.write(f'# {udn}\n## 1. Basic Info\n- UDN\t{udn}\n')
        for k, v in zip('firstname,lastname,gender,race,dob,alive,affect_final,uuid_case,phenotip'.split(','), 'FirstName,LastName,Gender,Race,DateBirth,Alive,Affect_state,UUID_case,Phenotip_ID'.split(',')):
            out.write(f'- {v}\t{proband[k]}\n')
        out.write('\n\n')

        # family members

        out.write('## 2. Family Members\n\n')
        out.write('| Relative | Name | UDN | Gender | Affected| Sequence_Uploaded | \n')
        out.write('| :----: | :----: | :----: | :----: | :----: | :----: |\n')
        for rel_to_proband, i in res.items():
            if rel_to_proband == 'Proband':
                continue

            # if gender not added
            out.write(f'| {rel_to_proband} | {i["firstname"]} {i["lastname"]} | {i["simpleid"]} | {i["gender"]} | {i["affect_final"]}| {i["seq_status"]} |\n')
        out.write('\n\n')

        # HPO terms
        out.write(f'## 3. HPO terms\n')

        for k, v in hpo:
            out.write(f'- {k}\t{v}\n')
        out.write('\n\n')
        # for k, v in hpo:
        #     out.write(f'- {v}\n')
        # out.write('\n\n')

        # symp,uuid_case
        out.write('## 4. Summary and diagnosis\n')
        for k, v in zip('similar_symp,comment,summary,evaluation'.split(','), 'Similar_Sympt,Comment,Summary,Evaluation'.split(',')):
            if proband[k]:
                out.write(f'\n\n### {v}\n\n{proband[k]}\n')
        out.write('\n\n')


    fn_pdf = f'{udn}.basic_info.pdf'
    # os.system(f"pandoc  -t pdf {udn}.basic_info.md --pdf-engine pdflatex -o {fn_pdf}")
    if not os.path.exists(fn_pdf):
        logger.info('converting basic info to pdf')
        os.system(f"pandoc  -t pdf {udn}.basic_info.md --pdf-engine xelatex -o {fn_pdf}")

    # build the HPO terms file

    with open(f'origin/{udn}.terms.txt', 'w') as out:
        for k, v in hpo:
            out.write(v+'\n')

    if not os.path.exists(f'pheno.keywords.txt'):
        os.system(f'cp origin/{udn}.terms.txt pheno.keywords.txt')


    with open(f'origin/{udn}.terms_hpo.txt', 'w') as out:
        for k, v in hpo:
            out.write(f'{v}\t{k}\n')

    with open(f'origin/{udn}.terms_pure_hpo.txt', 'w') as out:
        for k, v in hpo:
            out.write(k+'\n')

    # build the udn config file
    sex_convert = {'Female' : 2, 'Male': 1}

    cfg_info = {'father': '', 'mother':'', 'male':[], 'female':[], 'other': []}

    for rel_to_proband, irel in res.items():

        if rel_to_proband == 'Proband':
            continue
        try:
            rel_udn = irel['simpleid']
            rel_gender = irel['gender']
        except:
            tmp_set = set(['simpleid', 'affect_final', 'gender'])
            logger.info(f'{rel_to_proband}: key not found: {[_ for _ in tmp_set if _ not in irel]}, keys={irel.keys()}')
            continue

        if not irel.get('files') and sv_caller != 'pacbio':
            logger.warning(f'{rel_to_proband}: no files available, skip... ')
            continue

        try:
            rel_aff = affect_convert_to_01[irel['affect_final']]
        except:
            logger.warning(f"{rel_to_proband}: unkown affect state: {irel['affect_final']}")
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
    # brothers = '\n'.join([f'brother{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['brother'])])
    # sisters = '\n'.join([f'sister{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['sister'])])
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

    out = open(f'download.fastq.{udn}.txt', 'w')

    if 'bam' in update_aws_ft:
        out_bam = open(f'download.bam.{udn}.sh', 'w')
        html_bam = open(f'download.bam.{udn}.html', 'w')
        html_bam.write(f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>{udn} - bam </title>
    </head>
    <body>

        """)

    if 'cnv' in update_aws_ft:
        out_cnv = open(f'download.cnv.{udn}.sh', 'w')
        # out_other = open(f'download.other.{udn}.sh', 'w')
        out_cnv.write('mkdir origin 2>/dev/null\n')

    other_download = 0
    if len(set(update_aws_ft) - set(['bam', 'fastq', 'cnv'])) > 0:
        other_download = 1
        out_other = open(f'download.other.{udn}.sh', 'w')

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

    out_info.write('rel_to_proband\tfn\turl\tudn\tseq_type\tsize\tbuild\tmd5\turl_s3\n')


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
            md5 = ifl['md5']

            seq_type = ifl.get('seq_type')

            if re.match(r'.+\.bai$', fn) and 'bam' in update_aws_ft:
                out_bam.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
                print(f'<span><b>{rel_to_proband}:    </b></span><a href="{url}">{fn}</a></br></br>', file=html_bam)
            elif re.match(r'.+\.bam$', fn) and 'bam' in update_aws_ft:
                out_bam.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
                print(f'<span><b>{rel_to_proband}:    </b></span><a href="{url}">{fn}</a></br></br>', file=html_bam)
            elif re.match(r'.+\.cnv\.vcf(\.gz)?$', fn) and 'cnv' in update_aws_ft:
                if not os.path.exists(f'origin{fn}'):
                    if url == 'na':
                        logger.error(f'invalid CNV file url')
                    else:
                        out_cnv.write(f'size=$(stat -c %s origin/{fn} 2>/dev/null);if [[ "$size" -lt 1000 ]];then wget "{url}" -c -O "origin/{fn}";\nelse echo already exist: "origin/{fn}";fi\n')
            elif re.match(r'.+\.fastq.gz', fn) and 'fastq' in update_aws_ft:
                out.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
            elif other_download:
                out_other.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')

            out_info.write(f'{rel_to_proband}\t{fn}\t{url}\t{udn}\t{seq_type}\t{size}\t{build}\t{md5}\t{url_s3}\n')
            out_md5.write(f'{md5}  {fn}\n')

    try:
        html_bam.write("""</body>
    </html>""")
        html_bam.close()
    except:
        pass

    try:
        out.close()
    except:
        pass
    try:
        out_bam.close()
    except:
        pass
    out_md5.close()
    out_info.close()

    if 'bam' in update_aws_ft:
        logger.info('building IGV script')
        udn_igv.main('.', udn, logger)

    return driver

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
        udns = re.findall(r'UDN\d+', i)
        case = re.search(r'\(Case (\d+)\)', i)
        if not case:
            continue
        proband = f'{udns[0]}_case{case.group(1)}'
        res.append('@'.join([proband, *udns[1:]]))
    return ' '.join(res)


def upload_files(nodename, pw_accre_data, pw_accre_scratch, udn_raw, rename=False, udn_raw_old=None, initial=False):

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

        cmd = f"""sftp {nodename} <<< $'mkdir {pw_accre_data}/{udn_raw}/origin/\\nmkdir {pw_accre_scratch}/{udn_raw}/\\nput -r * {pw_accre_data}/{udn_raw}/\\n{fls}'  2>/dev/null >&2; touch {initial_upload_flag}"""
        os.system(cmd)
    else:
        fls = []  # file name can't contain space
        fls.append(f'mkdir {pw_accre_data}/{udn_raw}/')
        fls.append(f'mkdir {pw_accre_scratch}/{udn_raw}/')

        for ifl in upload_file_list:
            fls.append(f'put {ifl} {pw_accre_data}/{udn_raw}/')
            fls.append(f'put {ifl} {pw_accre_scratch}/{udn_raw}/')
        fls = '\\n'.join(fls)
        logger.info(f'update {upload_file_list} to  {nodename}: {udn_raw}')
        os.system(f"""sftp {nodename} <<< $'{fls}' >/dev/null 2>/dev/null""")


if __name__ == "__main__":

    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('udn', help="""one or multiple UDN ID, sep by space, could also contain the prefix, eg. 11_UDN123456.  if only selected family number are needed, append with @, eg. UDN123@UDN456_father@UDN_789@mother, in this case, other family member would not be resolved""", nargs='*')
    ps.add_argument(
        '-fn_cred', '-cred',
        help="""filename for the credential file, must include the UDN token, login email, login username(vunetID), login_password""")
    ps.add_argument('-create', '-encry', '-enc', help="""enter into create credential mode""", action='store_true')
    ps.add_argument('-sv_caller', '-svcaller', '-svsoftware', help="""the software used for sv calling, default=dragen, alt=pacbio""", default='dragen', choices=['pacbio', 'dragen'])
    ps.add_argument('-update_cred', '-uc',  help="""update the credential info""", action='store_true')
    ps.add_argument('-demo', help="""do not resolve the amazon link, just check the raw link. apply to both post and new query""", action='store_true')
    ps.add_argument('-ft', help="""could be multiple, if specified, would only update the amazon URL for these file types, if gzip only, add the .gz suffix. e.g. vcf.gz  would not match .vcf$""", nargs='*', default=None)
    ps.add_argument('-renew', '-new', '-update', help="""flag, update the amazon shared link. default is not renew""", action='store_true')
    ps.add_argument('-lite', help="""flag, donot download the cnv files and the bayler report""", action='store_true')
    ps.add_argument('-noupload', '-noup', help="""don't upload the files to ACCRE""", action='store_true')
    ps.add_argument('-force_upload', '-fu', help="""force upload all files to ACCRE""", action='store_true')
    ps.add_argument('-show', '-chrome', '-noheadless', help="""disable headless mode""", action='store_true')
    ps.add_argument('-showcred', help="""show the credential content and exit""", action='store_true')
    ps.add_argument('-reltable', '-rel', help="""force save the relative table, even in the parse_api_result mode""", action='store_true')
    ps.add_argument('-v', '-level', help="""set the logger level, default is info""", default='INFO', choices=['INFO', 'WARN', 'DEBUG', 'ERROR'])
    ps.add_argument('-nodename', '-node', '-remote', help="remote node name, default is va", choices=['va', 'pc'])
    args = ps.parse_args()
    print(args)

    udn_list = args.udn or [os.getcwd().rsplit('/')[-1]]

    force_upload = args.force_upload
    sv_caller = args.sv_caller

    save_rel_table = args.reltable
    headless = not args.show
    if platform == 'darwin':
        root = '/Users/files/work/cooperate/udn/cases'
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


    nodename = args.nodename
    pw_accre_data = '/data/cqs/chenh19/udn'
    pw_accre_scratch = '/fs0/members/chenh19/tmp/upload' if nodename == 'va' else '/scratch/h_vangard_1/chenh19/udn/upload'
    # pw_accre_upload = pw_accre_scratch

    upload_file_list = ['pheno.keywords.txt', 'download.*.txt', 'download.*.md5']  # upload these files to scratch and data of ACCRE

    gzip_only = set()
    update_aws_ft = set()

    if args.ft is None:
        update_aws_ft = ['cnv', 'fastq']
    else:
        err = 0
        args.ft = [_ for _ in args.ft if _.strip()]
        for i in args.ft:
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

    api_token_gateway = cred['token']
    api_token_fileservice = cred['fs_token']

    if args.showcred:
        logger.info(cred)
        sys.exit(1)


    udn_list = [_.strip().upper() for _ in udn_list if _.strip()]

    # validate udn
    p = re.compile(r'(?:.*)?(UDN\d+)')
    validated = []
    for i in udn_list:
        i = [_.strip() for _ in i.split('@')]
        udn_raw = i[0]

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
    driver = None
    cookie_token = None

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


    for udn_raw, udn, valid_family in validated:
        rename_remote = False
        udn_raw_old = ''
        if udn in udn_exist:
            udn_raw_old = udn_exist[udn].rsplit('/', 1)[-1]
            print(f'mv {root}/{udn_raw_old} {root}/{udn_raw}')
            if udn_raw != udn_raw_old:
                logger.warning(f'rename the folder, old = {udn_raw_old} , new = {udn_raw}')
                os.system(f'mv {root}/{udn_raw_old} {root}/{udn_raw}')
                rename_remote = True
        else:
            os.system(f'mkdir -p {root}/{udn_raw}/origin 2>/dev/null')


        pw_case = f'{root}/{udn_raw}'
        os.chdir(pw_case)

        logger = get_logger(f'{root}/{udn_raw}/{udn_raw}', level=args.v)
        logger.info(os.getcwd())

        fn_udn_api_pkl = f'{root}/{udn_raw}/{udn}.udn_api_query.pkl'

        if os.path.exists(fn_udn_api_pkl):
            logger.info('directly load API query result from pickle file')
            res = pickle.load(open(fn_udn_api_pkl, 'rb'))
            if demo:
                renew_amazon_link = False
            else:
                renew_amazon_link = args.renew

            driver = parse_api_res(res, cookie_token=cookie_token, update_aws_ft=update_aws_ft, renew_amazon_link=renew_amazon_link, demo=demo, udn_raw=udn_raw, valid_family=valid_family, gzip_only=gzip_only, driver=driver, sv_caller=sv_caller)

            with open(fn_udn_api_pkl, 'wb') as out:
                pickle.dump(res, out)
        else:
            res, driver = get_all_info(udn, cookie_token=cookie_token, demo=demo, get_aws_ft=update_aws_ft, udn_raw=udn_raw, valid_family=valid_family, lite_mode=lite, udn_proband=udn, gzip_only=gzip_only, driver=driver, sv_caller=sv_caller)
            with open(fn_udn_api_pkl, 'wb') as out:
                pickle.dump(res, out)

        if not lite:
            logger.info(f'downloading CNV vcf file')
            os.system(f'bash download.cnv.{udn}.sh')


        initial_upload_flag = 'origin/initial_uploaded'
        if upload and os.path.exists(fn_udn_api_pkl) and not force_upload and os.path.exists(initial_upload_flag):
            initail_upload_toggle = False
        else:
            initail_upload_toggle = True

        if upload:
            upload_files(nodename, pw_accre_data, pw_accre_scratch, udn_raw, rename=rename_remote, udn_raw_old=udn_raw_old, initial=initail_upload_toggle)

        print('\n########################\n\n')
