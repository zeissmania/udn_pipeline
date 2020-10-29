#! /usr/bin/env python3
"""
this is to get the files from UDN gateway
https://documenter.getpostman.com/view/1367615/RVu1JBWH#0ab3e9c3-792c-b292-6615-761a505df5bb

the input UDN should be the proband UDN, don't add the relative UDN ( they would be get based on the proband)
the output pickle would be
res
    proband
        udn_name
        gender
        age
        files
            file1
            file2
            file3
                file_att1
                file_att2
                file_att3
    rellative1
    relative2
    relative3




"""
import sys
import os
import time
import re
import pickle
import hashlib
from getpass import getpass
import base64
from cryptography.fernet import Fernet
import cryptography.fernet as fernet
import logging
import requests
from bs4 import BeautifulSoup as bs
# basic settings
base_url = 'https://gateway.undiagnosed.hms.harvard.edu/api'


class Credential():
    def __init__(self, password, fn_pickle=None):
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


def get_cookie():
    """
    get cookie using selenium
    """
    url = 'https://gateway.undiagnosed.hms.harvard.edu/login/dashboard/'
    timeout = 30

    from selenium import webdriver
    from selenium.webdriver.support.ui import WebDriverWait as wait
    path_chromedriver = '/Users/files/work/package/chromedriver'
    option = webdriver.ChromeOptions()
    option.add_argument('--headless')
    option.add_argument('--no-sandbox')
    driver = webdriver.Chrome(executable_path=path_chromedriver, options=option)

    driver.get(url)
    email = wait(driver, timeout).until(lambda x: x.find_element_by_xpath('//input[@type="email"]'))
    password = driver.find_element_by_xpath('//input[@type="password"]')
    button = driver.find_element_by_xpath('//button[@class="auth0-lock-submit"]')
    # driver.save_screenshot('udn_login_pre_sso.png')

    run_sso = False
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

    try:
        password.click()
        password.send_keys(cred['password'])
    except webdriver.remote.errorhandler.ElementNotVisibleException:
        run_sso = True

    button.click()

    if run_sso:
        username = wait(driver, timeout).until(lambda x: x.find_element_by_xpath('//input[@id="username"]'))
        password = driver.find_element_by_xpath('//input[@id="password"]')
        button = driver.find_element_by_xpath('//a[@class="ping-button normal allow"]')

        username.click()
        username.send_keys(cred['vunetID'])

        password.click()
        password.send_keys(cred['password'])

        button.click()

    for _ in range(10):
        if driver.current_url.find('gateway.undiagnosed.hms.harvard.edu') > -1:
            break
        time.sleep(1)
    else:
        logger.error('fail to login to UDN gateway')
        return 0

    cookies = driver.get_cookies()
    cookies_list = [f"{_['name']}={_['value']}" for _ in cookies]
    cookie = '; '.join(cookies_list)

    return cookie


def cleanup_memory():
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


def get_all_info(udn, cookie, rel_to_proband=None, res_all=None, info_passed_in=None):
    """
    rel_to_proband, if run as proband, this is None, otherwize, this func is called during geting the relatives of proband, and specified
    info_passed_in,  {'affected': rel_aff, 'seq_status': have_seq}   include the affect state, sequenced_state, only valid when called during geting the relatives of proband, and specified
    res_all,  when running specific relative, use this dict to accumulate among the family member
    """

    aff_convert = {'Unaffected': 'Unaffected', 'Affected': 'Affected', '3': 'Unaffected', '2': 'Affected', 'Unknown': 'Unknown'}
    res_all = res_all or {}
    rel_to_proband = rel_to_proband or 'Proband'
    info_passed_in = info_passed_in or {}

    res = {}
    res['rel_to_proband'] = rel_to_proband
    res['affect_from_proband_list'] = info_passed_in.get('affected') or 'Unknown'
    res['seq_status'] = info_passed_in.get('seq_status')

    cookie = cookie or get_cookie()
    header_cookie = get_header_cookie(cookie)
    if not cookie:
        logger.error('fail to get cookie, can\'t get the amazon s3 presigned URL, exit')
        return 0

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
        else:
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
        res['gender'] = json['patient']['gender']['name']
        res['firstname'] = json['patient']['patientfirst']
        res['lastname'] = json['patient']['patientlast']
        try:
            res['race'] = json['patient']['race'][0]['name']  # 1=male
        except:
            res['race'] = json['patient']['race']
        res['uuid'] = json['patient']['uuid']  #
        res['uuid_case'] = json['uuid']  #
        res['dob'] = json['patient']['dob']  # birthday

        affect_state_raw = json['patient']['affected'] # 2 = affected, 3=unaffected
        res['affect'] = affect_state_raw


        if res['affect_from_proband_list'] == 'affect_state_raw':
            res['affect_final'] = affect_state_raw
        else:
            res['affect_final'] = f'{affect_state_raw} / {res["affect_from_proband_list"]}'

        res['phenotip'] = json['patient']['phenotipsid']
        res['alive'] = not json['deceased']
        res['comment'] = json['comment']
        res['simpleid'] = json['patient']['simpleid']
        res['similar_symp'] = json['patient'].get('similarsymptomsexplain')

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
            response_relative = requests.request('GET', url_relative, headers=header_cookie)
        except:
            logger.error(f'fail to get relatives info: {udn}')
            return 0

        r = bs(response_relative.text, features='lxml')
        tb = r.find('tbody').find_all('tr')
        if len(tb) == 0:
            logger.error('Relative list is empty')
            return 0

        # for this part, even the sequence state is noted as "-"  not uploaded, when clicking the actual case, sometimes, the sequence file still exist
        for i in tb:
            try:
                seq_status = i.find(attrs={'class': 'sequploaded'})
                have_seq = seq_status.text
            except:
                logger.warning(f'fail to get the sequence uploaded info from relative: {i}')
                continue

            rel_to_proband_tmp = i.find(attrs={'class': 'relative'}).text
            rel_udn = i.find(attrs={'class': 'patientsimpleid'}).text
            rel_aff = i.find(attrs={'class': 'affected'}).text

            res_all = get_all_info(rel_udn, cookie, rel_to_proband=rel_to_proband_tmp, res_all=res_all, info_passed_in={'affected': rel_aff, 'seq_status': have_seq})
            logger.info(res_all.keys())

        # patient dignosis
        action = 'application'
        json = get_json(udn, action)

        d_json['doctor_review'] = json
        dump_json(d_json, fn_json_pk)

        if json:
            review = json[0]['application_review_set']
            if len(review) > 0:
                res['summary'] = review[0]['summary']
                res['evaluation'] = review[0]['evaluations']

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

    logger.info(f'now getting information of {rel_to_proband} : {udn}')

    if json:
        # json is a list
        if len(json) > 1:
            logger.warning(f'{len(json)}  records found for sequences API, would use the first one')
        json = json[0]
        res['bayler_report'] = [
            f'https://gateway.undiagnosed.hms.harvard.edu/patient/sequencing-report/{_}/'
            for _ in json['sequencereports']]

        files_raw = json['sequencingfiles']
        files = []

        for n_fl, ifl in enumerate(files_raw):
            if n_fl % 5 == 1 and n_fl > 1:
                time.sleep(5)
            completed = ifl['complete']
            fl_url = ifl['file_data']['locations'][0]['url']
            fl_size = ifl['file_data']['locations'][0]['filesize']
            fn = ifl['file_data']['filename']
            file_uuid = ifl['uuid']
            fl_expire = ifl['file_data']['expirationdate']
            fl_assembly = ifl['file_data']['metadata'].get('assembly')
            fl_md5 = ifl['file_data']['metadata'].get('md5sum')
            fl_type = ifl['file_data']['metadata'].get('description')

            if re.match(r'.+\.bai$', fn):
                pass
            elif re.match(r'.+\.bam$', fn):
                pass
            elif re.match(r'.+\.vcf', fn):
                pass
            elif re.match(r'.+\.fastq.gz$', fn):
                pass
            else:
                logger.warning(f'unkown filetype: {fn}')
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

            res_amazon_url = get_amazon_download_link(fn, file_uuid, header_cookie)
            try:
                header_cookie, download = res_amazon_url
            except:
                logger.error(f'error when getting_amazon_download_link return = {res_amazon_url} \nfn="{fn}"\nfile_uuid="{file_uuid}"\nheader_cookie="{header_cookie}"')
                sys.exit(1)

            files.append({'download': download, 'relative': rel_to_proband, 'file_uuid': file_uuid, 'fn': fn, 'url': fl_url, 'complete': completed, 'size': fl_size,
                          'expire': fl_expire, 'build': fl_assembly, 'assembly': fl_assembly, 'md5': fl_md5, 'type': fl_type})
        res['files'] = files
    else:
        logger.error(f'fail to get sequencing json')

    res_all[rel_to_proband] = res

    if rel_to_proband == 'Proband':
        try:
            parse_api_res(res_all, cookie=cookie)
        except:
            logger.error(f'fail to parse the result dict, try again')

    return res_all


def get_amazon_download_link(fn, file_uuid, header_cookie, n=0):
    url_download_info = f'https://gateway.undiagnosed.hms.harvard.edu/patient/downloadsequenceurl/{file_uuid}/'
    response_download_link = requests.request('GET', url_download_info, headers=header_cookie)
    if n > 20:
        logger.error(f'fail to get file download link: {fn}')
        return header_cookie, None
    if n > 5 and n % 5 == 1:
        logger.info(f'renew cookie: try: {n}')
        cleanup_memory()
        cookie = get_cookie()
        header_cookie = get_header_cookie(cookie)
    try:
        json = response_download_link.json()
    except Exception as e:
        if str(e.__class__).find('JSONDecodeError') > -1:
            # logger.info(f'not a json response,  url_download_info= "{url_download_info}", headers={header_cookie}')
            n += 1
            time.sleep(10)
            # logger.info(header_cookie)
            header_cookie, download_url = get_amazon_download_link(fn, file_uuid, header_cookie, n)
            if download_url:
                return header_cookie, download_url
        else:
            raise
    else:
        logger.info(f'got json: {fn}')
        try:
            return header_cookie, json['download_url']
        except:
            logger.error(f'{fn}: json does not contain the download link: {json}')
            return header_cookie, None

def get_header_cookie(cookie):
    return {'Connection': 'keep-alive', 'Pragma': 'no-cache', 'Cache-Control': 'no-cache', 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.162 Safari/537.36 Edg/80.0.361.109', 'Cookie': cookie}


def get_logger(prefix=None):
    prefix = prefix or 'udn_api_query'
    fn_log = f'{prefix}.log'
    fn_err = f'{prefix}.err'

    fmt = logging.Formatter(
        '%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    console.setLevel('INFO')

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


def parse_api_res(res, cookie=None, renew_amazon_link=False, pkl_fn=None):
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
    pw = os.getcwd().rsplit('/', 1)[-1]
    pw = f'/data/cqs/chenh19/udn/{pw}'
    if not os.path.isdir('origin'):
        os.system('mkdir origin 2>/dev/null')

    logger.info(f'keys for the result: {res.keys()}')
    proband_id = res['Proband']['simpleid']


    if renew_amazon_link:
        cookie = cookie or get_cookie()
        header_cookie = get_header_cookie(cookie)
        if not cookie:
            logger.error('fail to get cookie, can\'t get the amazon s3 presigned URL, exit')
            return 0


        for rel_to_proband, v in res.items():
            logger.info(f'renew amazon download link for {rel_to_proband}')

            if 'files' not in v:
                logger.info(f'{rel_to_proband}: no files found')
                continue
            for ifl in v['files']:
                fn = ifl['fn']
                file_uuid = ifl['file_uuid']

                res_amazon_url = get_amazon_download_link(fn, file_uuid, header_cookie)
                try:
                    header_cookie, download = res_amazon_url
                except:
                    logger.error(f'get_amazon_download_link return = {res_amazon_url} \nfn="{fn}"\nfile_uuid="{file_uuid}"\nheader_cookie="{header_cookie}"')
                    continue
                ifl['download'] = download

        pkl_fn = pkl_fn or f'{proband_id}.udn_api_query.pkl'
        if os.path.exists(pkl_fn) and os.path.getsize(pkl_fn) > 5000:
            os.system(f'mv {pkl_fn} {pkl_fn}.bk')
        pickle.dump(res, open(pkl_fn, 'wb'))

    # report for proband

    proband = res['Proband']
    udn = proband['simpleid']
    hpo = proband['symp']

    aff_convert = {'Affected': 1, 'Unaffected': 0, '3': 0, '2': 1}


    with open(f'{udn}.basic_info.md', 'w') as out:
        out.write(f'## 1. Basic Info\n- UDN\t{udn}\n')
        for k, v in zip('firstname,lastname,gender,race,dob,alive,affect,uuid_case,phenotip'.split(','), 'FirstName,LastName,Gender,Race,DateBirth,Alive,Affect_state,UUID_case,Phenotip_ID'.split(',')):
            out.write(f'- {v}\t{proband[k]}\n')
        out.write('\n\n')

        # family members

        out.write('## 2. Family Members\n\n')
        out.write('| Relative | Name | UDN | Gender | Affected| Sequence_Uploaded | \n')
        out.write('| :----: | :----: | :----: | :----: | :----: | :----: |\n')
        for rel_to_proband, i in res.items():
            if rel_to_proband == 'Proband':
                continue


            affect_state1 = aff_convert.get(i['affect_from_proband_list']) or i['affect_from_proband_list']
            affect_state2 = aff_convert.get(i['affect']) or i['affect']

            affect_state = affect_state1 if affect_state1 == affect_state2 else f'{affect_state1}/{affect_state2}'

            # if gender not added
            out.write(f'| {rel_to_proband} | {i["firstname"]} {i["lastname"]} | {i["simpleid"]} | {i["gender"]} | {affect_state}| {i["seq_status"]} |\n')
        out.write('\n\n')

        # HPO terms
        out.write(f'## 3. HPO terms\n')

        for k, v in hpo:
            out.write(f'- {k}\t{v}\n')
        out.write('\n\n')

        for k, v in hpo:
            out.write(f'- {v}\n')
        out.write('\n\n')


        # symp,uuid_case
        out.write('## 4. Summary and diagnosis\n')
        for k, v in zip('similar_symp,comment,summary,evaluation'.split(','), 'Similar_Sympt,Comment,Summary,Evaluation'.split(',')):
            out.write(f'### {v}\n{proband[k]}\n')
        out.write('\n\n')


    # build the HPO terms file
    with open(f'origin/{udn}.terms.txt', 'w') as out:
        for k, v in hpo:
            out.write(v+'\n')

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
            tmp_set = set(['simpleid', 'affect', 'gender'])
            logger.info(f'{rel_to_proband}: key not found: {[_ for _ in tmp_set if _ not in irel]}')
            continue

        if not irel.get('files'):
            logger.info(f'{rel_to_proband}: no files available, skip... ')
            continue

        try:
            affect_state1 = aff_convert.get(irel['affect_from_proband_list']) or irel['affect_from_proband_list']
            affect_state2 = aff_convert.get(irel['affect']) or irel['affect']
            rel_aff = affect_state1 if affect_state1 == affect_state2 else f'{affect_state1}/{affect_state2}'

        except:
            logger.error(f"{rel_to_proband}: unkown affect state: {irel['affect']} / {irel['affect_from_proband_list']}")
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

    # brothers = '\n'.join([f'brother{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['brother'])])
    # sisters = '\n'.join([f'sister{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['sister'])])
    male = '\n'.join([f'male{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['male'])])
    female = '\n'.join([f'female{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['female'])])
    other = '\n'.join([f'other{n+1}: {str(_)}' for n, _ in enumerate(cfg_info['other'])])


    cfg = f"""prj: {udn}  # the project name
path: {pw}
# if the file is under the path above, could use the filename only
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
    out = open(f'download.fastq.{udn}.txt', 'w')
    out_fastq_sh = open(f'download.fastq.{udn}.sh', 'w')
    out_bam = open(f'download.bam.{udn}.download_bam.sh', 'w')

    out_info = open(f'download.info.{udn}.txt', 'w')
    out_cnv = open(f'download.cnv.{udn}.sh', 'w')
    out_other = open(f'download.other.{udn}.sh', 'w')
    out_md5 = open(f'download.{udn}.md5', 'w')
    out_igv = open(f'{udn}.igv.files.txt', 'w')
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

    out_cnv.write('mkdir origin 2>/dev/null\n')

    pw_upload = f'/scratch/h_vangard_1/chenh19/udn/upload/{udn}'
    out_fastq_sh.write(f"""#! /usr/bin/env bash
mkdir -p {pw_upload}
cp download.fastq.{udn}.txt {pw_upload}
cp download.{udn}.md5 {pw_upload}
cd {pw_upload}
parallel '{{}}'  :::: {udn}.download_fastq.txt
cd -
    """)

    out_info.write('rel_to_proband\tfn\turl\tsize\tbuild\tmd5\turl_s3\n')

    for rel_to_proband, irel in res.items():
        if not irel.get('files'):
            continue
        for ifl in irel['files']:
            url = ifl['download']
            fn = ifl['fn']
            url_s3 = ifl['url']
            size = ifl['size']
            build = ifl['build']
            md5 = ifl['md5']

            if re.match(r'.+\.bai$', fn):
                out_igv.write(f'{url}\t{fn}\n')
                out_bam.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
                print(f'<span><b>{rel_to_proband}:    </b></span><a href="{url}">{fn}</a></br></br>', file=html_bam)
            elif re.match(r'.+\.bam$', fn):
                out_igv.write(f'{url}\t{fn}\n')
                out_bam.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
                print(f'<span><b>{rel_to_proband}:    </b></span><a href="{url}">{fn}</a></br></br>', file=html_bam)
            elif re.match(r'.+\.cnv\.vcf(\.gz)?$', fn):
                out_cnv.write(f'wget "{url}" -c -O "origin/{fn}"\n')
            elif re.match(r'.+\.fastq.gz', fn):
                out.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')
            else:
                out_other.write(f'nohup wget "{url}" -c -O "{fn}" > {fn}.download.log 2>&1 &\n')

            out_info.write(f'{rel_to_proband}\t{fn}\t{url}\t{size}\t{build}\t{md5}\t{url_s3}\n')
            out_md5.write(f'{md5}  {fn}\n')

    html_bam.write("""</body>
</html>""")

    html_bam.close()
    out.close()
    out_bam.close()
    out_md5.close()
    out_info.close()
    out_igv.close()
    out_fastq_sh.close()

if __name__ == "__main__":
    root = os.getcwd()
    logger = get_logger()

    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('udn', help="""one or multiple UDN ID, sep by space, could also contain the prefix, eg. 11_UDN123456""", nargs='*')
    ps.add_argument(
        '-fn_cred', '-cred',
        help="""filename for the credential file, must include the UDN token, login email, login username(vunetID), login_password""")
    ps.add_argument('-create', '-encry', '-enc', help="""enter into create credential mode""", action='store_true')
    ps.add_argument('-pw', help="""specify the output pw, default is create a folder same as UDN_ID""", default=None)
    ps.add_argument('-no_renew', '-norenew', help="""flag, if the pickle file already exist, renew the amazon link or not, default is renew""", action='store_true')
    args = ps.parse_args()
    fn_cred = args.fn_cred or 'udn.credential.pkl.encrypt'

    passcode = getpass('Input the passcode for the encrypted credential file: ')
    if not fn_cred:
        fn_cred = input('Specify the credential file name, if not exist, would create it: ')
    key = Credential(passcode, fn_cred)

    logger.info(fn_cred)


    if args.create:
        cred = create_cred()
    elif not os.path.exists(fn_cred):
        logger.warning(f'The credential file not found: {args.fn_cred}, would create it')
        cred = create_cred()
    else:
        cred = key.load()

    # unpack the cred
    update_cred = False
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


    cookie = get_cookie()
    header1 = {
        'Content-Type': 'application/json',
        'Authorization': f'Token {token}',
    }

    udn_list = args.udn
    udn_list = [_.strip().upper() for _ in udn_list if _.strip()]

    # validate udn
    p = re.compile(r'(?:.*)?(UDN\d+)$')
    validated = []
    for udn_raw in udn_list:
        m = re.match(p, udn_raw)
        if not m:
            logger.warning(f'Invalid UDN_ID: {udn_raw}')
            continue
        udn = m.group(1)
        validated.append([udn_raw, udn])

    for udn_raw, udn in validated:
        os.system(f'mkdir -p {root}/{udn_raw}/origin 2>/dev/null')
        os.chdir(f'{root}/{udn_raw}')

        logger = get_logger(f'{root}/{udn_raw}/{udn}.udn_api')
        logger.info(os.getcwd())

        fn_udn_api_pkl = f'{root}/{udn_raw}/{udn}.udn_api_query.pkl'

        if os.path.exists(fn_udn_api_pkl):
            logger.info('directly load API query result from pickle file')
            res = pickle.load(open(fn_udn_api_pkl, 'rb'))
            renew_amazon_link = not args.no_renew

            parse_api_res(res, cookie=cookie, renew_amazon_link=renew_amazon_link)
        else:
            res = get_all_info(udn, cookie=cookie)
            with open(fn_udn_api_pkl, 'wb') as out:
                pickle.dump(res, out)
