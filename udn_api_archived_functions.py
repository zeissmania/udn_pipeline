
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


def get_all_info_old(udn, cookie_token, rel_to_proband=None, res_all=None, info_passed_in=None, sequence_type_desired=None, demo=False, get_aws_ft=None, udn_raw=None, valid_family=None, lite_mode=None, udn_proband=None, gzip_only=None, driver=None, sv_caller='dragen', newname_prefix=None):
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
    fn_json_pk = f'intermed/{udn}.json.pkl'
    fn_json_pk_bk = f'intermed/{udn}.json.pkl.bk'
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

        with open(f'intermed/phenotip.{udn}.pdf', 'wb') as out:
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
                driver.save_screenshot(f'relatives_table.{udn}.png')
            except:
                logger.warning(f'{udn} : This case seems donnot have any relatives info !, check screenshot')
                driver.save_screenshot(f'relative_not.found.{udn}.png')
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

    # sequence_type_desired=None, demo=False, get_aws_ft=None, udn_raw=None, valid_family=None, lite_mode=None, udn_proband=None, gzip_only=None, driver=None, sv_caller='dragen', newname_prefix=None
            res_all, driver = get_all_info(rel_udn, cookie_token, rel_to_proband=rel_to_proband_tmp, res_all=res_all, info_passed_in={'affected': rel_aff, 'seq_status': have_seq}, sequence_type_desired=sequence_type_desired, demo=demo, get_aws_ft=get_aws_ft, udn_raw=udn_raw, valid_family=valid_family, lite_mode=lite_mode, udn_proband=udn_proband, gzip_only=gzip_only, driver=driver, sv_caller=sv_caller, newname_prefix=newname_prefix)
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

    # json_all is a list
    # each element is a dict, which is a sequencing request like
    # {'id': 1984,
    #   'patient': {'simpleid': 'UDN139529'},
    #   'sequencingfiles': [...]}
    seq_type_convert = {2: 'wes', 3: 'wgs', 4: 'rna', 5: 'reanalysis', 1: 'targeted variant'}
    seq_type_convert.update({str(k): v for k, v in seq_type_convert.items()})
    files = []
    res['bayler_report'] = []


    seq_type_available = [seq_type_convert.get(json['sequencing_type']) or json['sequencing_type'] for json in json_all]

    if len(json_all) > 1:
        logger.warning(f'multi seq type found ({len(seq_type_available)}): {seq_type_available}')
    elif not json or len(json_all) == 0:
        logger.error(f'\t{rel_to_proband}: fail to get sequencing json')

    use_all_seq_res = 0
    if len(set(sequence_type_desired) & set(seq_type_available)) == 0:
        use_all_seq_res = 1

    res['seq_json_all'] = json_all

    for json in json_all:
        seq_id = json["id"]
        res.setdefault('seq_id', []).append(seq_id)

        sequence_type_raw = json['sequencing_type']
        try:
            sequence_type = seq_type_convert[sequence_type_raw]
        except:
            logger.warning(f'{rel_to_proband} seq_id = {seq_id}: unkown sequence type: type code={sequence_type_raw}')
            continue

        if sequence_type in sequence_type_desired  or'all' in  sequence_type_desired or sequence_type == 'reanalysis':
            logger.info(f'\tdesired seq type found: {sequence_type}')
        elif use_all_seq_res:
            logger.warning(f'\tusing data from {sequence_type}')
        else:
            logger.warning(f'\tskipped due to unmatch sequence type, seq_id = {seq_id}, desired={sequence_type_desired} found={sequence_type}')
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
                fl_loc = ifl['file_data']['locations']
                fl_meta = ifl['file_data']['metadata']
            except:
                logger.error(f"wrong file json: {ifl}")
                raise
            else:
                try:
                    fl_loc = fl_loc[0]
                except:
                    fl_loc = 'NA'
                    logger.error(colored(f'file_locations is empty : {ifl}, \nfilename = {ifl.get("filename") or ""}\nurl= https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/{seq_id}/files/#\n\n', 'red'))


            fn = ifl['file_data']['filename']

            if re.match(r'.+\.bai$', fn):
                pass
            elif re.match(r'.+\.bam$', fn):
                pass
            elif re.match(r'.+\.g?vcf', fn):
                pass
            elif re.match(r'.+\.fastq.gz$', fn):
                pass
            else:
                logger.debug(f'unkown filetype: {sequence_type}: {fn}')
                continue

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
                # logger.debug(f'unexpected fl_info {fn}: {fl_loc}')
                fl_url = fl_size = fl_loc

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
            parse_api_res(res_all, cookie_token=cookie_token, udn_raw=udn_raw, valid_family=valid_family, driver=driver, sv_caller=sv_caller, newname_prefix=newname_prefix)
        except:
            logger.error(f'fail to parse the result dict, try again')
            raise
    return res_all, driver

def get_header_cookie(cookie_token):
    # csrftoken=2DIZXSDSd8ffFuwxfOcO3VKCC3W0a4v0
    # sessionid=0lkbrirzd0mz92o8h1ygpc2y3c7daf21
    # mod_auth_openidc_session=e36bf6ee-f612-4b5f-8e55-6f6ea7f3156c

    cookies_list = [f"{k}={v}" for k, v in cookie_token.items()]
    cookie = '; '.join(cookies_list)
    return {'Connection': 'keep-alive', 'Pragma': 'no-cache', 'Cache-Control': 'no-cache', 'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.162 Safari/537.36 Edg/80.0.361.109', 'Cookie': cookie}



def parse_seq_url(res):
    for rel_to_proband, v1 in res.items():
        print(rel_to_proband)
        try:
            urls = [f'https://gateway.undiagnosed.hms.harvard.edu/patient/sequence/{_}/files/' for _ in v1['seq_id']]
        except:
            logger.error(f'seq ID not found in API result, please check')
            break

        print('\t' + '\n\t'.join(urls))



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

        password.send_keys(Keys.RETURN)
        time.sleep(2)
        try:
            button.click()
        except:
            # logger.warning(f'fail to click SignOn button')
            pass

        # driver.save_screenshot('udn_login_sso.png')

    for _ in range(10):
        if driver.current_url.find('gateway.undiagnosed.hms.harvard.edu') > -1:
            break
        time.sleep(1)
    else:
        logger.error(f'fail to login to UDN gateway, url=\n{driver.current_url}')
        sys.exit(1)
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
