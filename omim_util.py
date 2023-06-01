
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


