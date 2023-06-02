import requests
import time
from bs4 import BeautifulSoup as bs
from selenium.webdriver.support.ui import WebDriverWait as wait
from chctool_lite import getlogger

global_flag = {'omim_blocked': 0}

logger_name = 'omim_utils'

pw_code = os.path.dirname(__file__)

def run_omim_scrapy(gn, gene_id, logger=None, driver=None):
    logger = logger or getlogger(logger_name)
    start = time.time()
    v = [gene_id, []]
    gene_id = str(int(float(gene_id)))
    driver, tmp = omim_gn_to_pheno(driver, gn, gene_id, logger)
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


def get_omim_map(logger=None):
    """
    input = genemap2.txt
    return = {pheno_id1: pheno_desc1}
    """
    logger = logger or getlogger(logger_name)

    if os.path.exists(f'{pw_code}/omim_map.pkl'):
        try:
            d_omim_map = pickle.load(open(f'{pw_code}/omim_map.pkl', 'rb'))
        except:
            pass
        else:
            return d_omim_map
    if os.path.exists(f'{pw_code}/genemap2.txt'):
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


def run_omim_scrapy_old(gn, gn_omim_id, logger, res_prev=None, driver=None):
    """
    get the omim disease discription
    return key=pheno_id, value=the description and clinicalfeatures text
    """
    headers = {'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/87.0.4280.67 Safari/537.36 Edg/87.0.664.52', 'sec-fetch-site': 'cross-site', 'origin': 'https://www.omim.org', 'referer': 'https://fonts.googleapis.com/', 'accept': '*/*', 'accept-encoding': 'gzip, deflate, br'}
    logger = logger or getlogger(logger_name)

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


def get_pheno_detail(driver, pheno_id, pheno_title, logger):
    """
    return
    None, there is error to get the description
    str,  the whole description
    """
    logger = logger or getlogger(logger_name)

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


def omim_gn_to_pheno(driver, gn, gene_id, logger):
    """
    return
    None - error, page not loaded
    'empty',  no linked phenotype for this gene
    list, element = [pheno_id, pheno_name, merged_title]
    """
    logger = logger or getlogger(logger_name)

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


def run_selenium(driver, url_omim, gn, gn_omim_id, logger):
    logger = logger or getlogger(logger_name)

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



def get_driver(driver, logger):
    logger = logger or getlogger(logger_name)

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

def parse_pheno_file(fn):
    logger 
    if not os.path.exists(fn):
        logger.error(f'pheno file not exist: {fn}')
        return None
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

    return l_pheno, l_pheno_raw
    
def omim_query(gene_list, fn_pheno, logger=None):
    """
    the input is the genelist  and the phenotype keyword file
    1. build the gene comment dict based on the previous manual input
    2. read the local OMIM disease discription, if not found, would run the online scrapy
    3. query the keyword file for each gene in the genelist, and rate(hit number) as reference for selecting

    the phenotype file, accept keyword in 3 mode
    word1+word2  this 2 word come together
    @word1 match exact word, would add \b at both side
    -word  do not include this word
    """
    logger = logger or getlogger(logger_name)
    driver=None

    tmp = parse_pheno_file(fn_pheno)
    if not tmp:
        return None
    pheno_list += tmp


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

