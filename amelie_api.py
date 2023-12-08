#! /usr/bin/env python3

import sys
import requests
import os
import gzip
import time
import pickle


url = 'https://amelie.stanford.edu/api'

BASEDIR = os.path.dirname(os.path.realpath(__file__))
# d_hpo = eval(open(f'{BASEDIR}/hpo_id_to_name.pdict').read())
d_hpo = pickle.load(open(f'{BASEDIR}/hpo_id_to_name.pkl', 'rb'))

def create_patient(taskname, genelist, pheno_hpo_list):
    r = requests.post(
        '%s/create_patient' % url,
        # verify=False,
        # headers=headers,
        # for the data,  null, false must all be string
        data={'vcf': False,
            'dominantAlfqCutoff': 0.1,
            'alfqCutoff': 0.5, # min. 0.1; max. 3.0 (percent)
            'filterByCount': 'false',
            'hmctCutoff': 1,
            'alctCutoff': 3,
            'vcfFile': 'null',
            'unaffectedDadVcf': 'null',
            'unaffectedMomVcf': 'null',
            'patientSex': 'null',
            'onlyPassVariants': 'false',
            'filterRelativesOnlyHom': 'false',
            'patientName': taskname,
            'phenotypesArea': "\n".join(pheno_hpo_list),
            'genotypesArea': '\n'.join(genelist),

            })

    return r.json()['patient_id']


def wait_for_processing(patient_id, timeout=120):
    interval = 10
    s = time.time()
    fail = 0
    for _ in range(int(timeout/interval)):
        if fail > 5:
            print(f'max fail time reached: {fail}, would skip AMELIE')
        s1 = time.time()
        print(f'\twait round {_}: {s1-s:.1f}')
        try:
            response = requests.get('%s/patient/%s' % (url, patient_id),)
                                # verify=False)
        except:
            print(f'fail to run requests')
            time.sleep(10)
            fail += 1
            continue
        time_response = time.time() - s1

        try:
            json = response.json()
        except:
            time_consumed = time.time() - s
            print(f'Fail to decode response json from AMELIE: time={time_consumed:.1f} s, response=\n\n{response.text}\n\n')
            fail += 1
            time.sleep(10)
            continue
        else:
            if json['processing'] != 1:
                return response
            time_consumed = time.time() - s
            print(f'\tamelie processing... retry-{_},  time_response={time_response:.1f}, elapsed time={time_consumed:.1f}s. processing={json["processing"]}')
            time.sleep(interval)

    return 0


def parse_json(prj, json_pdict, pw=None, pw_main=None, force=0):
    """
    extract the information from json
    r = json_pdict,   got from response.json(),  or read from the pdict file
    """
    pw = pw or os.getcwd()
    pw_main = pw_main or pw
    if isinstance(json_pdict, dict):
        json_pdict = [json_pdict]
    res_all = []
    for r in json_pdict:
        # useless key = ['created_on', 'email_subscription', 'has_dad', 'has_mom', 'has_user', 'has_variants', 'incidental_findings', 'is_vcf_patient', 'processing', 'swiper_cutoff',]
        # meta info keys = ['patient_id', 'patient_name', 'patient_sex']

        # keys in use:
        #  r ['not_added_genes', 'not_added_hpo_codes', 'gene_data', 'system_list']
        #  r['not_added_genes'] -> str, sep by comma. 'LOC646652, LOC100133077, FAM197Y2, LOC102724219, '
        #  r['not_added_hpo_codes'] -> str or None

        # r['gene_data'][0]
                # gene_name -> str
                # omim_link -> str
                # solutions -> list
                    # item[0]  -> dict  # a single paper
                        # 'hpo_ids'  -> list , all HPOs matched in this paper
                        # 'journal', 'pmid', 'pubmed_year', 'title'   -> str,  the information for the paper
                        # score -> float , the actual amelie score

        # r['system_list']  -> list, each item is a disease group. will group the input HPO into these groups
            # item[0]  -> dict, diseae category 1
                # phenotype -> list, each item is as hpo id-name match
                    # hpo_id -> str
                    # name  -> str
                # system -> str, the category name

        # write the error, not included HPO / gene
        missing_genes = r['not_added_genes']
        missing_hpo = r['not_added_hpo_codes']
        if missing_genes:
            with open(f'{pw_main}/err_amelie.{prj}.genes_not_added.txt', 'w') as out:
                print(missing_genes.replace(', ', '\n'), file=out)

        if missing_hpo:
            with open(f'{pw_main}/err_amelie.{prj}.HPO_not_added.txt', 'w') as out:
                print(missing_hpo.replace(', ', '\n'), file=out)

        # hpo_q = r['system_list'][0]['phenotypes']
        hpo_q_raw = [_['phenotypes'] for _ in r['system_list']]
        hpo_q = []
        for _ in hpo_q_raw:
            hpo_q.extend(_)

        hpo_q = [_['hpo_id']+ "\t" + _['name'] for _ in hpo_q]

        # res['query_hpo'] = hpo_q
        res = []
        genes = r['gene_data']
        for ign in genes:
            genename = ign['gene_name']
            ires = [genename]
            ires.append(ign.get('omim_link'))
            max_score = -1
            for paper in ign['solutions']:
                try:
                    del paper['id']
                    del paper['authors']
                except:
                    pass
                ihpo = []
                if paper['score'] > max_score:
                    max_score = paper['score']
                for _ in paper['hpo_ids']:
                    ihpo.append('')
                    ihpo.append(_)
                    ihpo.extend(d_hpo[_])
                paper['hpo'] = ihpo
                try:
                    del paper['hpo_ids']
                except:
                    pass

                ires.append(paper)
            ires.insert(0, max_score)
            res.append(ires)

        res = sorted(res, key=lambda x: x[0], reverse=True)
        res_all.extend(res)

    with open(f'{pw_main}/{prj}.amelie.lite.txt', 'w') as out:
        for ign in res_all:
            print(f'{ign[1]}\t{ign[0]:.2f}', file=out)

    # export the tsv file

    with gzip.open(f'{pw}/{prj}.amelie.parsed.raw.gz', 'wb') as out:
        for i in res_all:
            score, gn, omim, *papers = i
            out.write(f'{gn}\t{score}\n'.encode())
            for ipaper in papers:
                title = ipaper['title']
                pmid = ipaper['pmid']
                hpo = ipaper['hpo']
                out.write(f'\t{pmid}\t{title}\n'.encode())
                out.write(f'\t{hpo}\n\n'.encode())

    with open(f'{pw}/{prj}.amelie.parsed.pkl', 'wb') as out:
        pickle.dump(res_all, out)

    return res_all

def re_group_hpo(hpo_list):
    """
    currently, the hpo is in a same list, no hierarchy,
    like below. this func would return a dict, key = hpo
    ['',
 'HP:0005484',
 'Postnatal microcephaly',
 'Acquired microcephaly',
 'Deceleration of head growth',
 'Development of small head that was not present at birth',
 'Microcephaly, acquired',
 'Microcephaly, postnatal',
 'Postnatal deceleration of head circumference',
 'Secondary microcephaly',
 '',
 'HP:0012762',
 'Cerebral white matter atrophy',
 '',
 'HP:0002510',
 'Spastic tetraplegia',
 'Spastic quadriplegia',
 '',
 'HP:0001371',
 'Flexion contracture',
 'Contracture',
 'Contractures',
    """
    res = {}
    collect = 1
    for i in hpo_list:
        if i[:3] == 'HP:':
            hpoid = i
            res[i] = []
            collect = 1
        elif i and collect:
            res[hpoid].append(i)
        else:
            collect = 0
    return res

def query(prj, pheno, pdict, pw, pw_main):
    """
    query the phenotype, find the match
    pheno should be a list, put all the keywords in it
    pdict is the paresed amelie pdict
    """
    pheno = [_.lower() for _ in pheno]
    res = {}  # key = gene
    for ires in pdict:
        score, ign, omim = ires[:3]
        hit = {}  # key = phenotype
        for ipaper in ires[3:]:
            # ipaper is a dict, keys = ['journal', 'pmid', 'pubmed_year', 'score', 'title', 'hpo']
            pmid = ipaper['pmid']
            pmid_url = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
            hpo = re_group_hpo(ipaper['hpo'])  # a list

            for _, hpo_des in hpo.items():
                for ipheno in pheno:
                    for _ in hpo_des:
                        _ = _.lower()
                        if _.find(ipheno) > -1:
                            tmp = [pmid_url, ipaper['title'], f"{ipaper['journal']} ({ipaper['pubmed_year']})", f'PMID: {pmid}', hpo_des]
                            try:
                                hit[ipheno].append(tmp)
                            except:
                                hit[ipheno] = [tmp]
                            break
        if len(hit) > 0:
            res[ign] = [score, omim, hit]
        elif score > 10:
            print(f'\tWARNING: high amelie score, but phenotype not found in AMELIE description: gene={ign}, score={score}')


    # struc =
    # {gene1: 0=score, 1=omim, 2 = hits(below)
    #         pheno1: [
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #         ]
    #     }

    with open(f'{pw}/{prj}.amelie.matched_query.pkl', 'wb') as out:
        pickle.dump(res, out)


    # export into text file
    with open(f'{pw_main}/{prj}.amelie.matched_query.txt', 'w') as out:
        for ign, l1 in res.items():
            hits = l1[2]
            print('\n'.join([ign, str(l1[0]), l1[1]]), file=out)
            for ipheno, ihit in hits.items():
                for ipaper in ihit:
                    # pmid, title, journal, year, hpo_detail = ihit
                    hpo_detail = ipaper[4]
                    term_per_line = 10
                    hpo_detail_new = []
                    for _ in range(int(len(hpo_detail)/term_per_line)+1):
                        tmp = [str(i) for i in hpo_detail[_*term_per_line:(_ + 1)*term_per_line]]
                        hpo_detail_new.append('\t'.join(tmp))
                    hpo_detail = [f'{ign}\t{ipheno}\t{_}' for _ in hpo_detail_new]
                    print('\t' + '\n\t'.join(ipaper[:4]), file=out)
                    print('\t\t' + '\n\t\t'.join(hpo_detail), file=out)
                    print('\t' + '*' *30, file=out)
                    print('', file=out)
    # query(['intellect', 'palate', 'chin', 'hypotonia', 'seizure', 'speech'], d)
            print('#' * 30, file=out)
            print('\n\n', file=out)
    return res

def get_response(prj, genelist, pheno_hpo_list):
    timeout = 120
    print('\tCreating Patient ID')
    patient_id = create_patient(prj, genelist, pheno_hpo_list)
    print(f'\tPatient ID= {patient_id}')
    print("\tWait for response from AMELIE")
    response = wait_for_processing(patient_id, timeout)

    if response:
        res = response.json()
        # dump the result
        # with open(f'{prj}.amelie.pkl', 'wb') as out:
        #     pickle.dump(res, out)
        return res
    print(f'Fail to get response from AMELIE. timeout={timeout}s')
    return 0

def query_genelist_api(data):
    url = 'https://amelie.stanford.edu/api/gene_list_api/'
    r = requests.post(
        url,
        verify=False,
        data=data)
    try:
        return r.json()
    except:
        return 0


def main(prj, genelist, pheno_hpo_list, pw=None, force=False):
    # force toggle, run the parse and query process even if the file already exist
    print('\tnow running amelie')
    pw = pw or os.getcwd()

    # print(\tgenelist)
    if not isinstance(genelist, list):
        genelist = open(genelist).read().strip().split('\n')
        genelist = [_.strip() for _ in genelist if _.strip()]
    # print(f'\tlen genelist = {len(genelist)}')

    # print(\tpheno_hpo_list)
    if not isinstance(pheno_hpo_list, list):
        pheno_hpo_list = open(pheno_hpo_list).read().strip().split('\n')
        pheno_hpo_list = [_.strip() for _ in pheno_hpo_list if _.strip()]
    phenotypes = ','.join(pheno_hpo_list)

    n_genes = len(genelist)

    # print(f'\tlen HPO list = {len(pheno_hpo_list)}')
    json_pdict = f'{prj}.amelie.pkl'
    json_res = []  # is a list of list. each  element is a gene, like

    # ['CACNA1C',
    # [['26227324', 36.348135596969556],  # pid, score
    # ['28211989', 28.82416813203869],
    # ['27034553', 27.804954354849674],
    # ['22106044', 25.36226509548468],
    # ['33985586', 24.464163141829033]]]
    if not os.path.exists(json_pdict) or force:
        # deal with the situation that gene list larger than 1000
        n_gene_k = int(n_genes/1000)
        genelist_new = []
        for i in range(n_gene_k):
            genelist_new.append(genelist[1000 * i: 1000 * (i+1)])
        if n_genes % 1000 > 0:
            genelist_new.append(genelist[n_gene_k * 1000:])

        for i, igenelist in enumerate(genelist_new):
            data = {'patientName': f'{prj}_{i}',
                    'phenotypes': phenotypes,
                    'genes': ','.join(igenelist),
              }
            tmp = query_genelist_api(data)
            if tmp:
                json_res.extend(tmp)
            else:
                print('ERROR, fail to get AMELIE result')
        with open(json_pdict, 'wb') as o:
            pickle.dump(json_res, file=o)

    else:
        print('amelie.pkl exist: {json_pdict}')
        # json_res = eval(open(json_pdict).read())
        json_res = pickle.load(open(json_pdict, 'rb'))

    gene_score = []
    for i in json_res:
        gn = i[0]
        scores = sorted([_[1] for _ in i[1]], reverse=True)
        gene_score.append(f'{gn}\t{scores[0]:.2f}')

    fn_amelie_lite = f'{pw}/{prj}.amelie.lite.txt'
    print(fn_amelie_lite)
    with open(fn_amelie_lite, 'w') as o:
        print('\n'.join(gene_score), file=o)



def main_old(prj, genelist, pheno_hpo_list, pheno_for_match, pw=None, pw_main=None, force=False):
    # force toggle, run the parse and query process even if the file already exist
    print('\tnow running amelie')
    pw = pw or os.getcwd()
    pw_main = pw_main or pw  # the main files need to save to a separate folder
    os.chdir(pw)
    # print(\tgenelist)
    if not isinstance(genelist, list):
        genelist = open(genelist).read().strip().split('\n')
        genelist = [_.strip() for _ in genelist if _.strip()]
    # print(f'\tlen genelist = {len(genelist)}')

    # print(\tpheno_hpo_list)
    if not isinstance(pheno_hpo_list, list):
        pheno_hpo_list = open(pheno_hpo_list).read().strip().split('\n')
        pheno_hpo_list = [_.strip() for _ in pheno_hpo_list if _.strip()]

    if not isinstance(pheno_for_match, list):
        pheno_for_match = open(pheno_for_match).read().strip().split('\n')
        pheno_for_match = [_.split('#', 1)[0].strip().replace('+', ' ') for _ in pheno_for_match if _.strip()]
        pheno_for_match = [_ for _ in pheno_for_match if _]

    n_genes = len(genelist)

    # print(f'\tlen HPO list = {len(pheno_hpo_list)}')
    json_pdict = f'{prj}.amelie.pkl'
    json_res = []  # is a list of multiple amelie response
    if not os.path.exists(json_pdict):
        # deal with the situation that gene list larger than 1000
        n_gene_k = int(n_genes/1000)
        genelist_new = []
        for i in range(n_gene_k):
            genelist_new.append(genelist[1000 * i: 1000 * (i+1)])
        if n_genes % 1000 > 0:
            genelist_new.append(genelist[n_gene_k * 1000:])

        for igenelist in genelist_new:
            tmp = get_response(prj, igenelist, pheno_hpo_list)
            if tmp:
                json_res.append(tmp)
    else:
        print('amelie.pkl exist: {json_pdict}')
        # json_res = eval(open(json_pdict).read())
        json_res = pickle.load(open(json_pdict, 'rb'))

    # parse the hpo
    f_parsed_pdict = f'{prj}.amelie.parsed.pkl'
    if not os.path.exists(f_parsed_pdict) or force:
        if json_res:
            parsed_pdict = parse_json(prj, json_res, pw, pw_main)
        else:
            print('json result from AMILIE not available, quit')
            return 0
    else:
        parsed_pdict = pickle.load(open(f_parsed_pdict, 'rb'))

    # build the matched phenotpe terms in paper
    query(prj, pheno_for_match, parsed_pdict, pw, pw_main)



if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('prj', help="""task name""")
    ps.add_argument('genelist', help="""genelist, one gene per line""")
    ps.add_argument('file_pheno_hpo_list', help="""phenotype HPO ID list, one HPO ID per line""")
    ps.add_argument('file_pheno_for_match', help="""optional, default is same as pheno file. consice phenotype terms for check the match of phenotype with paper""", nargs='?', default=None)
    ps.add_argument('-pw', help="""default is current folder""", nargs='?', default=None)

    # file_pheno_for_match, the main different from file_pheno_hpo list is that, in this file, the terms don't need to be hpo exact term, it should be keywords, usually a single word per line, e.g. seizure,  intellect,  hypotonia

    args = ps.parse_args()

    genelist = args.genelist
    pheno_hpo_list = args.file_pheno_hpo_list

    file_pheno_for_match = args.file_pheno_for_match
    file_pheno_for_match = file_pheno_for_match or f'{prj}_terms.txt'

    pheno_for_match = [_.strip() for _ in open(args.file_pheno_for_match)]

    prj = args.prj
    pw = args.pw or os.getcwd()
    main(prj, genelist, pheno_hpo_list, pheno_for_match, pw)
