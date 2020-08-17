#! /usr/bin/env python3

#! /usr/bin/env python3

import sys, json
import requests
import os
import time


url = 'https://amelie.stanford.edu/api'

BASEDIR = os.path.dirname(os.path.realpath(__file__))
d_hpo = eval(open(f'{BASEDIR}/hpo_id_to_name.pdict').read())

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

            # 'genotypesArea': '\n'.join(genelist),
            # 'phenotypesArea': '\n'.join(hpo)
            })

    return r.json()['patient_id']


def wait_for_processing(patient_id, timeout=120):
    for _ in range(timeout):
        response = requests.get('%s/patient/%s' % (url, patient_id),)
                                # verify=False)
        if response.json()['processing'] != 1:
            return response
        time.sleep(1)
    return 0

def parse_json_old(prj, json_pdict):
    """
    extract the information from json
    r = json_pdict,   got from response.json(),  or read from the pdict file
    """
    r = json_pdict
    res = []
    # hpo_q = r['system_list'][0]['phenotypes']
    hpo_q_raw = [_['phenotypes'] for _ in r['system_list']]
    hpo_q = []
    for _ in hpo_q_raw:
        hpo_q.extend(_)

    hpo_q = [_['hpo_id']+ "\t" + _['name'] for _ in hpo_q]

    # res['query_hpo'] = hpo_q

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

    with open(f'{prj}.amelie_final_result.full.txt', 'w') as out:
        print('The Query phenotype is :\n\t', end='', file=out)
        print('\n\t'.join(hpo_q), file=out)

        print('The Genes are as follow', file=out)
        for ign in res:
            print(f'gene={ign[1]}\tscore={ign[0]:.2f}', file=out)
            print(f'OMIM: {ign[2]}', file=out)
            print(f'Papers', file=out)
            for ipaper in ign[3:]:
                print(f'\t{ipaper["title"]}', file=out)
                print(f'\tYear={ipaper["pubmed_year"]}\tPMID={ipaper["pmid"]}', file=out)
                print(f'\t\t', end='', file=out)
                print(f'\n\t\t'.join(ipaper["hpo"]), file=out)
                print(f'\t', end='', file=out)
                print('*' * 10 + '\n', file=out)

            print('#' * 30 + '\n\n', file=out)

    with open(f'{prj}.amelie.lite.txt', 'w') as out:
        for ign in res:
            print(f'{ign[1]}\t{ign[0]:.2f}', file=out)

    with open(f'{prj}.amelie.parsed.pdict', 'w') as out:
        print(res, file=out)

    # return res


def parse_json(prj, json_pdict):
    """
    extract the information from json
    r = json_pdict,   got from response.json(),  or read from the pdict file
    """
    r = json_pdict
    res = []
    # hpo_q = r['system_list'][0]['phenotypes']
    hpo_q_raw = [_['phenotypes'] for _ in r['system_list']]
    hpo_q = []
    for _ in hpo_q_raw:
        hpo_q.extend(_)

    hpo_q = [_['hpo_id']+ "\t" + _['name'] for _ in hpo_q]

    # res['query_hpo'] = hpo_q

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

    # with open(f'{prj}.amelie_final_result.full.txt', 'w') as out:
    #     print('The Query phenotype is :\n\t', end='', file=out)
    #     print('\n\t'.join(hpo_q), file=out)

    #     print('The Genes are as follow', file=out)
    #     for ign in res:
    #         # print('\t'.join([ign[1], '%.2f' % ign[0], ign[2]]))
    #         print(f'gene={ign[1]}\tscore={ign[0]:.2f}', file=out)
    #         print(f'OMIM: {ign[2]}', file=out)
    #         print(f'Papers', file=out)
    #         for ipaper in ign[3:]:
    #             print(f'\t{ipaper["title"]}', file=out)
    #             print(f'\tYear={ipaper["pubmed_year"]}\tPMID={ipaper["pmid"]}', file=out)
    #             print(f'\t\t', end='', file=out)
    #             print(f'\n\t\t'.join(ipaper["hpo"]), file=out)
    #             print(f'\t', end='', file=out)
    #             print('*' * 10 + '\n', file=out)

    #         print('#' * 30 + '\n\n', file=out)

    with open(f'{prj}.amelie.lite.txt', 'w') as out:
        for ign in res:
            print(f'{ign[1]}\t{ign[0]:.2f}', file=out)


    with open(f'{prj}.amelie.parsed.pdict', 'w') as out:
        print(res, file=out)

    return res

    # return res
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

def query(prj, pheno, pdict):
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
            pmid = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
            hpo = re_group_hpo(ipaper['hpo'])  # a list

            for ihpo, hpo_des in hpo.items():
                for ipheno in pheno:
                    for _ in hpo_des:
                        _ = _.lower()
                        if _.find(ipheno) > -1:
                            try:
                                hit[ipheno].append([pmid, ipaper['title'], ipaper['journal'], ipaper['pubmed_year'], hpo_des])
                            except:
                                hit[ipheno] = [[pmid, ipaper['title'], ipaper['journal'], ipaper['pubmed_year'], hpo_des]]
                            break
        if len(hit) > 0:
            res[ign] = [score, omim, hit]

    # struc =
    # {gene1: 0=score, 1=omim, 2 = hits(below)
    #         pheno1: [
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #         ]
    #     }

    with open(f'{prj}.amelie.matched_query.pdict', 'w') as out:
        print(res, file=out)


    # export into text file
    with open(f'{prj}.amelie.matched_query.txt', 'w') as out:
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
        with open(f'{prj}.amelie.pdict', 'w') as out:
            print(res, file=out)

        return res
    print(f'Fail to get response from AMELIE. timeout={timeout}s')
    return 0


def main(prj, genelist, pheno_hpo_list, pheno_for_match):
    print('\tnow running amelie')
    # print(\tgenelist)
    if not isinstance(genelist, list):
        genelist = open(genelist).read().strip().split('\n')
        genelist = [_.strip() for _ in genelist if _.strip()]
    # print(f'\tlen genelist = {len(genelist)}')

    # print(\tpheno_hpo_list)
    if not isinstance(pheno_hpo_list, list):
        pheno_hpo_list = open(pheno_hpo_list).read().strip().split('\n')
        pheno_hpo_list = [_.strip() for _ in pheno_hpo_list if _.strip()]

    # print(f'\tlen HPO list = {len(pheno_hpo_list)}')
    json_pdict = f'{prj}.amelie.pdict'
    if not os.path.exists(json_pdict):
        json_res = get_response(prj, genelist, pheno_hpo_list)
    else:
        json_res = eval(open(json_pdict).read())

    # parse the hpo
    f_parsed_pdict = f'{prj}.amelie.parsed.pdict'
    if not os.path.exists(f_parsed_pdict):
        if json_res:
            parsed_pdict = parse_json(prj, json_res)
        else:
            print('json result from AMILIE not available, quit')
            sys.exit(1)
    else:
        parsed_pdict = eval(open(f_parsed_pdict).read())

    # build the matched phenotpe terms in paper
    query(prj, pheno_for_match, parsed_pdict)


if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('prj', help="""task name""")
    ps.add_argument('genelist', help="""genelist, one gene per line""")
    ps.add_argument('file_pheno_hpo_list', help="""phenotype HPO ID list, one HPO ID per line""")
    ps.add_argument('file_pheno_for_match', help="""optional, default is same as pheno file. consice phenotype terms for check the match of phenotype with paper""", nargs='?', default=None)

    # file_pheno_for_match, the main different from file_pheno_hpo list is that, in this file, the terms don't need to be hpo exact term, it should be keywords, usually a single word per line, e.g. seizure,  intellect,  hypotonia

    args = ps.parse_args()

    genelist = args.genelist
    pheno_hpo_list = args.file_pheno_hpo_list

    file_pheno_for_match = args.file_pheno_for_match
    file_pheno_for_match = file_pheno_for_match or f'{prj}_terms.txt'

    pheno_for_match = [_.strip() for _ in open(args.file_pheno_for_match)]

    prj = args.prj
    main(prj, genelist, pheno_hpo_list, pheno_for_match)
