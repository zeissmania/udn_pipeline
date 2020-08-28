#! /usr/bin/env python3

#! /usr/bin/env python3

import sys
import requests
import os
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


def parse_json(prj, json_pdict, pw=None, pw_main=None, force=0):
    """
    extract the information from json
    r = json_pdict,   got from response.json(),  or read from the pdict file
    """
    pw = pw or os.getcwd()
    pw_main = pw_main or pw
    r = json_pdict

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

    with open(f'{pw_main}/{prj}.amelie.lite.txt', 'w') as out:
        for ign in res:
            print(f'{ign[1]}\t{ign[0]:.2f}', file=out)

    # export the tsv file

    with open(f'{pw}/{prj}.amelie.parsed.raw.txt', 'w') as out:
        for i in res:
            score, gn, omim, *papers = i
            print(f'{gn}\t{score}', file=out)
            for ipaper in papers:
                title = ipaper['title']
                pmid = ipaper['pmid']
                hpo = ipaper['hpo']
                print(f'\t{pmid}\t{title}', file=out)
                print(f'\t{hpo}', file=out)
                print(f'', file=out)


    with open(f'{pw}/{prj}.amelie.parsed.pkl', 'wb') as out:
        pickle.dump(res, out)

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
            pmid = f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'
            hpo = re_group_hpo(ipaper['hpo'])  # a list

            for _, hpo_des in hpo.items():
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
        elif score > 10:
            print(f'\tWARNING: amelie query no match with the input HPO terms: gene={ign}, score={score}')


    # struc =
    # {gene1: 0=score, 1=omim, 2 = hits(below)
    #         pheno1: [
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #                 (hpo_match1) [pmid, title, journal, year, hpo_detail_list],
    #         ]
    #     }

    with open(f'{pw}/{prj}.amelie.matched_query.pdict', 'w') as out:
        print(res, file=out)


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
        with open(f'{prj}.amelie.pkl', 'wb') as out:
            pickle.dump(res, out)
        return res
    print(f'Fail to get response from AMELIE. timeout={timeout}s')
    return 0



def main(prj, genelist, pheno_hpo_list, pheno_for_match, pw=None, pw_main=None, force=False):
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
        pheno_for_match = [_.strip() for _ in pheno_for_match if _.strip()]


    # print(f'\tlen HPO list = {len(pheno_hpo_list)}')
    json_pdict = f'{prj}.amelie.pkl'
    if not os.path.exists(json_pdict):
        json_res = get_response(prj, genelist, pheno_hpo_list)
    else:
        # json_res = eval(open(json_pdict).read())
        json_res = pickle.load(open(json_pdict, 'rb'))

    # parse the hpo
    f_parsed_pdict = f'{prj}.amelie.parsed.pkl'
    if not os.path.exists(f_parsed_pdict) or force:
        if json_res:
            parsed_pdict = parse_json(prj, json_res, pw, pw_main)
        else:
            print('json result from AMILIE not available, quit')
            sys.exit(1)
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
