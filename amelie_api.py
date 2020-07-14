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
        verify=False,
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


def wait_for_processing(patient_id):
    while True:
        response = requests.get('%s/patient/%s' % (url, patient_id),
                                verify=False)
        if response.json()['processing'] != 1:
            break
        time.sleep(1)


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


def get_response(prj, genelist, pheno_hpo_list):
    patient_id = create_patient(prj, genelist, pheno_hpo_list)
    print("Wait while vcfs are filtered and annotated and papers are ranked")
    wait_for_processing(patient_id)
    response = requests.get('%s/patient/%s' % (url, patient_id),
                            verify=False)
    res = response.json()
    # dump the result
    with open(f'{prj}.amelie.pdict', 'w') as out:
        print(res, file=out)

    return res


def main(prj, genelist, pheno_hpo_list):

    json_pdict = f'{prj}.amelie.pdict'
    if not os.path.exists(json_pdict):
        json_res = get_response(prj, genelist, pheno_hpo_list)
    else:
        json_res = eval(open(json_pdict).read())

    # parse the hpo
    # json_res = response.json()
    parse_json(prj, json_res)


if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('prj', help="""task name""")
    ps.add_argument('genelist', help="""genelist, one gene per line""")
    ps.add_argument('pheno_hpo_list', help="""phenotype HPO ID list, one HPO ID per line""")
    args = ps.parse_args()

    genelist = args.genelist
    genelist = open(genelist).read().strip().split('\n')
    pheno_hpo_list = args.pheno_hpo_list
    pheno_hpo_list = open(pheno_hpo_list).read().strip().split('\n')


    prj = args.prj
    main(prj, genelist, pheno_hpo_list)
