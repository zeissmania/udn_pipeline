#! /usr/bin/env python3

"""
input should be the project path
the folder should contains
0. the settings file (yaml), the config file
1. the VCF files
2. the phenotype/symptom txt file
3. query_terms, in this file, include the terms in short, to filter on the AMELIE result, each term per line


output:
anno_exp = f'{pw}/{sample_id}_{v[1]}.annotated.txt'
anno_filter = f'{pw}/{sample_id}_{v[1]}.filtered.txt
proband_genelist = f'{pw}/{proband}_{v[1]}.filtered.genelist
anno_extract = f'{pw}/{sample_id}_{v[1]}.extracted.txt
anno_tally = f'{pw}/{lb}.for_tally.pdict'

hpo = f'{proband}_terms_hpo.txt
hpo_pure = f'{proband}_terms_pure_hpo.txt
amelie_lite = f'{proband}.amelie.lite.txt
amelie_final = f'{proband}.amelie.matched_terms.txt

result_raw = f'{proband}.final_result.all.sorted.xlsx
result_selected = f'{proband}.final_result.selected.xlsx
report = f'{proband}.report.xlsx'

"""
import sys
import os
import glob
pw_code = os.path.dirname(os.path.realpath(f'{__file__}/..'))
sys.path.append(pw_code)
# print(pw_code)
from udn import udn_utils
# from udn import amelie_api
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('config', help="""the config file, yaml format""", nargs='?', default=None)
ps.add_argument('-omim', help="""force rerun the omim match report, no matter report exist or not""", action='store_true')
args = ps.parse_args()

config_file = args.config

if config_file is None:
    config_file = glob.glob('*.yaml')
    if len(config_file) == 0:
        print('ERROR: no yaml config file is found, exit')
        sys.exit(1)
    config_file = config_file[0]

case = udn_utils.UDN_case(config_file)
sv_caller = case.sv_caller
logger = case.logger
if case.done_phase3:
    # prj.report.xlsx
    case.omim_query(force_rerun=args.omim)
    logger.info('Report is ready')
elif case.done_phase2:
    # prj.selected.genes.txt
    case.omim_query(force_rerun=args.omim)
    logger.info('building the report...')
    import udn_report
    udn_report.main(case.prj, case.pw, sv_caller=sv_caller)

elif case.done_phase1:
    # prj.merged.sortedtsv
    # query the amelie API, build the amelie result files
    case.omim_query(force_rerun=args.omim)
    logger.info(f'FINAL STEP:  select the gene list manually, \n\t expected file = {case.pw}/{case.prj}.selected.genes.txt')
else:
    # filter /extract annotation
    # query the amelie API, build the amelie result files
    case.query_amelie(force=False)
    # get amelie dict
    case.amelie_dict = case.get_amelie_dict()

    if sv_caller == 'dragen':
        for sample_id in case.family:
            case.annotate(sample_id)
            case.anno_filter(sample_id)
            case.anno_extract(sample_id)
            case.family[sample_id]['sv_dict'] = case.group_sv_into_bin(sample_id)
            # match SV in proband with family
            case.run_proband_sv_match()
        # print(case.family[sample_id]['sv_dict'].keys())
    else:
        proband_id = case.proband_id
        case.annotate(proband_id)
        case.anno_filter(proband_id)
        case.anno_extract(proband_id)
        case.run_proband_sv_match_pacbio()

    # run OMIM match
    case.omim_query()

    logger.info('Phase 1 done, please select the top 10 gene list')
