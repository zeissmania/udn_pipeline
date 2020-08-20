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
from udn import udn_utils
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('config', help="""the config file, yaml format""")
args = ps.parse_args()

config_file = args.config

case = udn_utils.UDN_case(config_file)

# filter /extract annotation
for sample_id in case.family:
    case.annotate(sample_id)
    case.anno_filter(sample_id)
    case.anno_extract(sample_id)
    case.family[sample_id]['sv_dict'] = case.group_sv_into_bin(sample_id)

# query the amelie API, build the amelie result files
case.query_amelie()

# get amelie dict
case.amelie_dict = case.get_amelie_dict()

# match SV in proband with family
case.run_proband_sv_match()
