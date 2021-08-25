#! /usr/bin/env python3
"""
generate the IGV batch script
https://software.broadinstitute.org/software/igv/PortCommands
https://github.com/igvteam/igv/blob/master/src/main/resources/org/broad/igv/prefs/preferences.tab
input 1 = igv download file
input 2 = {udn}.selected.genes.txt. first column = mut_type

if input2 is not defined, would create the script just adding the bam tracks
"""

import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('udn', help="""udn ID, if not sepecified, would search the UDNxxxx.igv.files.txt""", nargs='?')
ps.add_argument('-pw', help="""path for the UDN project, default = pwd""", default=None)
args = ps.parse_args()

import os, sys
import re
pw = args.pw or os.getcwd()
udn = args.udn
if not udn:
    import glob
    tmp = glob.glob('UDN*.igv.files.txt')
    if len(tmp) == 0:
        print('Error, UDN not specified, UDNxxx.igv.files.txt not found, exit')
        sys.exit(1)
    elif len(tmp) > 1:
        print(f'multiple igv files found, would use first {tmp}')

    udn = tmp[0].split('.', 1)[0]

import logging


fmt = logging.Formatter('%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console = logging.StreamHandler(sys.stdout)
console.setFormatter(fmt)
console.setLevel('INFO')

logger = logging.getLogger(__file__)
logger.setLevel('DEBUG')
logger.addHandler(console)


pwigv = f'{pw}/igv'
os.system(f'mkdir -p {pwigv} 2>/dev/null')

# read-in the bam files
fn_in = f'{pw}/{udn}.igv.files.txt'
if not os.path.exists(fn_in):
    logger.error(f'{udn}: igv.files.txt file not found')
    sys.exit(1)

d_bam = {}
with open(fn_in) as fp:
    for i in fp:
        i = i.strip()
        if not i:
            continue
        url, fn = i.split('\t')
        ext = fn.rsplit('.', 1)[-1]
        lb = fn.replace('.bai', '')
        try:
            d_bam[ext][lb] = url
        except:
            d_bam[ext] = {lb: url}

def get_gene_region():
    # read-in the region/ genes
    for fn_in in [f'{pw}/{udn}.selected.genes.txt', f'{pw}/{udn}.selected.genes.xlsx', f'{pw}/selected.genes.txt', f'{pw}/selected.genes.xlsx']:
        if os.path.exists(fn_in):
            break
    else:
        logger.warning(f'selected.genes.txt/xlsx file not found, would only add bam tracks, no jumping and screenshot')
        return {}

    d_genes = {}

    ext = fn_in.rsplit('.', 1)[-1]
    if ext == 'txt':
        with open(fn_in) as fp:
            data = [i.strip() for i in fp if i.strip()]
            data = [i.split('\t') for i in data]
    elif ext == 'xlsx':
        import pandas as pd
        data = pd.read_excel(fn_in, keep_default_na=False)
        data = data.itertuples(False, None)


    for a in data:
        if len(a) < 9:
            print(f'error split: {a}')
            continue
        gn = a[9]
        try:
            s = a[4] + 0
            e = a[5] + 0
        except:
            try:
                s = int(re.sub(r"[,\s\"']+", '', a[4]))
                e = int(re.sub(r"[,\s\"']+", '', a[5]))
            except:

                print(f'fail to convert start/end to int: {a}')
                continue
        chr_ = a[3]
        try:
            chr_ = {'chrx': 'chrX', 'chry': 'chrY'}[chr_]
        except:
            pass
            # continue

        len_sv = e - s
        try:
            _ = d_genes[gn]
        except:
            d_genes[gn] = []   # do not add gn as the first element, because they are too big, no valid screenshot

        if len_sv > 20000:
            n = int(len_sv/10000)
            for _ in range(n-1):
                d_genes[gn].append(f'{chr_}:{s-2000 + _*10000}-{s + 2000 + (_ + 1)*10000}')
            d_genes[gn].append(f'{chr_}:{s-1000 + (n-1) *10000}-{e + 2000}')
        else:
            d_genes[gn].append(f'{chr_}:{s-2000}-{e+2000}')

    return d_genes

d_genes = get_gene_region()
# output igv batch script
with open(f'{pw}/igv.script.{udn}.txt', 'w') as out:
    out.write(f"""genome hg38
SAM.QUALITY_THRESHOLD 13
snapshotDirectory {pwigv}
maxPanelHeight 2000
squish
setSleepInterval 500
new

""")
    for lb, url_bam in d_bam['bam'].items():
        try:
            url_bai = d_bam['bai'][lb]
        except:
            logger.warning(f'bai file not found: {lb}')
            continue
        out.write(f'load "{url_bam}" index="{url_bai}"\n')

    out.write("""
load http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
collapse refGene.txt.gz

""")
    for gn, v in d_genes.items():
        for n, region in enumerate(v):
            # suffix = 'gene_overview' if n == 0 else n
            suffix = n + 1
            region_short = '.' + region.split(':')[-1]
            out.write(f'goto {region}\nsnapshot "{gn}.{suffix}{region_short}.png"\n\n')
