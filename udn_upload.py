#! /usr/bin/env python
"""
run the upload task
upload to emedgene / dropbox / shared drive
"""
import os, sys, re
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('info_file', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""")
ps.add_argument('dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db'])
ps.add_argument('remote_pw', help="""the remote path, start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345""")

ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = all""", nargs='*')
ps.add_argument('-pw', help="""download file output path, default is pwd""")

args = ps.parse_args()

info_file = args.info_file
dest = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene'}[args.dest]
ft_raw = args.ft or 'all'
ft_convert = {'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'vcf': 'vcf'}
remote_pw = args.remote_pw
remote_pw = re.sub('^/', '', remote_pw)
remote_pw = re.sub('/$', '', remote_pw)

pw = args.pw or os.getcwd()

import logging
prefix = 'download_upload_udn'
fn_log = f'{prefix}.log'
fn_err = f'{prefix}.err'

fmt = logging.Formatter('%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console = logging.StreamHandler(sys.stdout)
console.setFormatter(fmt)
console.setLevel('INFO')

fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
fh_err = logging.FileHandler(fn_err, mode='a', encoding='utf8')

fh_file.setLevel('DEBUG')
fh_file.setFormatter(fmt)

fh_err.setLevel('ERROR')
fh_err.setFormatter(fmt)

logger = logging.getLogger(__file__)
logger.setLevel('DEBUG')
logger.addHandler(console)
logger.addHandler(fh_file)
logger.addHandler(fh_err)

m = re.match(r'download\.info\.(UDN\d+)\.txt', info_file)
if m:
    udn = m.group(1)
else:
    udn = input(f'please input the UDN name for the proband: ')
    udn = udn.strip()

if not os.path.exists(info_file):
    logger.error(f'info file not found')
    sys.exit(1)

ft = []
wrong_ft = []
if ft_raw != 'all':
    for i in ft:
        try:
            ft.append(ft_convert[i])
        except:
            wrong_ft.append(i)
else:
    ft = 'all'

if len(wrong_ft) == 0:
    logger.error(f'unkown file type found: {wrong_ft}')


# check the remote folder
d_exist = {}
if dest == 'dropbox':
    res = os.popen(f'dbxcli ls /{remote_pw}').read().strip()
    if not res:
        logger.error(f'remote path for dropbox not exist!  exit...')
        sys.exit(1)

    fn_exist = f'{pw}/{udn}.dropbox.existing_files.txt'
    os.system(f'dbxcli ls -R /{remote_pw} -l > {fn_exist}')
    with open(fn_exist) as fp:
        fp.readline()
        for i in fp:
            i = i.strip()
            a = re.split('/', i, 1)
            fn = a[-1].rsplit('/', 1)[-1]
            d_exist[fn] = a[-1]

elif dest == 'emedgene':
    res = os.popen(f'aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/').read().strip()
    if not res:
        logger.info(f'remote path for emedgene not exist! exit... \naws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/')
        sys.exit(1)

    # 2020-09-21 17:46:59 5209980603 Vanderbilt/upload/GG251_PreUDN/BOB-00001_R1.fastq.gz
    fn_exist = f'{pw}/{udn}.emedgene.existing_files.txt'
    os.system(f'aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/ --recursive > {fn_exist}')
    with open(fn_exist) as fp:
        for i in fp:
            i = i.strip()
            a = re.split('Vanderbilt/upload/', i)
            fn = a[-1].rsplit('/', 1)[-1]
            d_exist[fn] = a[-1]



ft_convert = {'fastq.gz': 'fastq', 'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'bai':'bam', 'vcf.gz': 'vcf', 'vcf':'vcf'}
with open(info_file) as fp:
    for i in fp:
        a = i.split('\t')
        if a[0] == 'rel_to_proband':
            continue
        rel, fn, url = a[:3]

        ext = fn.rsplit('.', 1)[-1]
        ext2 = fn.rsplit('.', 2)[-2]

        ext = f'{ext2}.{ext}' if ext == 'gz' else ext
        i_ft = ft_convert.get(ext)

        if not i_ft:
            logger.warning(f'unkown file type: {fn}')
            continue

        if i_ft not in ft and ft != 'all':
            logger.info(f'skipped because of filetype filter: {fn}')
            continue

        with open(f'{pw}/{fn}.download_upload.{dest}.sh', 'w') as out:
            if not os.path.exists(f'{pw}/{fn}'):
                if fn in d_exist:
                    logger.info(f'File already uploaded! {fn}: {d_exist[fn]}')
                    continue
                print(f'wget "{url}" -c -O "{pw}/{fn}" > {pw}/{fn}.download.log 2>&1', file=out)
                if dest == 'dropbox':
                    print(f'dbxcli put {pw}/{fn} /{remote_pw}/{fn} > {pw}/fn.upload.{dest}.log 2>&1', file=out)
                elif dest == 'emedgene':
                    print(f'aws s3 cp {pw}/{fn} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn} > {pw}/fn.upload.{dest}.log 2>&1', file=out)
