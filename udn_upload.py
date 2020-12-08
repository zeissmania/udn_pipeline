#! /usr/bin/env python
"""
run the upload task
upload to emedgene / dropbox / shared drive
"""
import os, sys, re, glob
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db'])
ps.add_argument('remote_pw', help="""the remote path, if not exist, would create one. start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345""", default='', nargs='?')
ps.add_argument('info_file', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""", nargs='?')
ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = all""", nargs='*')
ps.add_argument('-pw', help="""download file output path, default is pwd""")

args = ps.parse_args()
print(args)

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


info_file = args.info_file
if not info_file:
    tmp = glob.glob('download.info.*.txt')
    if len(tmp) == 0:
        logger.error(f'no download info file found! exit...')
        sys.exit(1)
    elif len(tmp) > 1:
        logger.error(f'multiple info file found! \n{tmp}\nexit...')
        sys.exit(1)
    else:
        info_file = tmp[0]

dest = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene'}[args.dest]
ft_raw = args.ft or 'all'
ft_convert = {'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'vcf': 'vcf'}
remote_pw = args.remote_pw
remote_pw = re.sub('^/', '', remote_pw)
remote_pw = re.sub('/$', '', remote_pw)

pw = args.pw or os.getcwd()



import platform
node_name = platform.node()

if node_name.find('viccbiostat120') > -1:
    dock = 'singularity exec /mnt/d/dock/centos.sif '
elif node_name.find('vampire') > -1:
    dock = 'singularity exec /data/cqs/chenh19/dock/centos.sif '


m = re.match(r'download\.info\.(UDN\d+)\.txt', info_file)
if m:
    udn = m.group(1)
else:
    # udn = input(f'please input the UDN name for the proband: ')
    logger.error('wrong info file format, UDN id not found')
    sys.exit(1)
    # udn = udn.strip()

if not os.path.exists(info_file):
    logger.error(f'info file not found')
    sys.exit(1)

if not remote_pw:
    remote_pw = pw.rsplit('/', 1)[-1] if pw.find(udn) > -1 else udn
    logger.warning(f'remote path not specified, would create folder : {remote_pw}')

# create subfolder even the remote folder root is specified
if os.path.exists(f'{dest}.create_folder.log'):
    logger.info(f'subfolder already exist')
else:
    with open(info_file) as fp:
        fp.readline()
        sub_folders = set()
        for i in fp:
            a = i.split('\t')
            rel, iudn = re.sub(r'\W+', '_', a[0]), a[3]
            name = f'{rel}_{iudn}'
            if name not in sub_folders:
                sub_folders.add(name)
                logger.info(f'\tcreating subfolder: {remote_pw}/{name}')
                if dest == 'dropbox':
                    os.system(f'{dock} dbxcli mkdir {remote_pw}/{name}/ >>{dest}.create_folder.log 2>>{dest}.create_folder.log')
                elif dest == 'emedgene':
                    os.system(f'{dock} aws s3api put-object --bucket emg-auto-samples --key Vanderbilt/upload/{remote_pw}/{name}/ >>{dest}.create_folder.log 2>>{dest}.create_folder.log')

ft = []
wrong_ft = []
if ft_raw != 'all':
    for i in ft_raw:
        try:
            ft.append(ft_convert[i])
        except:
            wrong_ft.append(i)
else:
    ft = 'all'

if len(wrong_ft) > 0:
    logger.error(f'unkown file type found: {wrong_ft}')

# get the file expected size
d_file_size_exp = {}
with open(info_file) as fp:
    for i in fp:
        a = i.split('\t')
        rel, fn, url, fl_udn = a[:4]
        if a[0] == 'rel_to_proband':
            continue
        try:
            size_exp = int(a[5])
        except:
            logger.error(f'wrong filesize format: {a[5]}, line={i}')
            sys.exit(1)
        else:
            d_file_size_exp[fn] = size_exp


# check the remote folder
subfolders = []
d_exist = {}
if dest == 'dropbox':
    res = os.popen(f'{dock} dbxcli ls /{remote_pw}').read().strip()
    if not res:
        logger.error(f'remote path for dropbox not exist!  exit...')
        sys.exit(1)

    fn_exist = f'{pw}/{udn}.dropbox.existing_files.txt'
    os.system(f'{dock} dbxcli ls -R /{remote_pw} -l > {fn_exist}')
    with open(fn_exist) as fp:
        fp.readline()
        for i in fp:
            i = i.strip()
            a = re.split('/', i, 1)
            fn = a[-1].rsplit('/', 1)[-1]
            d_exist[fn] = a[-1]

            # extract the sub folders
            if a[0][0] == '-':
                sub_full_path = a[-1]
                sub_full_path = sub_full_path[1:] if sub_full_path[0] == '/' else sub_full_path
                sub_path = re.sub(f'^{remote_pw}', '', sub_full_path)
                if not sub_path:
                    logger.info(f'root remote path found: {sub_full_path}')
                    continue
                sub_path_end = sub_path.rsplit('/', 1)[-1]
                m = re.match(r'.*(UDN\d+)', sub_path_end)
                sub_path_udn = m.group(1) if m else ''

                subfolders.append([sub_path_udn, sub_path])
elif dest == 'emedgene':
    res = os.popen(f'{dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/').read().strip()
    if not res:
        logger.info(f'remote path for emedgene not exist! exit... \naws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/')
        sys.exit(1)

    # 2020-09-21 17:46:59 5209980603 Vanderbilt/upload/GG251_PreUDN/BOB-00001_R1.fastq.gz
    fn_exist = f'{pw}/{udn}.emedgene.existing_files.txt'
    os.system(f'{dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/ --recursive > {fn_exist}')
    with open(fn_exist) as fp:
        for i in fp:
            i = i.strip()
            a = re.split('Vanderbilt/upload/', i)
            fn = a[-1].rsplit('/', 1)[-1]

            file_size = re.split(r'\s+', i.strip())[2]
            try:
                file_size = int(file_size)
            except:
                logger.error(f'wrong remote emedgene file size: {file_size}')
                sys.exit(1)

            try:
                d_file_size_exp[fn]
            except:
                if fn and fn[-4:] != '.md5':
                    logger.warning(f'remote file not found locally line={i}, fn={fn}, file_size={file_size}')
            else:
                if file_size == d_file_size_exp[fn]:
                    d_exist[fn] = a[-1]

            # extract the sub folders
            if a[-1][-1] == '/':
                sub_full_path = a[-1]
                sub_full_path = sub_full_path[:-1] if sub_full_path[-1] == '/' else sub_full_path
                m = re.match(fr'^{remote_pw}(.*$)', sub_full_path)
                if m:
                    sub_path = m.group(1)
                    if not sub_path:
                        logger.info(f'root remote path found: {sub_full_path}')
                        continue
                    sub_path_end = sub_path.rsplit('/', 1)[-1]
                    m = re.match(r'.*(UDN\d+)', sub_path_end)
                    sub_path_udn = m.group(1) if m else None
                    if sub_path_udn:
                        subfolders.append([sub_path_udn, sub_path])
                    else:
                        logger.info(f'subfolder_without_UDN_ID: {sub_full_path}')
                else:
                    logger.warning(f'path not found under remote pw: {sub_full_path}')

logger.info(f'\n\tremote_dest={dest}\n\tremote_pw={remote_pw}\n\tfile_type_for_upload={ft}\n\tsubfolder={subfolders}')

os.system(f'mkdir {pw}/shell 2>/dev/null')
os.system(f'mkdir {pw}/shell_done 2>/dev/null')
os.system(f'mkdir {pw}/download 2>/dev/null')
os.system(f'mkdir {pw}/log 2>/dev/null')

ft_convert = {'fastq.gz': 'fastq', 'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'bai':'bam', 'vcf.gz': 'vcf', 'vcf':'vcf', 'gvcf.gz': 'vcf'}

# upload the md5 file
fn_md5 = f'download.{udn}.md5'
if os.path.exists(fn_md5) and fn_md5 not in d_exist:
    logger.info('md5sum file upload script building... ')
    with open(f'{pw}/shell/md5sum.upload.{dest}.sh', 'w') as out:
        if dest == 'dropbox':
            print(f'{dock} dbxcli put {pw}/{fn_md5} /{remote_pw}/{udn}_all_files.md5 > {pw}/log/upload.{dest}.{fn_md5}.log 2>&1', file=out)
        elif dest == 'emedgene':
            print(f'{dock} aws s3 cp {pw}/{fn_md5} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{udn}_all_files.md5 > {pw}/log/upload.{dest}.{fn}.log 2>&1', file=out)

with open(info_file) as fp:
    n_all = 0
    n_uploaded = 0
    n_downloaded = 0
    n_desired = 0
    for i in fp:
        a = i.split('\t')
        if a[0] == 'rel_to_proband':
            continue
        rel, fn, url, fl_udn = a[:4]
        size_exp = int(a[5])
        n_all += 1

        folder_extended = [_ for _ in subfolders if fl_udn == _[0]]
        if len(folder_extended) == 0:
            folder_extended = ''
        else:
            folder_extended = folder_extended[0][1]


        ext = fn.replace('.tbi', '').rsplit('.', 1)[-1]
        ext2 = fn.replace('.tbi', '').rsplit('.', 2)[-2]

        ext = f'{ext2}.{ext}' if ext == 'gz' else ext
        i_ft = ft_convert.get(ext)

        if not i_ft:
            logger.warning(f'unkown file type: {fn}')
            continue

        if i_ft not in ft and ft != 'all':
            logger.debug(f'skipped because of filetype filter: {fn}, i_ft={i_ft}, ft={ft}')
            continue

        n_desired += 1

        if not os.path.exists(f'{pw}/{fn}'):
            if fn in d_exist:
                n_uploaded += 1
                n_downloaded += 1
                logger.info(f'File already uploaded! {fn}: {d_exist[fn]}')
                os.system(f'mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh 2>/dev/null')
                continue

        fn_download = f'{pw}/download/{fn}'
        with open(f'{pw}/shell/{fn}.download_upload.{dest}.sh', 'w') as out:
            skip_download = 0
            if os.path.exists(fn_download):
                size = os.path.getsize(fn_download)
                if size == size_exp:
                    logger.info(f'file already downloaded, size match with expected {size_exp}: {fn}')
                    n_downloaded += 1
                    skip_download = 1
                else:
                    logger.warning(f'file already downloaded, but size not match: exp={size_exp} real={size} : {fn}')
                    # n_downloaded += 1
                    os.system(f'mv {fn_download} {fn_download}.bad 2>/dev/null')
            if not skip_download:
                print(f'wget "{url}" -c -O "{fn_download}" > {pw}/log/download.{fn}.log 2>&1', file=out)

            print(f"""local_size=$(stat {fn_download} -c %s 2>/dev/null)\nif [ "$local_size" -ne {size_exp} ];then\necho file size not match: {fn_download}: size=$local_size, expected={size_exp}\nmv {fn_download} {fn_download}.bad 2>/dev/null;\nexit;\nfi""", file=out)

            if dest == 'dropbox':
                print(f'{dock} dbxcli put {pw}/download/{fn} /{remote_pw}{folder_extended}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1', file=out)
            elif dest == 'emedgene':
                print(f'{dock} aws s3 cp {pw}/download/{fn} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}{folder_extended}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1', file=out)

    logger.info(f'total files = {n_all}, need_to_upload={n_desired}, already exist in server = {n_uploaded}, already downloaded = {n_downloaded}')
