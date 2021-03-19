#! /usr/bin/env python
"""
run the upload task
upload to emedgene / dropbox / shared drive
"""
import os, sys, re, glob
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db'], nargs='?', default='emedgene')
ps.add_argument('info_file', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""", nargs='?', default=None)
ps.add_argument('-remote', '-r', '-remote_pw', dest='remote_pw', help="""the remote path, if not exist, would create one. start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345""", default=None, nargs='?')
ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = fastq""", nargs='*')
ps.add_argument('-pw', help="""download file output path, default is pwd""")
ps.add_argument('-demo', help="""demo mode would not actually download or upload, or create remote folder""", action='store_true')
ps.add_argument('-rm', '-delete', help="""flag, if set, would delete the downloaded file from local disk when uploading is done""", action='store_true')
ps.add_argument('-asis', help="""use the remote_pw in the cmd line input, do not do the re-format, default is False""", action='store_true')

args = ps.parse_args()

demo = args.demo
asis = args.asis  # use the remote_pw in the cmd line input, do not do the re-format
rm = args.rm
pw = args.pw or os.getcwd()
# os.chdir(pw)
print(args)
print(f'current pw={pw}')
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
    tmp = glob.glob(f'{pw}/download.info.*.txt')
    if len(tmp) == 0:
        logger.error(f'no download info file found! exit...')
        sys.exit(1)
    elif len(tmp) > 1:
        logger.error(f'multiple info file found! \n{tmp}\nexit...')
        sys.exit(1)
    else:
        info_file = tmp[0]

dest = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene'}[args.dest]
ft_convert = {'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'vcf': 'vcf', 'bai': 'bam'}

ft = []
wrong_ft = []

if not args.ft:
    ft = ['fastq']
elif 'all' in args.ft:
    ft = 'all'
else:
    for i in args.ft:
        try:
            ft.append(ft_convert[i])
        except:
            wrong_ft.append(i)

if len(wrong_ft) > 0:
    logger.warning(f'unkown file type found: {wrong_ft}')

if len(ft) == 0:
    logger.error(f'No file type selected! exit')
    # sys.exit(1)

m = re.match(r'.*download\.info\.(UDN\d+)\.txt', info_file)
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


# remote path
pw = re.sub('/$', '', pw)
remote_pw = args.remote_pw

if not remote_pw:
    remote_pw = pw.rsplit('/', 1)[-1]
remote_pw = re.sub('^/', '', remote_pw)


if not asis:
    remote_pw = re.sub('upload_?', '', remote_pw, flags=re.I)
    remote_pw = re.sub(r'^\d+_', '', remote_pw, 1)
    m = re.match(r'(.*_)?(UDN\d+)(_.*)?', remote_pw)
    if not m:
        logger.warning(f'****** wrong remote path name: {remote_pw}, would use the UDN number: {udn}')
        remote_pw = udn
    else:
        m = m.groups()
        if m[0] is not None:
            p1 = re.sub('^_*', '', m[0])
            p1 = re.sub('_*$', '', p1)
        else:
            p1 = None

        if m[2] is not None:
            p2 = re.sub('^_*', '', m[2])
            p2 = re.sub('_*$', '', p2)
        else:
            p2 = None

        p_trailing = '_'.join([_ for _ in (p1, p2) if _])
        p_trailing = '_' + p_trailing if p_trailing else ''
        remote_pw = m[1] + p_trailing

# print(f'pw={pw}\nremote_pw={remote_pw}\ninfo={info_file}\nft={ft}\ndest={dest}')

# sys.exit(1)



import platform
node_name = platform.node()

if node_name.find('viccbiostat120') > -1:
    dock = 'singularity exec /mnt/d/dock/centos.sif '
elif node_name.find('vampire') > -1:
    dock = 'singularity exec /data/cqs/chenh19/dock/centos.sif '




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
                if dest == 'dropbox' and not demo:
                    os.system(f'{dock} dbxcli mkdir {remote_pw}/{name}/ >>{dest}.create_folder.log 2>>{dest}.create_folder.log')
                elif dest == 'emedgene' and not demo:
                    os.system(f'{dock} aws s3api put-object --bucket emg-auto-samples --key Vanderbilt/upload/{remote_pw}/{name}/ >>{dest}.create_folder.log 2>>{dest}.create_folder.log')


# get the file expected size
d_file_size_exp = {}
idx_size = 5
with open(info_file) as fp:

    for i in fp:
        a = i.split('\t')
        rel, fn, url, fl_udn = a[:4]
        if a[0].strip() == 'rel_to_proband':
            a = [_.strip() for _ in a]
            try:
                idx_size = a.index('size')
            except:
                logger.error('file size not found in the info file header! exit')
                sys.exit(1)
            continue
        try:
            size_exp = int(a[idx_size])
        except:
            line_repr = "\n".join(map(str, enumerate(a)))
            logger.error(f'wrong filesize format: {a[idx_size]}, line=\n{line_repr}')
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
                    logger.warning(f'****** file only found on remote server: line={i}, fn={fn}, file_size={file_size}')
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
                    logger.warning(f'****** path not found under remote pw: {sub_full_path}')

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
            print(f'{dock} aws s3 cp {pw}/{fn_md5} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{udn}_all_files.md5 > {pw}/log/upload.{dest}.{fn_md5}.log 2>&1', file=out)

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

        size_exp = int(a[idx_size])
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
            logger.warning(f'****** unkown file type: {fn}')
            continue

        if i_ft not in ft and ft != 'all':
            logger.debug(f'skipped because of filetype filter: {fn}, i_ft={i_ft}, ft={ft}')
            continue

        n_desired += 1

        if url == 'NA':
            logger.error(f'url not available: {fn}')
            continue

        if not os.path.exists(f'{pw}/{fn}'):
            if fn in d_exist:
                n_uploaded += 1
                n_downloaded += 1
                logger.info(f'File already uploaded! {fn}: {d_exist[fn]}, size={size_exp:,}')
                os.system(f'mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh 2>/dev/null')
                continue

        fn_download = f'{pw}/download/{fn}'
        with open(f'{pw}/shell/{fn}.download_upload.{dest}.sh', 'w') as out:
            skip_download = 0
            if os.path.exists(fn_download):
                size = os.path.getsize(fn_download)
                if size == size_exp:
                    logger.info(f'file already downloaded, size match with expected {size_exp:,}: {fn}')
                    n_downloaded += 1
                    skip_download = 1
                else:
                    logger.warning(f'****** file already downloaded, but size not match: exp={size_exp:,} real={size} : {fn}')
                    # n_downloaded += 1
                    # if not demo:
                    #     os.system(f'ln -s {fn_download} {fn_download}.bad 2>/dev/null')
            if not skip_download:
                print(f'wget "{url}" -c -O "{fn_download}" > {pw}/log/download.{fn}.log 2>&1', file=out)

            print(f"""local_size=$(stat {fn_download} -c %s 2>/dev/null)\nif [ -z "$local_size" ];then\n\techo fail to get file size {fn_download} >&2\nelif [ "$local_size" -ne {size_exp} ];then\n\techo file size not match: {fn_download}: size=$local_size, expected={size_exp:,} >&2\n\tln -s {fn_download} {fn_download}.bad 2>/dev/null;\n\texit;\nfi""", file=out)

            rm_cmd = ''

            if rm:
                rm_cmd = f'\trm {fn_download}\n\trm {pw}/log/upload.{dest}.{fn}.log\n\trm {pw}/log/download.{fn}.log'

            if dest == 'dropbox':
                print(f'{dock} dbxcli put {pw}/download/{fn} /{remote_pw}{folder_extended}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1', file=out)
            elif dest == 'emedgene':
                print(f'date +"%m-%d  %T"> {pw}/log/upload.{dest}.{fn}.log\n{dock} aws s3 cp {pw}/download/{fn} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}{folder_extended}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1\ndate +"%m-%d  %T"> {pw}/log/upload.{dest}.{fn}.log', file=out)
                print(f"""
sleep 10
remote_file_size=$({dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}{folder_extended}/{fn}|cut -d " " -f 3)
if [[ $remote_file_size -eq {size_exp} ]];then
    echo {fn} successfully uploaded to {remote_pw}{folder_extended}/  filesize={size_exp:,};
    mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh
    {rm_cmd}
else
    echo file size not match: actual=$remote_file_size, expected={size_exp:,}: {fn} >&2
fi""", file=out)

    n_need_upload = n_desired - n_uploaded
    logger.info(f'total files = {n_all}, need_to_upload={n_need_upload}, already exist in server = {n_uploaded}, already downloaded = {n_downloaded}')

    if n_need_upload > 0:
        print('\n', '!' * 50)
        print(f'{n_need_upload} files need to be uploaded')
    else:
        print('\n', '*' * 50)
        print('Done, all files already uploaded')
    print('\n\n\n')
