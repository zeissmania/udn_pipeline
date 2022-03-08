#! /usr/bin/env python
"""
run the upload task
upload to emedgene / dropbox / shared drive
the remote_pw  in info file would overwrite the -remote parameter

input = info file. columns
1=rel_to_proband, normal or 'na',
2=filename for upload
3=url for download, if don't need to download, use na
4=udn for this file, would put the file in the corresponsing subfolder. if na, would put the file under root remote folder
5=seq_type,  could be dropbox or wgs, wes or other.   if dropbox, would use dropbox for downloading the file. if other, would use regular wget
6=size,  int or 'na',  if na, would not check the downloaded and uploaded file size
7=build,  no effect
8=md5, no effect
9=url_s3, no effect
10=optional, remote pw

the following columns are must for the info file
'rel_to_proband', 'fn', 'url', 'udn', 'size'

optional = 'upload_type', 'download_type', 'remote_pw'

rule of remote pw

1. get relative(rel) and pure_udn(iudn) info from info.text
2. get the root remote_pw path, first check the  column remote_pw,
if not found, then used the remote_pw_in, the default folder

name = f'{rel}{iudn}'
remote_pw_final = f'{remote_pw}' if remote_flat else f'{remote_pw}/{name}'

"""
import os, sys, re, glob
import json
import platform
import time
from pprint import pprint as pp
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db'], nargs='?', default='emedgene')
ps.add_argument('info_file', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""", nargs='?', default=None)
ps.add_argument('-remote', '-r', '-remote_pw', dest='remote_pw', help="""the remote path, if not exist, would create one. start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345, if the remote folder is not in the root path, include all the path from root, e.g. cases_202107/34_269_AA_UDN616750""", default=None, nargs='?')
ps.add_argument('-remote_base', '-r_base', '-rbase', '-rb', help="remote path base, default is from the root of the remote folder")
ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = fastq""", nargs='*')
ps.add_argument('-pw', help="""the folder for each case, could be multiple, default = current folder, and would search for child folder for info.txt file""", nargs='*')
ps.add_argument('-profile', help="""aws profile, default=emedgene, other valid could be pacbio""", default='emedgene')
ps.add_argument('-demo', help="""demo mode, would not create remote folder, won't rename the local file with wrong file size""", action='store_true')
ps.add_argument('-rm', '-delete', help="""flag, if set, would delete the downloaded file from local disk when uploading is done""", action='store_true')
ps.add_argument('-noupload', '-noup', '-downonly', '-no_upload', '-download_only', '-no_up', '-down_only', help="""download only, don't upload the files""", action='store_true')
ps.add_argument('-asis', help="""use the remote_pw in the cmd line input, do not do the re-format, default is False""", action='store_true')
ps.add_argument('-lite', '-ck', help="""toggle, just check, donot build script""", action='store_true')
ps.add_argument('-force', '-f', help="""force create the shell file, even it already exist on the server""", action='store_true')
ps.add_argument('-skip_download', '-skip', '-nodown', help="""skip download the files""", action='store_true')
ps.add_argument('-allversion', '-all', help="""include all versions even the updated version exist""", action='store_true')
ps.add_argument('-remote_flat', '-no_remote_sub_folder', '-flat', help="""don't create remote subfolders""", action='store_true')
ps.add_argument('-force_upload', '-forceup', '-fu', help="""force upload the file to server, even the file size match with remote""", action='store_true')

args = ps.parse_args()


ft_convert = {'bai': 'bam', 'cnv.vcf': 'cnv', 'gvcf': 'gvcf', 'fq': 'fastq', 'vcf': 'vcf'}
ft_convert.update({_: _ for _ in ft_convert.values()})


# ps.add_argument('-prj', help="""project name for this run, if not set, would use the UDN infered from info file or current folder name""")

def getlogger():
    import logging
    prefix = 'download_upload_udn'
    fn_log = f'{prefix}.debug.log'
    fn_err = f'{prefix}.info.log'
    fmt = logging.Formatter('%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    console.setLevel('INFO')

    fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')
    fh_err = logging.FileHandler(fn_err, mode='w', encoding='utf8')

    fh_file.setLevel('DEBUG')
    fh_file.setFormatter(fmt)

    fh_err.setLevel('INFO')
    fh_err.setFormatter(fmt)

    logger = logging.getLogger(__file__)
    logger.setLevel('DEBUG')
    logger.addHandler(console)
    logger.addHandler(fh_file)
    logger.addHandler(fh_err)
    return logger


def build_remote_subfolder(name, dest):
    if demo:
        logger.info(f'demo mode, skip building remote subfolder')
        return 0
    if dest == 'dropbox':
        logger.info(f'\tcreating subfolder: {name}')
        os.system(f'{dock} dbxcli mkdir "{name}/" >>{dest}.create_folder.log 2>{dest}.create_folder.log')
    elif dest == 'emedgene':
        logger.info(f'building_subfolder Vanderbilt/upload/{name}/')
        # os.system(f'{dock} aws s3api put-object --profile {profile} --bucket emg-auto-samples --key Vanderbilt/upload/{name}/ >>{dest}.create_folder.log 2>{dest}.create_folder.log')
        os.system(f'{dock} aws s3api put-object --bucket emg-auto-samples --key "Vanderbilt/upload/{name}/" >>{dest}.create_folder.log 2>{dest}.create_folder.log')

def get_file_extension(fn):

    m = re.match(r'.*?\.(bam|bai|cnv|fastq|fq|gvcf|vcf)\b(\.gz)?', fn.lower())
    if m:
        try:
            return ft_convert[m.group(1)], m.group(2)
        except:
            return m.group(1), m.group(2)
    else:
        logger.warning(f'unclear file type: {fn}')
        return None, None


def get_remote_file_list(pw, dest, remote_pw):
    """
    return d_exist,  key = fn, value = [remote path, file size]
    """
    logger.info(f'remote_pw={remote_pw}, dest={dest}')
    remote_pw_plain = re.sub(r'\W+', '_', remote_pw)
    # fn_exist_complete = f'{pw}/remote.existing_files.{remote_pw_plain}.completed.txt'

    fn_exist = f'{pw}/remote.existing_files.{remote_pw_plain}.txt'
    # if not os.path.exists(fn_exist_complete):
    if dest == 'dropbox':
        os.system(f'{dock} dbxcli ls -R "/{remote_pw}" -l > {fn_exist}')
    elif dest == 'emedgene':
        os.system(f'{dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_pw}/" --recursive > {fn_exist}')

    d_exist = {}
    sub_folders = set()
    if dest == 'dropbox':

# Revision                        Size    Last modified Path
# -                               -       -             /CG2UDN/For Emedgene Upload
# 015be260c9d805e000000022be8dc00 9.5 KiB 1 day ago     /CG2UDN/cg2udn.key.xlsx
# 015be3846dafbcf000000022be8dc00 1.1 MiB 8 hours ago   /CG2UDN/Family VU097/Mom/20-097-002902_VUC2000652_D1_N1_LP200615017_UDI0067_HP200028_PL200021_SEQ_2006170017.292.final.vcf.gz

        with open(fn_exist) as fp:
            fp.readline()
            for i in fp:
                i = i.strip()
                a = re.split('/', i, 1)
                fn = a[-1].rsplit('/', 1)[-1]

                d_exist[fn] = [a[-1], 0]

                if i[0] == '-':
                    folder = a[-1].replace(remote_pw, '')
                    folder = re.sub(r'^/', '', folder)
                    folder = re.sub(r'/$', '', folder)
                    sub_folders.add(folder)
    elif dest == 'emedgene':
# 2021-02-22 15:45:30          0 Vanderbilt/upload/UDN782882_261_DP/Father_UDN098366/
# 2021-02-22 18:43:14 44042338086 Vanderbilt/upload/UDN782882_261_DP/Father_UDN098366/939022-UDN098366-D_R1.fastq.gz
# 2021-02-22 18:48:16 48156982037 Vanderbilt/upload/UDN782882_261_DP/Father_UDN098366/939022-UDN098366-D_R2.fastq.gz
# 2021-02-22 15:45:29          0 Vanderbilt/upload/UDN782882_261_DP/Mother_UDN954693/
# 2021-02-22 18:25:53 43856874123 Vanderbilt/upload/UDN782882_261_DP/Mother_UDN954693/939020-UDN954693-M_R1.fastq.gz
# 2021-02-22 18:53:31 47591653919 Vanderbilt/upload/UDN782882_261_DP/Mother_UDN954693/939020-UDN954693-M_R2.fastq.gz
# 2021-02-22 15:45:32          0 Vanderbilt/upload/UDN782882_261_DP/Proband_UDN782882/
# 2021-02-22 19:10:52 43542249072 Vanderbilt/upload/UDN782882_261_DP/Proband_UDN782882/939019-UDN782882-P_R1.fastq.gz
# 2021-02-22 19:49:40 47567355208 Vanderbilt/upload/UDN782882_261_DP/Proband_UDN782882/939019-UDN782882-P_R2.fastq.gz
# 2021-02-22 15:46:55       1829 Vanderbilt/upload/UDN782882_261_DP/UDN782882_all_files.md5
        with open(fn_exist) as fp:
            for i in fp:
                i = i.strip()
                a = re.split('Vanderbilt/upload/', i)
                fn = a[-1].rsplit('/', 1)[-1]
                if i[-1] == '/':
                    folder = a[-1].replace(remote_pw, '')
                    folder = re.sub(r'^/', '', folder)
                    folder = re.sub(r'/$', '', folder)
                    sub_folders.add(folder)


                file_size = re.split(r'\s+', i.strip())[2]
                try:
                    file_size = int(file_size)
                except:
                    logger.error(f'wrong remote emedgene file size: {file_size}')
                    sys.exit(1)

                if fn:
                    d_exist[fn] = [a[-1], file_size]

    return d_exist, sub_folders


def build_upload_script_simple(pw, fn, remote_pw, fn_new=None):
    fn_pure = fn_new or fn
    fn_pure = fn_pure.rsplit('/', 1)[-1]
    if not os.path.exists(fn):
        logger.warning(f'file not exist: {fn}')
        return 1

    if os.path.exists(f'{pw}/shell_done/{fn_pure}.upload.sh'):
        return 0

    fn_script = f'{pw}/shell/{fn_pure}.upload.sh'
    fn_status = f'{pw}/log/status.{fn_pure}.txt'
    with open(fn_script, 'w') as out:
        out.write(f"""#! /usr/bin/env bash
date +"%m-%d  %T">> {fn_status}
local_size=$(stat -c "%s" {fn} 2>/dev/null)
remote_file_size=$({dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn_pure}"|head -1|tr -s " "|cut -d " " -f 3)

if [[ "$local_size" -eq "$remote_file_size" ]];then
    echo already uploaded: {fn}
    exit 0
fi

{dock} aws s3 cp {fn} "s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn_pure}" > {pw}/log/upload.{dest}.{fn_pure}.log 2>&1
echo upload completed > {pw}/log/upload.{dest}.{fn_pure}.log
date +"%m-%d  %T">> {fn_status}

if [[ "$remote_file_size" -eq "$local_size" ]];then
    echo {fn} successfully uploaded to {remote_pw}/  filesize=$local_size >> {fn_status}
    mv "{fn_script}" {pw}/shell_done
fi
""")


def parse_info_file(pw, info_file, remote_pw_in=None, ft=None, gzip_only=None, updated_version_only=True, script_list=None, no_upload=False, remote_flat=False):
    """
    the info file is like
    rel_to_proband\tfn\turl\tudn\tseq_type\tsize\tbuild\tmd5\turl_s3\tdate_upload\tremote_pw\n
    the last column is optional
    return is a dict key = filename, value = dict,   url, remote_full_path, expected_size, download_type
    ft  - a list, the file type for keep, default is fastq
    remote_base_dir, default is empty, would search from the root of the emedgene server , which is emg-auto-samples/Vanderbilt/upload/
    gzip only, a set, if the ext is found in this set, then, the file must be gziped to get involved. e.g. if vcf in gzip_only, then,  the file named .vcf$ would not be included

    script_list,  in order to decide script run order across cases k = needdown, needup.
    """
    d = {}
    d_exist_all = {}


    script_list = script_list or {'needdown': [], 'needup': []}

    if not isinstance(gzip_only, set):
        gzip_only = set()

    ft = ft or ['fastq']
    # remote_base_pw = set()

    fn_all_folder = 'all_folders.txt'
    try:
        os.unlink(fn_all_folder)
    except:
        pass

    logger.info(f'remote_pw_in={remote_pw_in}')
    remote_pw_in = refine_remote_pw(remote_pw_in)


    logger.info(f'remote_pw afer refine={remote_pw_in}')
    try:
        d_exist, sub_folders = d_exist_all[remote_pw_in]
    except:
        d_exist, sub_folders = get_remote_file_list(pw, dest, remote_pw_in)
        d_exist_all[remote_pw_in] = (d_exist, sub_folders)

    # build script for md5 file
    if not no_upload:
        fn_md5 = glob.glob(f'{pw}/download.*.md5')
        if len(fn_md5) == 0:
            logger.warning(f'expected md5 file not found')
        else:
            if len(fn_md5) > 1:
                logger.warning(f'more than one md5 file found, only first would be uploaded: {fn_md5}')
            fn_md5 = fn_md5[0]
            fn_new = fn_md5.rsplit('/', 1)[-1].replace('download.', '')

            if fn_new in d_exist:
                logger.info(f'{fn_md5} already uploaded')
            else:
                # build the script
                # logger.info(f'upload md5 file')
                build_upload_script_simple(pw, fn_md5, remote_pw_in, fn_new=fn_new)

    logger.info('start parsing')
    with open(info_file) as fp:
        header = fp.readline()
        a = [_.strip() for _ in header.split('\t')]
        idx = {}
        for _ in ['upload_type', 'download_type', 'remote_pw', 'rel_to_proband', 'fn', 'url', 'udn', 'size', 'rename']:
        # rename, if set, will change both the filename and vcf header of this file to this (if no relative suffix found, will also add the relative)
        # rename is only valid for vcf files.

            try:
                idx[_] = a.index(_)
            except:
                idx[_] = 'NA'


        for i in fp:
            if not i.strip():
                continue
            a = i.split('\t')
            a = [_.strip() for _ in a]
            try:
                rel, fn, url, iudn, size_exp = [a[idx[_]] for _ in ['rel_to_proband', 'fn', 'url', 'udn', 'size']]
            except:
                logger.error(f'wrong info format: {a}')
                logger.error(f'idx = {idx}')
                raise

            download_type = 'wget'
            try:
                tmp = a[idx['download_type']].strip()
                if tmp:
                    download_type = tmp
            except:
                pass

            fatal = 0
            if download_type not in {'wget', 'dropbox'}:
                logger.error(f'invalid download type: {download_type}, valid = wget or dropbox')
                fatal = 1

            upload_type = 'emedgene'
            try:
                tmp = a[idx['upload_type']].strip()
                if tmp:
                    upload_type = tmp
            except:
                pass

            if upload_type not in {'emedgene', 'dropbox'}:
                logger.error(f'invalid upload_type: {upload_type}, valid = emedgene or dropbox')
                fatal = 1

            if fatal:
                sys.exit(1)


            # get the remote_pw
            try:
                remote_pw = a[idx['remote_pw']].strip()
                if not remote_pw:
                    if not remote_pw_in:
                        logger.error('must specify the remote_pw, in info file or by parameter, exit..')
                        return 1
                    remote_pw = remote_pw_in
                else:
                    remote_pw = refine_remote_pw(remote_pw)
            except:
                if not remote_pw_in:
                    logger.error('must specify the remote_pw, in info file or by parameter, exit..')
                    return 1
                remote_pw = remote_pw_in

            # remote_base_pw.add(remote_pw)


            try:
                d_exist, sub_folders = d_exist_all[remote_pw]
            except:
                d_exist, sub_folders = get_remote_file_list(pw, dest, remote_pw)
                d_exist_all[remote_pw] = (d_exist, sub_folders)

            # get the ext
            ext, gz = get_file_extension(fn)
            if 'all' not in ft and ext not in ft:
                logger.debug(f'file skipped due to file type not selected: {ext}  - {fn}')
                continue

            if ext in gzip_only and not gz:
                logger.info(f'file skipped due to gzip file only: {fn}')
                continue

            rel = re.sub(r'\W+', '_', rel)
            rel = '' if not rel or rel.lower() == 'na' else f'{rel}_'

            try:
                size_exp = int(size_exp)
            except:
                if size_exp.strip() and size_exp.lower() != 'na':
                    logger.warning(f'wrong filesize format: {size_exp}, line=@{a}@')
                size_exp = 'na'

            if iudn == 'na':
                iudn = ''
            name = f'{rel}{iudn}'
            name = re.sub(r'^[\W_]*(.*?)[\W_]*$', r'\g<1>', name)


            rename = None
            try:
                tmp = a[idx['rename']].strip()
                if tmp and ext == 'vcf':
                    m_tmp = re.match(r'.*\b(indel|snp)([a-z]*)[\W_]', fn, flags=re.I)
                    if m_tmp:
                        rename_suffix = '_' + ''.join(m_tmp.groups())
                    else:
                        rename_suffix = ''

                    rename = tmp
                    if not re.match('.*(father|dad|mother|mom|sister|brother|sibling|mate|pate|cousin|uncle|aunt|proband)', rename.lower()):
                        rename = f'{rename}_{rel[:-1]}{rename_suffix}'
            except:
                raise
                pass


            # if url.lower() == 'na':
            #     logger.warning(f'URL = NA : {fn}')

            fn_download = f'{pw}/download/{fn}'

            downloaded = 0
            if os.path.exists(fn_download):
                size_local = os.path.getsize(fn_download)
                if size_local == 0:
                    os.unlink(fn_download)
                    logger.warning('empty file: {fn_download}')
                elif size_exp == 'na':
                    downloaded = 1
                elif size_local == size_exp:
                    downloaded = 1
                else:
                    logger.warning(f'file downloaded, but size not match. exp = {size_exp}, local={size_local} : {fn_download}')
                    if not demo:
                        try:
                            os.symlink(fn_download, f'{fn_download}.bad')
                        except:
                            pass
            try:
                size_remote = d_exist[fn][1]
            except:
                uploaded = 0
                size_remote = 0
            else:
                if size_exp == 'na':
                    uploaded = 1
                elif size_exp == size_remote:
                    uploaded = 1
                else:
                    logger.warning(f'file uploaded, but size not match. exp = {size_exp}, remote={size_remote}')
                    uploaded = 0


            remote_pw_final = f'{remote_pw}' if remote_flat else f'{remote_pw}/{name}'
            remote_pw_final = re.sub(r'/$', '', remote_pw_final)

            v = {'size': size_exp, 'download_type': download_type, 'upload_type': upload_type, 'url': url, 'remote': remote_pw_final, 'downloaded': downloaded or uploaded, 'uploaded': uploaded, 'size_remote': f'{size_remote:,}', 'rename': rename}
            try:
                d[remote_pw][fn] = v
            except:
                d[remote_pw] = {fn: v}

    # check if all files have been uploaded
    # for remote_pw, v in d.items():
    #     remote_pw_plain = re.sub(r'\W+', '_', remote_pw)
    #     all_uploaded = 1 if all([_['uploaded'] for _ in v.values()]) else 0

    logger.debug(d)
    sub_folders_exist = set()
    for remote_pw, v in d_exist_all.items():
        i = set([f'{remote_pw}/{name}' for name in v[1]])
        sub_folders_exist |= i

    logger.debug(f'remote existing subfolders = {sub_folders_exist}')


    sub_folders_exp = set()
    for v1 in d.values():
        for v2 in v1.values():
            sub_folders_exp.add(v2['remote'])

    logger.debug(f'subfolder exp = {sub_folders_exp}')

    sub_folder_to_build = sub_folders_exp - sub_folders_exist
    if not demo and not no_upload:
        for i in sub_folder_to_build:
            build_remote_subfolder(i, dest=upload_type)
    if len(sub_folder_to_build) > 0:
        logger.info(f'demo build remote subfolder {sub_folder_to_build}')
    else:
        logger.info(f'all remote subfolder are ready')


    # check the status
    n_downloaded = 0
    n_uploaded = 0
    n_desired = 0
    need_upload = []
    invalid_url = []

    # check for the updated tag
    updated_tag = set()
    all_fn = set()

    for _, v1 in d.items():
        for fn, v in v1.items():
            all_fn.add(fn)
            if fn.lower().find('update') > 0:
                tag = fn.split('update')[0]
                tag = re.sub(r'\W*$', '', tag)
                updated_tag.add(tag)

    if len(updated_tag) > 0:
        logger.info(f'files with updated flag = {len(updated_tag)}')
        all_fn_new = {_ for _ in all_fn if _.split('.fastq')[0] not in updated_tag}

        excluded_file = all_fn - all_fn_new
        excluded_file_str = '\n\t'.join(sorted(excluded_file))
        len1 = len(excluded_file)
        if len1 > 0:
            if updated_version_only:
                logger.warning(f'{len1} files would be omitted due to old version:\n\t{excluded_file_str}')
            else:
                logger.warning(f'{len1} files have updated version, both version would be updated:\t{excluded_file_str}')
    else:
        excluded_file = set()

    for _ in list(d):
        v1 = d[_]
        for fn in list(v1):
            if fn in excluded_file:
                del d[_][fn]
                continue
            v = v1[fn]

            fn_script = f'{pw}/shell/{fn}.download_upload.emedgene.sh'
            n_desired += 1
            if v['url'].lower() == 'na':
                invalid_url.append(fn)
                continue
            if v['uploaded']:
                logger.info(f'file already uploaded: {fn}')
                n_uploaded += 1
                n_downloaded += 1
                continue

            if not v['downloaded']:
                need_upload.append(fn)
                script_list['needdown'].append(fn_script)
                # n_need_upload += 1
            else:
                n_downloaded += 1
                need_upload.append(fn)
                script_list['needup'].append(fn_script)

    n_need_upload = len(need_upload)
    logger.info(f'total files = {n_desired}, need_to_upload={n_need_upload}, already exist in server = {n_uploaded}, already downloaded = {n_downloaded}')

    n_invalid_url = len(invalid_url)
    if n_need_upload > 0 or n_invalid_url > 0:
        print('\n', '!' * 50)
        print(f'{n_downloaded}/{n_desired} files already downloaded')
        print(f'{n_need_upload}/{n_desired} files need to be uploaded')
        print('\t' + '\n\t'.join(need_upload))
        if n_invalid_url > 0:
            print(f'ERROR: {n_invalid_url} files with NA url:\n\t' + '\n\t'.join(invalid_url))

    else:
        print('\n', '*' * 50)
        print('Done, all files already uploaded')
    print('\n\n\n')

    ct = {'total': n_desired, 'downloaded': n_downloaded, 'uploaded': n_uploaded, 'need_to_upload': n_need_upload}

    return d, ct, script_list


def refine_remote_pw(remote_pw):
    if asis:
        return remote_pw

    remote_pw = re.sub('/$', '', remote_pw)
    remote_base_dir = ''
    tmp = remote_pw.rsplit('/', 1)

    if len(tmp) == 2:
        # the remote pw is like layer1/layer2 , e.g. not simple
        remote_base_dir = tmp[0] + '/'
        remote_pw = tmp[-1]

    # get the remote folder list
    fn_all_folder = 'all_folders.txt'
    if not os.path.exists(fn_all_folder):
        os.system(f'{dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_base_dir}" > {fn_all_folder} ')
    all_remote_folders = {}
    pattern = re.compile(r'.*(UDN\d+)(.*|$)')
    n = 0

    skip_words = ['Trio']  # if the path name contains these words, would not be considered as existed path

    for i in open(fn_all_folder):
        i = i.rsplit(' ', 1)[-1].strip()
        n += 1
        i = re.sub('/$', '', i)
        skip = 0
        for iword in skip_words:
            if i.find(iword) > -1:
                skip = 1
                break
        if skip:
            continue

        m = re.match(pattern, i)
        if not m:
            continue
        all_remote_folders[m.group(1)] = f'{remote_base_dir}{i}'

    if n == 0:
        logger.error(f'folder not exist: emg-auto-samples/Vanderbilt/upload/{remote_base_dir}')
        sys.exit(1)

    remote_pw = re.sub('^/', '', remote_pw)
    remote_pw = re.sub('upload_?', '', remote_pw, flags=re.I)
    remote_pw = re.sub(r'^\d+_', '', remote_pw, 1)
    m = re.match(r'(.*_)?(UDN\d+)(_.*)?', remote_pw)
    if not m:
        return remote_pw
    else:
        m = m.groups()
        udn_tmp = m[1]
        if udn_tmp in all_remote_folders:
            remote_pw = all_remote_folders[udn_tmp]
            logger.info(f'remote folder already exist: {remote_pw}')
            return remote_pw

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
        return remote_pw


def build_script(pw, d, info_file, no_upload=False):
    for _, v1 in d.items():
        for fn, v in v1.items():
            # {'size': size_exp, 'remote': f'{remote_pw}/{name}', 'url': url, 'downloaded': downloaded, 'uploaded': uploaded}
            url = v['url']
            size_exp = v['size']
            remote_pw = v['remote']
            download_type = v['download_type']
            upload_type = v['upload_type']
            ext = fn.rsplit('.', 1)[-1]
            rename = v['rename']

            fn_download = f'{pw}/download/{fn}'
            # fn_download_chunk = fn_download.replace('.gz', '')
            fn_script = f'{pw}/shell/{fn}.download_upload.{dest}.sh'
            fn_status = f'{pw}/log/status.{fn}.txt'

            if v['uploaded'] and not  args.force:
                os.system(f'mv  {fn_script} {pw}/shell_done/{fn}.download_upload.{dest}.sh 2>/dev/null')
                continue

            if url.lower() == 'na':
                continue

            with open(fn_script, 'w') as out:

                get_url = f"""url=$(awk -F "\\t" '$2=="{fn}" {{print $4}}' {info_file})"""
                print(get_url, file=out)
                rm_cmd = ''
                if rm:
                    rm_cmd = f'    rm {fn_download}\n    rm {pw}/log/upload.{dest}.{fn}.log\n    rm {pw}/log/download.{fn}.log'

                if download_type == 'dropbox':
                   download_cmd = f'dbxcli get "$url"  "{fn_download}" > {pw}/log/download.{fn}.log 2>&1'
                elif download_type == 'wget':
                    download_cmd = f'wget "$url" -c -O "{fn_download}" > {pw}/log/download.{fn}.log 2>&1'

                cmd = f"""
checklocal(){{
    local_size=$(stat -c "%s" {fn_download} 2>/dev/null)
    if [[ $local_size -eq 0 ]];then
        rm {fn_download} 2>/dev/null
        echo 1
        return
    fi

    # check local file
    if [[ "{size_exp}" -ne "na" ]] & [[ "$local_size" -ne "{size_exp}" ]];then
        echo local file not match with exp local_size=$local_size,  expected={size_exp}. skip uploading >> {fn_status}
        echo 1
        return
    else
        # the local size is ok
        echo 0
    fi

    # check md5sum
    if [[ ! -f {pw}/download/{fn}.md5 ]];then
        md5sum  {fn_download} >{pw}/download/{fn}.md5 2>{pw}/download/{fn}.md5.log
    fi

}}
"""
                if upload_type == 'emedgene':
                    cmd += f"""
checkremote(){{
    cd {pw}/download
    local_size=$(stat -c "%s" {fn_download} 2>/dev/null)
    # determine need to upload or not
    remote_file_size=$({dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn}"|grep "{fn}$" | tr -s " "|cut -d " " -f 3)
    if [[ "{size_exp}" != "na" ]] & [[ "$local_size" -ne "{size_exp}" ]];then
        echo "wrong_size"
    elif [[ "$remote_file_size" -eq "{size_exp}" ]];then
        echo {fn} successfully uploaded to {remote_pw}/  filesize={size_exp} >> {fn_status}
        echo already_uploaded
    elif [[ "{size_exp}" = "na" ]] & [[ "$local_size" -eq "$remote_file_size" ]];then
        echo size check is not specified, remote_file_size=$remote_file_size local_filesize=$local_size >> {fn_status}
        echo already_uploaded
    elif    [[ -z $remote_file_size ]];then
        echo not uploaded yet  >> {fn_status}
        echo upload
    else
        echo remote size not correct: remote = $remote_file_size, exp = {size_exp}
        echo upload
    fi

}}
"""
                else:
                    # there is no way for dropbox to get byte file size
                    cmd += f"""
checkremote(){{
    file_exist_from_remote=$({dock} dbxcli ls -l "/{remote_pw}/"|grep "{fn}")
    if [[ ! -z "$file_exist_from_remote" ]];then
        echo alrady_uploaded
    else
        echo upload
    fi

}}

                    """
                cmd += f"""
postupload(){{
    mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh
    {rm_cmd}
    exit 0
}}

download_status=$(checklocal)

"""
                if not no_upload:
                    cmd += f"""
upload_status=$(checkremote)

if [[ $upload_status = "already_uploaded" ]];then
    postupload
    exit 0
fi

"""
                cmd += f"""

if [[ $download_status -eq 1 ]];then
    {download_cmd}
    download_status=$(checklocal)
    if [[ $download_status -eq 1 ]];then
        echo download not completed >> {fn_status}
        exit 1
    fi
fi
"""
                # rename the vcf and bgzip
                if rename:
                    fn = f'{rename}.vcf.gz'
                    # prev command is
                    # bcftools view --no-version {fn_download}|bcftools reheader -s  <(echo "{rename}") |bgzip > {pw}/download/{fn}.tmp
                    # mv {pw}/download/{fn}.tmp {pw}/download/{fn}

                    cmd += f"""
# rename the fq file
vcf_deid {fn_download} -newname {rename} |bgzip > {pw}/download/{fn}
bgzip -r {pw}/download/{fn}
                    """

                if upload_type == 'dropbox':
                    upload_cmd = f'{dock} dbxcli put {pw}/download/{fn} "/{remote_pw}/{fn}" > {pw}/log/upload.{dest}.{fn}.log 2>&1'
                elif upload_type == 'emedgene':
                    upload_cmd = f'{dock} aws s3 cp "{pw}/download/{fn}" "s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn}" > {pw}/log/upload.{dest}.{fn}.log 2>&1\ndate +"%m-%d  %T">> {pw}/log/upload.{dest}.{fn}.log'

                if not no_upload:
                    cmd += f"""

# local size is ok, try upload
echo start uploading >> {fn_status}
date +"%m-%d  %T">> {fn_status}

{upload_cmd}
echo upload finish >> {fn_status}
date +"%m-%d  %T"  >> {fn_status}
upload_status=$(checkremote)
if [[ $upload_status = "already_uploaded" ]];then
    postupload
    exit 0
else
    echo "fail to upload"
    exit 1
fi
                """
                print(cmd, file=out)
            os.chmod(fn_script, 0o755)

def get_prj_name(pw):
    # get project name, not used
    if not args.prj:
        m = re.match(r'.*download\.info\.(UDN\d+)\.txt', info_file)
        if m:
            prj = m.group(1)
        else:
            prj = pw.rsplit('/', 1)[-1]
            prj = re.sub(r'\W+', '_', prj)
    else:
        prj = args.prj

    if prj in ['udn', 'upload', 'module', 'jb']:
        logger.error(f'invalid prj name specified, you may in the wrong working folder: {prj}, exit...')
        sys.exit(1)

def main(pw, script_list, info_file=None, remote_pw_in=None, updated_version_only=True, no_upload=False, remote_flat=False):
    if not info_file:
        tmp = glob.glob(f'{pw}/download.info.*.txt')
        if len(tmp) == 0:
            logger.error(f'no download info file found! exit...')
            return 1
        elif len(tmp) > 1:
            logger.error(f'multiple info file found! \n{tmp}\nexit...')
            return 1
        else:
            info_file = tmp[0]

    if not os.path.exists(info_file):
        logger.error(f'info file not found')
        sys.exit(1)

    # remote path
    pw = re.sub('/$', '', pw)
    remote_pw_in = remote_pw_in or pw.rsplit('/', 1)[-1]

    # print(f'pw={pw}, info={info_file} remote_pw_in={remote_pw_in}')
    # return 0

    # build folder
    for i in ['shell', 'shell_done', 'download', 'log']:
        os.makedirs(f'{pw}/{i}', exist_ok=True)
    # parse the info file
    d, ct, script_list = parse_info_file(pw, info_file, remote_pw_in=remote_pw_in, ft=ft, updated_version_only=updated_version_only, script_list=script_list, no_upload=no_upload, remote_flat=remote_flat)
    if not lite:
        build_script(pw, d, info_file, no_upload=no_upload)
    return ct, script_list


if __name__ == "__main__":

    logger = getlogger()

    convert1 = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene'}
    convert2 = {'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene', 'pacbio': 'pacbio', 'pac': 'pacbio'}
    dest = convert1[args.dest]
    demo = args.demo

    remote_flat = args.remote_flat
    lite = args.lite
    if lite:
        demo=True
    skip_download = args.skip_download
    force_upload = args.force_upload
    no_upload = args.noupload
    profile = args.profile
    profile = convert2[profile]
    asis = args.asis  # use the remote_pw in the cmd line input, do not do the re-format
    rm = args.rm
    remote_pw = args.remote_pw
    remote_base = args.remote_base
    print(args)

    node_name = platform.node()
    updated_version_only = not args.allversion

    if node_name.find('viccbiostat120') > -1:
        dock = 'singularity exec -B /mnt /mnt/d/dock/centos.sif '
    elif node_name.find('vampire') > -1:
        dock = 'singularity exec -B /fs0 /data/cqs/chenh19/dock/centos.sif '

    # get pw
    pw_raw = args.pw or [os.getcwd()]
    pw = []

    for ipw in pw_raw:
        tmp = ipw.rsplit('/', 1)[-1]
        info_file = glob.glob(f'{ipw}/download.info*.txt')
        len1 = len(info_file)
        if len1 == 0:
            sub_folders = glob.glob('*/')
            for isub in sub_folders:
                isub = f'{ipw}/{isub}'
                len2 = len(glob.glob(f'{isub}/download.info*.txt'))
                if len2 > 0:
                    pw.append(isub)
        else:
            pw.append(ipw)

    if len(pw) > 1 and remote_pw:
        logger.error('multiple local path found, but only 1 remote_pw specified, try specify one at each or specify remote_base, then use auto assigned name')
        sys.exit(1)

    logger.info(f'total cases = {len(pw)}\n{"*"*50}\n\n')
    pw = sorted(pw)

    # if len(pw) > 1 and not remote_base:
    #     logger.error(f'multiple local path found, you must specify a remote_base')
    #     sys.exit(1)

    # get the file type to upload
    gzip_only = set()
    ft = set()

    if not args.ft:
        ft = ['fastq']
    elif 'all' in args.ft:
        ft = ['all']
    else:
        err = 0
        for i in args.ft:
            ext = i.replace('.gz', '')
            if i[-2:] == 'gz':
                gzip_only.add(ext)
            try:
                ft.add(ft_convert[ext])
            except:
                logger.error(f'invalid filetype: {i}')
                err = 0
        if err:
            sys.exit(1)

    if len(ft) == 0:
        logger.error(f'No file type selected! exit')
        sys.exit(1)

    ct = {'total': 0, 'downloaded': 0, 'uploaded': 0, 'need_to_upload': 0}
    script_list = {'needup': [], 'needdown': []}

    for ipw in pw:
        ipw = re.sub(r'/$', '', ipw)
        pw_core = ipw.rsplit('/', 1)[-1]
        ipw = os.path.realpath(ipw)
        # logger.warning(f'path={ipw}')
        # continue

        fn_md5 = glob.glob(f'{ipw}/shell/*.md5.upload.sh')
        if len(fn_md5) > 0 and not no_upload:
            script_list['needup'].append(fn_md5[0])

        remote_base_prefix = f'{remote_base}/' if remote_base else ''
        remote_pw_in = remote_pw or f'{remote_base_prefix}{pw_core}'
        ct_prj, script_list = main(ipw, script_list, remote_pw_in=remote_pw_in, updated_version_only=updated_version_only, no_upload=no_upload, remote_flat=remote_flat)

        for k in ct:
            ct[k] += ct_prj[k]

        # ct = {'total': n_desired, 'downloaded': n_downloaded, 'uploaded': n_uploaded, 'need_to_upload': n_need_upload}

    print('\n' + '#' * 50)
    print(json.dumps(ct, indent=4))


    # build the script list
    n_needup = len(script_list['needup'])
    n_needdown = len(script_list['needdown'])

    print(f'total scripts = {n_needup + n_needdown}')

    with open(f'script_list.txt', 'w') as o:
        for k in ['needup', 'needdown']:
            if len(script_list[k]) > 0:
                print('\n'.join(script_list[k]), file=o)
