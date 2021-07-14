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
"""
import os, sys, re, glob
import platform
import time
from pprint import pprint as pp
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db'], nargs='?', default='emedgene')
ps.add_argument('info_file', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""", nargs='?', default=None)
ps.add_argument('-remote', '-r', '-remote_pw', dest='remote_pw', help="""the remote path, if not exist, would create one. start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345""", default=None, nargs='?')
ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = fastq""", nargs='*')
ps.add_argument('-pw', help="""download file output path, default is pwd""")
ps.add_argument('-profile', help="""aws profile, default=emedgene, other valid could be pacbio""", default='emedgene')
ps.add_argument('-demo', help="""demo mode would not actually download or upload, or create remote folder""", action='store_true')
ps.add_argument('-rm', '-delete', help="""flag, if set, would delete the downloaded file from local disk when uploading is done""", action='store_true')
ps.add_argument('-noupload', '-noup', '-downonly', '-no_upload', '-download_only', '-no_up', '-down_only', help="""download only, don't upload the files""", action='store_true')
ps.add_argument('-asis', help="""use the remote_pw in the cmd line input, do not do the re-format, default is False""", action='store_true')
ps.add_argument('-lite', '-ck', help="""toggle, just check, donot build script""", action='store_true')
ps.add_argument('-force', '-f', help="""force create the shell file, even it already exist on the server""", action='store_true')
ps.add_argument('-skip_download', '-skip', '-nodown', help="""skip download the files""", action='store_true')
ps.add_argument('-force_upload', '-forceup', '-fu', help="""force upload the file to server, even the file size match with remote""", action='store_true')

args = ps.parse_args()
args = ps.parse_args()


# ps.add_argument('-prj', help="""project name for this run, if not set, would use the UDN infered from info file or current folder name""")

def getlogger():
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
    return logger


def build_remote_subfolder(name):
    if dest == 'dropbox' and not demo:
        logger.info(f'\tcreating subfolder: {name}')
        os.system(f'{dock} dbxcli mkdir {name}/ >>{dest}.create_folder.log 2>{dest}.create_folder.log')
    elif dest == 'emedgene' and not demo:
        logger.info(f'building_subfolder Vanderbilt/upload/{name}/')
        # os.system(f'{dock} aws s3api put-object --profile {profile} --bucket emg-auto-samples --key Vanderbilt/upload/{name}/ >>{dest}.create_folder.log 2>{dest}.create_folder.log')
        os.system(f'{dock} aws s3api put-object --bucket emg-auto-samples --key Vanderbilt/upload/{name}/ >>{dest}.create_folder.log 2>{dest}.create_folder.log')

def get_file_extension(fn):
    m = re.match(r'.*?\.(bam|bai|cnv|fastq|fq|gvcf|vcf)\b', fn.lower())
    convert = {'bai': 'bam', 'fq': 'fastq', 'gvcf': 'vcf', 'cnv.vcf': 'cnv'}
    convert.update({_: _ for _ in convert.values()})
    if m:
        try:
            return convert[m.group(1)]
        except:
            return m.group(1)
    else:
        logger.warning(f'unclear file type: {fn}')
        return None

def get_remote_file_list(dest, remote_pw):
    """
    return d_exist,  key = fn, value = [remote path, file size]
    """
    logger.info(f'remote_pw={remote_pw}, dest={dest}')
    remote_pw_plain = re.sub(r'\W+', '_', remote_pw)
    # fn_exist_complete = f'{pw}/remote.existing_files.{remote_pw_plain}.completed.txt'

    fn_exist = f'{pw}/remote.existing_files.{remote_pw_plain}.txt'
    # if not os.path.exists(fn_exist_complete):
    if dest == 'dropbox':
        os.system(f'{dock} dbxcli ls -R /{remote_pw} -l > {fn_exist}')
    elif dest == 'emedgene':
        os.system(f'{dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/ --recursive > {fn_exist}')

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

def parse_info_file(info_file, remote_pw_in=None, ft=None):
    """
    the info file is like
    rel_to_proband\tfn\turl\tudn\tseq_type\tsize\tbuild\tmd5\turl_s3\tdate_upload\tremote_pw\n
    the last column is optional
    return is a dict key = filename, value = dict,   url, remote_full_path, expected_size, download_type
    ft  - a list, the file type for keep, default is fastq
    """
    d = {}
    d_exist_all = {}


    # get the remote folder list
    fn_all_folder = 'all_folders.txt'
    os.system(f'{dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/ > {fn_all_folder} ')
    all_remote_folders = {}
    pattern = re.compile(r'.*(UDN\d+)(.*|$)')
    for i in open(fn_all_folder):
        i = i.rsplit(' ', 1)[-1].strip()
        i = re.sub('/$', '', i)
        m = re.match(pattern, i)
        if not m:
            continue
        all_remote_folders[m.group(1)] = i

    ft = ft or ['fastq']
    # remote_base_pw = set()

    remote_pw_in = refine_remote_pw(remote_pw_in, all_remote_folders)
    logger.info(f'remote_pw_in={remote_pw_in}')

    with open(info_file) as fp:
        header = fp.readline()
        a = [_.strip() for _ in header.split('\t')]
        try:
            idx_size = a.index('size')
        except:
            logger.error('file size not found in the info file header! exit')
            sys.exit(1)

        for i in fp:
            if not i.strip():
                continue
            a = i.split('\t')
            a = [_.strip() for _ in a]
            try:
                rel, fn, url, iudn, download_type = a[:5]
                size_exp = a[idx_size]
            except:
                logger.error(f'wrong info format: {a}')
                continue
            download_type = download_type.lower()
            if download_type != 'dropbox':
                download_type = 'wget'

            # get the remote_pw
            try:
                remote_pw = a[9]
                remote_pw = refine_remote_pw(remote_pw, all_remote_folders)

            except:
                if not remote_pw_in:
                    logger.error('must specify the remote_pw, in info file or by parameter, exit..')
                    sys.exit(1)
                remote_pw = remote_pw_in

            # remote_base_pw.add(remote_pw)

            try:
                d_exist, sub_folders = d_exist_all[remote_pw]
            except:
                d_exist, sub_folders = get_remote_file_list(dest, remote_pw)
                d_exist_all[remote_pw] = (d_exist, sub_folders)

            # get the ext
            ext = get_file_extension(fn)
            if 'all' not in ft and ext not in ft:
                logger.debug(f'file skipped due to file type not selected: {ext}  - {fn}')
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
            name = re.sub(r'^[\W_]*(.*?)[\W_]*$', '\g<1>', name)

            # if url.lower() == 'na':
            #     logger.warning(f'URL = NA : {fn}')

            fn_download = f'{pw}/download/{fn}'

            downloaded = 0
            if os.path.exists(fn_download):
                size_local = os.path.getsize(fn_download)
                if size_exp == 'na':
                    downloaded = 1
                elif size_local == size_exp:
                    downloaded = 1
                else:
                    logger.warning(f'file downloaded, but size not match. exp = {size_exp}, local={size_local} : {fn_download}')
                    if size_local == 0:
                        os.remove(fn_download)
                    elif not demo:
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

            v = {'size': size_exp, 'download_type': download_type, 'url': url, 'remote': f'{remote_pw}/{name}', 'downloaded': downloaded or uploaded, 'uploaded': uploaded, 'size_remote': f'{size_remote:,}'}
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
    if not demo:
        for i in sub_folder_to_build:
            build_remote_subfolder(i)
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


    for _, v1 in d.items():
        for fn, v in v1.items():
            n_desired += 1
            if v['url'].lower() == 'na':
                invalid_url.append(fn)
            if v['uploaded']:
                logger.info(f'file already uploaded: {fn}')
                n_uploaded += 1
                n_downloaded += 1
                continue

            if not v['downloaded']:
                need_upload.append(fn)
                # n_need_upload += 1
            else:
                n_downloaded += 1
                need_upload.append(fn)
                # n_need_upload += 1

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

    return d

# udn

def refine_remote_pw(remote_pw, all_remote_folders):
    if asis:
        return remote_pw

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


def build_script(d):
    for _, v1 in d.items():
        for fn, v in v1.items():
            # {'size': size_exp, 'remote': f'{remote_pw}/{name}', 'url': url, 'downloaded': downloaded, 'uploaded': uploaded}
            url = v['url']
            size_exp = v['size']
            remote_pw = v['remote']
            download_type = v['download_type']

            fn_download = f'{pw}/download/{fn}'
            fn_download_chunk = fn_download.replace('.gz', '')
            fn_script = f'{pw}/shell/{fn}.download_upload.{dest}.sh'
            fn_status = f'{pw}/log/status.{fn}.txt'

            if v['uploaded'] and not  args.force:
                os.system(f'mv  {fn_script} {pw}/shell_done/{fn}.download_upload.{dest}.sh 2>/dev/null')
                continue

            with open(fn_script, 'w') as out:
                if not skip_download and (args.force or not v['downloaded']):
                    if download_type == 'dropbox':
                        print(f'dbxcli get "{url}"  "{fn_download}" > {pw}/log/download.{fn}.log 2>&1', file=out)
                    else:
                        print(f'wget "{url}" -c -O "{fn_download}" > {pw}/log/download.{fn}.log 2>&1', file=out)

                print(f"""echo -n "" > {fn_status}""", file=out)
                rm_cmd = ''
                if rm:
                    rm_cmd = f'    rm {fn_download}\n    rm {pw}/log/upload.{dest}.{fn}.log\n    rm {pw}/log/download.{fn}.log'

                if dest == 'dropbox':
                    print(f'{dock} dbxcli put {pw}/download/{fn} /{remote_pw}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1', file=out)
                elif dest == 'emedgene':
                    out.write(f"""
local_size=$(stat -c "%s" {fn_download} 2>/dev/null)
remote_file_size=$({dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn}|tr -s " "|cut -d " " -f 3)
# echo $remote_file_size $local_size  >> {fn_status}
# exit
# check local file
if [[ "{size_exp}" -ne "na" ]] & [[ "$local_size" -ne "{size_exp}" ]];then
    echo local file not match with exp actual remote size=$remote_file_size, local_size=$local_size,  expected={size_exp}. skip uploading >> {fn_status}
    exit
fi

# check md5sum
if [[ ! -f {fn_download_chunk}.md5 ]];then
    md5sum  {fn_download} >{fn_download_chunk}.md5 2>{fn_download_chunk}.md5.log
fi

# check file tail
# if [[ ! -f {fn_download_chunk}.gztool.tail.txt ]];then
#     gztool -t {fn_download} >{fn_download_chunk}.gztool.tail.txt 2>{fn_download_chunk}.gztool.log
# fi

if [[ ! -f {fn_download_chunk}.zcat.tail.txt ]];then
    zcat {fn_download} |tail -20 >{fn_download_chunk}.zcat.tail.txt 2>{fn_download_chunk}.zcat.log
fi

# build the report
udn_fq_report {fn_download}  -tofile

""")
                    if not no_upload:
                        cmd_tmp = '' if force_upload else f"""
    mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh
    {rm_cmd}
    exit
"""
                        out.write(f"""

# determine need to upload or not
remote_file_size=$({dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn}|tr -s " "|cut -d " " -f 3)
if [[ -z $remote_file_size ]];then
   echo not uploaded yet  >> {fn_status}
elif [[ "$remote_file_size" -eq "{size_exp}" ]];then
    echo {fn} successfully uploaded to {remote_pw}/  filesize={size_exp} >> {fn_status}
    {cmd_tmp}
elif [[ "{size_exp}" = "na" ]] & [[ "$local_size" -eq "$remote_file_size" ]];then
    echo size check is not specified, remote_file_size=$remote_file_size local_filesize=$local_size >> {fn_status}
    {cmd_tmp}
else
    echo before uploading, file size not match: actual remote size=$remote_file_size, local_size=$local_size,  expected={size_exp}: {fn} >> {fn_status}
    exit
fi


echo start uploading >> {fn_status}
date +"%m-%d  %T">> {fn_status}
{dock} aws s3 cp {pw}/download/{fn} s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn} > {pw}/log/upload.{dest}.{fn}.log 2>&1\ndate +"%m-%d  %T">> {pw}/log/upload.{dest}.{fn}.log
echo upload completed
date +"%m-%d  %T">> {fn_status}
remote_file_size=$({dock} aws s3 ls emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn}|head -1|tr -s " "|cut -d " " -f 3)
if [[ "$remote_file_size" -eq "{size_exp}" ]];then
    echo {fn} successfully uploaded to {remote_pw}/  filesize={size_exp} >> {fn_status}
    mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh
    {rm_cmd}
    exit
elif [[ "{size_exp}" = "na" ]];then
    if [[ "$local_size" -eq "$remote_file_size" ]];then
        echo size check is not specified, remote_file_size=$remote_file_size local_filesize=$local_size >> {fn_status}
        mv {pw}/shell/{fn}.download_upload.{dest}.sh {pw}/shell_done/{fn}.download_upload.{dest}.sh
        {rm_cmd}
        exit
    else
        echo ERROR! localsize - $local_size  not match with remote size - $remote_file_size >> {fn_status}
    fi
else
    echo after uploading, file size not match: local=$local_size , remote actual=$remote_file_size, expected={size_exp}: {fn} >> {fn_status}
fi""")

def get_prj_name():
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



if __name__ == "__main__":

    logger = getlogger()

    convert1 = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene'}
    convert2 = {'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene', 'pacbio': 'pacbio', 'pac': 'pacbio'}
    dest = convert1[args.dest]
    demo = args.demo
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
    pw = args.pw or os.getcwd()
    # os.chdir(pw)
    print(args)
    print(f'current pw={pw}')

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

    if not os.path.exists(info_file):
        logger.error(f'info file not found')
        sys.exit(1)

    # get the file type to upload
    ft = []
    wrong_ft = []
    ft_convert = {'fastq': 'fastq', 'fq': 'fastq', 'bam': 'bam', 'vcf': 'vcf', 'bai': 'bam'}

    if not args.ft:
        ft = ['fastq']
    elif 'all' in args.ft:
        ft = ['all']
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
        sys.exit(1)

    node_name = platform.node()

    if node_name.find('viccbiostat120') > -1:
        dock = 'singularity exec /mnt/d/dock/centos.sif '
    elif node_name.find('vampire') > -1:
        dock = 'singularity exec -B /fs0 /data/cqs/chenh19/dock/centos.sif '

    # remote path
    pw = re.sub('/$', '', pw)
    remote_pw_in = args.remote_pw or pw.rsplit('/', 1)[-1]

    # build folder
    for i in ['shell', 'shell_done', 'download', 'log']:
        os.makedirs(f'{pw}/{i}', exist_ok=True)

    # parse the info file
    d = parse_info_file(info_file, remote_pw_in=remote_pw_in, ft=ft)
    if not lite:
        build_script(d)
