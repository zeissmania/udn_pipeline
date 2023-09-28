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
from termcolor import colored
import time
from pprint import pprint as pp
import argparse as arg
from argparse import RawTextHelpFormatter
ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
ps.add_argument('local_files', help="""local files list, for simple mode only""", nargs='*')
ps.add_argument('-dest', help="""the destinition for the uploading, valid = dropbox/db, emedgene/em/ed, r/rdrive""", choices=['emedgene', 'em', 'ed', 'dropbox', 'db', 'rdrive', 'r'], nargs='?', default='emedgene')
ps.add_argument('-info_file', '-info', '-fn', help="""the file containing the amazon download link, usually like download.info.UDN945671.txt""", nargs='?', default=None)
ps.add_argument('-simple', help="""simple mode, for uploading local file to dropbox/emedgene only. no info file needed, just set the -remote  -dest, and put all the local file path as the positional arugments""", action='store_true')
ps.add_argument('-local_pw_layers', '-local_depth', '-d', '-localn', help="the local folder structure to keep, start from lower to upper. eg. a/b/c/1.txt  if set as 1, will keep c/1.txt, default = 0, if set as 2, will keep b/c/1.txt", type=int, default=0)
ps.add_argument('-remote', '-r', '-remote_pw', dest='remote_pw', help="""the remote path, if not exist, would create one. start from the root path, donot include the leading and ending slash. dropbox is like UDN12345/proband_UDN12345  for emedgene is like UDN12345, if the remote folder is not in the root path, include all the path from root, e.g. cases_202107/34_269_AA_UDN616750""", default=None, nargs='?')
ps.add_argument('-remote_base', '-r_base', '-rbase', '-rb', help="remote path base, default is from the root of the remote folder")
ps.add_argument('-ft', help="""the file type to upload, could be multiple types sep by space, such as fastq, fq, vcf, bam, default = fastq, if upload all files in the info file, use all""", nargs='*')
ps.add_argument('-pw', help="""the folder for each case, could be multiple, default = current folder, and would search for child folder for info.txt file""", nargs='*')
ps.add_argument('-profile', help="""aws profile, default=emedgene, other valid could be pacbio""", default='emedgene')
ps.add_argument('-demo', help="""demo mode, would not create remote folder, won't rename the local file with wrong file size""", action='store_true')
ps.add_argument('-rm', '-delete', help="""flag, if set, would delete the downloaded file from local disk when uploading is done""", action='store_true')
ps.add_argument('-noupload', '-noup', '-downonly', '-no_upload', '-download_only', '-no_up', '-down_only', help="""download only, don't upload the files""", action='store_true')
ps.add_argument('-asis', help="""use the remote_pw in the cmd line input, do not do the re-format, default is False""", action='store_true')
ps.add_argument('-lite', '-ck', help="""toggle, just check, donot build script""", action='store_true')
ps.add_argument('-force', '-f', help="""force create the shell file, even it already exist on the server""", action='store_true')
ps.add_argument('-force_clear', '-fc', help="""force clear the download folder without warning""", action='store_true')
ps.add_argument('-skip_download', '-skip', '-nodown', help="""skip download the files""", action='store_true')
ps.add_argument('-allversion', '-all', help="""include all versions even the updated version exist""", action='store_true')
ps.add_argument('-remote_flat', '-no_remote_sub_folder', '-flat', help="""don't create remote subfolders""", action='store_true')
ps.add_argument('-force_upload', '-forceup', '-fu', help="""force upload the file to server, even the file size match with remote""", action='store_true')
ps.add_argument('-forcedown', '-fd',  help="""force download the file even if the file already uploaded. but if the local file already exist, will skip""", action='store_true')
ps.add_argument('-clear', '-rmlog', help="""clear the log files for upload and download""", action='store_true')
ps.add_argument('-i', help="""ignore the files with url as NA""", action='store_true')
ps.add_argument('-suffix', '-s',  help="""add the suffix for all files in the info file, usually for reupload of previous uploaded file, usually need to use together with -asis""")
ps.add_argument('-nobuild', help="""do not build the remote pw""", action='store_true')

args = ps.parse_args()

nobuild = args.nobuild

ft_convert = {'bai': 'bam', 'cnv.vcf': 'cnv', 'gvcf': 'gvcf', 'fasta': 'fasta', 'fq': 'fastq', 'vcf': 'vcf',  'crai': 'cram', 'txt': 'txt'}
ft_convert.update({_: _ for _ in ft_convert.values()})

file_ct = {}

# ps.add_argument('-prj', help="""project name for this run, if not set, would use the UDN infered from info file or current folder name""")

def colored(text, color=None, on_color=None, attrs=None):
    # Copyright (c) 2008-2011 Volvox Development Team
    #
    # Permission is hereby granted, free of charge, to any person obtaining a copy
    # of this software and associated documentation files (the "Software"), to deal
    # in the Software without restriction, including without limitation the rights
    # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    # copies of the Software, and to permit persons to whom the Software is
    # furnished to do so, subject to the following conditions:
    #
    # The above copyright notice and this permission notice shall be included in
    # all copies or substantial portions of the Software.
    #
    # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    # THE SOFTWARE.
    #
    # Author: Konstantin Lepa <konstantin.lepa@gmail.com>

    ATTRIBUTES = dict(
            list(zip([
                'bold',
                'dark',
                '',
                'underline',
                'blink',
                '',
                'reverse',
                'concealed'
                ],
                list(range(1, 9))
                ))
            )
    del ATTRIBUTES['']


    HIGHLIGHTS = dict(
            list(zip([
                'on_grey',
                'on_red',
                'on_green',
                'on_yellow',
                'on_blue',
                'on_magenta',
                'on_cyan',
                'on_white'
                ],
                list(range(40, 48))
                ))
            )


    COLORS = dict(
            list(zip([
                'grey',
                'red',
                'green',
                'yellow',
                'blue',
                'magenta',
                'cyan',
                'white',
                ],
                list(range(30, 38))
                ))
            )
    COLORS['purple'] = 35


    RESET = '\u001b[0m'
    if os.getenv('ANSI_COLORS_DISABLED') is None:
        fmt_str = '\u001b[%dm%s'
        if color is not None:
            text = fmt_str % (COLORS[color], text)

        if on_color is not None:
            text = fmt_str % (HIGHLIGHTS[on_color], text)

        if attrs is not None:
            for attr in attrs:
                text = fmt_str % (ATTRIBUTES[attr], text)

        text += RESET
    return text

def red(s):
    return colored(str(s), 'red', attrs=['bold'])
def green(s):
    return colored(str(s), 'green', attrs=['bold'])


def getlogger():
    import logging
    class CustomFormatter(logging.Formatter):

        colors = {
            'black': '\u001b[30;20m',
            'red': '\u001b[31;20m',
            'r': '\u001b[31;20m',
            'bold_red': '\u001b[31;1m',
            'green': '\u001b[32;20m',
            'g': '\u001b[32;20m',
            'gb': '\u001b[32;1m',
            'yellow': '\u001b[33;20m',
            'blue': '\u001b[34;20m',
            'b': '\u001b[34;20m',
            'purple': '\u001b[35;1m',
            'p': '\u001b[35;1m',
            'grey': '\u001b[38;20m',
        }
        FORMATS = {
            logging.WARNING: colors['purple'],
            logging.ERROR: colors['bold_red'],
            logging.CRITICAL: colors['bold_red'],
        }

        def format(self, record):
            format_str = "%(asctime)s  %(levelname)-6s %(funcName)-20s  line: %(lineno)-5s  %(message)s"
            reset = "\u001b[0m"
            log_fmt = None
            record.msg = str(record.msg)
            if '@' in record.msg[:10]:
                try:
                    icolor, tmp = record.msg.split('@', 1)
                    log_fmt = self.colors.get(icolor)
                    if log_fmt:
                        record.msg = tmp
                except:
                    raise
                    pass
            else:
                log_fmt = self.FORMATS.get(record.levelno)
            if log_fmt:
                record.msg = log_fmt + record.msg + reset
            formatter = logging.Formatter(format_str, datefmt='%Y-%m-%d %H:%M:%S')
            return formatter.format(record)

    prefix = 'udn_upload'
    fn_log = f'{prefix}.log'

    fmt = logging.Formatter('%(asctime)s  %(levelname)-6s %(funcName)-20s  line: %(lineno)-5s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')

    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(CustomFormatter())
    console.setLevel('INFO')

    fh_file = logging.FileHandler(fn_log, mode='w', encoding='utf8')

    fh_file.setLevel('DEBUG')
    fh_file.setFormatter(fmt)

    try:
        logger = logging.getLogger(__file__)
    except:
        logger = logging.getLogger('terminal')
    logger.setLevel('DEBUG')
    logger.addHandler(console)
    logger.addHandler(fh_file)

    return logger

def build_remote_subfolder(name, dest):
    # if demo:
    #     logger.info(f'demo mode, skip building remote subfolder')
    #     return 0
    if dest == 'dropbox':
        cmd = f'{dock} dbxcli mkdir "{name}/" >>{dest}.create_folder.log 2>{dest}.create_folder.log'
    elif dest == 'emedgene':
        # os.system(f'{dock} aws s3api put-object --profile {profile} --bucket emg-auto-samples --key Vanderbilt/upload/{name}/ >>{dest}.create_folder.log 2>{dest}.create_folder.log')
        cmd = f'{dock} aws s3api put-object --bucket emg-auto-samples --key "Vanderbilt/upload/{name}/" >>{dest}.create_folder.log 2>{dest}.create_folder.log'
    elif dest == 'rdrive':
        cmd = f'{dock_smb} smbclient "//i10file.vumc.org/ped/"   -A /home/chenh19/cred/smbclient.conf -c \'mkdir "{name}" \' '
        # logger.info(cmd)
    
    return cmd

def get_file_extension(fn):

    m = re.match(r'.*?\.(bam|bai|cnv|fasta|fastq|fq|gvcf|vcf|cram|crai)\b(\.gz)?', fn.lower())
    if m:
        try:
            return ft_convert[m.group(1)], m.group(2)
        except:
            return m.group(1), m.group(2)
    else:
        logger.warning(f'unclear file type: {fn}')
        return None, None

def parse_smb(fn, remote_pw):
    res = {}
    d_exist = {}
    sub_folders = set()
    #   .                                   D        0  Mon Feb  4 10:04:18 2019
    #   ..                                  D        0  Thu Apr 27 12:58:22 2023
    #   UDN525928_WB                        D        0  Fri Nov 16 20:08:26 2018
    #   UDN525928_FB                        D        0  Fri Nov 16 19:58:31 2018
    #   .DS_Store                           A    10244  Thu Apr 27 12:57:36 2023

    # \GEN\UDN\RNAseq Data\Baylor_Blood\UDN525928_02\UDN525928_WB
    #   .                                   D        0  Fri Nov 16 20:08:26 2018
    #   ..                                  D        0  Mon Feb  4 10:04:18 2019
    #   LTG037_S8_L001_R1_001.bam           A 7340450542  Fri Nov 16 20:07:03 2018
    #   LTG037_S8_L001_R1_001.bam.bai       A  3981720  Fri Nov 16 20:07:06 2018

    # \GEN\UDN\RNAseq Data\Baylor_Blood\UDN525928_02\UDN525928_FB
    #   .                                   D        0  Fri Nov 16 19:58:31 2018
    #   ..                                  D        0  Mon Feb  4 10:04:18 2019
    #   LTG038_S3_L001_R1_001.bam.bai       A  3463032  Fri Nov 16 19:52:12 2018
    #   LTG038_S3_L001_R1_001.bam           A 7569750225  Fri Nov 16 19:52:10 2018
    pw = ''
    with open(fn) as f:
        for i in f:
            i = re.sub(r'\s+', ' ', i.strip())
            a = i.rsplit(' ')
            
            if a[0] == '.':
                # add the root folder
                sub_folders.add('')
            if i.startswith('\\GEN'):
                pw = i.strip().replace('\\', '/') + '/'
                folder = pw.replace(remote_pw, '')
                folder = re.sub(r'^/', '', folder)
                folder = re.sub(r'/$', '', folder)
                sub_folders.add(folder)
                continue
            
            
            
            if len(a) != 8 or a[0] in {'.', '..'} or a[1] == 'D':
                continue
            fn, size = a[0], a[2]
            try:
                size = int(size)
            except:
                continue
            fn_full = f'{pw}{fn}'
            d_exist[fn] = [fn_full, size]
            

    return d_exist, sub_folders

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
    elif dest == 'rdrive':
        cmd = f'{dock_smb} smbclient "//i10file.vumc.org/ped/"   -A /home/chenh19/cred/smbclient.conf -D "{remote_pw}" <<< $\'recurse on\\nls\' 2>/dev/null > {fn_exist}'
        # logger.info(cmd)
        os.system(cmd)

    d_exist = {}
    sub_folders = set()
    if dest == 'dropbox':
        # 5d66a231cc9838df7c751@@719714@@1 year ago@@/VUDP005/Proband Raw Files/Files.zip
        # 5d66a231cc9848df7c751@@14467565985@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L001_R2_001_701E.fastq.gz
        # 5d66a231cc9858df7c751@@13632995037@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L002_R1_001_701E.fastq.gz
        # 5d66a231cc9868df7c751@@9126032@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_markdup.bam.bai
        # 5d66a231cc9878df7c751@@7313019@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag.vcf
        # 5d66a231cc9888df7c751@@13677703607@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L001_R1_001_701E.fastq.gz
        # 5d66a231cc9898df7c751@@13597294972@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L003_R1_001_701E.fastq.gz
        # 5d66a231cc98a8df7c751@@14469268993@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L002_R2_001_701E.fastq.gz
        # 5d66a231cc98b8df7c751@@14008431602@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L004_R1_001_701E.fastq.gz
        # 5d66a231cc98c8df7c751@@11716344@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag.bed
        # 5d66a231cc98d8df7c751@@14727930576@@1 year ago@@/VUDP005/Proband Raw Files/2021-145-142-PGnome-TrioDiag_S212_L004_R2_001_701E.fastq.gz
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
    elif dest == 'rdrive':
        d_exist, sub_folders = parse_smb(fn_exist, remote_pw)

    return d_exist, sub_folders

def parse_md5(fn):
    # 8ba94de338d39c9dac9ad4fd8b9cff83  2094656-UDN675219-M_L001_R1_001.fastq.gz
    # 82d147b41d0ba93decb90452a50a10c5  2094656-UDN675219-M_L002_R1_001.fastq.gz
    # ed04a252512aeaccfe9f6ebf402cbe53  2094656-UDN675219-M_L003_R1_001.fastq.gz
    res = {}
    with open(fn) as f:
        for i in f:
            i = re.split(r'\s+', i.strip())
            
            try:
                md5_v = i[0]
                if md5_v != 'NA':
                    res[i[1]] = i[0]
            except:
                pass
    
    return res


def get_local_md5(fn):
    try:
        tmp = fn.replace('/download/', '/log/')
        with open(f'{tmp}.md5') as f:
            for i in f:
                v = re.split(r'\s+', i.strip())[0]
                return v
    except:
        return None

def add_fn_suffix(fn, fn_suffix):
    m = re.match(r'(.*?)(\.\w+)(.gz)?$', fn)
    tmp = m.groups()
    # logger.info(tmp)
    fn_pure, ext, gz = tmp
    gz = gz or ''
    if fn_suffix and fn_suffix.strip():
        fn_suffix = re.sub('^\W+', '', fn_suffix)
        
        return ''.join([fn_pure, f'.{fn_suffix}', ext, gz])
    return fn
    


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

    # download.info.UDN133971.txt
    udnmain = info_file.rsplit('/', 1)[-1].replace('download.info.', '').replace('.txt', '')
    ifile_ct = {}
    file_ct[udnmain] = ifile_ct

    if not isinstance(gzip_only, set):
        gzip_only = set()

    ft = ft or ['fastq']
    # remote_base_pw = set()

    fn_all_folder = 'all_folders.txt'
    try:
        os.unlink(fn_all_folder)
    except:
        pass
    
    remote_pw_in = remote_pw_in.replace('\\', '/')
    logger.info(f'remote_pw_in={remote_pw_in}')
    remote_pw_in_old = remote_pw_in
    remote_pw_in = refine_remote_pw(remote_pw_in, dest)

    if remote_pw_in != remote_pw_in_old:
        logger.warning(f'r@remote_pw changed after refine, current={remote_pw_in}, prev = {remote_pw_in_old}')
    else:
        logger.info(f'g@remote_pw unchanged: {remote_pw_in}')
    try:
        d_exist, sub_folders = d_exist_all[remote_pw_in]
    except:
        d_exist, sub_folders = get_remote_file_list(pw, dest, remote_pw_in)
        d_exist_all[remote_pw_in] = (d_exist, sub_folders)

    fn_md5 = glob.glob(f'{pw}/download.*.md5')
    if len(fn_md5) == 0:
        logger.warning(f'expected md5 file not found')
        exp_md5 = {}
    else:
        if len(fn_md5) > 1:
            logger.warning(f'more than one md5 file found, only first would be uploaded: {fn_md5}')
        fn_md5 = fn_md5[0]
    
        exp_md5 = parse_md5(fn_md5)

    logger.info('start parsing')
    remote_sub_pw_map = {}
    
    url_na = []
    with open(info_file) as fp:
        header = fp.readline()
        a = [_.strip() for _ in header.split('\t')]
        idx = {}
        for _ in ['upload_type', 'download_type', 'remote_pw', 'rel_to_proband', 'fn', 'url', 'udn', 'size', 'rename', 'sample_list']:
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
            reltmp = rel or 'NA'


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

            upload_type = None
            try:
                tmp = a[idx['upload_type']].strip()
                if tmp:
                    upload_type = tmp
            except:
                pass
            
            if upload_type is None:
                upload_type = dest

            if upload_type not in {'emedgene', 'dropbox', 'rdrive'}:
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
                    remote_pw = refine_remote_pw(remote_pw, dest)
            except:
                if not remote_pw_in:
                    logger.error('must specify the remote_pw, in info file or by parameter, exit..')
                    return 1
                remote_pw = remote_pw_in

            # remote_base_pw.add(remote_pw)
            remote_pw = remote_pw.replace('\\', '/')

            try:
                d_exist, sub_folders = d_exist_all[remote_pw]
            except:
                d_exist, sub_folders = get_remote_file_list(pw, dest, remote_pw)
                d_exist_all[remote_pw] = (d_exist, sub_folders)
            # get the ext
            ext, gz = get_file_extension(fn)
            fn_orig = fn
            if fn_suffix:
                fn = add_fn_suffix(fn, fn_suffix)
                logger.info(fn_suffix)
            # logger.info([remote_pw, fn, ext])
            
            if 'all' not in ft and ext not in ft:
                logger.debug(f'file skipped due to file type not selected: {ext}  - {fn}')
                continue

            if ext in gzip_only and not gz:
                logger.info(f'file skipped due to gzip file only: {fn}')
                continue

            if url == 'NA':
                url_na.append(fn)
                continue
            

            

            ifile_ct.setdefault(reltmp, {}).setdefault(ext, []).append(fn)

            rel = re.sub(r'\W+', '_', rel)
            rel = '' if not rel or rel.lower() == 'na' else f'{rel}_'

            try:
                size_exp = int(size_exp)
            except:
                if size_exp.strip() and size_exp.lower() != 'na':
                    if  size_exp.find('archived_file') > -1:
                        logger.warning(f'\t archived: {fn}')
                    else:
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
                    else:
                        # if rename_suffix:
                        #     with open('/fs0/members/chenh19/tmp/upload/upload_for_phillips/batch2/rename_vcf.txt', 'a') as o:
                        #         print(f'mv {pw}/download/{rename}.vcf.gz {pw}/download/{rename}{rename_suffix}.vcf.gz', file=o)

                        rename = f'{rename}{rename_suffix}'
                sample_list = a[idx['sample_list']].strip()
                sample_list = re.sub(r'\W+', ' ', sample_list)

            except:
                logger.error('rename column not found in header')
                logger.info(f'idx = {idx}')
                raise
            if rename:
                logger.info(f'fn = {fn}, rename = {rename}, suffix={rename_suffix}')

            # if url.lower() == 'na':
            #     logger.warning(f'URL = NA : {fn}')

            fn_download = f'{pw}/download/{fn}'

            downloaded = 0
            file_status_by_fn[fn] = []
            if os.path.exists(fn_download):
                size_local = os.path.getsize(fn_download)
                fn_pure = os.path.basename(fn_download)
                file_md5_exp = exp_md5.get(fn_pure) or 'NA'
                file_md5_local = get_local_md5(fn_download)
                
                md5_ok = 0
                if file_md5_exp is None or file_md5_exp == 'NA':
                    logger.debug(f'expected md5 not found for {fn_pure}')
                    md5_ok = 1
                elif file_md5_local is None:
                    logger.debug(red(f'fail to get local md5 for {fn_download}'))
                    md5_ok = 1
                elif file_md5_exp == file_md5_local:
                    md5_ok = 1
                else:
                    logger.warning(red(f'md5 not match for {fn_download}, exp = {file_md5_exp}, local = {file_md5_local}'))
                    file_status_by_fn[fn].append('md5_not_match')
                
                if size_local == 0:
                    os.unlink(fn_download)
                    logger.warning('empty file: {fn_download}')
                    file_status_by_fn[fn].append('empty_file')
                elif size_exp == 'na' and md5_ok:
                    downloaded = 1
                    file_status_by_fn[fn].append('downloaded')
                elif size_local == size_exp and md5_ok:
                    downloaded = 1
                    file_status_by_fn[fn].append('downloaded')
                else:
                    logger.warning(red(f'file downloaded, but size not match. exp = {size_exp}, local={size_local} : {fn_download}'))
                    file_status_by_fn[fn].append('size_unmatch')
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
                    file_status_by_fn[fn].append('uploaded')
                else:
                    logger.warning(f'file uploaded, but size not match. exp = {size_exp}, remote={size_remote}')
                    uploaded = 0
                    file_status_by_fn[fn].append(f'remote_size_unmatch')

            if not os.path.exists(fn_download) and uploaded == 0:
                file_status_by_fn[fn].append('not_started')
            

            remote_sub_pw = f'{remote_pw}' if remote_flat else f'{remote_pw}/{name}'
            
            if fn in d_exist:
                remote_sub_pw_exist = d_exist[fn][0].replace(fn, '')
                remote_sub_pw_exist = re.sub(r'\/+$', '', remote_sub_pw_exist)
                remote_sub_pw_map[remote_sub_pw] = remote_sub_pw_exist
            
            remote_pw_final = remote_sub_pw_map.get(remote_sub_pw) or remote_sub_pw

            remote_pw_final = re.sub(r'/+$', '', remote_pw_final)

            if not force_download:
                downloaded = downloaded

            v = {'size': size_exp, 'download_type': download_type, 'upload_type': upload_type, 'url': url, 'remote': remote_pw_final, 'downloaded': downloaded, 'uploaded': uploaded, 'size_remote': f'{size_remote:,}', 'rename': rename, 'sample_list': sample_list, 'fn_orig': fn_orig}
            try:
                d[remote_pw][fn] = v
            except:
                d[remote_pw] = {fn: v}

    if len(url_na) > 0:
        logger.warning(f'totally {len(url_na)} desired files has invalid NA url:\n\t')
        print('\n\t'.join(url_na))


    # logger.debug(d)
    sub_folders_exist = set()
        
    for remote_pw, v in d_exist_all.items():
        i = set([f'{remote_pw}/{name}' for name in v[1]])
        if len(v[0]) > 0 or len(v[1]) > 0:
            sub_folders_exist.add(remote_pw)
        sub_folders_exist |= i

    sub_folders_exp = set()
    for k1, v1 in d.items():
        sub_folders_exp.add(k1)
        for v2 in v1.values():
            sub_folders_exp.add(v2['remote'])

    logger.info(f'subfolder exp = {sub_folders_exp}')
    logger.info(f'subfolder existed = {sub_folders_exist}')

    sub_folder_to_build = sorted(sub_folders_exp - sub_folders_exist)
    if not demo and not no_upload:
        for i in sub_folder_to_build:
            build_remote_subfolder(i, dest=upload_type)
    if len(sub_folder_to_build) > 0:
        logger.info(f'demo build remote subfolder {sub_folder_to_build}')
    else:
        logger.info(f'g@all remote subfolder are ready')


    # check the status
    n_downloaded = 0
    n_uploaded = 0
    n_desired = 0
    need_upload = []
    
    # ok = uploaded

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
        
        
    
    file_status = {'ok': [], 'downloaded': [], 'raw': [], 'invalid_url': []}

    for _ in list(d):
        v1 = d[_]
        for fn in list(v1):
            if fn in excluded_file:
                del d[_][fn]
                continue
            v = v1[fn]

            fn_refine = re.sub(r'\s+', '_', fn)
            fn_script = f'{pw}/shell/{fn_refine}.download_upload.{dest}.sh'
            n_desired += 1

            if v['uploaded']:
                logger.info(colored(f'file already uploaded: {fn}', 'green'))
                n_uploaded += 1
                file_status['ok'].append(fn)

                if not force_download:
                    # n_downloaded += 1
                    continue
                
                
                
            if v['url'].lower() == 'na':
                file_status['invalid_url'].append(fn)
                continue
            
            if not v['downloaded']:
                if not v['uploaded']:
                    if not no_upload:
                        need_upload.append(fn)
                    script_list['needdown'].append(fn_script)
                    file_status['raw'].append(fn)
            else:
                n_downloaded += 1
                file_status['downloaded'].append(fn)
                if not no_upload:
                    need_upload.append(fn)
                    script_list['needup'].append(fn_script)


    need_upload = sorted(need_upload)
    n_need_upload = len(need_upload)
    logger.info(f'total files = {n_desired}, need_to_upload={n_need_upload}, already exist in server = {n_uploaded}, already downloaded = {n_downloaded}')

    invalid_url = file_status['invalid_url']
    n_invalid_url = len(invalid_url)
    n_need_download = n_desired - n_downloaded - n_uploaded
    if n_need_upload > 0 or n_invalid_url > 0 or n_need_download > 0:
        print('\n', '!' * 50)
        print(green(f'{n_downloaded}/{n_desired} files already downloaded'))
        if n_need_download > 0:
            print(red(f'{n_need_download}/{n_desired} files need to be downloaded'))
        if n_need_upload > 0:
            print(red(f'{n_need_upload}/{n_desired} files need to be uploaded'))
        
        str_tmp = [f'{_}  {red(file_status_by_fn.get(_))}' for _ in need_upload]
        
        print('\t' + '\n\t'.join(str_tmp))
        if n_invalid_url > 0:
            logger.info(f'\nERROR: {n_invalid_url} files with NA url')
            logger.debug('\n\t'.join(invalid_url))
    else:
        print('\n', '*' * 50)
        print(colored('Done, all files already uploaded', 'green', attrs=['bold']))
    print('\n\n\n')

    ct = {'total': n_desired, 'downloaded': n_downloaded, 'uploaded': n_uploaded, 'need_to_upload': n_need_upload}

    return d, ct, script_list, file_status


def refine_remote_pw(remote_pw, dest):
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
    if dest == "emedgene":
        os.system(f'{dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_base_dir}" > {fn_all_folder} ')
    elif dest == 'rdrive':
        cmd = f'{dock_smb} smbclient "//i10file.vumc.org/ped/"  -A "/home/chenh19/cred/smbclient.conf" -D "{remote_base_dir}" -c "ls" 2>/dev/null > {fn_all_folder}'
        os.system(cmd)

    all_remote_folders = {}
    pattern = re.compile(r'.*(UDN\d+)(.*|$)')
    n = 0

    skip_words = ['Trio']  # if the path name contains these words, would not be considered as existed path

    for i in open(fn_all_folder):
        if dest == 'rdrive':
            # UDN525928_02     D        0  Mon Feb  4 10:04:18 2019
            i = re.split(r'\s+D\s+', i.strip())[0].strip()
        else:
            # PRE UDN616750_AA/
            i = i.strip().rsplit(' ', 1)[-1].strip()
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
        logger.warning(f'base folder not exist: {remote_base_dir}')
        # sys.exit(1)

    
    remote_pw = re.sub('^/', '', remote_pw)
    remote_pw = re.sub('upload_?', '', remote_pw, flags=re.I)
    # remote_pw = re.sub(r'^\d+_', '', remote_pw, 1)
    m = re.match(r'(.*)?\s*(UDN\d+)\s*(.*)?', remote_pw)

    
    def extract_str(s):
        if s is None or not s.strip():
            return None
        
        s = re.sub('^_*', '', s.strip())
        s = re.sub('_*$', '', s)
        
        return s
    
    if m:
        m = m.groups()
        udn_tmp = m[1]
        if udn_tmp in all_remote_folders:
            remote_pw = all_remote_folders[udn_tmp]
            logger.info(f'g@remote folder already exist: {remote_pw}')
            return remote_pw

        p1 = extract_str(m[0])
        p2 = extract_str(m[2])
        
        if p1:
            m1 = re.match('(\d+)_(\d+.*)', p1)
            if m1:
                p1 = m1.group(2)


        p_trailing = '_'.join([_ for _ in (p1, p2) if _])
        p_trailing = '_' + p_trailing if p_trailing else ''
        remote_pw = m[1] + p_trailing
    
    return f'{remote_base_dir}{remote_pw}'

def get_fn_remote(fn_local, dir_local_files, local_pw_layers=0, fn_remote=None):
    
    fn_remote = fn_remote or fn_local
    fn_remote = os.path.basename(fn_remote)
    
    dir_local_files = os.path.abspath(dir_local_files)
    fn_local = os.path.abspath(fn_local)
    # logger.info(fn_remote)

    if local_pw_layers:
        tmp = fn_local.replace(dir_local_files, '').split('/')[:-1]
        try:
            pw_remote_addition = '/'.join(tmp[-local_pw_layers:])
        except:
            pass
        else:
            if not fn_remote or not fn_remote.startswith(pw_remote_addition):
                fn_remote = f'{pw_remote_addition}/{fn_remote}'
        # logger.info(pw_remote_addition)

    return fn_remote

def build_script_single(dest, remote_pw, fn_local, url_var=None, simple=False, fn_remote=None, fn_deid=None, sample_list=None, pw=None, download_type=None, size_exp=None, need_upload=True, local_file_size=None):
    # url_var is for the normal upload script with download step
    # the format is like  f"""url=$(awk -F "\\t" '$2=="{fn}" {{print ${idx_url} }}' {info_file})"""
    # if fn_deid is not None, then sample_list must be specified
    # fn_remote, if not specified, use the same as the local file

    fn = os.path.basename(fn_local)
    pw = pw or os.getcwd()
        
    local_file_size = local_file_size or {}

    fn_remote = fn_remote or fn
    # if not os.path.exists(fn_local):
    #     return None

    fn_refine = re.sub(r'\s+', '_', fn)
    fn_script = f'{pw}/shell/{fn_refine}.download_upload.{dest}.sh'
    fn_status = f'{pw}/log/status.{fn_refine}.txt'

    if simple == False and (url_var is None):
        logger.error(f'invalid url_var define: {url_var}, file= {fn_local}')
        return None
    cmd = []
    if url_var is not None:
        if url_var[:4] != 'url=':
            url_var = f'url="{url_var}"'

        cmd.append(url_var)


    if fn_deid is not None and sample_list is None:
        logger.error(f'when fn_deid is specified, the sample_list must also be specified')
        return None

    if simple:
        size_exp = local_file_size.get(fn_local) or os.path.getsize(fn_local)
        cmd_check_md5 = ''

    else:
        size_exp = size_exp or 'na'
        cmd_check_md5 = f"""# check md5sum

    md5_exp=$(grep -E -m1 "{fn}($|[^.])" {pw}/download.UDN*.md5 2>/dev/null|cut -d " " -f1)
    
    if [[ ! -f {pw}/log/{fn}.md5 && -n "$md5_exp" && "$md5_exp" != "NA" ]];then
    md5sum  "{fn_local}" > "{pw}/log/{fn}.md5" 2> "{pw}/log/{fn}.md5.log"
    fi

    md5_local=$(head -1 {pw}/log/{fn}.md5|cut -d " " -f1)

    if [[ -z $md5_exp || "$md5_exp" == "NA" ]];then
    echo "expected md5 not found : {fn}" >> "{fn_status}"
    echo 0
    elif [[ "$md5_exp" == "$md5_local" ]];then
    echo good, md5 match!  >> "{fn_status}"
    echo 0
    else
    echo 1
    echo "ERROR  md5 not match! {fn_local}; exp=@${{md5_exp}}@  local=@${{md5_local}}@"  >> "{fn_status}"

    echo $(date): "ERROR, md5 not match! {fn_local}; exp=@${{md5_exp}}@  local=@${{md5_local}}@"  >> /fs0/members/chenh19/tmp/upload/error.md5_not_match.files.txt
    fi
    """

    with open(fn_script, 'w') as out:
        rm_cmd = ''
        if rm:
            rm_cmd = f'    rm {fn_local}\n    rm {pw}/log/upload.{dest}.{fn}.log\n    rm {pw}/log/download.{fn}.log'
        if simple or download_type is None:
            download_cmd = ''
        if download_type == 'dropbox':
            download_cmd = f'dbxcli get "$url"  "{fn_local}" > {pw}/log/download.{fn}.log 2>&1'
        elif download_type == 'wget':
            download_cmd = f'wget -c -O "{fn_local}" "$url" > {pw}/log/download.{fn}.log 2>&1'



        cmd.append(f"""
    checklocal(){{
        local_size=$(stat -c "%s" {fn_local} 2>/dev/null)
        if [[ $local_size -eq 0 ]];then
        rm {fn_local} 2>/dev/null
        echo 1
        return
        fi

        # check local file
        if [[ "{size_exp}" -ne "na" ]] & [[ "$local_size" -ne "{size_exp}" ]];then
        echo local file not match with exp local_size=$local_size,  expected={size_exp}. skip uploading >> "{fn_status}"
        echo 1
        return
        fi
        {cmd_check_md5}
    }}
    """)
        # check remote
        fn_pure = os.path.basename(fn_remote)
        check_remote_logic = f"""
        local_size=$(stat -c "%s" "{fn_local}" 2>/dev/null)
        fnlog="{fn_status}"

        size_exp={size_exp}
        echo local_size=$local_size , remote size = $remote_file_size exp size = $size_exp >> "$fnlog"

        if [[ -z $remote_file_size ]];then
            # local file not exist
            if [[ -z $local_size ]];then
                echo "not downloaded yet" >> "$fnlog"
                echo "local_not_ready"
                return
            fi

            # exp size unkown, upload anyway
            if [[ $size_exp == "na" ]];then
                echo "size_exp not specified, upload anyway : {fn_local}" >>  "$fnlog"
                echo upload
                return
            fi

            if [[ $local_size -eq $size_exp ]];then
                echo "local size is correct, prepare to upload: {fn_local}" >>  "$fnlog"
                echo upload
                return
            else
                echo "local size not correct: local=$local_size , exp = $size_exp" >>  "$fnlog"
                echo "local_not_ready"
                return
            fi
            echo -e "unkown situation! \\n local=$local_size\\n remote=$remote_file_size \\n exp = $size_exp" >>  "$fnlog"

            exit 1
        elif [[ $size_exp == "na" ]];then
            if [[ "$local_size" -eq "$remote_file_size" ]];then
                echo size check is not specified, remote_file_size=$remote_file_size local_filesize=$local_size >>  "$fnlog"
                echo already_uploaded
                return
            else
                echo remote size not correct: remote = $remote_file_size, local = $local_size >>  "$fnlog"
                echo upload
                return
            fi
        elif [[ "$remote_file_size" -eq "{size_exp}" ]];then
                echo {fn} successfully uploaded to {remote_pw}/  filesize={size_exp} >>  "$fnlog"
                echo already_uploaded
                return
        elif [[ $size_exp -eq $local_size ]];then
            echo remote size not correct: remote = $remote_file_size, exp = {size_exp} >>  "$fnlog"
                echo upload
                return
        elif [[ $size_exp -ne $local_size ]];then
            echo {fn} local size not match exp, exp=$size_exp, local=$local_size >>  "$fnlog"
            echo "local_not_ready"
            return
        else
            echo -e "unkown situation! \\n local=$local_size\\n remote=$remote_file_size \\n exp = $size_exp" >>  "$fnlog"
            exit 1
        fi
        """
        
        
        if dest == 'emedgene':
            cmd.append(f"""
    checkremote(){{
        remote_file_size=$({dock} aws s3 ls "emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn_remote}"|grep -E -m1 "{fn_pure}$" | tr -s " "|cut -d " " -f 3)
        {check_remote_logic}

    }}
    """)
        elif dest == 'dropbox':
            # there is no way for dropbox to get byte file size
            cmd.append(f"""
    checkremote(){{
        file_exist_from_remote=$({dock} dbxcli ls -l "/{remote_pw}/"|grep -E -m1 "{fn_remote}($|[^.])")
        if [[ ! -z "$file_exist_from_remote" ]];then
            echo alrady_uploaded
        else
            echo upload
        fi
    }}

    """)
        elif dest == 'rdrive':
            tmp = fn_remote.rsplit('/', 1)
            fn_remote_pure = tmp[-1]
            if len(tmp) == 1:
                sub_pw = ''
            else:
                sub_pw = tmp[0]
            
            cmd.append(f"""
    checkremote(){{
        remote_file_size=$({dock_smb} smbclient "//i10file.vumc.org/ped/"   -A "/home/chenh19/cred/smbclient.conf" -D "{remote_pw}/" <<< $'recurse on\\nls "{sub_pw}" '|grep -E -m1 "{fn_remote_pure}($|[^.])"|head -1|tr -s " "|cut -d ' ' -f 4)
        {check_remote_logic}
    }}
            """)
        fn_refine = re.sub(r'\s+', '_', fn)
        cmd.append(f"""
    postupload(){{
        mv {pw}/shell/{fn_refine}.download_upload.{dest}.sh \\
            {pw}/shell_done/{fn_refine}.download_upload.{dest}.sh
        fntmp={pw}/log/{fn_refine}.tmplog

        tail -c 3000 {pw}/log/upload.{dest}.{fn_refine}.log  2>/dev/null|sed 's/\\r/\\n/g' > $fntmp;
        mv $fntmp {pw}/log/upload.{dest}.{fn_refine}.log
        tail -n 30 {pw}/log/download.{fn_refine}.log > $fntmp 2>/dev/null;
        mv $fntmp {pw}/log/download.{fn_refine}.log

        {rm_cmd}

        exit 0
    }}

    download_status=$(checklocal)

    """)
    
        # upload_status
        if not no_upload and not force_download:
            cmd.append(f"""
    upload_status=$(checkremote)
    if [[ $upload_status = "already_uploaded" ]];then
        postupload
        exit 0
    fi
    """)
        else:
            logger.info(f'g@{no_upload=}, {force_download=}')
        cmd.append(f"""

    if [[ $download_status -eq 1 ]];then
        {download_cmd}
        download_status=$(checklocal)
        if [[ $download_status -eq 1 ]];then
            echo download failed >> {fn_status}
            exit 1
        fi
    fi
    """)
        # rename the vcf and bgzip
        if fn_deid:
            fn = f'{fn_deid}.vcf.gz'
            # prev command is
            # bcftools view --no-version {fn_download}|bcftools reheader -s  <(echo "{rename}") |bgzip > {pw}/download/{fn}.tmp
            # mv {pw}/download/{fn}.tmp {pw}/download/{fn}
            # sample list example: 971146-  971147-

            cmd.append(f"""
    # rename the fq file
    vcf_deid {fn_local} -newname {fn_deid} -sample_list {sample_list}|bgzip > {pw}/download/{fn}
    bgzip -r {pw}/download/{fn}
    ln -sf {fn_local} {pw}/download/{fn_deid}.link

            """)

        if dest == 'dropbox':
            upload_cmd = f'{dock} dbxcli put "{fn_local}" "/{remote_pw}/{fn_remote}" > {pw}/log/upload.{dest}.{fn}.log 2>&1'
        elif dest == 'emedgene':
            upload_cmd = f'{dock} aws s3 cp "{fn_local}" "s3://emg-auto-samples/Vanderbilt/upload/{remote_pw}/{fn_remote}" > {pw}/log/upload.{dest}.{fn}.log 2>&1\ndate +"%m-%d  %T">> {pw}/log/upload.{dest}.{fn}.log'
        elif dest == 'rdrive':
            # upload_cmd = f"""{dock} smbclient  --socket-options='TCP_NODELAY IPTOS_LOWDELAY SO_KEEPALIVE SO_RCVBUF=16777216 SO_SNDBUF=16777216 SO_RCVTIMEO=120000 SO_SNDTIMEO=120000' "//i10file.vumc.org/ped/"   -A "/home/chenh19/cred/smbclient.conf"  <<< $'rm "{remote_pw}/{fn_remote}"\ntimeout 120\niosize 16384\nput "{fn_local}" "{remote_pw}/{fn_remote}" ' """
            
            fn_pure, ext = fn_remote.rsplit('.', 1)
            upload_cmd = f"""{dock_smb} smbclient "//i10file.vumc.org/ped/"   -A "/home/chenh19/cred/smbclient.conf"  <<< $'rm "{remote_pw}/{fn_remote}"\nput "{fn_local}" "{remote_pw}/{fn_remote}" ' """

        if not no_upload and need_upload:
            cmd.append(f"""
    upload_status=$(checkremote)

    if [[ $upload_status = "local_not_ready" ]];then
    echo local file note ready >> {fn_status}
    exit 1
    fi

    # local size is ok, try upload
    echo start uploading >> {fn_status}
    date +"%m-%d  %T">> {fn_status}

    {upload_cmd}
    echo upload finish >> {fn_status}
    date +"%m-%d  %T"  >> {fn_status}
    upload_status=$(checkremote)
    if [[ $upload_status = "already_uploaded" ]];then
    postupload
    echo "uploaded"
    exit 0
    else
    echo "fail to upload"
    exit 1
    fi
        """)
        print('\n'.join(cmd), file=out)
    os.chmod(fn_script, 0o755)
    return fn_script

def build_script(pw, d, info_file, no_upload=False):
    with open(info_file) as f:
        info_header = f.readline().strip().split('\t')
        idx_url = info_header.index('url') + 1

    for _, v1 in d.items():
        for fn, v in v1.items():
            # {'size': size_exp, 'remote': f'{remote_pw}/{name}', 'url': url, 'downloaded': downloaded, 'uploaded': uploaded}
            url = v['url']
            fn_orig = v['fn_orig']
            size_exp = v['size']
            remote_pw = v['remote']
            download_type = v['download_type']
            # ext = fn.rsplit('.', 1)[-1]
            fn_deid = v['rename']
            sample_list = v['sample_list']

            fn_download = f'{pw}/download/{fn}'
            # fn_download_chunk = fn_download.replace('.gz', '')

            need_upload = not v['uploaded']
            fn_refine = re.sub(r'\s+', '_', fn)
            fn_script = f'{pw}/shell/{fn_refine}.download_upload.{dest}.sh'
            if v['uploaded'] and not  args.force and not force_download:
                os.system(f'mv  {fn_script} {pw}/shell_done/{fn_refine}.download_upload.{dest}.sh 2>/dev/null')
                continue

            if url.lower() == 'na':
                continue
            url_var = f"""url=$(awk -F "\\t" '$2=="{fn_orig}" {{print ${idx_url} }}' {info_file})"""
            fn_script = build_script_single(dest, remote_pw, fn_download, url_var, simple=False, fn_deid=fn_deid, sample_list=sample_list, pw=pw, download_type=download_type, size_exp=size_exp, need_upload=need_upload)
            if fn_script is None:
                # logger.info(f'{dest=}, {remote_pw=}, {fn_download=}, {fn_deid=}, {pw=}, {download_type=}, {size_exp=}, {need_upload=}')
                logger.warning(f'fail to build script for {fn}')



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
    d, ct, script_list, file_status = parse_info_file(pw, info_file, remote_pw_in=remote_pw_in, ft=ft, updated_version_only=updated_version_only, script_list=script_list, no_upload=no_upload, remote_flat=remote_flat)
    if not lite:
        build_script(pw, d, info_file, no_upload=no_upload)
    return ct, script_list, file_status

def clearlog(force=False, args=None):
    # get folder list
    pw_list = os.popen(f'ls -d *UDN*/ vudp*/ VUDP*/ 2>/dev/null|sort|uniq').read().strip().split('\n')
    pw_list = [_[:-1] for _ in pw_list if _.strip()]
    
    if args.pw:
        pw_list = sorted(set(pw_list) & set(args.pw))
    
    n_pw = len(pw_list)
    
    logger.info(pw_list)
    
    if len(pw_list) == 0:
        info_file = glob.glob(f'UDN*.yaml')
        if len(info_file) > 0:
            pw_list = [os.path.abspath('.')]
        else:
            logger.info(f'gb@no UDN folder found under current folder')
            return 0
    logger.info(f'gb@there are {len(pw_list)} UDN folders found')
    log_fls = []

    for ipw in pw_list:
        fls = os.popen(f'find {ipw}/log -iname "*.log" -type f -size +20k 2>/dev/null').read().strip().split('\n')
        log_fls += [_.strip() for _ in fls if _.strip()]

    if len(log_fls) > 0:
        logger.warning(f'now trying to clear the logs, n = {len(log_fls)}')
    upload_fls = []
    download_fls = []

    for fn in log_fls:
        if fn.find('log/upload') > 0:
            upload_fls.append(fn)
        elif fn.find('log/download') > 0:
            download_fls.append(fn)

    cmd = ""
    if len(download_fls) > 0:
        download_fls = '\n'.join(download_fls)
        cmd += f"""
    while read fn;do
        tail -n 30 $fn  2>/dev/null> ${{fn}}.tmp;
        mv ${{fn}}.tmp $fn
    done <<EOF_
    {download_fls}\nEOF_
    """
    if len(upload_fls) > 0:
        upload_fls = '\n'.join(upload_fls)
        cmd += f"""
    while read fn;do
        tail -c 2000 $fn 2>/dev/null|sed 's/\\r/\\n/g' > ${{fn}}.tmp;
        mv ${{fn}}.tmp $fn
    done << EOF_
    {upload_fls}\nEOF_
    """
    

    if cmd:
        logger.warning(f'upload files = {len(upload_fls)}, download files = {len(download_fls)}')
        os.system(cmd)
    else:
        logger.info(f'gb@no log files need to be cleared')
        
    # deal with the download folder
    for ipw in pw_list:
        pw_download = f'{ipw}/download'
        if not os.path.exists(pw_download):
            logger.info(f'{ipw}, already removed')
            continue
        if not force:
            proceed = input(f'are you sure you need to remove folder\n  {pw_download}?  (y/n): ')
            if not proceed.lower().strip().startswith('y'):
                logger.warning(f'\tskip...')
                continue
        logger.info(f'gb@now removing {pw_download}')
        fn_flist = f'{ipw}/download.flist.txt'
        os.system(f'echo -e "\n\n$(date)\\n" >> {fn_flist}; ls {pw_download} -lha >> {fn_flist}; echo -e "***********\\n\\n" >> {fn_flist}; rm -rf {pw_download}')


def parse_ls_dump(fn):
    res = {}
    with open(fn) as f:
        f.readline()
        for i in f:
            a = re.split(r'\s+', i.strip())
            # -rw-r--r-- 1 chenh19 h_vangard_1    11716344 Jan 30 09:24 2022-344-011-PGnome-PatientOnlyDiag.bed
            try:
                size, ifl = a[4], a[8]
            except:
                logger.info(f'invalid line: {a}')
                sys.exit(1)
            if ifl in ['.', '..']:
                continue
            try:
                size = int(size)
            except:
                logger.info(f'invalid line in ls file list: {i[:-1]}')
                continue
            res[ifl] = size
    if len(res) == 0:
        return None, None
    
    return sorted(res), res


def get_local_files(args, ft):
    fls_raw = args.local_files
    if len(fls_raw) == 0:
        fls_raw = [os.path.abspath('.')]
    select_file = 0
    ls_dump = args.info_file
    file_size = {}
    if ls_dump:
        if not os.path.isfile(ls_dump):
            logger.error(f'file list not exist: {ls_dump}')
            sys.exit(1)
        flist, file_size = parse_ls_dump(ls_dump)
        if flist is None:
            logger.info(f'invalid file list file: {ls_dump}')
        select_file = 1
    elif len(fls_raw) == 1 and os.path.isdir(fls_raw[0]):
        # this is a folder
        pwtmp = os.path.abspath(fls_raw[0])
        logger.info(pwtmp)
        # -maxdepth 1
        flist = os.popen(f'find -L {pwtmp} -type f').read().split('\n')
        select_file = 1
    elif len(fls_raw) == 1 and not os.path.exists(fls_raw[0]):
        logger.error(f'file not exist: {fls_raw[0]}')
        sys.exit(1)
    else:
        logger.info(f'you directly specified the filelist')
        logger.info(fls_raw)
        flist = fls_raw
    
    if 'all' in ft:
        select_file = 0
    
    res = []
    for fn in flist:
        if not fn or fn.endswith('.log'):
            continue
        
        fn_refine = re.sub(r'\s+', '_', fn)
        if fn_refine != fn:
            logger.info(f'now renaming {fn}')
            os.system(f'mv "{fn}" "{fn_refine}" ')
            fn = fn_refine
        
        ext = fn.replace('.gz', '').rsplit('.', 1)[-1]
        if (select_file and ext in ft) or not select_file:
            res.append(fn)
        
    return res, file_size



if __name__ == "__main__":

    logger = getlogger()

    simple_mode = args.simple

    remote_flat = args.remote_flat
    lite = args.lite
    if lite:
        demo=True
    skip_download = args.skip_download
    force_upload = args.force_upload
    force_download = args.forcedown
    force_clear = args.force_clear
    no_upload = args.noupload
    profile = args.profile
    clearlog_flag = args.clear
    fn_suffix = args.suffix
    demo = args.demo
    ignore_na_url = args.i


    if fn_suffix:
        fn_suffix = re.sub(r'^\W+', '', fn_suffix)
        fn_suffix = re.sub(r'\W+$', '', fn_suffix)
    

    asis = args.asis  # use the remote_pw in the cmd line input, do not do the re-format
    rm = args.rm
    remote_pw = args.remote_pw
    remote_base = args.remote_base
    # print(args)

    node_name = platform.node()
    updated_version_only = not args.allversion

    if node_name.find('viccbiostat120') > -1:
        dock = 'singularity exec -B /mnt /mnt/d/dock/centos.sif '
        tmp = os.popen(f'which smbclient').read().strip()
        if tmp:
            dock_smb = ''
        else:
            dock_smb = dock
            
    elif node_name.find('vampire') > -1:
        if 'hanuman' in node_name:
            bind = '-B /fs0'
        else:
            bind = ''
        dock = f'singularity exec {bind} /data/cqs/chenh19/dock/centos.sif '
        dock_smb = dock


    convert1 = {'dropbox': 'dropbox', 'db': 'dropbox', 'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene', 'rdrive': 'rdrive', 'r': 'rdrive'}
    convert2 = {'emedgene': 'emedgene', 'em': 'emedgene', 'ed':'emedgene', 'pacbio': 'pacbio', 'pac': 'pacbio'}
    dest = convert1[args.dest]
    profile = convert2[profile]


    # get the file type to upload
    gzip_only = set()
    ft = set()
    
    file_status_by_fn = {}
    

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


    if simple_mode:
        remote_pw = args.remote_pw
        local_pw_layers = args.local_pw_layers  # keep the local parent folder
        # e.g. /a/b/c/d/1.fastq.gz   if layer = 0, will include no parent folder
        # if layer = 1, will upload to remote_pw/d/1.fastq.gz
        dir_local_files = args.local_files
        if not dir_local_files:
            logger.error(f'simple mode, must specify the root folder for the files need to be uploaded')
            sys.exit(1)
        
        
        
        if len(dir_local_files) > 1:
            logger.warning(f'multiple folders specified, will use the first')
        dir_local_files = dir_local_files[0]
        if not os.path.isdir(dir_local_files):
            logger.error('simple mode, first arg should be a folder, : {dir_local_files}')
            sys.exit(1)
        
        local_files, local_file_size = get_local_files(args, ft)
        logger.info(f'local files n = {len(local_files)}')
        # logger.info('\n' + '\n'.join(local_files)+ '\n\n')
        
        pw = os.getcwd()
        os.makedirs('log', exist_ok=True)
        os.makedirs('shell', exist_ok=True)
        os.makedirs('shell_done', exist_ok=True)

        local_files_ft_ct = {}
        for fn in local_files:
            ext = fn.replace('.gz', '').rsplit('.', 1)[-1]
            local_files_ft_ct.setdefault(ext, 0)
            local_files_ft_ct[ext] += 1

        if not remote_pw:
            logger.error('you must specify the remote path for uploading')
            sys.exit(1)
        
        if not demo:
            confirm_dest = input(f'You are uploading {len(local_files)} files: {local_files_ft_ct} to {dest}, is that correct? (yes/no):\n').lower()

            if confirm_dest not in {'yes', 'y'}:
                logger.error(f'you denied, please specify the correct -dest argument')
                sys.exit(1)
        
        remote_pw = re.sub(r'\\', '/', remote_pw)
        d_exist, sub_folders = get_remote_file_list(pw, dest, remote_pw)
        file_status = {'uploaded': [], 'need_to_upload': [], 'size_not_match': []}
        
        sub_folders = {f'{remote_pw}/{_}' if _ else remote_pw for _ in sub_folders}

        # logger.info('\n\t' + '\n\t'.join(sub_folders))
        
        
        # logger.info(d_exist)
        for fn in local_files:
            fn_pure = os.path.basename(fn)
            if fn_suffix:
                fn_pure = add_fn_suffix(fn_pure, fn_suffix)
            
            if fn_pure not in d_exist:
                file_status['need_to_upload'].append(fn)
            else:
                size = local_file_size.get(fn) or os.path.getsize(fn)
                # logger.info(f'{fn}: {size}')
                size_remote = d_exist[fn_pure][1]
                if size_remote == size:
                    file_status['uploaded'].append(fn)
                else:
                    file_status['need_to_upload'].append(fn)
                    file_status['size_not_match'].append(f'{fn}\tremote={size_remote}, local={size}, ratio remote/local = {size_remote/size:.3f}')
        size_not_match = file_status['size_not_match']
        fn_not_match = 'remote_size_not_match.txt'
        if len(size_not_match) > 0:
            logger.error(colored(f'{len(size_not_match)} files uploaded, but the size not match, please check: {fn_not_match}', 'red'))
            with open(fn_not_match, 'w') as o:
                print('\n'.join(size_not_match), file=o)
            print('\n'.join(size_not_match))
        else:
            try:
                os.unlink(fn_not_match)
            except:
                pass
        ntmp = len(file_status['uploaded'])
        if ntmp > 0:
            logger.info(f'g@{ntmp} files has been uploaded')
        #     # logger.info(colored(f'\nthe following files has been uploaded:\n\t' + '\n\t'.join(file_status['uploaded']), 'green'))
        #     logger.info(colored(f'\nthe following files has been uploaded:\n\t' + '\n\t'.join(file_status['uploaded']), 'green'))

        if len(file_status['need_to_upload']) == 0:
            print(colored('All done', 'green'))
            sys.exit(0)
        else:
            logger.info(f'now building scripts for {len(file_status["need_to_upload"])} files')


        scripts = []
        
        cmd_build_remote_pw = []

        remote_pw_list = [remote_pw]

        fls = sorted(file_status['need_to_upload'])
        purenames = [os.path.basename(_) for _ in fls]
        max_len = max([len(_) for _ in purenames])
        tmp = zip(fls, purenames)
        for fn, fn_pure in tmp:
            
            fn_remote = add_fn_suffix(fn, fn_suffix)
            fn_remote = get_fn_remote(fn, dir_local_files, local_pw_layers, fn_remote=fn_remote)
            # UDN620340_287/UDN763389_Father/9769-LR-0005_S1_L005_R1_001.fastq.gz 
            # if depth = 0, will be plain file name
            
            tmp = fn_remote.split('/')[:-1]
            if len(tmp) > 0:
                for i in range(len(tmp)):
                    tmp1 = '/'.join(tmp[:i+1])
                    tmp2 = f'{remote_pw}/{tmp1}'
                    if tmp2 not in remote_pw_list:
                        remote_pw_list.append(tmp2)
            
            fn_script = build_script_single(dest, remote_pw, fn, simple=True, need_upload=True, fn_remote=fn_remote, local_file_size=local_file_size)
            tmp = f'{remote_pw}/{fn_remote}'
            tmp = tmp.replace('/', red(' / '))
            print(fn_pure.ljust(max_len + 5) + f'remote:  {tmp}')
            
            if fn_script is None:
                logger.warning(f'fail to build script for {fn_pure}')
                continue
            scripts.append(fn_script)
        
        # logger.info(remote_pw_list)
        # logger.info(sub_folders)
        
        for pwtmp in remote_pw_list:
            if pwtmp  in sub_folders:
                continue
            cmd = build_remote_subfolder(pwtmp, dest)
            # if not demo:
            #     logger.info(f'building remote pw: {pwtmp}')
            #     os.system(cmd)
            # else:
            cmd_build_remote_pw.append(cmd)
        
        
        fn_script = f'shell/build_pw.sh'
        if len(cmd_build_remote_pw) > 0:
            with open(fn_script, 'w') as o:
                print('\n'.join(cmd_build_remote_pw), file=o)
            if not nobuild:
                logger.info(f'need to build remote folder: {fn_script}')
        else:
            try:
                os.unlink(fn_script)
            except:
                pass
            logger.info(f'g@all remote subfolders already exist: {remote_pw_list}')

        
        
        with open(f'script_list.txt', 'w') as o:
            print('\n'.join(scripts), file=o)
        sys.exit(0)

    if clearlog_flag:
        sys.exit(clearlog(force_clear, args))

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

    if args.pw:
        logger.info(f'input pw = {args.pw}')
        logger.info(f'parsed pw = {pw}')

    if len(pw) > 1 and remote_pw:
        logger.error('multiple local path found, but only 1 remote_pw specified, try specify one at each or specify remote_base, then use auto assigned name')
        sys.exit(1)

    logger.info(f'total cases = {len(pw)}\n{"*"*50}\n\n')
    pw = sorted(pw)

    # if len(pw) > 1 and not remote_base:
    #     logger.error(f'multiple local path found, you must specify a remote_base')
    #     sys.exit(1)

    ct = {'total': 0, 'downloaded': 0, 'uploaded': 0, 'need_to_upload': 0}
    script_list = {'needup': [], 'needdown': []}
    file_status_all = {}

    for ipw in pw:
        ipw = re.sub(r'/$', '', ipw)
        pw_core = ipw.rsplit('/', 1)[-1]
        ipw = os.path.realpath(ipw)
        # logger.warning(f'path={ipw}')
        # continue

        # fn_md5 = glob.glob(f'{ipw}/shell/*.md5.upload.sh')
        # if len(fn_md5) > 0 and not no_upload:
        #     script_list['needup'].append(fn_md5[0])

        remote_base_prefix = f'{remote_base}/' if remote_base else ''
        remote_pw_in = remote_pw or f'{remote_base_prefix}{pw_core}'
        
        remote_pw_in = remote_pw_in.replace('\\', '/')
        ct_prj, script_list, file_stat = main(ipw, script_list, remote_pw_in=remote_pw_in, updated_version_only=updated_version_only, no_upload=no_upload, remote_flat=remote_flat)
        
        for istat, v in file_stat.items():
            file_status_all.setdefault(istat, set()).update(set(v))

        for k in ct:
            ct[k] += ct_prj[k]

        # ct = {'total': n_desired, 'downloaded': n_downloaded, 'uploaded': n_uploaded, 'need_to_upload': n_need_upload}

    print('\n' + '#' * 50)
    print(json.dumps(ct, indent=4))

    # build the script list
    n_needup = len(script_list['needup'])
    n_needdown = len(script_list['needdown'])

    print(f'total scripts = {n_needup + n_needdown}')

    for udn, v1 in file_ct.items():
        print(udn)
        
        for rel, v2 in v1.items():
            line = f'    {rel}\n        ft,ok,down,raw,na_url'.replace(',', '\t')
            for ift, v3 in v2.items():
                n_invalid_url = n_ok = n_downloaded = n_raw = 0
                for ifn in v3:
                    if ifn in file_status_all['invalid_url']:
                        n_invalid_url += 1
                    elif ifn in file_status_all['ok']:
                        n_ok += 1
                    elif ifn in file_status_all['downloaded']:
                        n_downloaded += 1
                    elif ifn in file_status_all['raw']:
                        n_raw += 1

                n_ok = green(n_ok) if n_ok > 0 else 0
                n_downloaded = green(n_downloaded) if n_downloaded > 0 else 0
                n_raw = red(n_raw) if n_raw > 0 else green(0)
                n_invalid_url = red(n_invalid_url) if n_invalid_url > 0 else green(0)
                
                line += f'\n        {ift}\t{n_ok}\t{n_downloaded}\t{n_raw}\t{n_invalid_url}'
            print(line)

    with open(f'script_list.txt', 'w') as o:
        for k in ['needup', 'needdown']:
            if len(script_list[k]) > 0:
                print('\n'.join(script_list[k]), file=o)
