#! /usr/bin/env python3
"""
generate the IGV batch script
https://software.broadinstitute.org/software/igv/PortCommands
https://github.com/igvteam/igv/blob/master/src/main/resources/org/broad/igv/prefs/preferences.tab
input 1 = igv download file
input 2 = {udn}.selected.genes.txt. first column = mut_type

if input2 is not defined, would create the script just adding the bam tracks
"""
import os, sys
import re

def get_logger():
    import logging
    fmt = logging.Formatter('%(asctime)s  %(levelname)-9s   %(funcName)-10s   line: %(lineno)-5s   %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(fmt)
    console.setLevel('INFO')

    logger = logging.getLogger(__file__)
    logger.setLevel('DEBUG')
    logger.addHandler(console)
    return logger

class IGV():
    def __init__(self, pw, udn, logger=None):
        self.pw = pw
        self.udn = udn
        self.logger = logger or get_logger()

    def get_bam_list(self):
        # read-in the bam files
        logger = self.logger
        # fn_in = f'{self.pw}/{self.udn}.igv.files.txt'
        fn_in = f'{self.pw}/download.info.{self.udn}.txt'
        if not os.path.exists(fn_in):
            logger.error(f'{self.udn}: {fn_in} file not found')
            sys.exit(1)
        # logger.info(fn_in)

        d_bam = {}
        rank = {'Proband': 0, 'Father': 1, 'Mother': 2}
        sn_brother, sn_sister, sn_sibling, sn_other = 1, 1, 1, 1
        rel_rank = {}
        with open(fn_in) as fp:
            header = fp.readline().split('\t')
            # rel_to_proband	fn	url	udn	seq_type	size	build	md5	url_s3
            keys = ['rel_to_proband', 'fn', 'url']
            try:
                header_idx = {_: header.index(_) for _ in header if _ in keys}
            except:
                logger.error(f'some keys not found in header of {fn_in}:\nexpected={keys}\nheader={header}')
                sys.exit(1)
            for i in fp:
                i = i.strip()
                if not i:
                    continue
                a = i.split('\t')
                rel = a[header_idx['rel_to_proband']]
                if rel not in rel_rank:
                    try:
                        rel_rank[rel] = rank[rel]
                    except:
                        if rel.lower().startswith('brother'):
                            sn_brother += 1  # should be less than 90
                            rel_rank[rel] = 10 + sn_brother
                        elif rel.lower().startswith('bister'):
                            sn_sister += 1 # should be less than 100
                            rel_rank[rel] = 100 + sn_sister
                        elif rel.lower().startswith('sibling'):
                            sn_sibling += 1  # should be less than 100
                            rel_rank[rel] = 200 + sn_sibling
                        else:
                            sn_other += 1
                            rel_rank[rel] = 300 + sn_other # no limit


                fn = a[header_idx['fn']]
                url = a[header_idx['url']]
                ext = fn.rsplit('.', 1)[-1]
                lb = fn.replace('.bai', '')
                if ext not in set(['bam', 'bai']):
                    continue
                try:
                    d_bam[lb]
                except:
                    d_bam[lb] = [rel_rank[rel], 'error', 'error']
                if ext == 'bam':
                    d_bam[lb][1] = url
                elif ext == 'bai':
                    d_bam[lb][2] = url

        # sort by pre-defined relative order
        bam_list = sorted(d_bam.values(), key=lambda _: _[0])
        return bam_list

    def get_gene_region(self):
        # read-in the region/ genes
        logger = self.logger
        pw = self.pw
        udn = self.udn
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

def main(pw, udn, logger):
    # output igv batch script
    pw = os.path.realpath(pw)
    pwigv = f'{pw}/igv'
    os.system(f'mkdir -p {pwigv} 2>/dev/null')
    igv_obj = IGV(pw, udn, logger)
    bam_list = igv_obj.get_bam_list()
    if len(bam_list) == 0:
        logger.error(f'no bam file found in {pw}/{udn}.igv.files.txt')

    d_genes = igv_obj.get_gene_region()

    with open(f'{pw}/igv.script.{udn}.txt', 'w') as out:
        out.write(f"""genome hg38
    SAM.QUALITY_THRESHOLD 13
    snapshotDirectory {pwigv}
    maxPanelHeight 2000
    squish
    setSleepInterval 500
    new

    """)
        for _rank, url_bam, url_bai in bam_list:
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

if __name__ == "__main__":

    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('udn', help="""udn ID, if not sepecified, would search the UDNxxxx.igv.files.txt""", nargs='?')
    ps.add_argument('-pw', help="""path for the UDN project, default = pwd""", default=None)
    args = ps.parse_args()


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
    logger = get_logger()
    main(pw, udn, logger)
