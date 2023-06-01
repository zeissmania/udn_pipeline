
def check_amelie_missing(self):
    pw = self.pw
    prj = self.prj
    logger = self.logger

    fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
    fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'

    total_genes = line_count(fn_genelist)
    total_hpo = line_count(fn_hop_pure_id)

    if total_genes > 1000:
        logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), max input is 1000')

    if os.path.exists(f'{pw}/err_amelie.{prj}.genes_not_added.txt'):
        bad_genes = line_count(f'{pw}/err_amelie.{prj}.genes_not_added.txt')
        logger.warning(f'AMELIE missing genes: {bad_genes}/{total_genes} not included in AMELIE query')


def query_amelie(self, force=False) -> 'build all the amelie result files':
    """
    run the amelie API
    """
    pw = self.pw
    prj = self.prj
    logger = self.logger

    fn = f'{pw}/{prj}.amelie.lite.txt'

    if os.path.exists(fn) and not force:
        try:
            self.amelie_dict = self.get_amelie_dict()
        except:
            logger.info(f'rerun amelie query')
            os.system(f'rm {fn} 2>/dev/null')
            force=True
        else:
            # self.check_amelie_missing()
            return 0
    logger.info('getting information from AMELIE')

    fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
    fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'
    fn_pheno = f'{pw}/pheno.keywords.txt'

    total_genes = line_count(fn_genelist)
    # total_hpo = line_count(fn_hop_pure_id)

    if total_genes > 1000:
        logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), would split the genelist')

    # main(prj, genelist, pheno_hpo_list, pw=None, pw_main=None, force=False)
    # amelie_api.main(prj, fn_genelist, fn_hop_pure_id, pw=pw, force=force)
    amelie_api.main_old(prj, fn_genelist, fn_hop_pure_id, pheno_for_match=fn_pheno,pw=pw, force=force)
    logger.info('amelie done')
    # self.check_amelie_missing()

    # # link the final amelie files
    # os.system(f'ln {pw}/{intermediate_folder}/{prj}.amelie.matched_query.txt {pw} 2>/dev/null')
    # os.system(f'ln {pw}/{intermediate_folder}/{prj}.amelie.lite.txt {pw} 2>/dev/null')

def get_amelie_dict(self):
    """
    after querying the amelie API
    convert {prj}.amelie.lite.txt(only the final score of gene) to dict
    """
    pw = self.pw
    prj = self.prj
    logger = self.logger
    fn = f'{pw}/{prj}.amelie.lite.txt'

    try:
        d_amelie = [_.strip().split('\t') for _ in open(fn)]
        d_amelie = {k: v for k, v in d_amelie}
        logger.info(f'AMELIE count = {len(d_amelie)}')
    except:
        if not os.path.exists(fn):
            logger.error(f'amelie result not found: {fn}, would use dummy amelie score')
            fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
            return {gn.strip(): 0 for gn in open(fn_genelist) if gn.strip()}

    return d_amelie
