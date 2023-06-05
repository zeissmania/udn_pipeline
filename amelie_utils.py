from . import amelie_api
from chctool_lite import getlogger, line_count

logger_name = 'query_amelie'
intermediate_folder = 'intermediate'
def query_amelie(pw, prj, logger=None, force=False) -> 'build all the amelie result files':
    """
    run the amelie API
    """
    logger = logger or getlogger(logger_name)

    fn = f'{pw}/{prj}.amelie.lite.txt'

    if os.path.exists(fn) and not force:
        try:
            amelie_dict = get_amelie_dict(pw, prj, logger)
        except:
            logger.info(f'rerun amelie query')
            os.system(f'rm {fn} 2>/dev/null')
            force=True
        else:
            return amelie_dict
    logger.info('getting information from AMELIE')

    fn_genelist = f'{pw}/{intermediate_folder}/{prj}.genelist'
    fn_hop_pure_id = f'{pw}/{intermediate_folder}/{prj}.terms_pure_hpo.txt'
    fn_pheno = f'{pw}/pheno.keywords.txt'

    total_genes = line_count(fn_genelist)
    if total_genes > 1000:
        logger.warning(f'the gene list for AMELIE is larger than 1000, ({total_genes}), would split the genelist')
        
    amelie_api.main_old(prj, fn_genelist, fn_hop_pure_id, pheno_for_match=fn_pheno, pw=pw, force=force)
    logger.info('amelie done')
    return get_amelie_dict(pw, prj, logger)

def get_amelie_dict(pw, prj, logger=None):
    """
    after querying the amelie API
    convert {prj}.amelie.lite.txt(only the final score of gene) to dict
    """
    fn = f'{pw}/{prj}.amelie.lite.txt'
    logger = logger or getlogger(logger_name)

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
