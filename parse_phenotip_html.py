"""
from UDN gateway -> phenotip, there is a part which include all the phenotype of this proband
this script is to parse the html into phenotype text
"""

def parse(html) -> 'str, all the phenotypes':
    """
    html is the html text
    include the elements like
     <span class="search-term symptom" title="Click to disable">Abnormality of the cerebral white matter</span>
     <span class="search-term not_symptom" title="Click to disable">Tracheoesophageal fistula</span>
    """

    from bs4 import BeautifulSoup as bs
    r = bs(html)
    x = r.find_all('span', attrs={'class': 'search-term symptom'})
    return '\n'.join([_.text for _ in x])


if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('html_text', help="""html_text""")
    args = ps.parse_args()

    html = args.html_text
    print(parse(html))
