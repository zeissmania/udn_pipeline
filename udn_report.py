"""

"""
import sys, os, re
from openpyxl.styles import PatternFill, Border, Side, Alignment, Font, NamedStyle, Color
import openpyxl
from openpyxl.utils import get_column_letter
import pandas as pd

import xlsxwriter


# # def s():
# #     fno = 'test.xlsx'
# #     wb.save(fno)
# #     !open {fno}
# wb = openpyxl.Workbook()
# ws = wb.active
# ws.title = 'Original'


def main(prj, pw=None, fn_selected_genes=None):
    """
    expected file = selected.genes.xlsx
    prj = UDNxxxx
    pw, default is current folder
    fn_selected_genes , default is '{prj}.selected.genes.txt'
    first column should be the sv_type, valid = 'denovo, heter, xlink'
    other columns are just as the final tsv file generated by the pipeline

    # family info, is a list or tuple, like
    # ['Proand 123456', 'Mother (unaff) 22222', 'Father(unaff) 33333']
    """


    pw = pw or os.getcwd()
    for fn_tmp in [fn_selected_genes, f'{pw}/{prj}.selected.genes.xlsx', f'{pw}/selected.genes.xlsx']:
        if fn_tmp and os.path.exists(fn_tmp):
            fn_selected_genes = fn_tmp
            break
    else:
        print(f'ERROR, selected.genes.xlsx file not found!')
        return 1

    # group the selected SV into groups
    sv_selected = {}
    sv_type_conversion = {'denovo': 'De Novo',
                        'de novo': 'De Novo',
                        'd': 'De Novo',
                        'de': 'De Novo',
                        'heter': 'Heterozygous',
                        'he': 'Heterozygous',
                        'h': 'Heterozygous',
                        'hetero': 'Heterozygous',
                        'het': 'Heterozygous',
                        'xlink': 'X-linked',
                        'x-link': 'X-linked',
                        'x': 'X-linked',
                        'xl': 'X-linked',
                        }


    wb = xlsxwriter.Workbook(f'{pw}/{prj}.report.xlsx')
    ws = wb.add_worksheet('Original')


    # read the raw data
    d = pd.read_excel(fn_selected_genes, keep_default_na=False)
    header = d.columns
    prev_columns = list(header[:21])

    # sort the family
    family_info_raw = list(header[21:])
    family_info_priority = []
    d_family_priority = {'proband': 1, 'mother': 2, 'father':3, 'sister': 4, 'brother': 5}

    for i in family_info_raw:
        rel_to_proband = i.split(' ', 1)[0].lower()
        rel_to_proband = re.split('\W+', rel_to_proband)[0]
        try:
            family_info_priority.append(d_family_priority[rel_to_proband])
        except:
            family_info_priority.append(999)

    family_info = sorted(zip(family_info_raw, family_info_priority), key=lambda x: x[1])
    family_info = [_[0] for _ in family_info]

    header = prev_columns + family_info
    d = d[header]

    # data = [_.strip().split('\t') for _ in fh]

    bad = 0
    for i in d.itertuples(False, None):
        sv_type = i[0]
        if sv_type.lower() not in sv_type_conversion:
            if not sv_type.strip():
                print(f'mut type not specified: {i}')
                bad = 1
            try:
                if i[1] != "AMELIE":
                    print(f'wrong SV catergory found! svtype={sv_type}')
            except:
                print(f'wrong line: {i}')
            continue
        sv_type_long = sv_type_conversion[sv_type.lower()]
        try:
            sv_selected[sv_type_long].append(i[1:])
        except:
            sv_selected[sv_type_long] = [i[1:]]

    if bad:
        print('Fail to build report, add mut type of each SV in selected.genes first')
        return 1
    n_family = len(family_info)
    total_col = 11 + n_family


    # define formats
    formats = {}
    header1 = wb.add_format({'font_size': 16, 'bold': True, })
    header2 = wb.add_format({'font_size': 14, 'bold': True, })
    formats['header_sv'] = wb.add_format({'font_size': 14, 'bold': True, 'bg_color': '#d9d9d9', 'text_wrap': True})
    formats['header_main'] = wb.add_format({'font_size': 11, 'bold': True, 'bg_color': '#dce6f1', 'border': 2, 'border_color': '#d9d9d9', 'text_wrap': True})
    formats['header_main'].set_align('center')
    formats['header_main'].set_align('vcenter')

    formats['fmt_txt'] = wb.add_format({'font_size': 11, 'bold': False, 'border': 2, 'border_color': '#d9d9d9'})
    formats['fmt_txt'].set_align('center')
    formats['fmt_txt'].set_align('vcenter')

    formats['fmt_comment'] = wb.add_format({'font_size': 11, 'bold': False, 'border': 2, 'border_color': '#d9d9d9', 'text_wrap': True})
    formats['fmt_comment'].set_text_wrap()
    formats['fmt_comment'].set_align('vcenter')
    formats['fmt_comment'].set_align('left')

    formats['bold'] = wb.add_format({'font_size': 11, 'bold': True, 'text_wrap': True})
    formats['bold'].set_text_wrap()
    formats['bold'].set_align('vcenter')

    wrap = wb.add_format({'text_wrap': True, 'align': 'vcenter'})
    # set column width

    # ignore errors
    ws.ignore_errors({'number_stored_as_text': 'A1:V200'})


    width = [17, 16, 19.5, 10] + [10] * n_family + [12, 12, 12, 8, 9.5, 8, 13]
    for n, v in enumerate(width):
        ws.set_column(n, n, v)

    ws.set_column(11 + n_family, 11 + n_family, 90, cell_format=wrap)

    # first row
    row = 0
    ws.merge_range(row, 0, row, total_col, prj, cell_format=header1)

    note = 'Research Variants (need Sanger confirmation unless indicated under Comments)'
    row += 1
    ws.merge_range(row, 0, row, total_col, note, cell_format=header2)

    # note = 'new row'
    # row += 1
    # ws.merge_range(row, 0, row, total_col, note, cell_format=header2)

    # sort the rank with in he same sv type
    sv_selected_sorted = [(k, sorted(v, key=lambda _: int(_[0]))) for k, v in sv_selected.items()]

    # sort the order of sv type
    sv_selected_sorted = sorted(sv_selected_sorted, key=lambda _: int(_[1][0][0]))

    row += 1
    for sv_type_long, info in sv_selected_sorted:
        # sv_type_long = sv_type_conversion[sv_type]
        row = add_header(ws, row, sv_type_long, family_info, formats)

        for info1 in info:

            res = add_data(ws, row, info1, n_family, formats)
            if res is not None:
                row = res
    wb.close()



def add_header(ws, row, sv_type_long, family_info, formats):
    n_family = len(family_info)
    # sv_type line
    ws.merge_range(row, 0, row, 11 + n_family, f'Structural Variants({sv_type_long})', cell_format=formats['header_sv'])

    row += 1
    merged_cell_col_idx = [0, 2, 3] + [_ + 4 for _ in range(n_family)] + [_ + 7 + n_family for _ in [0, 1, 2, 3, 4]]
    merged_cell_value = ['Gene', 'Change', 'Effect'] + family_info +['Baylor WGS', 'Emedgene', 'Yu Shyr', 'PreUDN Panel', 'Comments']
    for col, v in zip(merged_cell_col_idx, merged_cell_value):
        ws.merge_range(row, col, row + 3, col, v, cell_format=formats['header_main'])


    # pos
    col_idx = [1] + [4 + n_family + _ for _ in [0, 1, 2]]
    values = [['Chr', 'Position', 'rs#', ''], ['Quality', 'GQ', 'Coverage', ''], ['GnomAD', 'DGV(dup/del)', 'GERP', 'CADD'], ['Missense Z', 'LoF pLI', '', '']]

    for col, v in zip(col_idx, values):
        for i_row, i in enumerate(v):
            ws.write_string(row + i_row, col, i, cell_format=formats['header_main'])
    return row + 4

def add_data(ws, row, data, n_family, formats):
    try:
        rank, _, chr_, s, e, sv_type, qual, exon_span_tag, gn, sv_len, exon_span, dgv_gain, dgv_loss, gnomad, _, ddd_disease, _, phenotype, inher = data[:19]
    except:
        print('wrong line format! ')
        print(len(data[:19]))
        print(data[:19])
    # exon_span_tag = '' if exon_span_tag == np.nan else exon_span_tag
    cn_proband = data[20]
    cn_family = list(data[21:])

    if len(cn_family) != n_family - 1:
        print(f'error, the family number count in file({len(cn_family)}) and Family info({n_family}) list are differnet !')
        sys.exit(1)

    pos = f'{s}-{e}'
    dgv = f'{dgv_gain}/{dgv_loss}'

    comment = phenotype

    # write the pos
    col = 1
    try:
        ws.write_string(row, col, chr_, cell_format=formats['fmt_txt'])
    except:
        print('wrong row, data = {data}, chr={chr_}')
        return None
    ws.merge_range(row + 1, col, row + 2, col, pos, cell_format=formats['fmt_txt'])
    ws.write_string(row + 3, col, sv_len, cell_format=formats['fmt_txt'])

    # write merged cells
    col_idx = [0, 3] + [4 + _ for _ in range(n_family)] + [4 + n_family + _ for _ in [3, 4, 6] ]
    values = [gn, sv_type, cn_proband] + cn_family + ['', '', '']
    for col, v in zip(col_idx, values):
         ws.merge_range(row, col, row + 3, col, v, cell_format=formats['fmt_txt'])


    # write the split cols
    col_idx = [2] + [4 + n_family + _ for _ in [0, 1, 2, 5]]
    values = [
            ['', exon_span, exon_span_tag, ''],
            [qual, '', '', ''],
            [dgv, gnomad, '', ''],
            ['', '', '', ''],
            [rank, '✔', '', '']
    ]
    for col, v in zip(col_idx, values):
       for n, iv in enumerate(v):
           iv = str(iv)
           ws.write_string(row + n, col, iv, cell_format=formats['fmt_txt'])

    # write the comment
    col = 11 + n_family
    ws.merge_range(row, col, row + 3, col, '', cell_format=formats['fmt_comment'])
    tmp = comment.split('**')
    comment_list = []
    bold = formats['bold']
    for n, _ in enumerate(tmp):
        if n % 2:
            comment_list.append(bold)
            comment_list.append(_)
        elif not _:
            continue
        else:
            comment_list.append(formats['fmt_comment'])
            comment_list.append(_)
    # print(comment_list)
    if len(comment_list) == 1:
        ws.write(row, col, comment_list[0])
    else:
        ws.write_rich_string(row, col, *comment_list)

    # prepare to set the row height
    # each line is about 100 char, height = 15 px
    tmp = comment.split('\n')
    tmp = [int(len(_)/105) + 1 for _ in tmp]
    height = sum(tmp) - 3 # exclude the 3 merged row height
    # print(f"expected_lines = {height}")
    ws.set_row(row + 3, height * 14)

    return row + 4



if __name__ == "__main__":
    import argparse as arg
    from argparse import RawTextHelpFormatter
    ps = arg.ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    ps.add_argument('prj', help="""case name/title""")
    ps.add_argument('fn', help="""selected SV list, 1st column = SV type(denovo/heter/xlink)
    2nd column is the rank
    last 2/3 column (col 21-end) must be the copy number of proband and family members
    """, nargs='?')
    # ps.add_argument('-family', help="""family info, put in the header line. sep by comma
    # like  'Proand 123456', 'Mother (unaff) 22222', 'Father(unaff) 33333'""", nargs='+')
    args = ps.parse_args()

    fn = args.fn
    prj = args.prj
    # family_info = ' '.join(args.family)
    # family_info = family_info.split(',')
    # family_info = [_.strip() for _ in family_info if _.strip()]
    pw = os.getcwd()
    if fn and not os.path.exists(fn):
        print(f'error, selected SV list file not exist ! {fn}')
        sys.exit(1)

    main(prj, pw, fn)
