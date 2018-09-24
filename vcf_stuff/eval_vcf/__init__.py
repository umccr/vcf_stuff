from os.path import dirname, abspath

import pandas as pd


def package_path():
    return dirname(abspath(__file__))


def vcf_stats_to_df(stat_by_sname):
    idx = pd.MultiIndex.from_arrays([
        ['Sample', 'SNP', 'SNP', 'SNP', 'SNP'   , 'SNP',  'SNP', 'SNP', 'SNP', 'INDEL', 'INDEL', 'INDEL' , 'INDEL' , 'INDEL', 'INDEL', 'INDEL', 'SNP',],
        [''      , 'TP' , 'FP' , 'FN' , 'Recall', 'Prec', 'F1' , 'F2' , 'F3' , 'TP'   , 'FP'   , 'FN'    , 'Recall', 'Prec' , 'F1'   , 'F2'   , 'F3' ,]
        ],
        names=['1', '2'])

    data = []
    s_truth = i_truth = None
    for sname, stats in stat_by_sname.items():
        s_truth, s_tp, s_fp, s_fn, s_rec, s_prec, s_f1, s_f2, s_f3, \
        i_truth, i_tp, i_fp, i_fn, i_rec, i_prec, i_f1, i_f2, i_f3 = stats
        data.append({
            ('Sample', ''): sname,
            ('SNP', 'TP'): int(s_tp),
            ('SNP', 'FP'): int(s_fp),
            ('SNP', 'FN'): int(s_fn),
            ('SNP', 'Recall'): float(s_rec),  # TP/truth - grows with TP grow
            ('SNP', 'Prec'): float(s_prec),  # TP/called - falls with FP grow
            ('SNP', 'F1'): float(s_f1),
            ('SNP', 'F2'): float(s_f2),  # Value recall twice higher than precision
            ('SNP', 'F3'): float(s_f3),
            ('INDEL', 'TP'): int(i_tp),
            ('INDEL', 'FP'): int(i_fp),
            ('INDEL', 'FN'): int(i_fn),
            ('INDEL', 'Recall'): float(i_rec),
            ('INDEL', 'Prec'): float(i_prec),
            ('INDEL', 'F1'): float(i_f1),
            ('INDEL', 'F2'): float(i_f2),
            ('INDEL', 'F3'): float(i_f3),
        })
    # data.append({
    #     ('Sample', ''): 'Truth',
    #     ('SNP', 'TP'): int(s_truth),
    #     ('SNP', 'FP'): 0,
    #     ('SNP', 'FN'): int(s_truth),
    #     ('SNP', 'Prec'): 1.0,
    #     ('SNP', 'Recall'): 1.0,
    #     ('INDEL', 'TP'): int(i_truth),
    #     ('INDEL', 'FP'): 0,
    #     ('INDEL', 'FN'): int(i_truth),
    #     ('INDEL', 'Prec'): 1.0,
    #     ('INDEL', 'Recall'): 1.0
    # })
    return pd.DataFrame(data, columns=idx)
