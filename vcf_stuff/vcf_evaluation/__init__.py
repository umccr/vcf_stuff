import pandas as pd


def stats_to_df(stat_by_sname):
    idx = pd.MultiIndex.from_arrays([
        ['Sample', 'SNP', 'SNP', 'SNP', 'SNP',  'SNP'   , 'INDEL', 'INDEL', 'INDEL' , 'INDEL', 'INDEL'],
        [''      , 'TP' , 'FP' , 'FN' , 'Prec', 'Recall', 'TP'   , 'FP'   , 'FN'    , 'Prec' , 'Recall']
        ],
        names=['1', '2'])

    data = []
    s_truth = i_truth = None
    for sname, stats in stat_by_sname.items():
        s_truth, s_tp, s_fp, s_fn, s_prec, s_rec, i_truth, i_tp, i_fp, i_fn, i_prec, i_rec = stats
        data.append({
            ('Sample', ''): sname,
            ('SNP', 'TP'): int(s_tp),
            ('SNP', 'FP'): int(s_fp),
            ('SNP', 'FN'): int(s_fn),
            ('SNP', 'Prec'): float(s_prec),  # TP/called
            ('SNP', 'Recall'): float(s_rec),  # TP/truth
            ('INDEL', 'TP'): int(i_tp),
            ('INDEL', 'FP'): int(i_fp),
            ('INDEL', 'FN'): int(i_fn),
            ('INDEL', 'Prec'): float(i_prec),
            ('INDEL', 'Recall'): float(i_rec)
        })
    data.append({
        ('Sample', ''): 'Truth',
        ('SNP', 'TP'): int(s_truth),
        ('SNP', 'FP'): 0,
        ('SNP', 'FN'): int(s_truth),
        ('SNP', 'Prec'): 1.0,
        ('SNP', 'Recall'): 1.0,
        ('INDEL', 'TP'): int(i_truth),
        ('INDEL', 'FP'): 0,
        ('INDEL', 'FN'): int(i_truth),
        ('INDEL', 'Prec'): 1.0,
        ('INDEL', 'Recall'): 1.0
    })
    return pd.DataFrame(data, columns=idx)


def dislay_stats_df(df):
    """ Pretty-printing the table on the screen
    """
    with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.width', None,
            'display.float_format', lambda v: '{:,.2f}%'.format(100.0*v)
            ):
        print(df.to_string(index=True))
