# Calculating evaluation statistics (recall, precision, F1, etc) to display in pretty table

import seaborn as sns
import warnings
import itertools

def _f_measure(b, prec, recall):
    #assert (b**2 * prec + recall) != 0, (b, prec, recall)
    return (1 + b**2) * prec * recall / (b**2 * prec + recall)

def _cnt_stats(cdf_c, aln):
    snp = dict()
    ind = dict()
    
    snp['TP'] = len(cdf_c.query(f'{aln}_t == "tp" & is_snp'))
    ind['TP'] = len(cdf_c.query(f'{aln}_t == "tp" & not is_snp'))
    snp['FN'] = len(cdf_c.query(f'{aln}_t == "fn" & is_snp'))
    ind['FN'] = len(cdf_c.query(f'{aln}_t == "fn" & not is_snp'))
    snp['FP'] = len(cdf_c.query(f'{aln}_t == "fp" & is_snp'))
    ind['FP'] = len(cdf_c.query(f'{aln}_t == "fp" & not is_snp'))
    
    snp['P']      = snp['TP'] + snp['FP']
    snp['Prec']   = snp['TP'] / snp['P'] if snp['P'] else 0
    snp['T']      = snp['TP'] + snp['FN']
    snp['Recall'] = snp['TP'] / snp['T'] if snp['T'] else 0

    ind['P']      = ind['TP'] + ind['FP']
    ind['Prec']   = ind['TP'] / ind['P'] if ind['P'] else 0
    ind['T']      = ind['TP'] + ind['FN']
    ind['Recall'] = ind['TP'] / ind['T'] if ind['T'] else 0

    snp['F1'] = _f_measure(1, snp['Prec'], snp['Recall'])
    snp['F2'] = _f_measure(2, snp['Prec'], snp['Recall'])
    ind['F1'] = _f_measure(1, ind['Prec'], ind['Recall'])
    ind['F2'] = _f_measure(2, ind['Prec'], ind['Recall'])
    return snp, ind
    
def make_stats(bn, query_fmt=None, recalc_ensemble=False, add_truth=False):
    stats_by_sample = dict()
    #bns = ([bn] if isinstance(bn, str) else bn) if bn else benchmarks_names
    print(f'Benchmark: {bn}' + (f', filtered: {query_fmt}' if query_fmt else ''))
    b = benchmark_d[bn]
    df = b['df']
    for aln in b['aligners']:
        if query_fmt:
            df_c = df.copy()
            rej = df_c.query('not ' + query_fmt.format(**locals()))
            df_c = df_c.drop(rej.query(f'{aln}_t == "fp"').index)
            for i in rej.query(f'{aln}_t == "tp"').index:
                df_c.loc[i,f'{aln}_t'] = 'fn'
        else:
            df_c = df

        lbls_called_by = dict() 
        for clr in b['callers']:
            cdf_c = df_c.query(f'caller == "{clr}"')
            stats_by_sample[f'{aln} {clr}'] = _cnt_stats(cdf_c, aln)                
            if clr != 'ensemble':
                lbls_called_by[clr] = set(cdf_c.query(f'{aln}_t == "tp" | {aln}_t == "fp"').index)

        if recalc_ensemble:
            ensemble_lbls = set()
            for c1, c2 in itertools.combinations(lbls_called_by, 2):
                ensemble_lbls |= (lbls_called_by[c1] & lbls_called_by[c2])
            ensemble_rej = set().union(*lbls_called_by.values()) - ensemble_lbls
            edf = df_c.query(f'caller == "ensemble"').copy()
            for i in (ensemble_rej & set(edf.index)):
                if edf.loc[i,f'{aln}_t'] == 'tp':
                    edf.loc[i,f'{aln}_t'] = 'fn'
                if edf.loc[i,f'{aln}_t'] == 'fp':
                    edf.loc[i,f'{aln}_t'] = 'tn'
            stats_by_sample[f'ens_re {aln}'] = _cnt_stats(edf, aln)  
    
    if add_truth:
        some_aln = b["aligners"][0]
        some_clr = b["callers"][0]
        truth_snp = len(df.query(f'({some_aln}_t in ["fp", "fn"]) and caller=="{some_clr}" and is_snp'))
        truth_ind = len(df.query(f'({some_aln}_t in ["fp", "fn"]) and caller=="{some_clr}" and not is_snp'))
        stats_by_sample[f'{bn} truth'] = {'T': truth_snp}, {'T': truth_ind}
    return stats_by_sample







# Printing statistics in a pretty table

import plotly.plotly as py
import plotly.graph_objs as go
import colorlover as cl
import seaborn as sns

max_and_min = {
    ('SNP', 'Prec'):     (0.6, 1),
    ('INDEL', 'Prec'):   (0.6, 1),
    ('SNP', 'Recall'):   (0.8, 1),
    ('INDEL', 'Recall'): (0.8, 1),
    ('SNP', 'F2'):       (0,   1),
    ('INDEL', 'F2'):     (0,   1),
}

def stats_to_df(stat_by_sname):
    metircs = ['FP' , 'FN' , 'Recall', 'Prec', 'F2']
    idx = pd.MultiIndex.from_arrays([
            ['Sample'] + ['SNP']*len(metircs) + ['INDEL']*len(metircs), 
            [''] + metircs + metircs],
        names=['', 'Sample'])
    data = []
    s_truth = i_truth = None
    for sname, (snp, ind) in stat_by_sname.items():
        d = {('Sample', ''): sname}
        for lbl, st in [('SNP', snp), ('INDEL', ind)]:
            for m in metircs:
                d[(lbl, m)] = st.get(m, m['T'] if m == 'TP' else (0 if m in ['FP', 'FN'] else 1.0))
        data.append(d)
    return pd.DataFrame(data, columns=idx)

def _print_stats_txt(stats_by_sample):
    df = stats_to_df(stats_by_sample)
    df = df.set_index(df[('Sample', '')].values).drop([('Sample', '')], axis=1) 
    with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.width', None,
            'display.float_format', lambda v: '{:,.2f}%'.format(100.0*v)
        ):
        print(df.to_string(index=True))

def display_stats(stats_by_sample):
    def _highlight_max(vals):
        BLUE_HUE = 211
        RED_HUE = 1
        MIN_BLUE, MAX_BLUE = 50, 100  # colors for minv..maxv
        MIN_RED, MAX_RED = 40, 100    # colors for 0..minv
        styles = []
        for v in vals:
            style = ''
            if vals.name in max_and_min:
                minv, maxv = max_and_min[vals.name]
                if v > minv:  
                    k = (v - minv) / (maxv - minv) 
                    l = MAX_BLUE - k * (MAX_BLUE - MIN_BLUE)
                    bg_style = f'background-color: hsl({BLUE_HUE}, 100%, {l}%)'
                    text_clr = 'white' if l < 70 else 'black'
                else:
                    assert minv != 0, minv, k
                    k = (v - 0) / (minv - 0) 
                    l = MIN_RED + k * (MAX_RED - MIN_RED)
                    bg_style = f'background-color: hsl({RED_HUE}, 100%, {l}%)'
                    text_clr = 'white' if l < 70 else 'black'
                style = f'{bg_style}; color: {text_clr}; '
            styles.append(style)
        return styles

    df = stats_to_df(stats_by_sample)
    df = df.set_index(df[('Sample', '')].values).drop([('Sample', '')], axis=1) 
    return df.style\
        .format(lambda v: ('{:,.2f}%' if v !=1 else '{}%').format(100.0*v) if isinstance(v, float) else v)\
        .apply(_highlight_max)\
        .set_table_styles([dict(selector='td', props=[('padding', '1px 4px 1px 4px')]),
                           dict(selector='th', props=[('text-align', 'left')]), 
                           dict(selector='tr', props=[('line-height', '10px')]),
                           dict(selector='td', props=[('font-size', '14px')]),
                          ])