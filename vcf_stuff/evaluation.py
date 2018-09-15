from os.path import dirname, abspath, join, basename, isfile
import pandas as pd
import numpy as np
import math


###################
##### Metrics #####
###################

def f_measure(b, prec, recall):
    return (1 + b**2) * prec * recall / (b**2 * prec + recall) if prec + recall > 0 else 0


#########################
##### TXT rendering #####
#########################

def dislay_stats_df(df):
    """ Pretty-printing the table on the screen
    """
    with pd.option_context(
            'display.max_rows', None,
            'display.max_columns', None,
            'display.width', None,
            'display.float_format', lambda v: '{:,.2f}%'.format(100.0*v),
            ):
        print(df.to_string(index=True, na_rep='.'))


#############################
##### Jupyter rendering #####
#############################

def df_to_html(df):
    max_and_min = {
        'Prec':   (0.6, 1),
        'Recall': (0.8, 1),
        'F2':     (0,   1),
    }

    def _highlight_max(vals):
        BLUE_HUE = 211
        RED_HUE = 1
        MIN_BLUE, MAX_BLUE = 50, 100  # colors for minv..maxv
        MIN_RED, MAX_RED = 40, 100    # colors for 0..minv
        styles = []
        for v in vals:
            style = ''
            if not math.isnan(v):
                if len(vals.name) == 2 and vals.name[1] in max_and_min:
                    minv, maxv = max_and_min[vals.name[1]]
                    if v > minv:
                        k = (v - minv) / (maxv - minv)
                        l = MAX_BLUE - k * (MAX_BLUE - MIN_BLUE)
                        bg_style = f'background-color: hsl({BLUE_HUE}, 100%, {l}%)'
                        text_clr = 'white' if l < 70 else 'black'
                    else:
                        assert minv != 0, (minv, v)
                        k = (v - 0) / (minv - 0)
                        l = MIN_RED + k * (MAX_RED - MIN_RED)
                        bg_style = f'background-color: hsl({RED_HUE}, 100%, {l}%)'
                        text_clr = 'white' if l < 70 else 'black'
                    style = f'{bg_style}; color: {text_clr}; '
            styles.append(style)
        return styles

    df = df.set_index(df[('Sample')].values).drop([('Sample')], axis=1)
    return df.style\
        .format(lambda v: ('{:,.2f}%' if v !=1 else '{}%').format(100.0*v) if isinstance(v, float) else v)\
        .apply(_highlight_max)\
        .set_table_styles([dict(selector='td', props=[('padding', '1px 4px 1px 4px')]),
                           dict(selector='th', props=[('text-align', 'left')]),
                           dict(selector='tr', props=[('line-height', '10px')]),
                           dict(selector='td', props=[('font-size', '14px')]),
                          ])\
        .render()
