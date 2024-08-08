import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import xiutilities.bi_fdr as fdr


def fdr_cutoff_plot(df, xlim=0.05, score_col='match_score', rescore_col='xiMLScore', fdr_steps=0.001,
                    top_col='top_ranking', top_col_rescored='xiMLScore_top_ranking',
                    native_fdr_col=None, rescore_fdr_col=None, mismatch_col=None, ylabel='CSMs',
                    print_points=[]):
    # Calculate FDR if none provided
    if native_fdr_col is None:
        df.drop(['fdr_native', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_native'] = fdr.calculate_bi_fdr(df, score_col=score_col)
    else:
        df.loc[:, 'fdr_native'] = df.loc[:, native_fdr_col]

    if rescore_fdr_col is None:
        df.drop(['fdr_native', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_rescored'] = fdr.calculate_bi_fdr(df, score_col=rescore_col)
    else:
        df.loc[:, 'fdr_rescored'] = df.loc[:, rescore_fdr_col]
    
    df = df[
        ((df['fdr_native'] <= xlim) | (df['fdr_rescored'] <= xlim)) & df['isTT']
    ]

    fdr_x = np.arange(0, xlim+fdr_steps, fdr_steps)

    samples_native_cumsum = pd.DataFrame(index=fdr_x)
    samples_rescored_cumsum = pd.DataFrame(index=fdr_x)
    samples_gained_cumsum = pd.DataFrame(index=fdr_x)
    samples_lost_cumsum = pd.DataFrame(index=fdr_x)
    rescore_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    native_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    
    # Calculate sample number for FDR cutoff
    native_below = df[
        (df['fdr_native'] <= xlim) &
        df['isTT']
    ]
    samples_native_cumsum['count'] = samples_native_cumsum.index.to_series().apply(
        lambda val: len(native_below[native_below['fdr_native'] <= val])
    )
    
    rescored_below = df[
        (df['fdr_rescored'] <= xlim) &
        df['isTT']
    ]
    samples_rescored_cumsum['count'] = samples_rescored_cumsum.index.to_series().apply(
        lambda val: len(rescored_below[rescored_below['fdr_rescored'] <= val])
    )
    
    # Calculate gain
    gain_cands = df[
        (df['fdr_rescored'] <= xlim) &
        df['isTT']
    ]
    samples_gained_cumsum['count'] = samples_gained_cumsum.index.to_series().apply(
        lambda val: len(gain_cands[
            (gain_cands['fdr_native'] > val) &
            (gain_cands['fdr_rescored'] <= val)
        ]) + len(rescored_below[
                (rescored_below['fdr_rescored'] <= val) &
                (~rescored_below[top_col]) & rescored_below[top_col_rescored]
        ])
    )

    # Calculate loss
    loss_cands = df[
        (df['fdr_native'] <= xlim) &
        df['isTT']
    ]
    samples_lost_cumsum['count'] = samples_lost_cumsum.index.to_series().apply(
        lambda val: len(loss_cands[
            (
                (loss_cands['fdr_native'] <= val) &
                (loss_cands['fdr_rescored'] > val)
            )
        ]) + len(native_below[
                (native_below['fdr_native'] <= val) &
                native_below[top_col] & (~native_below[top_col_rescored])
        ])
    )

    if mismatch_col is not None:
        # Calculate rescore mismatch 
        rescore_mismatch_cand = df[
            (df['fdr_rescored'] <= xlim) &
            df['isTT'] & df[mismatch_col]
        ]
        rescore_mismatch_cumsum['count'] = rescore_mismatch_cumsum.index.to_series().apply(
            lambda val: len(rescore_mismatch_cand[
                (rescore_mismatch_cand['fdr_rescored'] <= val)
            ])
        )

        # Calculate native mismatch
        native_mismatch_cand = df[
            (df['fdr_native'] <= xlim) &
            df['isTT'] & df[mismatch_col]
        ]
        native_mismatch_cumsum['count'] = native_mismatch_cumsum.index.to_series().apply(
            lambda val: len(native_mismatch_cand[
                (
                    (native_mismatch_cand['fdr_native'] <= val)
                )
            ])
        )
    
    fig, ax = plt.subplots(figsize=(12, 6), dpi=100)

    data = {
        'native': samples_native_cumsum['count'],
        'rescored': samples_rescored_cumsum['count'],
        'gained': samples_gained_cumsum['count'],
        'lost': samples_lost_cumsum['count']
    }
    
    sns.lineplot(
        ax=ax,
        data=data
    )
    
    plt.xlabel('FDR cutoff')
    plt.ylabel(ylabel)
    
    ax.set_xlim(0, xlim)
    xfilter = samples_native_cumsum.loc[:xlim].index
    
    ylim = max([
        samples_native_cumsum.loc[xfilter, 'count'].max(),
        samples_rescored_cumsum.loc[xfilter, 'count'].max()
    ])*1.1
    
    ax.set_ylim(
        0,
        ylim
    )
    
    for x in print_points:
        if x > fdr_x[-1]:
            continue
        x_idx = 0
        while x > fdr_x[x_idx] or x_idx >= len(fdr_x):
            x_idx += 1
        y_native = samples_native_cumsum['count'].iloc[x_idx]
        y_rescored = samples_rescored_cumsum['count'].iloc[x_idx]
        y = max([y_native, y_rescored])
        text = f"{y_native}→{y_rescored}"
        ax.plot(
            x,
            y,
            'xk'
        )
        plt.text(
            x=x,
            y=y+ylim*0.01,
            s=text
        )
    
    if mismatch_col is not None:
        data2 = {
            'native mismatch': native_mismatch_cumsum['count']/data['native']*100,
            'rescore mismatch': rescore_mismatch_cumsum['count']/data['rescored']*100,
            'fdr': fdr_x*100,
        }
        ax2 = ax.twinx()
        
        sns.lineplot(
            ax=ax2,
            data=data2
        )
    
    return fig


def fdr_mismatch_plot(df, xlim=0.05, score_col='match_score', rescore_col='xiMLScore', fdr_steps=0.001,
                      native_fdr_col=None, rescore_fdr_col=None, mismatch_col='mismatch', ylabel='proportion'):
    # Calculate FDR if none provided
    if native_fdr_col is None:
        df.drop(['fdr_native', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_native'] = fdr.calculate_bi_fdr(df, score_col=score_col)
    else:
        df.loc[:, 'fdr_native'] = df.loc[:, native_fdr_col]

    if rescore_fdr_col is None:
        df.drop(['fdr_native', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_rescored'] = fdr.calculate_bi_fdr(df, score_col=rescore_col)
    else:
        df.loc[:, 'fdr_rescored'] = df.loc[:, rescore_fdr_col]
    
    df = df.loc[
        (df['fdr_native'] <= xlim) | (df['fdr_rescored'] <= xlim)
    ]

    fdr_x = np.arange(0, xlim+fdr_steps, fdr_steps)

    fdr_series = pd.DataFrame(index=fdr_x, data={'fdr': fdr_x})
    samples_native_cumsum = pd.DataFrame(index=fdr_x)
    samples_rescored_cumsum = pd.DataFrame(index=fdr_x)
    rescore_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    native_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    
    # Calculate sample number for FDR cutoff
    native_below = df.loc[df['fdr_native'] <= xlim]
    samples_native_cumsum['count'] = samples_native_cumsum.index.to_series().apply(
        lambda val: len(native_below[native_below['fdr_native'] <= val])
    )
    
    rescored_below = df.loc[df['fdr_rescored'] <= xlim]
    samples_rescored_cumsum['count'] = samples_rescored_cumsum.index.to_series().apply(
        lambda val: len(rescored_below.loc[rescored_below['fdr_rescored'] <= val])
    )
    
    # Calculate rescore mismatch 
    rescore_mismatch_cand = df.loc[
        (df['fdr_rescored'] <= xlim) &
        df['isTT'] & df[mismatch_col]
    ]
    rescore_mismatch_cumsum.loc[:, 'count'] = rescore_mismatch_cumsum.index.to_series().apply(
        lambda val: len(rescore_mismatch_cand[
            (rescore_mismatch_cand['fdr_rescored'] <= val)
        ])
    )

    # Calculate native mismatch
    native_mismatch_cand = df.loc[
        (df['fdr_native'] <= xlim) &
        df['isTT'] & df[mismatch_col]
    ]
    native_mismatch_cumsum.loc[:, 'count'] = native_mismatch_cumsum.index.to_series().apply(
        lambda val: len(native_mismatch_cand[
            (
                (native_mismatch_cand['fdr_native'] <= val)
            )
        ])
    )
    
    fig, ax = plt.subplots(figsize=(12, 6), dpi=100)

    data = {
        'native mismatch': native_mismatch_cumsum['count']/samples_native_cumsum['count'],
        'rescore mismatch': rescore_mismatch_cumsum['count']/samples_rescored_cumsum['count'],
        'fdr': fdr_series['fdr'],
    }
    
    sns.lineplot(
        ax=ax,
        data=data
    )
    
    plt.xlabel('FDR cutoff')
    plt.ylabel(ylabel)
    
    ax.set_xlim(0, xlim)
    
    return fig, data
