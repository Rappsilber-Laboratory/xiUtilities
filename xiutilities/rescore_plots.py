import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import xiutilities.bi_fdr as fdr


def fdr_cutoff_plot(df,
                    xlim=0.05,
                    score_col='match_score',
                    rescore_col='xiMLScore',
                    fdr_steps=0.001,
                    top_col='top_ranking',
                    top_col_rescored='xiMLScore_top_ranking',
                    standard_fdr_col=None,
                    rescore_fdr_col=None,
                    mismatch_col=None,
                    target_col='isTT',
                    ylabel='CSMs',
                    print_points=[],
                    figsize=(12, 6),
                    dpi=100):
    """
    Plot how many CSMs fall below certain FDR thresholds (x-axis) for a standard and a rescored FDR.

    :param df: DataFrame of CSMs
    :param xlim: Up to which FDR cutoff to plot
    :param score_col: Name of the standard score column
    :param rescore_col: Name of the rescored score column
    :param fdr_steps: Step size of FDR cutoff points
    :param top_col: Name of the standard top ranking column
    :param top_col_rescored: Name of the rescored top ranking column
    :param standard_fdr_col: Name of the standard FDR column
    :param rescore_fdr_col: Name of the rescored FDR column
    :param mismatch_col: Name of the mismatch indication column (if not provided, no mismatch counts are plotted)
    :param target_col: Name of the target indication column
    :param ylabel: Label for the y-axis
    :param figsize: Figure size
    :param dpi: Figure DPI
    :param print_points: List of points where to print the relation of standard and rescored matches
    :return: A tuple of the figure and the DataFrame it originates from
    """
    # Calculate FDR if none provided
    if standard_fdr_col is None:
        df.drop(['fdr_standard', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_standard'] = fdr.calculate_bi_fdr(df, score_col=score_col)
    else:
        df.loc[:, 'fdr_standard'] = df.loc[:, standard_fdr_col]

    if rescore_fdr_col is None:
        df.drop(['fdr_standard', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_rescored'] = fdr.calculate_bi_fdr(df, score_col=rescore_col)
    else:
        df.loc[:, 'fdr_rescored'] = df.loc[:, rescore_fdr_col]
    
    df = df[
        ((df['fdr_standard'] <= xlim) | (df['fdr_rescored'] <= xlim)) & df[target_col]
    ]

    fdr_x = np.arange(0, xlim+fdr_steps, fdr_steps)

    samples_standard_cumsum = pd.DataFrame(index=fdr_x)
    samples_rescored_cumsum = pd.DataFrame(index=fdr_x)
    samples_gained_cumsum = pd.DataFrame(index=fdr_x)
    samples_lost_cumsum = pd.DataFrame(index=fdr_x)
    rescore_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    standard_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    
    # Calculate sample number for FDR cutoff
    standard_below = df[
        (df['fdr_standard'] <= xlim) &
        df[target_col]
    ]
    samples_standard_cumsum['count'] = samples_standard_cumsum.index.to_series().apply(
        lambda val: len(standard_below[standard_below['fdr_standard'] <= val])
    )
    
    rescored_below = df[
        (df['fdr_rescored'] <= xlim) &
        df[target_col]
    ]
    samples_rescored_cumsum['count'] = samples_rescored_cumsum.index.to_series().apply(
        lambda val: len(rescored_below[rescored_below['fdr_rescored'] <= val])
    )
    
    # Calculate gain
    gain_cands = df[
        (df['fdr_rescored'] <= xlim) &
        df[target_col]
    ]
    samples_gained_cumsum['count'] = samples_gained_cumsum.index.to_series().apply(
        lambda val: len(gain_cands[
            (gain_cands['fdr_standard'] > val) &
            (gain_cands['fdr_rescored'] <= val)
        ]) + len(rescored_below[
                (rescored_below['fdr_rescored'] <= val) &
                (~rescored_below[top_col]) & rescored_below[top_col_rescored]
        ])
    )

    # Calculate loss
    loss_cands = df[
        (df['fdr_standard'] <= xlim) &
        df[target_col]
    ]
    samples_lost_cumsum['count'] = samples_lost_cumsum.index.to_series().apply(
        lambda val: len(loss_cands[
            (
                (loss_cands['fdr_standard'] <= val) &
                (loss_cands['fdr_rescored'] > val)
            )
        ]) + len(standard_below[
                (standard_below['fdr_standard'] <= val) &
                standard_below[top_col] & (~standard_below[top_col_rescored])
        ])
    )

    if mismatch_col is not None:
        # Calculate rescore mismatch 
        rescore_mismatch_cand = df[
            (df['fdr_rescored'] <= xlim) &
            df[target_col] & df[mismatch_col]
        ]
        rescore_mismatch_cumsum['count'] = rescore_mismatch_cumsum.index.to_series().apply(
            lambda val: len(rescore_mismatch_cand[
                (rescore_mismatch_cand['fdr_rescored'] <= val)
            ])
        )

        # Calculate standard mismatch
        standard_mismatch_cand = df[
            (df['fdr_standard'] <= xlim) &
            df[target_col] & df[mismatch_col]
        ]
        standard_mismatch_cumsum['count'] = standard_mismatch_cumsum.index.to_series().apply(
            lambda val: len(standard_mismatch_cand[
                (
                    (standard_mismatch_cand['fdr_standard'] <= val)
                )
            ])
        )
    
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    data = {
        'standard': samples_standard_cumsum['count'],
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
    xfilter = samples_standard_cumsum.loc[:xlim].index
    
    ylim = max([
        samples_standard_cumsum.loc[xfilter, 'count'].max(),
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
        y_standard = samples_standard_cumsum['count'].iloc[x_idx]
        y_rescored = samples_rescored_cumsum['count'].iloc[x_idx]
        y = max([y_standard, y_rescored])
        text = f"{y_standard}→{y_rescored}"
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
            'standard mismatch': standard_mismatch_cumsum['count']/data['standard']*100,
            'rescore mismatch': rescore_mismatch_cumsum['count']/data['rescored']*100,
            'fdr': fdr_x*100,
        }
        ax2 = ax.twinx()
        
        sns.lineplot(
            ax=ax2,
            data=data2
        )
    
    return fig


def fdr_mismatch_plot(df,
                      xlim=0.05,
                      plot_fdr=True,
                      score_col='match_score',
                      rescore_col='xiMLScore',
                      fdr_steps=0.001,
                      target_col='isTT',
                      standard_fdr_col=None,
                      rescore_fdr_col=None,
                      mismatch_col='mismatch',
                      ylabel='proportion',
                      figsize=(12, 6),
                      dpi=100):
    """
    Plot mismatch proportion over FDR cutoff

    :param df: Input DF
    :param xlim: X plot limit
    :param plot_fdr: Optional FDR plotting
    :param score_col: standard score column
    :param rescore_col: Rescored score column
    :param fdr_steps: FDR step resolution
    :param target_col: Column indicating target-target CSMs
    :param standard_fdr_col: standard FDR
    :param rescore_fdr_col: Rescored FDR
    :param mismatch_col: Column indicating mismatched CSMs
    :param ylabel: Label for the y-axis
    :param figsize: Figure size
    :param dpi: Figure DPI
    :return: Tuple of the resulting figure and the underlying data
    """
    # Calculate FDR if none provided
    if standard_fdr_col is None:
        df.drop(['fdr_standard', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_standard'] = fdr.calculate_bi_fdr(df, score_col=score_col)
    else:
        df.loc[:, 'fdr_standard'] = df.loc[:, standard_fdr_col]

    if rescore_fdr_col is None:
        df.drop(['fdr_standard', 'fdr_rescored'], axis=1, inplace=True, errors='ignore')
        df.loc[:, 'fdr_rescored'] = fdr.calculate_bi_fdr(df, score_col=rescore_col)
    else:
        df.loc[:, 'fdr_rescored'] = df.loc[:, rescore_fdr_col]
    
    df = df.loc[
        (df['fdr_standard'] <= xlim) | (df['fdr_rescored'] <= xlim)
    ]

    fdr_x = np.arange(0, xlim+fdr_steps, fdr_steps)

    fdr_series = pd.DataFrame(index=fdr_x, data={'fdr': fdr_x})
    samples_standard_cumsum = pd.DataFrame(index=fdr_x)
    samples_rescored_cumsum = pd.DataFrame(index=fdr_x)
    rescore_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    standard_mismatch_cumsum = pd.DataFrame(index=fdr_x)
    
    # Calculate sample number for FDR cutoff
    standard_below = df.loc[df['fdr_standard'] <= xlim]
    samples_standard_cumsum['count'] = samples_standard_cumsum.index.to_series().apply(
        lambda val: len(standard_below[standard_below['fdr_standard'] <= val])
    )
    
    rescored_below = df.loc[df['fdr_rescored'] <= xlim]
    samples_rescored_cumsum['count'] = samples_rescored_cumsum.index.to_series().apply(
        lambda val: len(rescored_below.loc[rescored_below['fdr_rescored'] <= val])
    )
    
    # Calculate rescore mismatch 
    rescore_mismatch_cand = df.loc[
        (df['fdr_rescored'] <= xlim) &
        df[target_col] & df[mismatch_col]
    ]
    rescore_mismatch_cumsum.loc[:, 'count'] = rescore_mismatch_cumsum.index.to_series().apply(
        lambda val: len(rescore_mismatch_cand[
            (rescore_mismatch_cand['fdr_rescored'] <= val)
        ])
    )

    # Calculate standard mismatch
    standard_mismatch_cand = df.loc[
        (df['fdr_standard'] <= xlim) &
        df[target_col] & df[mismatch_col]
    ]
    standard_mismatch_cumsum.loc[:, 'count'] = standard_mismatch_cumsum.index.to_series().apply(
        lambda val: len(standard_mismatch_cand[
            (
                (standard_mismatch_cand['fdr_standard'] <= val)
            )
        ])
    )
    
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    data = {
        'standard mismatch': standard_mismatch_cumsum['count']/samples_standard_cumsum['count'],
        'rescore mismatch': rescore_mismatch_cumsum['count']/samples_rescored_cumsum['count'],
    }
    if plot_fdr:
        data['fdr'] = fdr_series['fdr']
    
    sns.lineplot(
        ax=ax,
        data=data
    )
    
    plt.xlabel('FDR cutoff')
    plt.ylabel(ylabel)
    
    ax.set_xlim(0, xlim)
    
    return fig, data
