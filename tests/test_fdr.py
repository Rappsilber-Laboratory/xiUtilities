import xiutilities.bi_fdr
import pandas as pd
import numpy as np
import random
import string


def test_mono_fdr():
    tt_score_offset = 1_000
    df = pd.DataFrame(columns=['decoy_class', 'match_score'])
    dd_scores = np.array(range(100_000))
    tt_scores = dd_scores.copy()
    td_scores = np.concatenate([dd_scores, dd_scores])
    tt_scores += tt_score_offset
    df = pd.concat([
        df,
        pd.DataFrame({
            'match_score': dd_scores
        }).assign(
            decoy_class='DD'
        )
    ])
    df = pd.concat([
        df,
        pd.DataFrame({
            'match_score': tt_scores
        }).assign(
            decoy_class='TT'
        )
    ])
    df = pd.concat([
        df,
        pd.DataFrame({
            'match_score': td_scores
        }).assign(
            decoy_class='TD'
        )
    ]).reset_index(drop=True)
    df = df.sort_values(['decoy_class']).copy()  # Sort by decoy_class to avoid one-off error
    df['fdr'] = xiutilities.bi_fdr.calculate_fdr(
        df,
        decoy_class='decoy_class'
    )

    assert len(
        df[
            df['fdr'] == 0
        ]
    ) == tt_score_offset
