import xiutilities.bi_fdr
from xiutilities import pandas_utils
import pandas as pd
import numpy as np
import random
import string


def test_async_apply():
    random.seed(0)
    df: pd.DataFrame = pd.DataFrame(index=range(100_000))
    df['random_string1'] = df.index.to_series().apply(
        lambda _: ''.join(random.choices(string.ascii_letters, k=20))
    )
    df['random_string2'] = df.index.to_series().apply(
        lambda _: ''.join(random.choices(string.ascii_letters, k=20))
    )
    df['random_string_lower'] = pandas_utils.async_apply(
        df=df,
        f=lambda r: r['random_string1'].lower() + r['random_string2'].lower(),
        axis=1
    )
    assert all(
        np.equal(
            (df['random_string1'].str.len() + df['random_string2'].str.len()),
            df['random_string_lower'].str.len()
        )
    )
    assert all(
        np.equal(
            df['random_string_lower'].str.lower(),
            df['random_string_lower']
        )
    )


def test_swap_columns():
    random.seed(0)
    df: pd.DataFrame = pd.DataFrame(index=range(100_000))
    df['random_string1'] = df.index.to_series().apply(
        lambda _: ''.join(random.choices(string.ascii_lowercase, k=20))
    )
    df['random_string2'] = df.index.to_series().apply(
        lambda _: ''.join(random.choices(string.ascii_lowercase, k=20))
    )
    df['swap_mask'] = np.greater(
        df['random_string1'],
        df['random_string2']
    )
    df = pandas_utils.swap_columns(
        df=df,
        mask=df['swap_mask'],
        column_pairs=[['random_string1', 'random_string2']]
    )
    assert all(
        np.less_equal(
            df['random_string1'],
            df['random_string2']
        )
    )


def test_mono_fdr():
    tt_score_offset = 1_000
    df = pd.DataFrame(columns=['decoy_class', 'match_score'])
    dd_scores = np.array(range(100_000))
    tt_scores = dd_scores + tt_score_offset
    td_scores = np.concatenate([dd_scores, dd_scores])
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
    ]).reset_index(drop=True).copy()
    df['fdr'] = xiutilities.bi_fdr.calculate_fdr(
        df,
        decoy_class='decoy_class'
    )

    assert len(
        df[
            df['fdr'] == 0
        ]
    ) >= tt_score_offset  # FIXME
