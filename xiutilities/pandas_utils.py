import uuid
from collections.abc import Iterable
from math import ceil

import numpy as np
import pandas as pd
import multiprocess as mp


def async_apply(df, f, target_size=10_000, max_cpu=1000, **kwargs):
    """
    Use mutltiprocessing to apply ``f()`` to DataFrame ``df``.

    :param df: Input DataFrame
    :param f: Function to use in the ``apply()``
    :param target_size: Target size for slices handled in parallel
    :param max_cpu: Maximum number of cores to use
    :param kwargs: Keyword arguments forwarded to the ``df.apply()``
    :return: Result of the ``apply()``
    """
    n_slices_target = ceil(len(df) / target_size)
    n_slices = min([
        mp.cpu_count(),
        n_slices_target,
        max_cpu
    ])
    slice_size = ceil(len(df) / n_slices)
    slices = [
        df.iloc[i * slice_size: (i + 1) * slice_size]
        for i in range(n_slices)
    ]

    with mp.Pool(n_slices) as pool:
        def job(x):
            return x.apply(f, **kwargs)

        results = pool.map(job, slices)

    return pd.concat(results).copy()


def swap_columns(df, mask, column_pairs: list[tuple]) -> pd.DataFrame:
    """
    Swap column pairs where ``mask`` is ``True``.

    :param df: Input DataFrame
    :type df: DataFrame
    :param mask: Mask where to swap the columns
    :type mask: Series|ndarray|list
    :param column_pairs: List of column pairs
    :type column_pairs: list[tuple]
    :return: Input DataFrame with swapped columns
    :rtype: DataFrame
    """
    for c1, c2 in column_pairs:
        df_swap = df.rename({
            c1: c2,
            c2: c1,
        }, axis=1)
        df.loc[mask, [c1, c2]] = df_swap.loc[mask, [c1, c2]]
    return df.copy()


def serialize_columns(df):
    nan_filter = None
    for c in df.columns:
        if nan_filter is None:
            nan_filter = ~df[c].isna()
        else:
            nan_filter &= ~df[c].isna()

    # Try to find row without NaNs
    type_row = df[nan_filter]

    if len(type_row) > 0:
        # If some found take the first one
        type_row = type_row.iloc[:1]
    else:
        # If none found take the first row od the DF
        type_row = df.iloc[:1]

    for c in df.columns:
        col_row = type_row
        # Check if type_row has a valid value for column
        if col_row.iloc[0][c] is None:
            # If not search for a row with a non-NaN value
            col_row = df[~df[c].isna()]
            if len(col_row) == 0:
                # Only NaN values for column
                continue
        # Take the non-NaN value found
        col_row = col_row.iloc[0]
        if type(col_row[c]) is uuid.UUID:
            df[df[c].isna()] = ''
            df[c] = df[c].astype(str, errors='ignore')
        elif type(col_row[c]) is str:
            pass
        elif isinstance(col_row[c], Iterable):
            df[df[c].isna()] = ''
            df[c] = df[c].apply(
                lambda x: None if x is None else ';'.join(np.array(x).astype(str))
            )
        elif np.issubdtype(type(col_row[c]), np.generic):
            pass
        else:
            df[df[c].isna()] = ''
            df[c] = df[c].astype(str, errors='ignore')
    df.columns = df.columns.map(str)
    return df
