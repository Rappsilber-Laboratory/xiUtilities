import pandas as pd
import multiprocess as mp
from math import ceil


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
