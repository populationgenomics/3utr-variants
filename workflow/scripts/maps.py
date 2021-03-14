"""
Compute MAPS score
    Input:
        TSV file with variant counts, singleton counts and expected singletons
    Output:
        MAPS table as TSV
"""
from typing import Iterable
import pandas as pd


def aggregate_maps(df: pd.DataFrame, grouping: Iterable[str]) -> pd.DataFrame:
    """
    Pandas implementation of aggregated MAPS score
    :param df: pandas Dataframe containing 'singleton_count', 'variant_count',
        'expected_singletons' and group columns
        Converted output of utr3variants.maps.count_for_maps()
    :param grouping: columns from count_df to group by
    """
    df = (
        df.groupby(grouping)
        .agg(
            dict(
                variant_count=sum,
                singleton_count=sum,
                expected_singletons=sum,
            )
        )
        .reset_index()
    )

    df['ps'] = df.singleton_count / df.variant_count
    df['maps'] = (df.singleton_count - df.expected_singletons) / df.variant_count
    df['maps_sem'] = (df.ps * (1 - df.ps) / df.variant_count) ** 0.5
    return df


if __name__ == '__main__':
    count_df = pd.read_table(snakemake.input.counts.__str__())
    maps_df = aggregate_maps(count_df, grouping=snakemake.wildcards.variant_subset)
    maps_df.to_csv(snakemake.output.maps)
