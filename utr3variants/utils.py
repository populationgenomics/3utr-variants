"""
Basic functions
"""
import pybedtools
import pandas as pd


def get_most_expressed(
    df,
    aggregate_column,
    expression_column,
    interval_columns=None,
    new_column='most_expressed',
):
    """
    Determine and annotate column for most expressed interval feature
    :param df: interval dataframe
    :param aggregate_column: column to aggregate by (e.g. gene)
    :param expression_column: expression column to maximise over
    :param interval_columns: columns used to define interval
    :param new_column: column for truth values of most expressed features
    """
    if interval_columns is None:
        interval_columns = ['chrom', 'start', 'end', 'strand']

    # remove duplicate intervals, keep most expressed
    df.sort_values(
        by=interval_columns + [aggregate_column, expression_column], inplace=True
    )
    df.drop_duplicates(subset=interval_columns, keep='first', inplace=True)

    # infer most-used intervals
    df[new_column] = False
    df.loc[
        df.groupby([aggregate_column])[expression_column].idxmax(), new_column
    ] = True

    return df


def encode_annotations(
    df, annotation_columns, encode_column, sep='|', fillna='', verbose=True
):
    """
    Encode annotation columns into a single column of concatenated strings
    """
    if verbose:
        print('Encode annotations...')
    df[encode_column] = (
        df[annotation_columns].fillna(fillna).astype(str).agg(sep.join, axis=1)
    )
    return df


def extract_annotations(
    df: pd.DataFrame,
    annotation_string: str,
    annotation_columns: list,
    sep: str = '|',
    database: str = None,
    feature: str = None,
    verbose: bool = True,
):
    """
    Extract annotations and annotate database and feature
    :param df: DataFrame of features to extract annotations to
    :param annotation_string: column in df of concatenated annotation strings
    :param annotation_columns: list of annotation columns to extract to, number of
        elements must match number of entries df[name_column]
    :param sep: separator used in df[name_column] for splitting annotation string to
        annotation columns
    :param database: database of features
    :param feature: feature type e.g. PAS, hexamer
    :return: df with new columns for each annotation, database and feature
    """
    if verbose:
        print('Extract annotations...')
    df[annotation_columns] = df[annotation_string].str.split(sep, expand=True)
    if database is not None and database not in df.columns:
        df['database'] = database
    if feature is not None and feature not in df.columns:
        df['feature'] = feature
    return df


def convert_chromosome_bed(
    interval: pybedtools.Interval, style: str
) -> pybedtools.Interval:
    """
    Convert chromosome style of interval to match chr_style
    """
    interval.chrom = convert_chromosome(interval.chrom, style)
    return interval


def convert_chromosome(chromosome, style) -> str:
    """
    Convert chromosome style of interval to match chr_style
    """
    if style == '' and 'chr' in chromosome:
        chromosome = chromosome.replace('chr', '')
    elif style == 'chr' and 'chr' not in chromosome:
        chromosome = f'chr{chromosome}'
    return chromosome
