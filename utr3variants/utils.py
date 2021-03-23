"""
Basic functions
"""
import pybedtools
import pandas as pd


def extract_annotations(
    df: pd.DataFrame,
    annotation_string: str,
    annotations_columns: list,
    sep: str = '|',
    database: str = None,
    feature: str = None,
):
    """
    Extract annotations and annotate database and feature
    :param df: DataFrame of features to extract annotations to
    :param annotation_string: column in df of concatenated annotation strings
    :param annotations_columns: list of annotation columns to extract to, number of
        elements must match number of entries df[name_column]
    :param sep: separator used in df[name_column] for splitting annotation string to
        annotation columns
    :param database: database of features
    :param feature: feature type e.g. PAS, hexamer
    :return: df with new columns for each annotation, database and feature
    """
    df[annotations_columns] = df[annotation_string].str.split(sep, expand=True)
    df['database'] = database
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
