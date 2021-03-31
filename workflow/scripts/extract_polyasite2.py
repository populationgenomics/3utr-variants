"""
Extract PolyA sites and hexamers from PolyASite 2.0

Input:
    PolyASite 2.0 database file
Output:
    TSV file of PAS and hexamer intervals with all annotations as columns
"""

from collections import OrderedDict
import numpy as np
import pandas as pd
import pybedtools
from utr3variants.utils import (
    encode_annotations,
    extract_annotations,
    convert_chromosome,
    get_most_expressed,
)

COLUMN_TYPES = OrderedDict(
    chrom=str,
    start=int,
    end=int,
    cluster_id=str,
    score=float,
    strand=str,
    percent_expressed=float,
    n_protocols=int,
    expression=float,
    cluster_annotation=str,
    hexamer_motif=str,
)

MERGED_COLUMN_TYPES = OrderedDict(
    chrom=str,
    start=int,
    end=int,
    name=str,
    score=float,
    strand=str,
    gchrom=str,
    gstart=int,
    gend=int,
    gname=str,
    gscore=str,
    gstrand=str,
    g13=str,
    g14=str,
    g15=str,
    g16=str,
    g17=str,
    g18=str,
    g19=str,
)

CANONICAL_CHROMOSOMES = list(range(1, 22)) + ['X', 'Y']

if __name__ == '__main__':
    genome = snakemake.config['assembly_ucsc']
    annotations = snakemake.params.annotations
    pas_db = snakemake.input.db.__str__()
    genes_file = snakemake.input.genes.__str__()
    output = snakemake.output

    df = pd.read_csv(pas_db, sep='\t', names=COLUMN_TYPES.keys(), dtype=COLUMN_TYPES)

    # match chromosome style to genome annotation
    df = df[df['chrom'].isin(CANONICAL_CHROMOSOMES)]
    df['chrom'] = df.chrom.apply(lambda x: convert_chromosome(x, style='chr'))

    # put annotation into name
    df = encode_annotations(df, annotations, 'name')

    print('Create BedTools object')
    bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
    bed = pybedtools.BedTool.from_dataframe(df[bed_cols])

    # add genes
    print('Annotate overlapping genes')
    genes_bed = pybedtools.BedTool(genes_file).sort()
    df_genes = (
        bed.sort()
        .closest(genes_bed, s=True, fu=True, D='ref')
        .to_dataframe(names=MERGED_COLUMN_TYPES.keys(), dtype=MERGED_COLUMN_TYPES)
    )

    # annotate most expressed
    print('Annotate most expressed PAS')
    df_genes = get_most_expressed(
        df_genes,
        aggregate_column='gname',
        expression_column='score',
        interval_columns=['chrom', 'start', 'end', 'strand'],
        new_column='most_expressed',
    )
    annotations.append('most_expressed')
    df_genes['name'] = df_genes['name'] + '|' + df_genes['most_expressed'].map(str)

    # drop unnecessary columns
    df_genes.drop(
        columns=[x for x in df_genes.columns if x.startswith('g')], inplace=True
    )

    intervals_df = extract_annotations(
        df_genes,
        annotation_string='name',
        annotation_columns=annotations,
        database='PolyASite2',
        feature='PAS',
    )

    if 'hexamer_motif' in annotations:
        print('Extract hexamers')
        hexamers_df = intervals_df.copy().replace('nan', np.nan, regex=True)
        hexamers_df.dropna(subset=['hexamer_motif'], inplace=True)

        # split by different signals per PAS by ;
        hexamers_df['hexamer_motif_split'] = hexamers_df['hexamer_motif'].str.split(';')
        hexamers_df = hexamers_df.explode('hexamer_motif_split')

        # split each signal by @ into different columns
        hexamers_df[['hexamer_motif', 'rel_pos', 'hex_start']] = hexamers_df[
            'hexamer_motif_split'
        ].str.split('@', expand=True)

        # retrieve hexamer interval
        hexamers_df.hex_start = hexamers_df.hex_start.astype(int)
        hexamers_df.start = (hexamers_df.hex_start - 1).where(
            hexamers_df.strand == '+', hexamers_df.hex_start - 6
        )
        hexamers_df.end = hexamers_df.hex_start.where(
            hexamers_df.strand == '-', hexamers_df.hex_start - 1 + 6
        )

        # remove unnecessary columns
        hexamers_df.drop(
            columns=['hexamer_motif_split', 'rel_pos', 'hex_start'], inplace=True
        )

        # append hexamer regions translated into BED coordinates
        hexamers_df['database'] = 'PolyASite2'
        hexamers_df['feature'] = 'hexamer'
        intervals_df['hexamer_motif'] = np.nan
        intervals_df = pd.concat([intervals_df, hexamers_df])

    print('save...')
    # get stats
    intervals_df['feature'].value_counts().to_csv(output.stats, sep='\t')
    intervals_df.to_csv(output.intervals, index=False, sep='\t')
