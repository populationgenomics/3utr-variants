"""
Extract PolyA sites and hexamers from PolyASite 2.0

Input:
    PolyASite 2.0 database file
Output:
    TSV file of PAS and hexamer intervals with all annotations as columns
"""

from collections import OrderedDict
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


if __name__ == '__main__':
    genome = snakemake.config['assembly_ucsc']
    annotations = snakemake.params.annotations
    pas_db = snakemake.input.db.__str__()
    genes_file = snakemake.input.genes.__str__()
    output = snakemake.output

    df_db = pd.read_csv(pas_db, sep='\t', names=COLUMN_TYPES.keys(), dtype=COLUMN_TYPES)

    # match chromosome style to genome annotation
    df_db['chrom'] = df_db.chrom.apply(lambda x: convert_chromosome(x, style='chr'))

    if 'most_expressed' not in annotations:
        df_db = encode_annotations(df_db, annotations, 'name')
    else:
        df_db = encode_annotations(
            df_db,
            annotation_columns=[x for x in annotations if x != 'most_expressed'],
            encode_column='name',
        )

        bed_cols = ['chrom', 'start', 'end', 'name', 'score', 'strand']
        bed = pybedtools.BedTool.from_dataframe(df_db[bed_cols])

        print('Annotate overlapping genes')
        genes_bed = pybedtools.BedTool(genes_file).sort()
        df_db = (
            bed.sort()
            .closest(genes_bed, s=True, fu=True, D='ref')
            .to_dataframe(names=MERGED_COLUMN_TYPES.keys(), dtype=MERGED_COLUMN_TYPES)
        )

        print('Annotate most expressed PAS')
        df_db = get_most_expressed(
            df_db,
            aggregate_column='gname',
            expression_column='score',
            interval_columns=['chrom', 'start', 'end', 'strand'],
            new_column='most_expressed',
        )
        df_db['name'] = df_db['name'] + '|' + df_db['most_expressed'].map(str)

        # drop unnecessary columns
        df_db.drop(
            columns=[x for x in df_db.columns if x.startswith('g')], inplace=True
        )

    intervals_df = extract_annotations(
        df_db,
        annotation_string='name',
        annotation_columns=annotations,
        database='PolyASite2',
        feature='PAS',
    )

    if 'hexamer_motif' in annotations:
        print('Extract hexamers')
        hexamers_df = intervals_df[intervals_df['hexamer_motif'] != ''].copy()

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
        intervals_df['hexamer_motif'] = ''
        intervals_df = pd.concat([intervals_df, hexamers_df])

    print('save...')
    intervals_df['feature'].value_counts().to_csv(output.stats, sep='\t')  # get stats
    intervals_df.to_csv(output.intervals, index=False, sep='\t')
