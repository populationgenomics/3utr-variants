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
from utr3variants.utils import extract_annotations, convert_chromosome

column_types = OrderedDict(
    chrom=str,
    start=int,
    end=int,
    cluster_id=str,
    score=float,
    strand=str,
    percent_expressed=float,
    n_protocols=int,
    avg_TPM=float,
    cluster_annotation=str,
    hexamer_motif=str,
)

if __name__ == '__main__':
    genome = snakemake.config['assembly_ucsc']
    annotations = snakemake.params.annotations
    pas_db = snakemake.input.db.__str__()
    output = snakemake.output

    print('read database')
    df = pd.read_csv(pas_db, sep='\t', names=column_types.keys(), dtype=column_types)

    # match chromosome style to genome annotation
    df['chrom'] = df.chrom.apply(lambda x: convert_chromosome(x, style='chr'))

    # put annotation into name
    df['name'] = df[annotations].astype(str).agg('|'.join, axis=1)

    # create bedtools object/table
    print('Create BedTools object')
    # BedTool object needed for hexamer extraction
    bed_cols = [
        'chrom',  # chrom
        'start',  # start
        'end',  # end
        'name',  # name
        'score',  # score
        'strand',  # strand
    ]
    bed = pybedtools.BedTool.from_dataframe(df[bed_cols])

    # TODO: add genes & most expressed by gene
    # bed.closest(genes_bed, s=True, fu=True, D=True)

    # convert back to dataframe with annotations
    intervals_df = extract_annotations(
        bed.to_dataframe(),
        annotation_string='name',
        annotations_columns=annotations,
        database='PolyASite2',
        feature='PAS',
    )

    if 'hexamer_motif' in annotations:
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
