"""
Merge different UTR intervals and annotate
"""
import tempfile
from functools import reduce
import pandas as pd
from pandas.api.types import CategoricalDtype
import pybedtools
from utr3variants.utils import (
    encode_annotations,
    extract_annotations,
    convert_chromosome,
)

BED_COLS = ['chrom', 'start', 'end', 'name', 'score', 'strand']
CANONICAL_CHROMOSOMES = [str(x) for x in list(range(1, 23))] + ['X', 'Y']
CANONICAL_CHROMOSOMES_CHR = [f'chr{chr}' for chr in CANONICAL_CHROMOSOMES]

chr_dtype = CategoricalDtype(CANONICAL_CHROMOSOMES, ordered=True)


def rename_interval(
    interval: pybedtools.Interval,
    prefix: str = None,
    name: str = None,
    sep='|',
    keep_col=-1,
) -> pybedtools.Interval:
    """
    Change name of interval by adding a prefix and/or overwriting the old name

    :param interval:
    :param prefix: prefix to add to existing name if specified
    :param name: value to overwrite old name, depending on value of keep_col
    :param sep: character to separate prefix from rest of name
    :param keep_col: column index of annotation column to keep in interval.name list
        separated by sep. Overwrite all annotations with name if keep_col=-1 (default)
        e.g. for interval.name = 'AAUAAA|No' and keep_col = 0: keep 'AAUAAA', drop 'No'
    """
    if name is None:
        name = interval.name
        if keep_col >= 0:
            name = name.split(sep)[keep_col]
    if prefix is not None:
        name = sep.join(filter(None, [prefix, name]))
    interval.name = name
    return interval


if __name__ == '__main__':
    chr_style = snakemake.params.chr_style_gnomAD
    annotations = snakemake.params.annotations
    chainfile = snakemake.input.chainfile.__str__()

    utrs = pd.read_table(
        snakemake.input.utrs.__str__(), sep='\t', dtype=str, names=BED_COLS
    )  # pybedtools.BedTool(snakemake.input.utrs.__str__())
    polya_db = pd.read_table(snakemake.input.PolyA_DB.__str__(), sep='\t', dtype=str)
    polya_site = pd.read_table(
        snakemake.input.PolyASite2.__str__(), sep='\t', dtype=str
    )

    # liftover PolyASite2 to hg19
    pasite_anno = [x for x in annotations if x in polya_site.columns]
    polya_site = encode_annotations(polya_site, pasite_anno, 'name')
    polya_site_bed = pybedtools.BedTool.from_dataframe(polya_site[BED_COLS]).liftover(
        chainfile
    )
    polya_site = extract_annotations(
        polya_site_bed.to_dataframe(dtype=str), 'name', pasite_anno
    )
    polya_site = polya_site[polya_site.chrom.isin(CANONICAL_CHROMOSOMES_CHR)]

    # convert chromosome format to match gnomAD dataset
    print('Convert chromosome styles')
    utrs.chrom = utrs.chrom.apply(lambda x: convert_chromosome(x, chr_style))
    polya_db.chrom = polya_db.chrom.apply(lambda x: convert_chromosome(x, chr_style))
    polya_site.chrom = polya_site.chrom.apply(
        lambda x: convert_chromosome(x, chr_style)
    )

    # subtract from 3'UTR
    utrs_bed = pybedtools.BedTool.from_dataframe(utrs)
    polya_db_bed = pybedtools.BedTool.from_dataframe(polya_db)
    polya_site_bed = pybedtools.BedTool.from_dataframe(polya_site)
    with tempfile.TemporaryDirectory() as tmp_dir:
        other_utr = (
            utrs_bed.subtract(polya_db_bed)  # pylint: disable=too-many-function-args
            .subtract(polya_site_bed)
            .saveas(f'{tmp_dir}/other_utr.bed')
            .to_dataframe(dtype=str)
        )
        other_utr['database'] = 'GENCODE'
        other_utr['feature'] = '3UTR'

    # merge all intervals via pandas
    df = reduce(
        lambda df_left, df_right: pd.merge(
            df_left,
            df_right,
            how='outer',
            on=list(set(df_left.columns) & set(df_right.columns)),
        ),
        [other_utr, polya_db, polya_site],
    )

    # sort intervals
    df['chrom'] = df['chrom'].astype(str).astype(chr_dtype)
    df.sort_values(by=['chrom', 'start', 'end', 'strand'], inplace=True)

    print('Create hail-parsable locus intervals')
    df['locus_interval'] = (
        '('
        + df.chrom.astype(str)
        + ':'
        + df.start.map(str)
        + '-'
        + df.end.map(str)
        + ']'
    )

    # Encode annotations
    df = encode_annotations(df, annotations, 'name', verbose=True)

    # save
    print('save...')
    df[['locus_interval', 'strand'] + annotations].to_csv(
        snakemake.output.intervals, sep='\t', index=False
    )
    df[['chrom', 'start', 'end', 'name', 'score', 'strand']].to_csv(
        snakemake.output.bed, sep='\t', index=False, header=False
    )
