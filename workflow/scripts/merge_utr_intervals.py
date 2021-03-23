"""
Merge different UTR intervals and annotate
"""
import tempfile
from functools import reduce
import pandas as pd
import pybedtools
from utr3variants.utils import convert_chromosome_bed, convert_chromosome


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
    utrs = pybedtools.BedTool(snakemake.input.utrs.__str__())
    polya_db = pd.read_table(snakemake.input.PolyA_DB.__str__(), sep='\t')

    with tempfile.TemporaryDirectory() as tmp_dir:
        # convert chromosome format to match gnomAD dataset
        utrs = utrs.each(convert_chromosome_bed, chr_style)
        polya_db.chrom = polya_db.chrom.apply(
            lambda x: convert_chromosome(x, chr_style)
        )

        # subtract from 3'UTR
        polya_db_bed = pybedtools.BedTool.from_dataframe(polya_db)
        other_utr = (
            utrs.subtract(  # pylint: disable=too-many-function-args
                polya_db_bed
            ).saveas(f'{tmp_dir}/other_utr.bed')
            # convert other 3'UTR to dataframe
            .to_dataframe()
        )
        other_utr['database'] = 'GENCODE'
        other_utr['feature'] = '3UTR'

    # merge all intervals via pandas
    common_cols = list(set(other_utr.columns) & set(polya_db.columns))

    df = reduce(
        lambda df_left, df_right: pd.merge(
            df_left, df_right, how='outer', on=common_cols
        ),
        [other_utr, polya_db],
    )

    # convert to hail-parsable locus interval
    df['locus_interval'] = (
        '(' + df.chrom.map(str) + ':' + df.start.map(str) + '-' + df.end.map(str) + ']'
    )

    # remove redundant columns
    df = df[['locus_interval', 'strand', 'database', 'feature'] + annotations]

    # save
    df.to_csv(snakemake.output.intervals, sep='\t', index=False)
