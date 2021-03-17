"""
Merge different UTR intervals and annotate
"""
import tempfile
import pybedtools
from utr3variants.utils import convert_chr_style

# map annotation to column in |-separated annotation entry
ANNOTATION_MAP = {'overall': -1, 'hexamer': 0, 'conserved': 1}


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
    utrs = pybedtools.BedTool(snakemake.input.utrs.__str__())
    hexamers = pybedtools.BedTool(snakemake.input.hexamers.__str__())
    pas = pybedtools.BedTool(snakemake.input.pas.__str__())

    with tempfile.TemporaryDirectory() as tmp_dir:
        other_utr = (
            utrs.subtract(hexamers)  # pylint: disable=too-many-function-args
            .subtract(pas)
            .each(rename_interval, name='other_UTR|other_UTR|other_UTR')
            .saveas(f'{tmp_dir}/other_utr.bed')
        )
        hexamers = hexamers.each(
            rename_interval,
            prefix='hexamer',
        ).saveas(f'{tmp_dir}/hexamer.bed')
        pas = pas.each(
            rename_interval,
            prefix='PAS',
        ).saveas(f'{tmp_dir}/pas.bed')

        # concatenate regions
        bed = other_utr.cat(hexamers, pas, postmerge=False).sort()

        # convert chromosome format to match gnomAD dataset
        bed = bed.each(
            convert_chr_style, chr_style=snakemake.params['chr_style_gnomAD']
        )

        df = bed.saveas(f'{tmp_dir}/merged.bed').to_dataframe()

    # split annotation columns
    annotations = snakemake.params.annotations
    df[annotations] = df.name.str.split('|', expand=True)

    # convert to hail-parsable locus interval
    df['locus_interval'] = (
        '(' + df.chrom.map(str) + ':' + df.start.map(str) + '-' + df.end.map(str) + ']'
    )

    # remove redundant columns
    df = df[['locus_interval', 'strand'] + annotations]

    # save
    print(df.head())
    df.to_csv(snakemake.output.intervals, sep='\t', index=False)
