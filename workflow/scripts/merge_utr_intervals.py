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
    name: str = '',
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
    if keep_col >= 0:
        name = interval.name.split(sep)[keep_col]
    if prefix is not None:
        name = sep.join(filter(None, [prefix, name]))
    interval.name = name
    return interval


if __name__ == '__main__':
    utrs = pybedtools.BedTool(snakemake.input.utrs.__str__())
    hexamers = pybedtools.BedTool(snakemake.input.hexamers.__str__())
    pas = pybedtools.BedTool(snakemake.input.pas.__str__())

    annotation_type = snakemake.wildcards.annotation

    with tempfile.TemporaryDirectory() as tmp_dir:
        hexamers = (
            hexamers.subtract(pas)  # pylint: disable=too-many-function-args
            .each(
                rename_interval,
                name='hexamer',
                keep_col=ANNOTATION_MAP[annotation_type],
            )
            .saveas(f'{tmp_dir}/hexamer.bed')
        )
        pas = pas.each(
            rename_interval,
            name='PAS',
            keep_col=ANNOTATION_MAP[annotation_type],
        ).saveas(f'{tmp_dir}/pas.bed')
        # concatenate regions
        bed = pas.cat(hexamers, postmerge=False).sort()

        if annotation_type == 'overall':  # include other UTRs
            other_utr = (
                utrs.subtract(hexamers)  # pylint: disable=too-many-function-args
                .subtract(pas)
                .each(rename_interval, name='other_UTR')
                .saveas(f'{tmp_dir}/other_utr.bed')
            )
            bed = bed.cat(other_utr, postmerge=False).sort()

    # convert chromosome format to match gnomAD dataset
    bed = bed.each(convert_chr_style, chr_style=snakemake.params['chr_style_gnomAD'])
    bed.saveas(snakemake.output.intervals)
