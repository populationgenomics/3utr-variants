"""
Basic functions
"""
import pybedtools
import hail as hl


def convert_chr_style(
    interval: pybedtools.Interval, chr_style: str
) -> pybedtools.Interval:
    """
    Convert chromosome style of interval to match chr_style
    """
    if chr_style == '' and 'chr' in interval.chrom:
        interval.chrom = interval.chrom.replace('chr', '')
    elif chr_style == 'chr' and 'chr' not in interval.chrom:
        interval.chrom = f'chr{interval.chrom}'
    return interval


def import_interval_table(paths, interval_field, **kwargs):
    """
    Import table containing a locus interval field
    :param paths: paths to be passed to hail.import_table
    :param interval_field: name of interval field to be parsed by
        hail.parse_locus_interval and keyed by
    """
    ht = hl.import_table(paths, **kwargs)
    return ht.transmute(
        **{interval_field: hl.parse_locus_interval(ht[interval_field])}
    ).key_by(interval_field)
