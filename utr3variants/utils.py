"""
Basic functions
"""
import pybedtools


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
