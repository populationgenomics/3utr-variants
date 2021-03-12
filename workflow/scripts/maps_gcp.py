"""
Compute MAPS on UTR variants of gnomAD hail table in GCP
Functions imported from scripts/prepare_gnomad.py and scripts/maps_score.py

    Input:
        UTR interval bed file
        gnomAD hail table
        context hail table
        mutation hail table
    Output:
        MAPS hail table

Note:
The files utr3variants/annotate_gnomad.py and utr3variants/maps.py need to be copied to
the dataproc session when running this script.
See rules/evaluation.smk for a generic hailctl command.
"""

import hail as hl

# pylint: disable=E0401
from annotate_gnomad import filter_gnomad, annotate_for_maps, annotate_by_intervals
from maps import maps  # pylint: disable=E0401


def main(args):
    """
    Subset and annotate gnomAD hail table variants by intervals and compute MAPS score
    """

    hl.init(default_reference=args.genome_assembly)

    intervals = hl.import_bed(args.intervals)
    mutation_ht = hl.read_table(args.mutation_ht)
    context_ht = hl.read_table(args.context_ht)
    ht = hl.read_table(args.gnomAD_ht)

    print('Filter')
    subset_interval = hl.parse_locus_interval(args.chr_subset)
    ht = filter_gnomad(ht, [subset_interval])

    print('Annotate')
    ht = annotate_for_maps(ht, context_ht, mutation_ht)
    ht = annotate_by_intervals(ht, intervals, new_column='UTR_group')

    print('MAPS score')
    maps_ht = maps(ht, grouping=['UTR_group'])
    if args.verbose:
        maps_ht.show()

    print('save...')
    maps_ht.export(args.output)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute MAPS score')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    parser.add_argument(
        '--intervals', required=True, help='Interval text file for subsetting gnomAD'
    )
    parser.add_argument(
        '--gnomAD_ht',
        required=True,
        help='Google cloud URL or directory path to gnomAD hail table',
    )
    parser.add_argument(
        '--context_ht',
        required=True,
        help='Google cloud URL or directory path to gnomAD context hail table.',
    )
    parser.add_argument(
        '--mutation_ht',
        required=True,
        help='Google cloud URL or directory path to gnomAD mutation hail table',
    )
    parser.add_argument(
        '--genome_assembly',
        required=True,
        help='Genome assembly identifier e.g. GRCh37',
    )
    parser.add_argument(
        '--chr_subset',
        required=True,
        help='Chromosome region to subset variants to',
    )
    parser.add_argument(
        '--verbose',
        required=False,
        action='store_true',
        help='Show more output',
    )
    main(args=parser.parse_args())
