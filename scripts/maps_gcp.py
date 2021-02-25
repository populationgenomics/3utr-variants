"""
Compute MAPS on UTR variants of gnomAD hail table in GCP
Functions imported from scripts/prepare_gnomad.py and scripts/maps_score.py

When run as script, subset and annotate gnomAD hail table and compute MAPS score
    Input:
        UTR interval bed file
        gnomAD hail table
        context hail table
        mutation hail table
    Output:
        MAPS hail table
The files scripts/prepare_gnomad.py and scripts/maps_score.py need to be copied to
dataproc when running this script.
See rules/evaluation.smk for a generic hailctl command.
"""

import hail as hl

from annotate_gnomad import filter_gnomad, annotate_for_maps  # pylint: disable=E0401
from maps import maps  # pylint: disable=E0401


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
        '--log',
        required=True,
        help='GCP link for log output in bucket',
    )
    args = parser.parse_args()

    hl.init(default_reference=args.genome_assembly)

    # TODO: include annotation
    # intervals = hl.import_bed(bed_file, reference_genome=reference_genome')
    with open(args.intervals, 'r') as f:
        intervals = f.readlines()
    intervals = [
        hl.parse_locus_interval(x, reference_genome=args.genome_assembly)
        for x in intervals
    ]

    mutation_ht = hl.read_table(args.mutation_ht)
    context_ht = hl.read_table(args.context_ht)

    ht = hl.read_table(args.gnomAD_ht)
    ht = filter_gnomad(ht, intervals)
    ht = annotate_for_maps(ht, context_ht)

    print('MAPS score')
    # TODO: use different grouping
    maps_ht = maps(ht, mutation_ht, additional_grouping=['protein_coding'])
    maps_ht.show()

    print('save...')
    maps_ht.export(args.output)
    hl.copy_log(args.log)
