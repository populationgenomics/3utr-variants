"""GCP script for MAPS from complete gnomAD set"""

import hail as hl
import gnomad.utils.vep

from prepare_gnomAD import get_worst_consequence_with_non_coding, prepare_ht
from maps_score import maps

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute MAPS score')
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('--intervals', required=True)
    parser.add_argument('--gnomAD_ht', required=True)
    parser.add_argument('--context_ht', required=True)
    parser.add_argument('--mutation_ht', required=True)
    parser.add_argument('--genome_assembly', required=True)
    args = parser.parse_args()

    hl.init(default_reference=args.genome_assembly)

    # intervals = hl.import_bed(bed_file, reference_genome=reference_genome')
    with open(args.intervals, 'r') as f:
        intervals = f.readlines()
    # snp_ht = snp_ht.filter(hl.is_defined(intervals[snp_ht.locus]))

    mutation_ht = hl.read_table(args.mutation_ht)

    # subset context table
    context_ht = hl.read_table(args.context_ht)
    context_ht = hl.filter_intervals(context_ht, [intervals])

    # prepare gnomAD table
    ht = hl.read_table(args.gnomAD_ht)
    ht = hl.filter_intervals(ht, [intervals])
    ht = ht.filter(hl.len(ht.filters) == 0)
    ht = gnomad.utils.vep.filter_vep_to_canonical_transcripts(ht)
    print(f'entries in gnomAD after filtering: {ht.count()}')

    ht = get_worst_consequence_with_non_coding(ht)
    context = context_ht[ht.key]
    snp_ht = prepare_ht(
        ht=ht.annotate(context=context.context, methylation=context.methylation),
        trimer=True,
    )

    print('MAPS score')
    maps_ht = maps(snp_ht, mutation_ht, additional_grouping=['protein_coding'])

    print('save...')
    maps_ht.write(args.output)
