"""
Compute MAPS score (local run)

Compute MAPS score on local gnomAD hail table annotated in scripts/prepare_gnomad.py
    Input:
        UTR interval bed file
        prepared gnomAD hail table directory
    Output:
        MAPS hail table directory
"""
import hail as hl
from utr3variants.maps import maps


if __name__ == '__main__':
    ref = snakemake.config['genome_assembly']
    interval_file = snakemake.input['intervals'].__str__()
    gnomad_path = snakemake.input['gnomAD'].__str__()
    mutation_path = snakemake.config['gnomAD']['mutation_rate_ht'].__str__()
    maps_out = snakemake.output['maps']

    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log['hail'],
        default_reference=ref,
    )

    # intervals = hl.import_bed(bed_file, reference_genome=reference_genome')
    with open(interval_file, 'r') as f:
        intervals = f.readlines()
    intervals = [hl.parse_locus_interval(x, reference_genome=ref) for x in intervals]

    mutation_ht = hl.read_table(mutation_path)
    snp_ht = hl.read_table(gnomad_path)

    print('filter intervals')
    snp_ht = hl.filter_intervals(snp_ht, intervals)
    # snp_ht = snp_ht.filter(hl.is_defined(intervals[snp_ht.locus]))

    print('MAPS score')
    maps_ht = maps(snp_ht, mutation_ht, additional_grouping=['protein_coding'])
    maps_ht.show()

    print('save...')
    maps_ht.export(maps_out)
