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
from utr3variants.annotate_gnomad import annotate_by_intervals
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

    intervals = hl.import_bed(interval_file)
    mutation_ht = hl.read_table(mutation_path)
    ht = hl.read_table(gnomad_path)

    print('Annotate by intervals')
    ht = annotate_by_intervals(ht, intervals, new_column='UTR_group', repartition=False)

    print('MAPS score')
    maps_ht = maps(ht, mutation_ht, additional_grouping=['UTR_group'])
    maps_ht.show()

    print('save...')
    maps_ht.export(maps_out)
