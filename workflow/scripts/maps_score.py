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
import pandas as pd
from utr3variants.annotate_gnomad import annotate_for_maps, annotate_by_intervals
from utr3variants.maps import maps_precomputed


if __name__ == '__main__':
    ref = snakemake.config['genome_assembly']
    interval_file = snakemake.input['intervals'].__str__()
    gnomad_path = snakemake.input['gnomAD'].__str__()
    mutation_path = snakemake.config['gnomAD']['mutation_rate_ht']
    context_path = snakemake.config['gnomAD']['context_ht']
    mut_fit_path = snakemake.input['mut_fit'].__str__()

    maps_out = snakemake.output['maps']

    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log['hail'],
        default_reference=ref,
    )

    with open(mut_fit_path, 'r') as f:
        intercept, slope = [float(line.rstrip('\n')) for line in f]

    intervals = hl.import_bed(interval_file)
    mutation_ht = hl.read_table(mutation_path)
    ht = hl.read_table(gnomad_path)
    context_ht = hl.read_table(context_path)
    print('Done reading objects.')

    ht = annotate_for_maps(ht, context_ht)
    ht = annotate_by_intervals(ht, intervals, new_column='UTR_group')

    maps_ht = maps_precomputed(
        ht, mutation_ht, mutability_fit=(intercept, slope), grouping=['UTR_group']
    )
    # maps_ht = maps(ht, mutation_ht, additional_grouping=['UTR_group'])

    print('save...')
    maps_ht.export(maps_out)
    maps_df = pd.read_table(maps_out)
    print(maps_df.head())
