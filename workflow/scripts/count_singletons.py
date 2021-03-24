"""
Count variants, singletons and determine expected singletons for MAPS score (local run)

    Input:
        UTR interval bed file
        gnomAD hail table directory subset and annotated by scripts/prepare_gnomad.py
    Output:
        Exported table in TSV
"""
import hail as hl
import pandas as pd
from utr3variants.annotate_gnomad import annotate_by_intervals, import_interval_table
from utr3variants.maps import count_for_maps


if __name__ == '__main__':
    ref = snakemake.config['genome_assembly']
    interval_file = snakemake.input['intervals'].__str__()
    gnomad_path = snakemake.input['gnomAD'].__str__()
    mutation_path = snakemake.config['gnomAD']['mutation_rate_ht']

    counts_out = snakemake.output.counts

    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log['hail'],
        default_reference=ref,
    )
    mutation_ht = hl.read_table(mutation_path)
    ht = hl.read_table(gnomad_path)
    intervals = import_interval_table(interval_file, 'locus_interval')

    annotations = snakemake.params.annotations
    for annotation in annotations:
        ht = annotate_by_intervals(ht, intervals, annotation_column=annotation)

    count_ht = count_for_maps(
        ht,
        mutation_ht,
        additional_grouping=annotations,
        skip_mut_check=snakemake.config['gnomAD']['skip_checks'],
    )

    print('save...')
    count_ht.export(counts_out)

    # Preview of counts table
    count_df = pd.read_table(counts_out)
    print(count_df.head())
