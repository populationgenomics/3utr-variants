"""
Subset and annotate gnomAD hail table and save locally
    Input:
        gnomAD hail table URL (specified in config)
        context hail table URL (specified in config)
    Output:
        directory for subset and annotated gnomAD hail table
"""
import hail as hl
from utr3variants.annotate_gnomad import filter_gnomad, annotate_for_maps


if __name__ == '__main__':
    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log['hail'],
        default_reference=snakemake.config['genome_assembly'],
    )

    context_ht = hl.read_table(snakemake.config['gnomAD']['context_ht'])
    ht = hl.read_table(snakemake.config['gnomAD']['gnomAD_ht'])

    subset_interval = hl.parse_locus_interval(snakemake.config['gnomAD']['subset'])
    ht = filter_gnomad(ht, [subset_interval])
    ht = annotate_for_maps(ht, context_ht)

    print('save...')
    ht.write(snakemake.output['gnomAD_ht'])
