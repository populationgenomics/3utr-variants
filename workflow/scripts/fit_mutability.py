"""
Fit mutability (local run)

Compute linear model on mutation rates
    Input:
        gnomAD hail table annotated in scripts/prepare_gnomad.py
    Output:
        linear model fit paramters
"""
import hail as hl
from utr3variants.annotate_gnomad import annotate_for_maps
from utr3variants.maps import fit_mutability


if __name__ == '__main__':
    ref = snakemake.config['genome_assembly']
    gnomad_path = snakemake.input['gnomAD'].__str__()
    mutation_path = snakemake.config['gnomAD']['mutation_rate_ht']
    context_path = snakemake.config['gnomAD']['context_ht']

    mut_out = snakemake.output['fit']

    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log['hail'],
        default_reference=ref,
    )

    ht = hl.read_table(gnomad_path)
    mutation_ht = hl.read_table(mutation_path)
    context_ht = hl.read_table(context_path)

    ht = annotate_for_maps(ht, context_ht)
    intercept, slope = fit_mutability(
        ht, mutation_ht, mut_check=False
    )  # TODO: turn on check

    print('save...')
    with open(mut_out, 'w') as f:
        f.write(f'{intercept}\n{slope}\n')
