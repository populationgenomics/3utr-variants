"""Compute MAPS score"""

from typing import List, Union
import pybedtools
import hail as hl
import sys


def bed_to_hail_interval(file: str, ref: str) -> List[hl.IntervalExpression]:
    """
    Merge overlapping intervals to reduce number of hail interval queries

    :param file: this is a first param
    :param ref: name of reference genome 'GRCh37' or 'GRChr38'
    :returns: list of converted hail.IntervalExpression objects
    """
    bed_intervals = pybedtools.BedTool(file)
    print(f'number of bed intervals: {len(bed_intervals)}')

    # shape bed intervals for parsing hail locus interval
    chr_style = '' if ref == 'GRCh37' else 'chr'

    def format_interval(x: pybedtools.Interval):
        chromosome = str(x.chrom)
        if chr_style == '' and 'chr' in chromosome:
            chromosome = chromosome.replace('chr', '')
        elif chr_style == 'chr' and 'chr' not in chromosome:
            chromosome = f'chr{chromosome}'
        return f'{chromosome}:{x.start + 1}-{x.stop + 1}'

    locus_interval = map(format_interval, bed_intervals)
    return [hl.parse_locus_interval(x, reference_genome=ref) for x in locus_interval]


def downsampling_counts_expr(
    ht: Union[hl.Table, hl.MatrixTable],
    pop: str = 'global',
    variant_quality: str = 'adj',
    singleton: bool = False,
    impose_high_af_cutoff: bool = False,
) -> hl.expr.ArrayExpression:
    indices = hl.zip_with_index(ht.freq_meta).filter(
        lambda f: (f[1].size() == 3)
        & (f[1].get('group') == variant_quality)
        & (f[1].get('pop') == pop)
        & f[1].contains('downsampling')
    )
    sorted_indices = hl.sorted(indices, key=lambda f: hl.int(f[1]['downsampling'])).map(
        lambda x: x[0]
    )

    # TODO: this likely needs to be fixed for aggregations that return missing

    def get_criteria(i):
        if singleton:
            return hl.int(ht.freq[i].AC == 1)
        elif impose_high_af_cutoff:
            return hl.int((ht.freq[i].AC > 0) & (ht.freq[i].AF <= 0.001))
        else:
            return hl.int(ht.freq[i].AC > 0)

    return hl.agg.array_sum(hl.map(get_criteria, sorted_indices))


def count_variants(
    ht: hl.Table,
    count_singletons: bool = False,
    count_downsamplings=(),
    additional_grouping=(),
    partition_hint: int = 100,
    omit_methylation: bool = False,
    return_type_only: bool = False,
    force_grouping: bool = False,
    singleton_expression: hl.expr.BooleanExpression = None,
    impose_high_af_cutoff_here: bool = False,
):
    """
    Count variants by context, ref, alt, methylation_level
    """

    grouping = hl.struct(context=ht.context, ref=ht.ref, alt=ht.alt)
    if not omit_methylation:
        grouping = grouping.annotate(methylation_level=ht.methylation_level)
    for group in additional_grouping:
        grouping = grouping.annotate(**{group: ht[group]})

    if count_singletons:
        # singleton = hl.any(
        #    lambda f: (f.meta.size() == 1)
        #    & (f.meta.get('group') == 'adj')
        #    & (f.AC[1] == 1),
        #    ht.freq,
        # )
        if singleton_expression is None:
            singleton_expression = ht.freq[0].AC == 1

    if count_downsamplings or force_grouping:
        # Slower, but more flexible (allows for downsampling agg's)
        output = {
            'variant_count': hl.agg.count_where(ht.freq[0].AF <= 0.001)
            if impose_high_af_cutoff_here
            else hl.agg.count()
        }
        for pop in count_downsamplings:
            output[f'downsampling_counts_{pop}'] = downsampling_counts_expr(
                ht, pop, impose_high_af_cutoff=impose_high_af_cutoff_here
            )
        if count_singletons:
            output['singleton_count'] = hl.agg.count_where(singleton_expression)
            for pop in count_downsamplings:
                output[
                    f'singleton_downsampling_counts_{pop}'
                ] = downsampling_counts_expr(ht, pop, singleton=True)
        return (
            ht.group_by(**grouping).partition_hint(partition_hint).aggregate(**output)
        )
    else:
        agg = {'variant_count': hl.agg.counter(grouping)}
        if count_singletons:
            agg['singleton_count'] = hl.agg.counter(
                hl.agg.filter(singleton_expression, grouping)
            )

        if return_type_only:
            return agg['variant_count'].dtype
        else:
            return ht.aggregate(hl.struct(**agg))


def maps(
    ht: hl.Table,
    mutation_ht: hl.Table,
    additional_grouping: List[str] = [],
    singleton_expression: hl.expr.BooleanExpression = None,
    skip_worst_csq: bool = False,
) -> hl.Table:
    if not skip_worst_csq:
        additional_grouping.insert(0, 'worst_csq')
    ht = count_variants(
        ht,
        count_singletons=True,
        additional_grouping=additional_grouping,
        force_grouping=True,
        singleton_expression=singleton_expression,
    )
    ht = ht.annotate(
        mu=mutation_ht[
            hl.struct(
                context=ht.context,
                ref=ht.ref,
                alt=ht.alt,
                methylation_level=ht.methylation_level,
            )
        ].mu_snp,
        ps=ht.singleton_count / ht.variant_count,
    )
    if not ht.all(hl.is_defined(ht.mu)):
        print('Some mu were not found...')
        print(
            ht.aggregate(hl.agg.filter(hl.is_missing(ht.mu), hl.agg.take(ht.row, 1)[0]))
        )
        sys.exit(1)
    syn_ps_ht = ht.filter(ht.worst_csq == 'synonymous_variant')
    syn_ps_ht = syn_ps_ht.group_by(syn_ps_ht.mu).aggregate(
        singleton_count=hl.agg.sum(syn_ps_ht.singleton_count),
        variant_count=hl.agg.sum(syn_ps_ht.variant_count),
    )
    syn_ps_ht = syn_ps_ht.annotate(
        ps=syn_ps_ht.singleton_count / syn_ps_ht.variant_count
    )

    lm = syn_ps_ht.aggregate(
        hl.agg.linreg(
            syn_ps_ht.ps, [1, syn_ps_ht.mu], weight=syn_ps_ht.variant_count
        ).beta
    )
    print(f'Got MAPS calibration model of: slope: {lm[1]}, intercept: {lm[0]}')
    ht = ht.annotate(expected_singletons=(ht.mu * lm[1] + lm[0]) * ht.variant_count)

    agg_ht = ht.group_by(*additional_grouping).aggregate(
        singleton_count=hl.agg.sum(ht.singleton_count),
        expected_singletons=hl.agg.sum(ht.expected_singletons),
        variant_count=hl.agg.sum(ht.variant_count),
    )
    agg_ht = agg_ht.annotate(
        ps=agg_ht.singleton_count / agg_ht.variant_count,
        maps=(agg_ht.singleton_count - agg_ht.expected_singletons)
        / agg_ht.variant_count,
    )
    agg_ht = agg_ht.annotate(
        maps_sem=(agg_ht.ps * (1 - agg_ht.ps) / agg_ht.variant_count) ** 0.5
    )
    return agg_ht


if __name__ == '__main__':
    reference_genome = snakemake.config['genome_assembly']
    bed_file = snakemake.input['bed'].__str__()
    gnomad_path = snakemake.input['gnomAD'].__str__()
    mutation_path = snakemake.config['gnomAD']['mutation_rate_ht'].__str__()
    out_dir = snakemake.output['maps']

    hl.init(
        local=f'local[{snakemake.threads}]',
        log=snakemake.log[0],
        default_reference=reference_genome,
    )

    intervals = bed_to_hail_interval(bed_file, ref=reference_genome)
    # intervals = hl.import_bed(bed_file, reference_genome=reference_genome')
    mutation_ht = hl.read_table(mutation_path)
    snp_ht = hl.read_table(gnomad_path)

    print('filter intervals')
    snp_ht = hl.filter_intervals(snp_ht, intervals)
    # snp_ht = snp_ht.filter(hl.is_defined(intervals[snp_ht.locus]))

    print('MAPS score')
    maps_ht = maps(snp_ht, mutation_ht, additional_grouping=['protein_coding'])

    print('save...')
    maps_ht.write(out_dir)
