"""
MAPS score functions from https://github.com/macarthur-lab/gnomad_lof
"""
import sys
from typing import Union, Iterable
import hail as hl


def downsampling_counts_expr(
    ht: Union[hl.Table, hl.MatrixTable],
    pop: str = 'global',
    variant_quality: str = 'adj',
    singleton: bool = False,
    impose_high_af_cutoff: bool = False,
) -> hl.expr.ArrayExpression:
    """
    unmodified from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L49
    """
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
        if singleton:  # pylint: disable=R1705
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
    unmodified from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L68
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

    if count_downsamplings or force_grouping:  # pylint: disable=R1705
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

        if return_type_only:  # pylint: disable=R1705
            return agg['variant_count'].dtype
        else:
            return ht.aggregate(hl.struct(**agg))


def add_mutation_rate(ht: hl.Table, mutation_ht: hl.Table, mut_check: bool = True):
    """
    Add mutation rates indexed by context, ref/alt alleles & methylation context
    :param ht: hail table to be annotated
    :param mutation_ht: mutation rates indexed by context, ref/alt & methylation_context
    :param mut_check: whether to check for missing mutation rates. This check can be
        skipped to save time e.g. if it has already been done before
    """
    ht = ht.annotate(
        mu=mutation_ht[
            hl.struct(
                context=ht.context,
                ref=ht.ref,
                alt=ht.alt,
                methylation_level=ht.methylation_level,
            )
        ].mu_snp
    )

    # edit: force skip check if certain mutation ht is fine
    if mut_check:
        print('Check mutation rates...')
        if not ht.all(hl.is_defined(ht.mu)):
            print('Some mu were not found...')
            print(
                ht.aggregate(
                    hl.agg.filter(hl.is_missing(ht.mu), hl.agg.take(ht.row, 1)[0])
                )
            )
            sys.exit(1)

    return ht


def fit_mutability(count_ht: hl.Table, worst_csq='synonymous_variant'):
    """
    Subset to worst consequence and fit model to mutation rates
    :param count_ht: must contain 'worst_csq', 'variant_counts' and 'singleton_counts'
    :param worst_csq: name of worst consequence to normalise by
    :return: (slope, intercept) of mutability fit
    """
    ht = count_ht.filter(count_ht.worst_csq == worst_csq)
    ht = ht.group_by(ht.mu).aggregate(
        singleton_count=hl.agg.sum(ht.singleton_count),
        variant_count=hl.agg.sum(ht.variant_count),
    )
    ht = ht.annotate(ps=ht.singleton_count / ht.variant_count)
    print('Fit linear model...')
    intercept, slope = ht.aggregate(
        hl.agg.linreg(ht.ps, [1, ht.mu], weight=ht.variant_count).beta
    )
    print(f'Got MAPS calibration model of: slope: {slope}, intercept: {intercept}')
    return slope, intercept


def count_for_maps(
    ht: hl.Table,
    mutation_ht: hl.Table,  # pylint: disable=W0621
    additional_grouping=None,  # change from mutable default to None
    singleton_expression: hl.expr.BooleanExpression = None,
    skip_mut_check: bool = False,  # added
) -> hl.Table:
    """
    Count singletons aggregated by provided grouping & provide expected singleton counts
    adapted from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L267
    :param ht: gnomAD variant hail table, must contain 'worst_csq'
    :param mutation_ht: mutation rates table
    :param additional_grouping: additional grouping to aggregate by
    :param singleton_expression: singleton definition
    :param skip_mut_check: skip checking for mutation rate ht
    """
    if additional_grouping is None:
        additional_grouping = []
    additional_grouping.insert(0, 'worst_csq')
    ht = count_variants(
        ht,
        count_singletons=True,
        additional_grouping=additional_grouping,
        force_grouping=True,
        singleton_expression=singleton_expression,
    )
    ht = add_mutation_rate(ht, mutation_ht, mut_check=not skip_mut_check)
    slope, intercept = fit_mutability(ht)
    return ht.annotate(
        expected_singletons=(ht.mu * slope + intercept) * ht.variant_count
    )


def aggregate_maps_hail(count_ht: hl.Table, grouping: Iterable[str]) -> hl.Table:
    """
    Aggregate MAPS scores of variant count hail table from count_for_maps()
    adapted from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L267

    :param count_ht: contains variant counts, singleton counts and expected singletons
    :param grouping: names of annotation columns in count_ht to aggregate by
    :return: hail table with aggregated MAPS statistics
    """
    agg_ht = count_ht.group_by(*grouping).aggregate(
        singleton_count=hl.agg.sum(count_ht.singleton_count),
        expected_singletons=hl.agg.sum(count_ht.expected_singletons),
        variant_count=hl.agg.sum(count_ht.variant_count),
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
