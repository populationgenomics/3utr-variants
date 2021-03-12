"""
MAPS score functions from https://github.com/macarthur-lab/gnomad_lof
"""
import sys
from typing import Union, Tuple
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


def fit_mutability(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mut_check: bool = True,
    singleton_expression: hl.expr.BooleanExpression = None,
    worst_csq: str = 'synonymous_variant',
) -> Tuple[float, float]:
    """
    Fit linear model to mutation rates against worst consequence class
    :param ht: gnomAD hail table, needs to have annotated 'worst_csq'
    :param mutation_ht: mutation rates hail table
        indexed by context, ref, alt and methylation_level
    :param mut_check: whether to check for missing mutation rates
    :param singleton_expression: singleton definition
    :param worst_csq: name of worst consequence class to calibrate against,
        default: synonymous_variant
    :return: linear model fit slope & intercept
    """
    ht = count_variants(
        ht,
        count_singletons=True,
        additional_grouping=('worst_csq',),
        force_grouping=True,
        singleton_expression=singleton_expression,
    )
    ht = add_mutation_rate(ht, mutation_ht, mut_check=mut_check)

    # Subset to worst consequence class variants
    syn_ps_ht = ht.filter(ht.worst_csq == worst_csq)
    syn_ps_ht = syn_ps_ht.group_by(syn_ps_ht.mu).aggregate(
        singleton_count=hl.agg.sum(syn_ps_ht.singleton_count),
        variant_count=hl.agg.sum(syn_ps_ht.variant_count),
    )
    syn_ps_ht = syn_ps_ht.annotate(
        ps=syn_ps_ht.singleton_count / syn_ps_ht.variant_count
    )

    print('Linear regression...')
    intercept, slope = syn_ps_ht.aggregate(
        hl.agg.linreg(
            syn_ps_ht.ps, [1, syn_ps_ht.mu], weight=syn_ps_ht.variant_count
        ).beta
    )
    print(f'Got MAPS calibration model of: slope: {slope}, intercept: {intercept}')
    return intercept, slope


def maps_precomputed(
    ht: hl.Table,
    mutation_ht: hl.Table,
    mutability_fit: Tuple[float, float],
    grouping=None,
    singleton_expression: hl.expr.BooleanExpression = None,
    mut_check: bool = False,
) -> hl.Table:
    """
    Compute MAPS based on precomputed mutability fit
    :param ht: gnomAD hail table, annotated with mutation rates
    :param mutation_ht: mutation rates for mutability adjustment
    :param mutability_fit: model fit of mutation rates as tuple of (intersect, slope)
    :param grouping: tuple of group names to aggregate MAPS by
    :param singleton_expression: singleton definition
    :param mut_check: whether to check for missing mutation rates
    """
    ht = count_variants(
        ht,
        count_singletons=True,
        additional_grouping=grouping,
        force_grouping=True,
        singleton_expression=singleton_expression,
    )
    ht = add_mutation_rate(ht, mutation_ht, mut_check=mut_check)

    intercept, slope = mutability_fit
    ht = ht.annotate(expected_singletons=(ht.mu * slope + intercept) * ht.variant_count)

    agg_ht = ht.group_by(*grouping).aggregate(
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


def maps(
    ht: hl.Table,
    mutation_ht: hl.Table,  # pylint: disable=W0621
    additional_grouping=None,  # change from mutable default to None
    singleton_expression: hl.expr.BooleanExpression = None,
    skip_worst_csq: bool = False,
    skip_mut_check: bool = False,  # added
) -> hl.Table:
    """
    adapted from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L267

    :param skip_mut_check: skip checking for mutation rate ht
    """
    if additional_grouping is None:
        additional_grouping = []  # use mutable default
    if not skip_worst_csq:
        additional_grouping.insert(0, 'worst_csq')

    print('Count variants')
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

    # edit: force skip check if certain mutation ht is fine
    if skip_mut_check:
        print('Check mutation rates...')
        if not ht.all(hl.is_defined(ht.mu)):
            print('Some mu were not found...')
            print(
                ht.aggregate(
                    hl.agg.filter(hl.is_missing(ht.mu), hl.agg.take(ht.row, 1)[0])
                )
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

    print('Fit linear model...')
    lm = syn_ps_ht.aggregate(
        hl.agg.linreg(
            syn_ps_ht.ps, [1, syn_ps_ht.mu], weight=syn_ps_ht.variant_count
        ).beta
    )
    print(f'Got MAPS calibration model of: slope: {lm[1]}, intercept: {lm[0]}')
    ht = ht.annotate(expected_singletons=(ht.mu * lm[1] + lm[0]) * ht.variant_count)

    print('Annotate MAPS statistics...')
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
