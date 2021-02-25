"""
Annotation functions for gnomAD hail table
"""
from typing import Union, Any
import gnomad.utils.vep
import hail as hl


def annotate_by_intervals(
    ht, intervals_ht, annotation_column='target', new_column=None
):
    """
    Annotate a locus-keyed hail table by interval-level annotation

    :param  ht: hail table to be annotated
    :param intervals_ht: hail table with annotation on interval-basis
    :param annotation_column: name of annotation column from intervals_ht
        default: 'target' annotation column of imported UCSC BED hail table
    """
    if new_column is None:
        new_column = annotation_column

    ht = ht.filter(hl.is_defined(intervals_ht[ht.locus]))
    anno_values = intervals_ht.aggregate(
        hl.agg.collect_as_set(intervals_ht[annotation_column])
    )

    # build annotation expression
    expr = hl.case()
    for value in anno_values:
        print(value)
        interval_sub = intervals_ht.filter(intervals_ht[annotation_column] == value)
        expr = expr.when(hl.is_defined(interval_sub[ht.locus]), value)
    expr = expr.default('na')

    return ht.annotate(**{new_column: expr})


def filter_gnomad(ht: hl.Table, intervals: Any = None, verbose: bool = True):
    """
    Annotate gnomAD hail table with necessary

    :params ht: gnomAD hail table to be annotated
    :params intervals: list of hl.IntervalExpression objects (optional)
    :params verbose: whether to print a message after filtering
    """
    if intervals is not None:
        ht = hl.filter_intervals(ht, intervals)
    ht = ht.filter(hl.len(ht.filters) == 0)
    ht = gnomad.utils.vep.filter_vep_to_canonical_transcripts(ht)
    if verbose:
        print(f'entries in gnomAD after filtering: {ht.count()}')
    return ht


def annotate_for_maps(ht, context_ht):
    """
    Include annotations for MAPS analogue to
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint/summary_statistics.py#L96  # noqa: E501

    :param ht: gnomAD hail table to be annotated
    :param context_ht: context data to include as annotation
    """
    ht = get_worst_consequence_with_non_coding(ht)
    context_ht = context_ht[ht.key]
    return prepare_ht(
        ht=ht.annotate(context=context_ht.context, methylation=context_ht.methylation),
        trimer=True,
    )


########################################################################################
# Functions from https://github.com/macarthur-lab/gnomad_lof
########################################################################################


def get_worst_consequence_with_non_coding(ht):
    """
    copied and adapted from gnomad_lof
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint/summary_statistics.py#L18 # noqa: E501
    removed commented statements
    """

    def get_worst_csq(
        csq_list: hl.expr.ArrayExpression, protein_coding: bool
    ) -> hl.struct:
        lof = hl.null(hl.tstr)
        no_lof_flags = hl.null(hl.tbool)
        # lof_filters = hl.null(hl.tstr)
        # lof_flags = hl.null(hl.tstr)
        if protein_coding:
            all_lofs = csq_list.map(lambda x: x.lof)
            lof = hl.literal(['HC', 'OS', 'LC']).find(all_lofs.contains)
            csq_list = hl.if_else(  # changed from cond to ifelse
                hl.is_defined(lof), csq_list.filter(lambda x: x.lof == lof), csq_list
            )
            no_lof_flags = hl.or_missing(
                hl.is_defined(lof),
                csq_list.any(lambda x: (x.lof == lof) & hl.is_missing(x.lof_flags)),
            )
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        worst_csq = hl.literal(gnomad.utils.vep.CSQ_ORDER).find(all_csq_terms.contains)
        return hl.struct(
            worst_csq=worst_csq,
            protein_coding=protein_coding,
            lof=lof,
            no_lof_flags=no_lof_flags,
        )

    protein_coding = ht.vep.transcript_consequences.filter(
        lambda x: x.biotype == 'protein_coding'
    )
    return ht.annotate(
        **hl.case(missing_false=True)
        .when(hl.len(protein_coding) > 0, get_worst_csq(protein_coding, True))
        .when(
            hl.len(ht.vep.transcript_consequences) > 0,
            get_worst_csq(ht.vep.transcript_consequences, False),
        )
        .when(
            hl.len(ht.vep.regulatory_feature_consequences) > 0,
            get_worst_csq(ht.vep.regulatory_feature_consequences, False),
        )
        .when(
            hl.len(ht.vep.motif_feature_consequences) > 0,
            get_worst_csq(ht.vep.motif_feature_consequences, False),
        )
        .default(get_worst_csq(ht.vep.intergenic_consequences, False))
    )


def prepare_ht(ht, trimer: bool = False):
    """
    adapted from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/constraint_basics.py#L199 # noqa: E501
    removed annotate_coverage parameter
    include as inner function: trimer_from_heptamer
        https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L162 # noqa: E501
    """
    if trimer:

        def trimer_from_heptamer(
            t: Union[hl.MatrixTable, hl.Table]
        ) -> Union[hl.MatrixTable, hl.Table]:
            trimer_expr = hl.if_else(hl.len(t.context) == 7, t.context[2:5], t.context)
            return (
                t.annotate_rows(context=trimer_expr)
                if isinstance(t, hl.MatrixTable)
                else t.annotate(context=trimer_expr)
            )

        ht = trimer_from_heptamer(ht)

    str_len = 3 if trimer else 7
    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f'[ATCG]{{{str_len}}}')
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows(
            (hl.len(ht.ref) == 1)
            & (hl.len(ht.alt) == 1)
            & ht.context.matches(f'[ATCG]{{{str_len}}}')
        )
        ht = annotate_variant_types(collapse_strand(ht), not trimer)

    annotation = {
        'methylation_level': hl.case()
        .when(ht.cpg & (ht.methylation.MEAN > 0.6), 2)
        .when(ht.cpg & (ht.methylation.MEAN > 0.2), 1)
        .default(0)
    }

    if isinstance(ht, hl.Table):
        ht = ht.annotate(**annotation)
    else:
        ht = ht.annotate_rows(**annotation)

    return ht


def annotate_variant_types(
    t: Union[hl.MatrixTable, hl.Table], heptamers: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L110 # noqa: E501
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (
        ((t.ref == 'A') & (t.alt == 'G'))
        | ((t.ref == 'G') & (t.alt == 'A'))
        | ((t.ref == 'T') & (t.alt == 'C'))
        | ((t.ref == 'C') & (t.alt == 'T'))
    )
    cpg_expr = (
        (t.ref == 'G') & (t.alt == 'A') & (t.context[mid_index - 1 : mid_index] == 'C')
    ) | (
        (t.ref == 'C')
        & (t.alt == 'T')
        & (t.context[mid_index + 1 : mid_index + 2] == 'G')
    )
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (
        hl.case()
        .when(t.cpg, 'CpG')
        .when(t.transition, 'non-CpG transition')
        .default('transversion')
    )
    variant_type_model_expr = hl.if_else(t.cpg, t.context, 'non-CpG')
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(
            variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
        )
    return t.annotate(
        variant_type=variant_type_expr, variant_type_model=variant_type_model_expr
    )


def collapse_strand(
    ht: Union[hl.Table, hl.MatrixTable]
) -> Union[hl.Table, hl.MatrixTable]:
    """
    from:
    https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L36 # noqa E501
    include as inner function reverse_complement_bases
        https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L22 # noqa E501
    include as inner function flip_base
        https://github.com/macarthur-lab/gnomad_lof/blob/master/constraint_utils/generic.py#L27
    """

    def reverse_complement_bases(
        bases: hl.expr.StringExpression,
    ) -> hl.expr.StringExpression:
        def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
            return (
                hl.switch(base)
                .when('A', 'T')
                .when('T', 'A')
                .when('G', 'C')
                .when('C', 'G')
                .default(base)
            )

        return hl.delimit(
            hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), ''
        )

    collapse_expr = {
        'ref': hl.if_else(
            ((ht.ref == 'G') | (ht.ref == 'T')),
            reverse_complement_bases(ht.ref),
            ht.ref,
        ),
        'alt': hl.if_else(
            ((ht.ref == 'G') | (ht.ref == 'T')),
            reverse_complement_bases(ht.alt),
            ht.alt,
        ),
        'context': hl.if_else(
            ((ht.ref == 'G') | (ht.ref == 'T')),
            reverse_complement_bases(ht.context),
            ht.context,
        ),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T'),
    }
    return (
        ht.annotate(**collapse_expr)
        if isinstance(ht, hl.Table)
        else ht.annotate_rows(**collapse_expr)
    )
