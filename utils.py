"""
Basic functions
"""
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
