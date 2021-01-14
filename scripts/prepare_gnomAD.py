from gnomad.utils.vep import *


def get_worst_consequence_with_non_coding(ht):
    """
    copied and adapted from gnomad_lof
    """

    def get_worst_csq(csq_list: hl.expr.ArrayExpression, protein_coding: bool) -> hl.struct:
        lof = hl.null(hl.tstr)
        no_lof_flags = hl.null(hl.tbool)
        # lof_filters = hl.null(hl.tstr)
        # lof_flags = hl.null(hl.tstr)
        if protein_coding:
            all_lofs = csq_list.map(lambda x: x.lof)
            lof = hl.literal(['HC', 'OS', 'LC']).find(lambda x: all_lofs.contains(x))
            csq_list = hl.cond(hl.is_defined(lof), csq_list.filter(lambda x: x.lof == lof), csq_list)
            no_lof_flags = hl.or_missing(hl.is_defined(lof),
                                         csq_list.any(lambda x: (x.lof == lof) & hl.is_missing(x.lof_flags)))
            # lof_filters = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_filter).filter(lambda x: hl.is_defined(x))), '|')
            # lof_flags = hl.delimit(hl.set(csq_list.map(lambda x: x.lof_flags).filter(lambda x: hl.is_defined(x))), '|')
        all_csq_terms = csq_list.flatmap(lambda x: x.consequence_terms)
        worst_csq = hl.literal(CSQ_ORDER).find(lambda x: all_csq_terms.contains(x))
        return hl.struct(
            worst_csq=worst_csq,
            protein_coding=protein_coding,
            lof=lof,
            no_lof_flags=no_lof_flags,
            # lof_filters=lof_filters,
            # lof_flags=lof_flags
        )

    protein_coding = ht.vep.transcript_consequences.filter(lambda x: x.biotype == 'protein_coding')
    return ht.annotate(
        **hl.case(missing_false=True)
            .when(hl.len(protein_coding) > 0, get_worst_csq(protein_coding, True))
            .when(hl.len(ht.vep.transcript_consequences) > 0,
                  get_worst_csq(ht.vep.transcript_consequences, False))
            .when(hl.len(ht.vep.regulatory_feature_consequences) > 0,
                  get_worst_csq(ht.vep.regulatory_feature_consequences, False))
            .when(hl.len(ht.vep.motif_feature_consequences) > 0,
                  get_worst_csq(ht.vep.motif_feature_consequences, False))
            .default(get_worst_csq(ht.vep.intergenic_consequences, False)))


def prepare_ht(ht, trimer: bool = False, annotate_coverage: bool = True):
    if trimer:
        def trimer_from_heptamer(t: Union[hl.MatrixTable, hl.Table]) -> Union[hl.MatrixTable, hl.Table]:
            trimer_expr = hl.cond(hl.len(t.context) == 7, t.context[2:5], t.context)
            return t.annotate_rows(context=trimer_expr) if isinstance(t, hl.MatrixTable) else t.annotate(
                context=trimer_expr)

        ht = trimer_from_heptamer(ht)

    str_len = 3 if trimer else 7
    if isinstance(ht, hl.Table):
        ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)
    else:
        ht = ht.annotate_rows(ref=ht.alleles[0], alt=ht.alleles[1])
        ht = ht.filter_rows((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))
        ht = annotate_variant_types(collapse_strand(ht), not trimer)

    annotation = {
        'methylation_level': hl.case().when(
            ht.cpg & (ht.methylation.MEAN > 0.6), 2
        ).when(
            ht.cpg & (ht.methylation.MEAN > 0.2), 1
        ).default(0)
    }

    if annotate_coverage:
        annotation['exome_coverage'] = ht.coverage.exomes.median

    return ht.annotate(**annotation) if isinstance(ht, hl.Table) \
        else ht.annotate_rows(**annotation)


def annotate_variant_types(t: Union[hl.MatrixTable, hl.Table], heptamers: bool = False) \
        -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (((t.ref == "A") & (t.alt == "G")) | ((t.ref == "G") & (t.alt == "A")) |
                       ((t.ref == "T") & (t.alt == "C")) | ((t.ref == "C") & (t.alt == "T")))
    cpg_expr = (((t.ref == "G") & (t.alt == "A") & (t.context[mid_index - 1:mid_index] == 'C')) |
                ((t.ref == "C") & (t.alt == "T") & (t.context[mid_index + 1:mid_index + 2] == 'G')))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = (hl.case()
                         .when(t.cpg, 'CpG')
                         .when(t.transition, 'non-CpG transition')
                         .default('transversion'))
    variant_type_model_expr = hl.cond(t.cpg, t.context, "non-CpG")
    if isinstance(t, hl.MatrixTable):
        return t.annotate_rows(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)
    else:
        return t.annotate(variant_type=variant_type_expr, variant_type_model=variant_type_model_expr)


def collapse_strand(ht: Union[hl.Table, hl.MatrixTable]) -> Union[hl.Table, hl.MatrixTable]:
    def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
        def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
            return (hl.switch(base)
                    .when('A', 'T')
                    .when('T', 'A')
                    .when('G', 'C')
                    .when('C', 'G')
                    .default(base))

        return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')

    collapse_expr = {
        'ref': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.ref), ht.ref),
        'alt': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                       reverse_complement_bases(ht.alt), ht.alt),
        'context': hl.cond(((ht.ref == 'G') | (ht.ref == 'T')),
                           reverse_complement_bases(ht.context), ht.context),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T')
    }
    return ht.annotate(**collapse_expr) if isinstance(ht, hl.Table) else ht.annotate_rows(**collapse_expr)


if __name__ == '__main__':
    hl.init()

    gnomAD_out_ht = snakemake.output[0]

    context_ht = hl.read_table(snakemake.config["gnomAD"]["context_ht"])
    coverage_ht = hl.read_table(snakemake.config["gnomAD"]["coverage_ht"])
    ht = hl.read_table(snakemake.config["gnomAD"]["gnomAD_ht"])

    ht = ht.annotate(coverage=coverage_ht[ht.locus].median)
    ht = ht.filter((hl.len(ht.filters) == 0))  # & get_an_adj_criteria(ht, sample_sizes))
    ht = filter_vep_to_canonical_transcripts(ht)

    ht = get_worst_consequence_with_non_coding(ht)

    context = context_ht[ht.key]
    snp_ht = prepare_ht(
        ht=ht.annotate(context=context.context, methylation=context.methylation),
        trimer=True,
        annotate_coverage=False)
    snp_ht = snp_ht.persist()

    snp_ht.write(gnomAD_out_ht)
