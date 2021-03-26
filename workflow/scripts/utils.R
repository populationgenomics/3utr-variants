# Title     : MAPS for consequence
# Objective : Compute reference values for worst_csq
# Created by: mumichae
# Created on: 22/3/21

library(tools)
library(magrittr)

## CONSTANTS
# from: https://github.com/macarthur-lab/gnomad_lof/blob/master/R/constants.R#L103
lof_like <- c(
  'frameshift_variant', 'essential_splice', 'stop_gained', 'splice_donor_variant',
  'splice_acceptor_variant'
)
mis_like <- c(
  'missense_variant', 'inframe_indel', 'stop_lost', 'mature_miRNA_variant', 'start_lost'
)
syn_like <- c(
  'synonymous_variant', '3_prime_UTR_variant', '5_prime_UTR_variant',
  'splice_region_variant', 'extended_splice', 'stop_retained_variant',
  'non_coding_transcript_exon_variant', 'upstream_gene_variant',
  'downstream_gene_variant', 'intron_variant', 'intergenic_variant',
  'regulatory_region_variant'
)


format_vep_category <- function(category_list) {
  # from: https://github.com/macarthur-lab/gnomad_lof/blob/master/R/constants.R#L111
  return(category_list %>%
           gsub("_", " ", .) %>%
           gsub('stop gained', 'nonsense', .) %>%
           gsub("inframe indel", "in-frame indel", .) %>%
           gsub("initiator codon", "start lost", .) %>%
           gsub(" variant", "", .) %>%
           gsub("transcript exon", "transcript", .) %>%
           gsub(" prime ", "'", .) %>%
           gsub("probably damaging", "prob damaging", .) %>%
           gsub("possibly damaging", "poss damaging", .))
}


maps <- function(count_dt, grouping) {
  grouping <- grouping[!is.na(grouping)]
  maps_dt <- count_dt[, .(
    variant_count = sum(variant_count),
    singleton_count = sum(singleton_count),
    expected_singletons = sum(expected_singletons)
  ), by = grouping]

  maps_dt[, ps := singleton_count / variant_count]
  maps_dt[, maps := (singleton_count - expected_singletons) / variant_count]
  maps_dt[, maps_sem := sqrt(ps * (1 - ps) / variant_count)]
  maps_dt
}


maps_reference <- function(count_dt) {
  # get worst_csq MAPS for reference lines
  dt_csq <- maps(count_dt, grouping = 'worst_csq')
  dt_csq[worst_csq %in% lof_like, consequence := 'LoF']
  dt_csq[worst_csq %in% mis_like, consequence := 'missense']
  #dt_csq[worst_csq %in% syn_like, consequence := 'other synonymous-like']
  dt_csq[worst_csq == 'synonymous_variant', consequence := 'synonymous']
  dt_csq <- maps(dt_csq, grouping = 'consequence')
  dt_csq <- na.omit(dt_csq)
  dt_csq
}
