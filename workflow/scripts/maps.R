# Title     : MAPS Score
# Objective : Compute MAPS score and save plot
# Created by: mumichae
# Created on: 18/3/21

library(tools)
library(data.table)
library(ggplot2)
library(magrittr)

## CONSTANTS TODO: move to separate script for worst_csq
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


# TODO: move to separate file
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


count_dt <- fread(snakemake@input$counts)

# TODO: move worst_csq to separate script
# get worst_csq MAPS for reference lines
dt_csq <- maps(count_dt, grouping = 'worst_csq')
dt_csq[worst_csq %in% lof_like, consequence := 'LoF']
dt_csq[worst_csq %in% mis_like, consequence := 'missense']
#dt_csq[worst_csq %in% syn_like, consequence := 'other synonymous-like']
dt_csq[worst_csq == 'synonymous_variant', consequence := 'synonymous']
dt_csq <- maps(dt_csq, grouping = 'consequence')

# get aggregation MAPS
aggregations <- strsplit(snakemake@wildcards$aggregation, '-')[[1]]
dt <- maps(count_dt, grouping = aggregations)
fwrite(dt, snakemake@output$tsv, sep = '\t')


# Plot only top 100 variants
dt <- dt[variant_count > 100]
setorderv(dt, cols = aggregations)

new_columns <- c('anno_x', 'shape', 'facet')[seq_along(aggregations)]
setnames(dt, old = aggregations, new = new_columns)
# dt[, anno_x := paste0(anno_x, '\n(', variant_count, ')')]

# Flatten for certain x-values
if ('shape' %in% names(dt) & aggregations[1] %in% c('hexamer', 'conserved')) {
  dt_lines <- dt[shape %in% anno_x]
  dt <- dt[!shape %in% anno_x]
}

title <- paste0('MAPS on ', snakemake@wildcards$chr_subset,
                ' (', snakemake@params$chr_subset, ')')
dodge_width <- 0.8

if ('shape' %in% names(dt)) {
  p <- ggplot(dt, aes(reorder(anno_x, variant_count), maps, shape = shape, group = shape))
} else {
  p <- ggplot(dt, aes(reorder(anno_x, variant_count), maps))
}

if ('facet' %in% names(dt)) {
  p <- p +
    facet_grid('facet~.') +
    ggtitle(title, subtitle = paste('Facet:', toTitleCase(aggregations[3])))
}

if (exists('dt_lines')) {
  p <- p + geom_hline(
      aes(yintercept = maps, color = shape),
      data = dt_lines[, .(maps, shape)]
    )
}

p <- p +
  geom_hline(yintercept = 0, colour = 'grey50', linetype = 'dotted') +
  geom_hline(  # add lines for reference
    aes(yintercept = maps, color = consequence),
    data = dt_csq[, .(maps, consequence)],
    linetype = 'dashed'
  ) +
  geom_errorbar(
    aes(ymin = maps - maps_sem, ymax = maps + maps_sem),
    width = 0.15,
    position = position_dodge(width = dodge_width)
  ) +
  geom_point(position = position_dodge(width = dodge_width)) +
  labs(
    title = title,
    x = toTitleCase(aggregations[1]),
    y = 'MAPS',
    shape = toTitleCase(aggregations[2]),
    color = 'Consequence'
  ) +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks = element_blank(),
  )

ggsave(snakemake@output$png, width = 8, height = 6)
