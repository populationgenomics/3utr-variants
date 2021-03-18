# Title     : MAPS Score
# Objective : Compute MAPS score and save plot
# Created by: mumichae
# Created on: 18/3/21

library(ggplot2)
library(data.table)

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


dt <- fread(snakemake@input$counts)

chr_subset <- snakemake@params$chr_subset
aggregations <- strsplit(snakemake@wildcards$aggregation, '-')[[1]]

dt <- maps(dt, grouping = aggregations)
fwrite(dt, snakemake@output$tsv, sep = '\t')


# Plot only top 100 variants
dt <- dt[variant_count > 100]
setorderv(dt, cols = aggregations)

new_columns <- c('anno_x', 'color', 'facet')[seq_along(aggregations)]
setnames(dt, old = aggregations, new = new_columns)
# dt[, anno_x := paste0(anno_x, '\n(', variant_count, ')')]

title <- paste0('MAPS on: ', snakemake@wildcards$chr_subset, ' (', chr_subset, ')')
dodge_width <- 0.8

if('color' %in% names(dt)) {
  p <- ggplot(dt, aes(reorder(anno_x, variant_count), maps, color=color, group=color))
} else {
  p <- ggplot(dt, aes(reorder(anno_x, variant_count), maps))
}

if ('facet' %in% names(dt)) {
  p <- p + facet_grid('facet~.')
}

p <- p +
  geom_hline(yintercept = 0, color = 'grey') +
  geom_errorbar(
      aes(ymin = maps - maps_sem, ymax = maps + maps_sem),
      width = 0.2,
      position = position_dodge(width = dodge_width)
  ) +
  geom_point(position = position_dodge(width = dodge_width)) +
  labs(title = title, x = aggregations[1], y = 'MAPS') +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks = element_blank(),
  )

ggsave(snakemake@output$png, width = 8, height = 6)
