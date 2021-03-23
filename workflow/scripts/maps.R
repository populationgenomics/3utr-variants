# Title     : MAPS Score
# Objective : Compute MAPS score and save plot
# Created by: mumichae
# Created on: 18/3/21

library(data.table)
library(ggplot2)
source(snakemake@input$utils)

count_dt <- fread(snakemake@input$counts)

# get aggregation MAPS
aggregations <- strsplit(snakemake@wildcards$aggregation, '-')[[1]]
dt <- maps(count_dt, grouping = aggregations)
fwrite(dt, snakemake@output$tsv, sep = '\t')

# get reference MAPS
dt_csq <- maps_reference(count_dt)

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
