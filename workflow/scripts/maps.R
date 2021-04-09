# Title     : MAPS Score
# Objective : Compute MAPS score and save plot
# Created by: mumichae
# Created on: 18/3/21

library(data.table)
library(ggplot2)
source(snakemake@input$utils)

count_dt <- fread(snakemake@input$counts)

# get reference MAPS
dt_csq <- maps_reference(count_dt)

# get aggregation MAPS
aggregations <- strsplit(snakemake@wildcards$aggregation, '-')[[1]]
dt <- maps(count_dt, grouping = aggregations)
if (aggregations[1] %in% c('expression', 'percent_expressed')) {
  setorderv(dt, aggregations[1], 1)
} else {
  setorder(dt, -maps)
}
fwrite(dt, snakemake@output$tsv, sep = '\t')

# Plot only top n variants
dt <- dt[variant_count > snakemake@params$variant_count_min]

new_columns <- c('anno_x', 'shape', 'facet')[seq_along(aggregations)]

# Flatten for 3'UTR and other variants
if (length(aggregations) == 3) {
  if ('database' %in% aggregations) {
    flatten_values <- c('GENCODE', 'gnomAD')
    #dt_lines <- dt[database %in% flatten_values]
    dt <- dt[!database %in% flatten_values]
  } else if ('feature' %in% aggregations) {
    flatten_values <- c('3UTR', 'other variant')
    #dt_lines <- dt[feature %in% flatten_values]
    dt <- dt[!feature %in% flatten_values]
  }
  #setnames(dt_lines, old = aggregations, new = new_columns)
}

setnames(dt, old = aggregations, new = new_columns)
dt[, anno_x := factor(anno_x, levels = unique(anno_x), ordered = TRUE)]

title <- paste0('MAPS on ', snakemake@wildcards$chr_subset,
                ' (', snakemake@params$chr_subset, ')')
dodge_width <- 0.5

if ('shape' %in% names(dt)) {
  p <- ggplot(dt, aes(anno_x, maps, color = shape, group = shape))
} else {
  p <- ggplot(dt, aes(anno_x, maps))
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
  geom_point(position = position_dodge(width = dodge_width)) +
  labs(
    title = title,
    x = toTitleCase(aggregations[1]),
    y = 'MAPS',
    color = toTitleCase(aggregations[2])
  ) +
  geom_hline(  # add lines for reference
    aes(yintercept = maps, group = consequence),
    data = dt_csq[, .(maps, consequence)],
    linetype = 'dashed',
    color = 'grey50'
  ) +
  geom_text(
    aes(x = 0, y = maps, label = consequence, vjust = 1.3),
    hjust = -0.05,
    data = dt_csq,
    inherit.aes = FALSE,
    size = 2.5
  ) +
  geom_text(
    aes(y = min(maps - maps_sem) - 0.03, label = variant_count),
    size = 2.5,
    angle = 30,
    position = position_dodge(dodge_width)
  ) +
  geom_errorbar(
    aes(ymin = maps - maps_sem, ymax = maps + maps_sem),
    width = 0.15,
    position = position_dodge(width = dodge_width)
  ) +
  #scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks = element_blank(),
  )

ggsave(snakemake@output$png, width = 8, height = 6)
