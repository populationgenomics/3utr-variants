"""
Plots for a single variant subset
    Input: TSV file containing MAPS scores
    Output: Plot of MAPS for all consequences
"""
# flake8: noqa

import pandas as pd
from plotnine import *  # pylint: disable=W0614,W0401

df = pd.read_table(snakemake.input[0], sep='\t')

df['UTR_group'] = df['UTR_group'] + '|'
df[['interval_set', 'annotation']] = df['UTR_group'].str.split('|', expand=True)[[0, 1]]

df = df[df['variant_count'] >= 10]
break_scale = round((df['maps'] + df['maps_sem']).max() / 5, 2)
breaks = [round(x * break_scale, 2) for x in range(-5, 6)]

chr_subset = snakemake.config['gnomAD']['subset']
title = f'{chr_subset} - {snakemake.wildcards.variant_subset}'
point_position = position_dodge(1)

maps_plot = (
    ggplot(df)
    + aes('worst_csq', 'maps', color='annotation')
    + geom_point(position=point_position)
    + geom_errorbar(
        aes(ymin='maps-maps_sem', ymax='maps+maps_sem'),
        width=0.1,
        position=point_position,
    )
    + geom_text(aes(label='variant_count'), size=12, position=point_position)
    + geom_hline(yintercept=0, color='grey')
    + facet_grid('interval_set~.')
    + ggtitle(title=title)
    + scale_y_continuous(breaks=breaks)
    + scale_color_brewer(type="qual", palette='Dark2')
    + theme_bw(base_size=16)
    + theme(
        legend_position='top',
        axis_text_x=element_text(angle=90, hjust=1),
        axis_ticks=element_blank(),
    )
)

ggsave(maps_plot, snakemake.output.maps, width=15, height=15, dpi=150)
