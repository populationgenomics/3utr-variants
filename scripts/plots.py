"""
Plots comparing all variant subsets
    Input: TSV file containing MAPS scores
    Output: Plot of MAPS for all variant subsets
"""
# flake8: noqa

import pandas as pd
from plotnine import *  # pylint: disable=W0614,W0401

df = pd.read_table(snakemake.input[0], sep='\t')

df = df[df.worst_csq == '3_prime_UTR_variant']
df[['interval_set', 'annotation']] = df['UTR_group'].str.split('|', expand=True)

chr_subset = snakemake.config['gnomAD']['subset']
title = f'{chr_subset} worst_csq: 3\'UTR variant'
point_position = position_dodge(1)

maps_plot = (
    ggplot(df)
    + aes('interval_set', 'maps', color='annotation')
    + geom_point(position=point_position)
    + geom_errorbar(
        aes(ymin='maps-maps_sem', ymax='maps+maps_sem'),
        width=0.1,
        position=point_position,
    )
    + geom_text(aes(label='variant_count'), position=point_position)
    + geom_hline(yintercept=0, color='grey')
    + facet_grid('variant_subset~.')
    + ggtitle(title=title)
    + scale_color_brewer(type="qual", palette='Dark2')
    + theme_classic()
    + theme(
        legend_position='right',
        # axis_text_x=element_text(angle=90, hjust=1),
        axis_ticks=element_blank(),
    )
)

ggsave(maps_plot, snakemake.output['maps'], width=8, height=9)
