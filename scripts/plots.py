"""
Plots comparing all variant subsets
"""
# flake8: noqa

import pandas as pd
from plotnine import *  # pylint: disable=W0614,W0401

df = pd.read_csv(snakemake.input[0], sep='\t')

maps_plot = (
    ggplot(df)
    + aes('worst_csq', 'maps', color='protein_coding')
    + geom_point()
    + geom_hline(yintercept=0)
    + facet_grid('variant_subset~.')
    + ggtitle(title=snakemake.config['gnomAD']['subset'])
    + theme_classic()
    + scale_color_brewer(type="qual", palette='Set1')
    + theme(
        legend_position='top',
        axis_text_x=element_text(angle=45, hjust=1),
        axis_ticks=element_blank(),
    )
)

ggsave(maps_plot, snakemake.output['maps'], width=8, height=9)
