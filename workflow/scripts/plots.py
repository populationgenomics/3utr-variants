"""
Plots comparing all variant subsets
    Input: TSV file containing MAPS scores
    Output: Plot of MAPS for all variant subsets
"""
# flake8: noqa

import pandas as pd
from plotnine import *  # pylint: disable=W0614,W0401

annotations = snakemake.params.annotations

df = pd.read_table(snakemake.input[0], sep='\t')
df['annotation_value'] = df[annotations].apply(lambda x: '|'.join(x.dropna()), 1)
df['dummy'] = 'dummy'

maps_syn = df[df.worst_csq == 'synonymous_variant']['maps']
maps_miss = df[df.worst_csq == 'missense_variant']['maps']
maps_stop = df[df.worst_csq == 'stop_gained']['maps']

# sdf[df.worst_csq.isin('synonymous_variant', 'missense_variant', 'stop_gained')]['annotation_value'] = 'worst_csq'
worst_csq_ref = ['synonymous_variant', 'missense_variant', 'stop_gained']
# df.annotation_value = 'worst_csq' if df.worst_csq.isin(worst_csq_ref) else df.annotation_value
# df.annotation_value = df.annotation_value.where(df.worst_csq.isin(worst_csq_ref), df.annotation_value)
df.loc[df.worst_csq.isin(worst_csq_ref), 'annotation_value'] = 'worst_csq'

df = df[df.annotation_value != '']

chr_subset = snakemake.params['chr_subset']
title = f'{chr_subset}'
point_position = position_dodge(1)

maps_plot = (
    ggplot(df)
    + aes('dummy', 'maps', color='annotation_value')
    + geom_point(position=point_position)
    + geom_errorbar(
        aes(ymin='maps-maps_sem', ymax='maps+maps_sem'),
        width=0.1,
        position=point_position,
    )
    # + geom_text(aes(label='variant_count'), position=point_position)
    + geom_hline(yintercept=0, color='grey')
    + geom_hline(yintercept=maps_syn, color='grey')
    + geom_hline(yintercept=maps_miss, color='orange')
    + geom_hline(yintercept=maps_stop, color='red')
    + facet_grid('~annotation')
    + ggtitle(title=title)
    # + scale_color_brewer(type="qual", palette='Dark2')
    + theme_bw()
    + theme(
        legend_position='bottom',
        # axis_text_x=element_text(angle=90, hjust=1),
        axis_ticks=element_blank(),
    )
)

ggsave(maps_plot, snakemake.output['maps'], width=15, height=5)
