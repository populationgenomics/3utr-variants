"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

ALL_ANNOTATIONS = ['database', 'feature'] + \
        config['PolyA_DB']['annotation_columns'] + \
        config['PolyASite2']['annotation_columns']
GS = GSRemoteProvider()


rule count_singletons_GCP:
    # Count singeltons on Google Cloud
    input:
        intervals=rules.merge_UTR_intervals.output.intervals,
    output:
        maps=config["bucket"] + '/{chr_subset}/variant_count.tsv',
    params:
        script='workflow/scripts/count_singletons_gcp.py',
        annotate_gnomad='utr3variants/annotate_gnomad.py',
        maps='utr3variants/maps.py',
        chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset],
        annotations=ALL_ANNOTATIONS
    shell:
        """
        INTERVAL_PATH="gs://{config[bucket]}/$(basename {input})"
        gsutil cp {input} $INTERVAL_PATH
        hailctl dataproc submit {config[cluster]} \
            --pyfiles {params.annotate_gnomad},{params.maps} \
            {params.script}  \
                -o gs://{output.maps} \
                --intervals $INTERVAL_PATH \
                --annotations {params.annotations} \
                --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
                --context_ht {config[gnomAD][context_ht]} \
                --mutation_ht {config[gnomAD][mutation_rate_ht]} \
                --genome_assembly {config[genome_assembly]} \
                --chr_subset {params.chr_subset} \
                --skip_checks {config[gnomAD][skip_checks]}
        gsutil rm $INTERVAL_PATH
        """


rule prepare_gnomAD:
    # Prepare gnomAD hail table (for local run)
    threads: 10
    output:
        gnomAD_ht=directory(output_root / '{chr_subset}/gnomAD.ht')
    params:
        chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset]
    log: hail=str(output_root / 'logs/prepare_gnomAD_hail_{chr_subset}.log')
    script: '../scripts/prepare_gnomad.py'


rule count_variants_local:
    # Count variants locally
    input:
        intervals=rules.merge_UTR_intervals.output.intervals,
        gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
        counts=output_root / '{chr_subset}/local/variant_counts.tsv',
    params:
        annotations=ALL_ANNOTATIONS,
    log: hail=str(output_root / 'logs/local_{chr_subset}_counts.log')
    script: '../scripts/count_singletons.py'


def count_file(wildcards):
    """
    Get variant counts file based on wildcards
    Handles local or remote files (depending on config)
    """
    return rules.count_variants_local.output.counts \
        if wildcards.run_location == 'local' \
        else GS.remote(rules.count_singletons_GCP.output.maps,keep_local=True)


rule clean_variant_counts:
    input: counts=count_file
    output: counts=output_root / '{chr_subset}/{run_location}/variant_counts_cleaned.tsv'
    run:
        import pandas as pd

        df = pd.read_table(input.counts,sep='\t')
        if 'hexamer_motif' in ALL_ANNOTATIONS:
            df = df[df.hexamer_motif != 'Arich']  # remove A-rich hexamers
        if 'percent_expressed' in ALL_ANNOTATIONS:
            df['percent_expressed'] = pd.cut(df.percent_expressed.astype(float), 10).astype(str)
        # fill in blank values
        df['database'].fillna('gnomAD',inplace=True)
        df['feature'].fillna('other variant',inplace=True)
        for anno in [x for x in ALL_ANNOTATIONS if x not in ['feature', 'database']]:
            df[anno].fillna(df.feature,inplace=True)
        df.to_csv(output.counts,sep='\t',index=False)


rule MAPS:
    # Compute MAPS according to aggregations in wildcard and plot
    input:
        counts=rules.clean_variant_counts.output,
        utils='workflow/scripts/utils.R'
    output:
        tsv=output_root / '{chr_subset}/{run_location}/MAPS/{aggregation}.tsv',
        png=output_root / '{chr_subset}/{run_location}/MAPS/{aggregation}.png',
    params:
        chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset],
        variant_count_min=10,
    conda:
        "../envs/utr-variants-r.yml"
    script: '../scripts/maps.R'


def gather_files(wildcards, target, **kwargs):
    """
    Expand target expression, depending on wildcards.run_location
    """
    if wildcards.run_location == 'local':
        return expand(target,**kwargs)
    return GS.remote(expand(target.__str__(),**kwargs),keep_local=True)


# rule gather_MAPS:
#     # Merge different MAPS results into single table
#     input:
#         lambda wildcards: expand(
#             rules.MAPS.output.maps,
#             annotation=interval_annotations+['worst_csq'],
#             allow_missing=True
#         )
#     output: output_root / '{chr_subset}/{run_location}/MAPS.tsv'
#     run:
#         import pandas as pd
#
#         df_list = []
#         for annotation, file in zip(interval_annotations+['worst_csq'],input):
#             df = pd.read_table(file.__str__(),sep='\t')
#             df['annotation'] = annotation
#             df_list.append(df)
#         df_all = pd.concat(df_list)
#         df_all.to_csv(output[0],index=False,sep='\t')
#
#
# rule plots:
#     input: rules.MAPS.output.maps
#     output:
#         maps=output_root / '{chr_subset}/{run_location}/MAPS/{annotation}.png'
#     params:
#         chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset]
#     script: '../scripts/plots_single.py'
#
#
# rule plots_all:
#     # Plot all MAPS results in single plot
#     input:
#         maps=rules.gather_MAPS.output
#     output:
#         maps=output_root / '{chr_subset}/{run_location}/MAPS.png'
#     params:
#       chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset],
#       annotations=interval_annotations
#     script: '../scripts/plots.py'
