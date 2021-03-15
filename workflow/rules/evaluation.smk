"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

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
        annotations=interval_annotations
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
                --chr_subset {params.chr_subset}
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
        counts=output_root / '{chr_subset}/MAPS/local/variant_counts.tsv',
    params:
        annotations=interval_annotations
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


rule MAPS:
    # Compute MAPS locally
    input:
        counts=count_file
    output:
        maps=output_root / '{chr_subset}/MAPS/{run_location}/MAPS_{variant_subset}.tsv',
    log:
        hail=str(
            output_root / 'logs/{run_location}_{chr_subset}_{variant_subset}_MAPS.log'
        )
    script: '../scripts/maps.py'


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
#         lambda wildcards: gather_files(
#             wildcards,
#             rules.MAPS.output.maps,
#             variant_subset=variant_subsets.keys(),
#             allow_missing=True
#         )
#     output: output_root / '{chr_subset}/MAPS/all_{run_location}.tsv'
#     run:
#         import pandas as pd
#
#         df_list = []
#         for variant_subset, file in zip(variant_subsets.keys(),input):
#             df = pd.read_table(file.__str__(),sep='\t')
#             df['variant_subset'] = variant_subset
#             df_list.append(df)
#         df_all = pd.concat(df_list)
#         df_all.to_csv(output[0],index=False,sep='\t')


rule plots:
    input: rules.MAPS.output.maps
    output:
        maps=output_root / '{chr_subset}/plots/{run_location}/MAPS_{variant_subset}.png'
    params:
        chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset]
    script: '../scripts/plots_single.py'


# rule plots:
#     # Plot all MAPS results in single plot
#     input:
#         #maps=rules.gather_MAPS.output,
#         single_plots=expand(
#             rules.plot_single.output,
#             variant_subset=interval_annotations,
#             allow_missing=True
#         )
#     #output:
#         #maps=output_root / '{chr_subset}/plots/{run_location}/MAPS_all.png'
#     #params:
#     #   chr_subset=lambda wildcards: config['chr_subsets'][wildcards.chr_subset]
#     script: '../scripts/plots.py'
