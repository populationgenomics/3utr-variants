"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

GS = GSRemoteProvider()

# wildcards for feature extraction rules
variant_subsets = {
    'PolyADB-overall': expand(
        rules.merge_UTR_intervals.output.intervals,
        annotation='overall'
    ),
    'PolyADB-hexamer': expand(
        rules.merge_UTR_intervals.output.intervals,
        annotation='hexamer'
    ),
    'PolyADB-conserved': expand(
        rules.merge_UTR_intervals.output.intervals,
        annotation='conserved'
    )
}


rule MAPS_GCP:
    # Compute MAPS on Google Cloud
    input: lambda wildcards: variant_subsets[wildcards.variant_subset]
    output:
        maps=GS.remote(
            f'{config["bucket"]}/MAPS_{{variant_subset}}.tsv',
            keep_local=True
        ),
    params:
        gnomad_prepare='utr3variants/annotate_gnomad.py',
        maps_score='utr3variants/maps.py'
    shell:
        """
        INTERVAL_PATH="gs://{config[bucket]}/$(basename {input})"
        gsutil cp {input} $INTERVAL_PATH
        hailctl dataproc submit {config[cluster]} \
            --pyfiles {params.gnomad_prepare},{params.maps_score} \
            scripts/maps_gcp.py  \
                -o gs://{output.maps} \
                --intervals $INTERVAL_PATH \
                --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
                --context_ht {config[gnomAD][context_ht]} \
                --mutation_ht {config[gnomAD][mutation_rate_ht]} \
                --genome_assembly {config[genome_assembly]} \
                --chr_subset {config[gnomAD][subset]}
        gsutil rm $INTERVAL_PATH
        """


rule prepare_gnomAD:
    # Prepare gnomAD hail table (for local run)
    threads: 10
    output:
        gnomAD_ht=directory(output_root / 'gnomAD.ht')
    log: hail=str(output_root / 'logs/prepare_gnomAD_hail.log')
    script: '../scripts/prepare_gnomad.py'


rule MAPS_local:
    # Compute MAPS locally
    input:
        intervals=lambda wildcards: variant_subsets[wildcards.variant_subset],
        gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
        maps=output_root / 'MAPS/local/{variant_subset}.tsv',
    log: hail=str(output_root / 'logs/MAPS_local_{variant_subset}_hail.log')
    threads: 10
    script: '../scripts/maps_score.py'


def maps_file(wildcards):
    """
    Get MAPS file based on wildcards defined in variant_subsets
    Handles local or remote files (depending on config)
    """
    return rules.MAPS_local.output.maps \
        if wildcards.run_location == 'local' \
        else rules.MAPS_GCP.output.maps


def gather_files(wildcards, target, **kwargs):
    """
    Expand target expression, depending on wildcards.run_location
    """
    subsets = variant_subsets.keys()
    if wildcards.run_location == 'local':
        return expand(target,**kwargs)
    return GS.remote(expand(target.__str__(),**kwargs))


rule gather_MAPS:
    # Merge different MAPS results into single table
    input: lambda w: gather_files(w,maps_file(w),variant_subset=variant_subsets.keys())
    output: output_root / 'MAPS/all_{run_location}.tsv'
    run:
        import pandas as pd

        df_list = []
        for variant_subset, file in zip(variant_subsets.keys(),input):
            df = pd.read_table(file.__str__(),sep='\t')
            df['variant_subset'] = variant_subset
            df_list.append(df)
        df_all = pd.concat(df_list)
        df_all.to_csv(output[0],index=False,sep='\t')


rule plot_single:
    input: maps_file
    output:
        maps=output_root / 'plots/{run_location}/MAPS_{variant_subset}.png'
    script: '../scripts/plots_single.py'


rule plots:
    # Plot all MAPS results in single plot
    input:
        maps=rules.gather_MAPS.output,
        single_plots=expand(
            rules.plot_single.output,
            variant_subset=variant_subsets.keys(),
            allow_missing=True
        )
    output:
        maps=output_root / 'plots/MAPS_{run_location}.png'
    script: '../scripts/plots.py'
