"""
Evaluate different 3' UTR subsets with MAPS score
Includes preparation of interval formats, gnomAD hail table, local & CGP runs
"""

rule convert_bedfile:
    input: lambda wildcards: variant_subsets[wildcards.variant_subset]
    output: output_root/'intervals/{variant_subset}_hail.txt'
    params:
        chr_style_hail = '' if config['genome_assembly'] == 'GRCh37' else 'chr'
    script: '../scripts/convert_bed.py'

rule MAPS_GCP:
    input: rules.convert_bedfile.output
    output:
        maps=output_root/'MAPS/{variant_subset}.csv'
    params:
        gnomad_prepare='scripts/prepare_gnomad.py',
        maps_score='scripts/maps_score.py'
    shell:
        """
        # TODO: manage GCP files with snakemake.remote.GS remote provider
        OUTPUT_MAPS="gs://{config[bucket]}/MAPS_{wildcards.variant_subset}.csv"
        LOG_GCP="gs://{config[bucket]}/MAPS_{wildcards.variant_subset}.log"
        hailctl dataproc submit {config[cluster]} \
            --files {input} \
            --pyfiles {params.gnomad_prepare},{params.maps_score} \
            scripts/maps_gcp.py  \
                -o $OUTPUT_MAPS \
                --intervals $(basename {input}) \
                --gnomAD_ht {config[gnomAD][gnomAD_ht]} \
                --context_ht {config[gnomAD][context_ht]} \
                --mutation_ht {config[gnomAD][mutation_rate_ht]} \
                --genome_assembly {config[genome_assembly]} \
                --log $LOG_GCP

        echo download from gs://{config[bucket]}/${{OUTPUT_GCP}}
        gsutil cp -r "gs://{config[bucket]}/${{OUTPUT_GCP}}" {output}
        """

rule prepare_gnomAD:
    threads: 10
    output:
        gnomAD_ht=directory(output_root/'gnomAD.ht')
    log: hail=str(output_root/'logs/prepare_gnomAD_hail.log')
    script: '../scripts/prepare_gnomad.py'

rule MAPS_local:
    input:
         intervals=rules.convert_bedfile.output,
         gnomAD=rules.prepare_gnomAD.output.gnomAD_ht
    output:
        maps=output_root/'MAPS/{variant_subset}_local.csv',
    log: hail=str(output_root/'logs/MAPS_local_{variant_subset}_hail.log')
    threads: 10
    script: '../scripts/maps_score.py'
